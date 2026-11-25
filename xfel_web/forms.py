"""Turn an HTML form into the dictionaries the map pipelines expect."""

from __future__ import annotations

import math
import re
import subprocess
from pathlib import Path

from werkzeug.datastructures import FileStorage
from werkzeug.utils import secure_filename

from xfel_web.ccp4 import fft_grid_resolution

ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.+\-]{0,80}$")
LABEL_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]{0,20}$")
# PDB atom names, including primes written as * or ' and charged names such as N1+.
_ATOM_NAME_RE = re.compile(r"^[A-Za-z0-9*'+-]{1,4}$")

PHASE_SOURCES = ("light_validation", "dark_validation", "light_refmac", "dark_refmac")


class FormError(ValueError):
    pass


def _text(form, name: str) -> str:
    return (form.get(name) or "").strip()


def _flag(form, name: str) -> bool:
    return _text(form, name).lower() in {"1", "true", "on", "yes"}


def _float(form, name: str, label: str) -> float:
    raw = _text(form, name)
    try:
        return float(raw)
    except ValueError as exc:
        raise FormError(f"{label} must be a number.") from exc


def _int(form, name: str, label: str) -> int:
    raw = _text(form, name)
    try:
        return int(raw)
    except ValueError as exc:
        raise FormError(f"{label} must be an integer.") from exc


def _identifier(form, name: str, label: str) -> str:
    value = _text(form, name)
    if not ID_RE.match(value):
        raise FormError(f"{label} should be a short name such as 6w1v or 6w1v_omitO6.")
    return value


def _pdb_coordinates(path: Path, label: str) -> Path:
    if path.suffix.lower() not in {".pdb", ".ent"}:
        raise FormError(f"{label} must be a PDB coordinate file.")
    return path


def _harmonic_sigma(form) -> float:
    raw = _text(form, "harmonic_sigma")
    sigma = 0.1 if not raw else _float(form, "harmonic_sigma", "Harmonic sigma")
    if (_flag(form, "harmonic_added") or _flag(form, "use_harmonic_atoms")) and sigma <= 0:
        raise FormError("Harmonic sigma must be positive.")
    return sigma


def _harmonic_protein_sigma(form) -> float:
    raw = _text(form, "harmonic_protein_sigma")
    sigma = 0.1 if not raw else _float(form, "harmonic_protein_sigma", "Protein harmonic sigma")
    if sigma <= 0:
        raise FormError("Protein harmonic sigma must be positive.")
    return sigma


def _protein_harmonic(form, refine: bool) -> tuple[bool, float]:
    """Whether protein atoms outside the site are restrained, and the sigma.

    The switch is read only while REFMAC refinement is on. Off keeps the
    default sigma and applies nothing.
    """
    if not refine or not _flag(form, "harmonic_protein"):
        return False, 0.1
    return True, _harmonic_protein_sigma(form)


def _candidate_atoms(form) -> str:
    """Which reference atoms may be placed on an omit-map peak."""
    value = _text(form, "candidate_atoms") or "missing"
    if value not in {"missing", "all"}:
        raise FormError(
            "Atoms eligible for placement must be the heavy atoms absent from the omit model, "
            "or every heavy atom in the reference model."
        )
    return value


def _odd_neighborhood(form) -> int:
    size = _int(form, "neighborhood_size", "Neighborhood size")
    if size < 1 or size % 2 == 0:
        raise FormError("Neighborhood size must be a positive odd integer.")
    return size


def _map_center_fields(
    form,
    mode_key: str,
    x_key: str,
    y_key: str,
    z_key: str,
    label: str,
) -> tuple[str, tuple[float, float, float] | None]:
    """How a clip is centered: the site residue, or an explicit point."""
    mode = _text(form, mode_key) or "residue"
    if mode not in {"residue", "cartesian"}:
        raise FormError(f"{label} must be the site residue or a Cartesian point.")
    if mode == "residue":
        return mode, None
    point = (
        _float(form, x_key, f"{label} X"),
        _float(form, y_key, f"{label} Y"),
        _float(form, z_key, f"{label} Z"),
    )
    if not all(math.isfinite(value) for value in point):
        raise FormError(f"{label} coordinates must be finite numbers, in ångströms.")
    return mode, point


def _map_center(form) -> tuple[str, tuple[float, float, float] | None]:
    """How the omit-map clip is centered: the site residue, or an explicit point."""
    return _map_center_fields(form, "map_center", "center_x", "center_y", "center_z", "Map center")


def _levels(text: str, label: str) -> list[float]:
    parts = text.replace(",", " ").split()
    if not parts:
        raise FormError(f"{label} needs at least one number.")
    try:
        return [float(part) for part in parts]
    except ValueError as exc:
        raise FormError(f"{label} must be a list of numbers.") from exc


def viewer_sigma_bounds(minimum: float, maximum: float) -> tuple[float, float]:
    """Inclusive slider limits for the one isosurface in the interactive viewer."""
    if not math.isfinite(minimum) or not math.isfinite(maximum):
        raise FormError("Viewer sigma limits must be finite numbers.")
    if not minimum < maximum:
        raise FormError("Viewer sigma minimum must be lower than the maximum.")
    return minimum, maximum


def parse_viewer_sigma(form) -> tuple[float, float]:
    return viewer_sigma_bounds(
        _float(form, "viewer_sigma_min", "Viewer sigma minimum"),
        _float(form, "viewer_sigma_max", "Viewer sigma maximum"),
    )


def _optional_levels(text: str, fallback: list[float], label: str) -> list[float]:
    if not text.strip():
        return list(fallback)
    levels = _levels(text, label)
    if any(level <= 0 for level in levels):
        raise FormError(f"{label} uses positive contour levels. The negative contour is drawn automatically.")
    return sorted(levels)


def parse_extra_residues(text: str, *, what: str = "Extra residue") -> list[dict]:
    extras = []
    for line_no, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.replace(",", " ").split()
        if len(parts) < 2:
            raise FormError(f"{what} line {line_no} needs a chain and a residue number.")
        chain, resi_text = parts[0], parts[1]
        if len(chain) != 1:
            raise FormError(f"{what} line {line_no}: chain must be one character.")
        try:
            resi = int(resi_text)
        except ValueError as exc:
            raise FormError(f"{what} line {line_no}: residue number is not an integer.") from exc
        label = parts[2] if len(parts) > 2 else f"{chain}{resi}"
        if not LABEL_RE.match(label):
            raise FormError(
                f"{what} label {label!r} must start with a letter and use only letters, digits, and underscores."
            )
        extras.append({"chain": chain, "resi": resi, "label": label})
    return extras


def _site_chain_value(text: str, label: str) -> str:
    chain = text.strip()
    if len(chain) != 1 or not chain.isalnum():
        raise FormError(f"{label} must be a single letter or digit, as stored in the PDB file.")
    return chain


def _site_label_value(text: str, label: str) -> str:
    value = text.strip() or "site"
    if not LABEL_RE.match(value):
        raise FormError(f"{label} must start with a letter and use only letters, digits, and underscores.")
    return value


def additional_site_from_form(form) -> dict | None:
    """The optional second site, or None when its switch is off."""
    if not _flag(form, "use_additional_site"):
        return None
    chain = _site_chain_value(_text(form, "alt_site_chain"), "Additional site chain")
    resi = _int(form, "alt_site_resi", "Additional site residue")
    label = _site_label_value(_text(form, "alt_site_label"), "Additional site label")
    if "alt_map_center" in form:
        mode, point = _map_center_fields(
            form,
            "alt_map_center",
            "alt_center_x",
            "alt_center_y",
            "alt_center_z",
            "Additional map center",
        )
    else:
        mode, point = "residue", None
    return {
        "chain": chain,
        "resi": resi,
        "label": label,
        "map_center": mode,
        "map_center_xyz": point,
        "extra_residues": parse_extra_residues(
            _text(form, "alt_extra_residues"),
            what="Additional neighbouring residue",
        ),
    }


def _values(form, name: str) -> list[str]:
    """Every submitted value for one field, including repeated row inputs."""
    if hasattr(form, "getlist"):
        values = form.getlist(name)
    else:
        raw = form.get(name, [])
        if raw is None or raw == "":
            values = []
        elif isinstance(raw, (list, tuple)):
            values = list(raw)
        else:
            values = [raw]
    return [str(item).strip() for item in values]


def _form_rows(form, names: tuple[str, ...]) -> list[tuple[str, ...]]:
    columns = [_values(form, name) for name in names]
    width = max((len(column) for column in columns), default=0)
    rows = []
    for index in range(width):
        rows.append(tuple(column[index] if index < len(column) else "" for column in columns))
    return rows


def _atom_field(text: str, label: str) -> dict | None:
    """``chain residue atom``, or None when the field is blank."""
    if not text:
        return None
    parts = text.replace(",", " ").split()
    if len(parts) != 3:
        raise FormError(f"{label} should be 'chain residue atom', for example 'A 602 MN1'.")
    chain, resi_text, name = parts
    if len(chain) != 1 or not chain.isalnum():
        raise FormError(f"{label}: chain must be a single letter or digit.")
    try:
        resi = int(resi_text)
    except ValueError as exc:
        raise FormError(f"{label}: residue number is not an integer.") from exc
    if not _ATOM_NAME_RE.match(name):
        raise FormError(f"{label}: atom {name!r} is not a PDB atom name.")
    return {"chain": chain, "resi": resi, "name": name.upper()}


def _optional_number(text: str, label: str):
    if not text:
        return None
    try:
        return float(text)
    except ValueError as exc:
        raise FormError(f"{label} must be a number.") from exc


def parse_distance_restraints(rows: list[tuple[str, ...]]) -> list[dict]:
    """Rows of two atoms, an optional distance, and an optional sigma."""
    specs = []
    for index, row in enumerate(rows, start=1):
        atom1_text, atom2_text, value_text, sigma_text = (list(row) + ["", "", "", ""])[:4]
        label = f"Distance row {index}"
        atom1 = _atom_field(atom1_text, f"{label}, first atom")
        atom2 = _atom_field(atom2_text, f"{label}, second atom")
        if atom1 is None and atom2 is None and not value_text and not sigma_text:
            continue
        if atom1 is None or atom2 is None:
            raise FormError(f"{label} needs two atoms, each as 'chain residue atom'.")
        specs.append(
            {
                "atom1": atom1,
                "atom2": atom2,
                "distance": _optional_number(value_text, f"{label} distance"),
                "sigma": _optional_number(sigma_text, f"{label} sigma"),
            }
        )
    return specs


def parse_angle_restraints(rows: list[tuple[str, ...]]) -> list[dict]:
    """Rows of three atoms, an optional angle, and an optional sigma."""
    specs = []
    for index, row in enumerate(rows, start=1):
        texts = (list(row) + ["", "", "", "", ""])[:5]
        label = f"Angle row {index}"
        atoms = [
            _atom_field(texts[0], f"{label}, first atom"),
            _atom_field(texts[1], f"{label}, vertex"),
            _atom_field(texts[2], f"{label}, third atom"),
        ]
        if all(atom is None for atom in atoms) and not texts[3] and not texts[4]:
            continue
        if any(atom is None for atom in atoms):
            raise FormError(f"{label} needs three atoms, each as 'chain residue atom'.")
        specs.append(
            {
                "atom1": atoms[0],
                "vertex": atoms[1],
                "atom3": atoms[2],
                "angle_deg": _optional_number(texts[3], f"{label} angle"),
                "sigma_deg": _optional_number(texts[4], f"{label} sigma"),
            }
        )
    return specs


def parse_dihedral_restraints(rows: list[tuple[str, ...]]) -> list[dict]:
    """Rows of four atoms, an optional angle, and an optional sigma."""
    specs = []
    for index, row in enumerate(rows, start=1):
        texts = (list(row) + ["", "", "", "", "", ""])[:6]
        label = f"Dihedral row {index}"
        atoms = [_atom_field(texts[offset], f"{label}, atom {offset + 1}") for offset in range(4)]
        if all(atom is None for atom in atoms) and not texts[4] and not texts[5]:
            continue
        if any(atom is None for atom in atoms):
            raise FormError(f"{label} needs four atoms, each as 'chain residue atom'.")
        specs.append(
            {
                "atoms": atoms,
                "angle_deg": _optional_number(texts[4], f"{label} angle"),
                "sigma_deg": _optional_number(texts[5], f"{label} sigma"),
            }
        )
    return specs


def parse_exclude_restraints(text: str) -> list[tuple[str, int, tuple[str, ...]]]:
    """Lines of ``chain residue atoms``. ``*`` means every atom in the residue."""
    specs = []
    for line_no, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.replace(",", " ").split()
        if len(parts) < 3:
            raise FormError(
                f"Exclude-restraint line {line_no} should be 'chain residue atoms', "
                "for example 'A 416 *' or 'A 189 CA CB CG'."
            )
        chain, resi_text, *atoms = parts
        if len(chain) != 1 or not chain.isalnum():
            raise FormError(
                f"Exclude-restraint line {line_no}: chain must be a single letter or digit."
            )
        try:
            resi = int(resi_text)
        except ValueError as exc:
            raise FormError(
                f"Exclude-restraint line {line_no}: residue number is not an integer."
            ) from exc
        if any(atom == "*" for atom in atoms) and atoms != ["*"]:
            raise FormError(
                f"Exclude-restraint line {line_no}: '*' means every atom and must be the only atom token."
            )
        normalized = []
        for atom in atoms:
            if atom != "*" and not _ATOM_NAME_RE.match(atom):
                raise FormError(
                    f"Exclude-restraint line {line_no}: atom {atom!r} is not a PDB atom name."
                )
            normalized.append("*" if atom == "*" else atom.upper())
        specs.append((chain, resi, tuple(normalized)))
    return specs


def parse_harmonic_atoms(text: str) -> list[dict]:
    """Lines of ``chain residue atom``, with an optional altloc."""
    atoms = []
    seen: set[tuple] = set()
    for line_no, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.replace(",", " ").split()
        if len(parts) not in {3, 4}:
            raise FormError(
                f"Harmonic-atom line {line_no} should be 'chain residue atom', "
                "for example 'A 416 MN1' or 'A 416 O5 B'."
            )
        chain, resi_text, name, *rest = parts
        if len(chain) != 1 or not chain.isalnum():
            raise FormError(
                f"Harmonic-atom line {line_no}: chain must be a single letter or digit."
            )
        try:
            resi = int(resi_text)
        except ValueError as exc:
            raise FormError(
                f"Harmonic-atom line {line_no}: residue number is not an integer."
            ) from exc
        if not _ATOM_NAME_RE.match(name):
            raise FormError(
                f"Harmonic-atom line {line_no}: atom {name!r} is not a PDB atom name."
            )
        alt = rest[0] if rest else ""
        if alt and (len(alt) != 1 or not alt.isalnum()):
            raise FormError(
                f"Harmonic-atom line {line_no}: altloc must be a single letter or digit."
            )
        atom = {
            "chain": chain,
            "resi": resi,
            "name": name.upper(),
            "alt": alt.upper(),
            "icode": "",
        }
        key = (atom["chain"], atom["resi"], atom["name"], atom["alt"])
        if key in seen:
            continue
        seen.add(key)
        atoms.append(atom)
    return atoms


def _specified_harmonic(form, refine: bool) -> tuple[list[dict], float]:
    """Atoms the omit-model form holds with an external harmonic restraint."""
    if not refine or not _flag(form, "use_harmonic_atoms"):
        return [], 0.1
    atoms = parse_harmonic_atoms(_text(form, "harmonic_atom_lines"))
    if not atoms:
        raise FormError("Harmonic restraints need at least one atom, such as 'A 416 MN1'.")
    return atoms, _harmonic_sigma(form)


def _restraint_block(form, flag: str, parser, field_names: tuple[str, ...], sigma_name: str, sigma_label: str, default_sigma: float):
    """Restraints for one switch. Off keeps the default sigma and applies nothing."""
    if not _flag(form, flag):
        return [], default_sigma
    try:
        specs = parser(_form_rows(form, field_names))
    except ValueError as exc:
        raise FormError(str(exc)) from exc
    sigma = _float(form, sigma_name, sigma_label)
    if sigma <= 0:
        raise FormError(f"{sigma_label} must be positive.")
    return specs, sigma


def _active_restraints(form, refine: bool):
    """External restraints are parsed only while REFMAC refinement is on."""
    if not refine:
        return [], 0.5, [], 25.0, [], 10.0
    distances, distance_sigma = _restraint_block(
        form, "use_distance_restraints", parse_distance_restraints,
        ("distance_atom_1", "distance_atom_2", "distance_value", "distance_row_sigma"),
        "distance_sigma", "Distance sigma", 0.5,
    )
    angles, angle_sigma = _restraint_block(
        form, "use_angle_restraints", parse_angle_restraints,
        ("angle_atom_1", "angle_atom_2", "angle_atom_3", "angle_value", "angle_row_sigma"),
        "angle_sigma", "Angle sigma", 25.0,
    )
    dihedrals, dihedral_sigma = _restraint_block(
        form, "use_dihedral_restraints", parse_dihedral_restraints,
        (
            "dihedral_atom_1", "dihedral_atom_2", "dihedral_atom_3", "dihedral_atom_4",
            "dihedral_value", "dihedral_row_sigma",
        ),
        "dihedral_sigma", "Dihedral sigma", 10.0,
    )
    return distances, distance_sigma, angles, angle_sigma, dihedrals, dihedral_sigma


def _exclude_restraints(form, refine: bool) -> list[tuple[str, int, tuple[str, ...]]]:
    """Residues left out of REFMAC restraints. The switch is read only while REFMAC is on."""
    if not refine or not _flag(form, "use_exclude_restraints"):
        return []
    return parse_exclude_restraints(_text(form, "exclude_restraints"))


def _link_record_file(form, files, output_dir: Path, project: Path) -> Path | None:
    """External LINK records for reconstruction. None uses the starting omit model.

    The file is read only while the switch is on. Off ignores a path or upload.
    """
    if not _flag(form, "restore_links"):
        return None
    path = resolve_optional_input_file(
        form, files, "link_records", "LINK record file", output_dir, project
    )
    if path is None:
        return None
    from xfel_web.reconstruction import read_link_records

    if not read_link_records(path):
        raise FormError(
            "The LINK record file contains no LINK records. "
            "Leave the field empty to use the LINK records in the starting omit model."
        )
    return path


def _custom_dictionary(form, files, output_dir: Path, project: Path):
    if not _flag(form, "use_custom_dictionary"):
        return None
    return resolve_optional_input_file(
        form, files, "custom_dictionary", "custom dictionary", output_dir, project
    )


def resolve_output_dir(text: str, project: Path) -> Path:
    if not text.strip():
        raise FormError("Choose an output directory.")
    if any(char in text for char in "\n\r\0"):
        raise FormError("The output directory contains a newline.")
    path = Path(text.strip()).expanduser()
    if not path.is_absolute():
        path = project / path
    path = path.resolve()
    if '"' in str(path) or "\\" in str(path):
        raise FormError("Paths cannot contain quotes or backslashes.")
    if path == Path(path.anchor) or len(path.parts) < 3:
        raise FormError("Choose a specific output directory, not the top of the disk.")
    blocked = {
        project.resolve(),
        (project / "data").resolve(),
        (project / "web-application").resolve(),
    }
    if path in blocked or (project / "data").resolve() in path.parents:
        raise FormError("Put the output in its own folder. The project directory and the data folder are not used as output.")
    if path.exists() and not path.is_dir():
        raise FormError("The output path exists and is not a directory.")
    return path


def _has_gz_extension(path: Path) -> bool:
    return path.name.lower().endswith(".gz")


def decompress_gzip(path: Path, dest_dir: Path, label: str, dest_name: str | None = None) -> Path:
    """Return the file unchanged, or a gunzip copy when the name ends in `.gz`."""
    if not _has_gz_extension(path):
        return path
    inner = (dest_name or path.name[:-3]).strip()
    name = secure_filename(inner) or "input"
    if name.lower().endswith(".gz"):
        name = name[:-3] or "input"
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / name
    if dest.resolve() == path.resolve():
        dest = dest_dir / f"plain_{name}"
    owned_copy = path.resolve().parent == dest_dir.resolve() and path.resolve() != dest.resolve()
    try:
        with dest.open("wb") as handle:
            subprocess.run(
                ["gunzip", "-c", str(path)],
                stdout=handle,
                stderr=subprocess.PIPE,
                check=True,
            )
    except FileNotFoundError as exc:
        dest.unlink(missing_ok=True)
        raise FormError("gunzip is not available, so a .gz input cannot be decompressed.") from exc
    except subprocess.CalledProcessError as exc:
        dest.unlink(missing_ok=True)
        detail = (exc.stderr or b"").decode("utf-8", "replace").strip()
        message = f"{label} could not be decompressed with gunzip."
        if detail:
            message = f"{message} {detail}"
        raise FormError(message) from exc
    if not dest.is_file() or dest.stat().st_size == 0:
        dest.unlink(missing_ok=True)
        raise FormError(f"{label} decompressed to an empty file.")
    if owned_copy:
        path.unlink(missing_ok=True)
    return dest


def resolve_optional_input_file(form, files, key: str, label: str, output_dir: Path, project: Path) -> Path | None:
    """A file the form may leave blank. A filled path or upload must exist."""
    upload = files.get(f"file_{key}")
    has_upload = isinstance(upload, FileStorage) and bool(upload.filename)
    if not has_upload and not _text(form, f"path_{key}"):
        return None
    path = resolve_input_file(form, files, key, label, output_dir, project)
    if not path.is_file() or path.stat().st_size == 0:
        raise FormError(f"{label} was not found or is empty.")
    return path


def resolve_input_file(form, files, key: str, label: str, output_dir: Path, project: Path) -> Path:
    folder = output_dir / "inputs"
    upload = files.get(f"file_{key}")
    if isinstance(upload, FileStorage) and upload.filename:
        folder.mkdir(parents=True, exist_ok=True)
        filename = secure_filename(upload.filename) or key
        dest = folder / f"{key}_{filename}"
        upload.save(dest)
        if dest.stat().st_size == 0:
            raise FormError(f"{label} was uploaded empty.")
        if _has_gz_extension(dest):
            return decompress_gzip(dest, folder, label)
        return dest
    raw = _text(form, f"path_{key}")
    if not raw:
        raise FormError(f"Choose a file or enter a path for {label}.")
    path = Path(raw).expanduser()
    if not path.is_absolute():
        path = project / path
    path = path.resolve()
    if '"' in str(path) or "\\" in str(path):
        raise FormError(f"The path for {label} cannot contain quotes or backslashes.")
    if _has_gz_extension(path):
        return decompress_gzip(path, folder, label, dest_name=f"{key}_{path.name}")
    return path


_BREF_TOKENS = {"ISOT", "ANIS", "MIXE", "OVER"}
_TYPE_TOKENS = {"REST", "UNRE"}


def _resolution_limits(form, source: str) -> tuple[float | None, float | None]:
    """Explicit cutoffs are required only when that source is selected."""
    if source != "explicit":
        return None, None
    lowres = _float(form, "lowres", "Low-resolution cutoff")
    dmin = _float(form, "dmin", "High-resolution cutoff")
    if not (lowres > dmin > 0):
        raise FormError(
            "The low-resolution cutoff must be a larger number of ångströms than the high-resolution cutoff."
        )
    return lowres, dmin


def _fft_grid(form, source: str, lowres: float | None, dmin: float | None) -> float | None:
    """A custom FFT grid is read only while its switch is on."""
    if not _flag(form, "use_fft_grid"):
        return None
    value = _float(form, "fft_dmin", "FFT grid resolution for the final map")
    if source == "explicit":
        try:
            fft_grid_resolution(lowres, dmin, value)
        except ValueError as exc:
            raise FormError(str(exc)) from exc
    return value


def _refmac_settings(form) -> tuple[int, str, str, bool]:
    """Cycles, ``REFI BREF``, ``REFI TYPE``, and whether the section is on.

    The section off means zero cycles. The keyword defaults are still stored
    so a later zero-cycle job has a complete protocol.
    """
    if not _flag(form, "use_refmac"):
        return 0, "ISOT", "REST", False
    cycles = _int(form, "refmac_ncyc", "REFMAC cycles")
    if cycles < 0:
        raise FormError("REFMAC cycles cannot be negative. Use 0 to calculate the map without moving atoms.")
    bref = (_text(form, "refmac_bref") or "ISOT").upper()
    if bref not in _BREF_TOKENS:
        raise FormError("B-factor refinement must be ISOTropic, ANISotropic, MIXEd, or OVERall.")
    kind = (_text(form, "refmac_type") or "REST").upper()
    if kind not in _TYPE_TOKENS:
        raise FormError("Refinement type must be RESTrained or UNREstrained.")
    return cycles, bref, kind, True


def _initial_refmac_cycles(form, refine: bool, step_cycles: int) -> int:
    """Cycles for the starting omit model. A blank field follows the later steps.

    The field is read only while REFMAC refinement is on. Off means zero cycles.
    """
    if not refine:
        return 0
    if not _text(form, "initial_refmac_ncyc"):
        return step_cycles
    cycles = _int(
        form,
        "initial_refmac_ncyc",
        "Number of refinement cycles for the initial omit model",
    )
    if cycles < 0:
        raise FormError(
            "Number of refinement cycles for the initial omit model cannot be negative. "
            "Use 0 to calculate that map without moving atoms."
        )
    return cycles


def _common(form, project: Path) -> dict:
    output_dir = resolve_output_dir(_text(form, "output_dir"), project)
    ccp4_root = Path(_text(form, "ccp4_root") or "/Applications/ccp4-9").expanduser()
    resolution_source = _text(form, "resolution_source") or "sf"
    if resolution_source not in {"sf", "explicit"}:
        raise FormError("Resolution source must be the structure-factor file or explicit cutoffs.")
    lowres, dmin = _resolution_limits(form, resolution_source)
    fft_dmin = _fft_grid(form, resolution_source, lowres, dmin)
    chain = _site_chain_value(_text(form, "site_chain"), "Site chain")
    resi = _int(form, "site_resi", "Site residue")
    label = _site_label_value(_text(form, "site_label"), "Site label")
    contour_mode = _text(form, "contour_mode") or "sigma"
    surface = _text(form, "surface") or "mesh"
    if contour_mode not in {"sigma", "absolute"}:
        raise FormError("Contour mode must be sigma or absolute.")
    if surface not in {"mesh", "solid"}:
        raise FormError("Surface representation must be mesh or solid.")
    viewer_sigma_min, viewer_sigma_max = parse_viewer_sigma(form)
    png_width = _int(form, "png_width", "Figure width")
    png_height = _int(form, "png_height", "Figure height")
    if not (200 <= png_width <= 4000 and 200 <= png_height <= 4000):
        raise FormError("Figure width and height must be between 200 and 4000 pixels.")
    return {
        "output_dir": output_dir,
        "ccp4_root": ccp4_root,
        "resolution_source": resolution_source,
        "lowres": lowres,
        "dmin": dmin,
        "fft_dmin": fft_dmin,
        "map_border": _float(form, "map_border", "Extended map border"),
        "site_chain": chain,
        "site_resi": resi,
        "site_label": label,
        "site_copies": [(chain, resi)],
        "additional_site": additional_site_from_form(form),
        "extra_residues": parse_extra_residues(_text(form, "extra_residues")),
        "export_radius": _float(form, "export_radius", "Export radius"),
        "carve": _float(form, "carve", "Carve radius"),
        "contour_mode": contour_mode,
        "surface": surface,
        "viewer_sigma_min": viewer_sigma_min,
        "viewer_sigma_max": viewer_sigma_max,
        "show_oec_bonds": _flag(form, "show_oec_bonds"),
        "draw_distances": _flag(form, "draw_distances"),
        "png_width": png_width,
        "png_height": png_height,
        "ray_trace": _flag(form, "ray_trace"),
        "recompute_map": _flag(form, "recompute_map"),
        "recompute_fragment": _flag(form, "recompute_fragment"),
    }


def _difference_levels(form) -> tuple[list[float], list[float]]:
    levels = _levels(_text(form, "sigma_levels"), "Sigma levels")
    if any(level <= 0 for level in levels):
        raise FormError("Sigma levels for a difference map are positive. Negative density is drawn at the same magnitudes.")
    positive = _optional_levels(_text(form, "sigma_positive"), levels, "Positive sigma levels")
    negative = _optional_levels(_text(form, "sigma_negative"), levels, "Negative sigma levels")
    return sorted(positive), sorted(negative)


def parse_job(workflow: str, form, files, project: Path) -> dict:
    cfg = _common(form, project)
    output_dir = cfg["output_dir"]

    def take(key: str, label: str) -> Path:
        return resolve_input_file(form, files, key, label, output_dir, project)

    if workflow == "fofo":
        positive, negative = _difference_levels(form)
        phase = _text(form, "phase_source") or "light_validation"
        if phase not in PHASE_SOURCES:
            raise FormError("Unknown phase source.")
        cfg.update(
            {
                "dark_id": _identifier(form, "dark_id", "Dark dataset name"),
                "light_id": _identifier(form, "light_id", "Light dataset name"),
                "phase_source": phase,
                "positive_levels": positive,
                "negative_levels": negative,
                "files": {
                    "dark_pdb": take("dark_pdb", "dark coordinates"),
                    "light_pdb": take("light_pdb", "light coordinates"),
                    "dark_sf": take("dark_sf", "dark structure factors"),
                    "light_sf": take("light_sf", "light structure factors"),
                    "dark_validation": take("dark_validation", "dark validation coefficients"),
                    "light_validation": take("light_validation", "light validation coefficients"),
                },
            }
        )
    elif workflow == "refinement":
        positive, negative = _difference_levels(form)
        cycles, bref, refi_type, refine = _refmac_settings(form)
        distances, distance_sigma, angles, angle_sigma, dihedrals, dihedral_sigma = _active_restraints(form, refine)
        center_mode, center_xyz = _map_center(form)
        custom_dictionary = _custom_dictionary(form, files, output_dir, project) if refine else None
        harmonic_protein, harmonic_protein_sigma = _protein_harmonic(form, refine)
        harmonic_atoms, harmonic_sigma = _specified_harmonic(form, refine)
        cfg.update(
            {
                "dark_id": _identifier(form, "dark_id", "Omit-model name"),
                "light_id": _identifier(form, "light_id", "Observed-data name"),
                "refmac_ncyc": cycles,
                "refmac_bref": bref,
                "refmac_type": refi_type,
                "exclude_restraints": _exclude_restraints(form, refine),
                "harmonic_atoms": harmonic_atoms,
                "harmonic_sigma": harmonic_sigma,
                "harmonic_protein": harmonic_protein,
                "harmonic_protein_sigma": harmonic_protein_sigma,
                "map_center": center_mode,
                "map_center_xyz": center_xyz,
                "distance_restraints": distances,
                "distance_sigma": distance_sigma,
                "angle_restraints": angles,
                "angle_sigma": angle_sigma,
                "dihedral_restraints": dihedrals,
                "dihedral_sigma": dihedral_sigma,
                "custom_dictionary": custom_dictionary,
                "positive_levels": positive,
                "negative_levels": negative,
                "files": {
                    "dark_pdb": take("dark_pdb", "omit-model coordinates"),
                    "light_pdb": take("light_pdb", "light coordinates"),
                    "light_sf": take("light_sf", "observed structure factors"),
                },
            }
        )
    elif workflow == "reconstruct":
        positive, negative = _difference_levels(form)
        cycles, bref, refi_type, refine = _refmac_settings(form)
        n_atoms = _int(form, "n_atoms", "Number of atoms to reconstruct")
        if n_atoms < 1:
            raise FormError("Number of atoms to reconstruct must be at least 1.")
        density_cutoff = _float(form, "density_cutoff", "Density cutoff")
        if not math.isfinite(density_cutoff):
            raise FormError("Density cutoff must be a finite number.")
        n_check = _int(form, "n_check_map", "Map check interval")
        if n_check < 0:
            raise FormError("Map check interval cannot be negative. Use 0 to show only the final model.")
        spline_seed = _flag(form, "spline_seed")
        if spline_seed:
            spline_order = _int(form, "spline_order", "Spline order")
            if not 0 <= spline_order <= 5:
                raise FormError("Spline order must be an integer from 0 through 5.")
        else:
            spline_order = 3
        distances, distance_sigma, angles, angle_sigma, dihedrals, dihedral_sigma = _active_restraints(form, refine)
        center_mode, center_xyz = _map_center(form)
        custom_dictionary = _custom_dictionary(form, files, output_dir, project) if refine else None
        harmonic_added = refine and _flag(form, "harmonic_added")
        harmonic_protein, harmonic_protein_sigma = _protein_harmonic(form, refine)
        exclude_restraints = _exclude_restraints(form, refine)
        restore_links = _flag(form, "restore_links")
        cfg["show_oec_bonds"] = False
        cfg["ray_trace"] = True
        cfg["recompute_map"] = True
        cfg["recompute_fragment"] = True
        cfg.update(
            {
                "dark_id": _identifier(form, "dark_id", "Omit-model name"),
                "light_id": _identifier(form, "light_id", "Observed-data name"),
                "refmac_ncyc": cycles,
                "initial_refmac_ncyc": _initial_refmac_cycles(form, refine, cycles),
                "refmac_bref": bref,
                "refmac_type": refi_type,
                "map_center": center_mode,
                "map_center_xyz": center_xyz,
                "distance_restraints": distances,
                "distance_sigma": distance_sigma,
                "angle_restraints": angles,
                "angle_sigma": angle_sigma,
                "dihedral_restraints": dihedrals,
                "dihedral_sigma": dihedral_sigma,
                "exclude_restraints": exclude_restraints,
                "custom_dictionary": custom_dictionary,
                "positive_levels": positive,
                "negative_levels": negative,
                "n_atoms": n_atoms,
                "restore_links": restore_links,
                "link_records": _link_record_file(form, files, output_dir, project),
                "density_cutoff": density_cutoff,
                "n_check_map": n_check,
                "candidate_atoms": _candidate_atoms(form),
                "harmonic_added": harmonic_added,
                "harmonic_sigma": _harmonic_sigma(form) if harmonic_added else 0.1,
                "harmonic_protein": harmonic_protein,
                "harmonic_protein_sigma": harmonic_protein_sigma,
                "neighborhood_size": _odd_neighborhood(form),
                "spline_seed": spline_seed,
                "spline_order": spline_order,
                "geometry_bonds": True,
                "viewer_light_label": "Reference site",
                "files": {
                    "dark_pdb": _pdb_coordinates(take("dark_pdb", "starting omit-model coordinates"), "Starting omit-model coordinates"),
                    "light_pdb": _pdb_coordinates(take("light_pdb", "reference model"), "Reference model"),
                    "light_sf": take("light_sf", "observed structure factors"),
                },
            }
        )
    elif workflow == "twofofc":
        levels = _levels(_text(form, "sigma_levels"), "Sigma levels")
        cfg.update(
            {
                "pdb_id": _identifier(form, "pdb_id", "Structure name"),
                "delta_b": _float(form, "delta_b", "Sharpening ΔB"),
                "sigma_levels": levels,
                "files": {
                    "pdb": take("pdb", "coordinates"),
                    "sf": take("sf", "structure factors"),
                    "validation": take("validation", "validation 2Fo−Fc coefficients"),
                },
            }
        )
    else:
        raise FormError(f"Unknown workflow {workflow!r}.")
    return cfg
