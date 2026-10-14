"""Turn an HTML form into the dictionaries the three pipelines expect."""

from __future__ import annotations

import re
from pathlib import Path

from werkzeug.datastructures import FileStorage
from werkzeug.utils import secure_filename

ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.+\-]{0,80}$")
LABEL_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]{0,20}$")

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


def _chain(form) -> str:
    chain = _text(form, "site_chain")
    if len(chain) != 1 or not chain.isalnum():
        raise FormError("Site chain must be a single letter or digit, as stored in the PDB file.")
    return chain


def _levels(text: str, label: str) -> list[float]:
    parts = text.replace(",", " ").split()
    if not parts:
        raise FormError(f"{label} needs at least one number.")
    try:
        return [float(part) for part in parts]
    except ValueError as exc:
        raise FormError(f"{label} must be a list of numbers.") from exc


def _optional_levels(text: str, fallback: list[float], label: str) -> list[float]:
    if not text.strip():
        return list(fallback)
    levels = _levels(text, label)
    if any(level <= 0 for level in levels):
        raise FormError(f"{label} uses positive contour levels. The negative contour is drawn automatically.")
    return sorted(levels)


def parse_extra_residues(text: str) -> list[dict]:
    extras = []
    for line_no, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.replace(",", " ").split()
        if len(parts) < 2:
            raise FormError(f"Extra residue line {line_no} needs a chain and a residue number.")
        chain, resi_text = parts[0], parts[1]
        if len(chain) != 1:
            raise FormError(f"Extra residue line {line_no}: chain must be one character.")
        try:
            resi = int(resi_text)
        except ValueError as exc:
            raise FormError(f"Extra residue line {line_no}: residue number is not an integer.") from exc
        label = parts[2] if len(parts) > 2 else f"{chain}{resi}"
        if not LABEL_RE.match(label):
            raise FormError(
                f"Extra residue label {label!r} must start with a letter and use only letters, digits, and underscores."
            )
        extras.append({"chain": chain, "resi": resi, "label": label})
    return extras


def parse_copies(chain: str, resi: int, text: str) -> list[tuple[str, int]]:
    copies = [(chain, resi)]
    for line_no, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.replace(",", " ").split()
        if len(parts) < 2:
            raise FormError(f"Additional copy line {line_no} needs a chain and a residue number.")
        copy_chain = parts[0]
        if len(copy_chain) != 1:
            raise FormError(f"Additional copy line {line_no}: chain must be one character.")
        try:
            copy_resi = int(parts[1])
        except ValueError as exc:
            raise FormError(f"Additional copy line {line_no}: residue number is not an integer.") from exc
        pair = (copy_chain, copy_resi)
        if pair not in copies:
            copies.append(pair)
    return copies


def parse_distance_restraints(text: str) -> list[tuple]:
    specs = []
    for line_no, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.replace(",", " ").split()
        if len(parts) not in {2, 3, 4}:
            raise FormError(
                f"Distance line {line_no} should be 'ATOM1 ATOM2', "
                "'ATOM1 ATOM2 distance', or 'ATOM1 ATOM2 distance sigma'."
            )
        atom_a, atom_b = parts[0].upper(), parts[1].upper()
        try:
            if len(parts) == 2:
                specs.append((atom_a, atom_b))
            elif len(parts) == 3:
                specs.append((atom_a, atom_b, float(parts[2])))
            else:
                specs.append((atom_a, atom_b, float(parts[2]), float(parts[3])))
        except ValueError as exc:
            raise FormError(f"Distance line {line_no}: distance and sigma must be numbers.") from exc
    return specs


def parse_angle_restraints(text: str) -> list[tuple]:
    specs = []
    for line_no, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.replace(",", " ").split()
        if len(parts) not in {3, 4, 5}:
            raise FormError(
                f"Angle line {line_no} should be 'ATOM1 VERTEX ATOM3', "
                "optionally followed by the angle and a sigma, both in degrees."
            )
        atom_a, vertex, atom_c = (part.upper() for part in parts[:3])
        try:
            if len(parts) == 3:
                specs.append((atom_a, vertex, atom_c))
            elif len(parts) == 4:
                specs.append((atom_a, vertex, atom_c, float(parts[3])))
            else:
                specs.append((atom_a, vertex, atom_c, float(parts[3]), float(parts[4])))
        except ValueError as exc:
            raise FormError(f"Angle line {line_no}: angle and sigma must be numbers, in degrees.") from exc
    return specs


def parse_monomers(text: str) -> list[str] | None:
    names = [part.strip().upper() for part in text.replace(",", " ").split() if part.strip()]
    return names or None


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


def resolve_input_file(form, files, key: str, label: str, output_dir: Path, project: Path) -> Path:
    upload = files.get(f"file_{key}")
    if isinstance(upload, FileStorage) and upload.filename:
        folder = output_dir / "inputs"
        folder.mkdir(parents=True, exist_ok=True)
        filename = secure_filename(upload.filename) or key
        dest = folder / f"{key}_{filename}"
        upload.save(dest)
        if dest.stat().st_size == 0:
            raise FormError(f"{label} was uploaded empty.")
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
    return path


def _common(form, project: Path) -> dict:
    output_dir = resolve_output_dir(_text(form, "output_dir"), project)
    ccp4_root = Path(_text(form, "ccp4_root") or "/Applications/ccp4-9").expanduser()
    resolution_source = _text(form, "resolution_source") or "sf"
    if resolution_source not in {"sf", "explicit"}:
        raise FormError("Resolution source must be the structure-factor file or explicit cutoffs.")
    lowres = _float(form, "lowres", "Low-resolution cutoff")
    dmin = _float(form, "dmin", "High-resolution cutoff")
    if resolution_source == "explicit" and not (lowres > dmin > 0):
        raise FormError("The low-resolution cutoff must be a larger number of ångströms than the high-resolution cutoff.")
    chain = _chain(form)
    resi = _int(form, "site_resi", "Site residue")
    label = _text(form, "site_label") or "site"
    if not LABEL_RE.match(label):
        raise FormError("Site label must start with a letter and use only letters, digits, and underscores.")
    contour_mode = _text(form, "contour_mode") or "sigma"
    surface = _text(form, "surface") or "mesh"
    if contour_mode not in {"sigma", "absolute"}:
        raise FormError("Contour mode must be sigma or absolute.")
    if surface not in {"mesh", "solid"}:
        raise FormError("Surface representation must be mesh or solid.")
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
        "map_border": _float(form, "map_border", "Map border"),
        "site_chain": chain,
        "site_resi": resi,
        "site_label": label,
        "site_copies": parse_copies(chain, resi, _text(form, "extra_copies")),
        "extra_residues": parse_extra_residues(_text(form, "extra_residues")),
        "export_radius": _float(form, "export_radius", "Export radius"),
        "carve": _float(form, "carve", "Carve radius"),
        "contour_mode": contour_mode,
        "surface": surface,
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
        cycles = _int(form, "refmac_ncyc", "REFMAC cycles")
        if cycles < 0:
            raise FormError("REFMAC cycles cannot be negative. Use 0 to calculate the map without moving atoms.")
        try:
            distances = parse_distance_restraints(_text(form, "distance_restraints"))
            angles = parse_angle_restraints(_text(form, "angle_restraints"))
        except ValueError as exc:
            raise FormError(str(exc)) from exc
        cfg.update(
            {
                "dark_id": _identifier(form, "dark_id", "Omit-model name"),
                "light_id": _identifier(form, "light_id", "Observed-data name"),
                "refmac_ncyc": cycles,
                "distance_restraints": distances,
                "distance_sigma": _float(form, "distance_sigma", "Distance sigma"),
                "angle_restraints": angles,
                "angle_sigma": _float(form, "angle_sigma", "Angle sigma"),
                "use_libcheck": _flag(form, "use_libcheck"),
                "lib_monomers": parse_monomers(_text(form, "lib_monomers")),
                "positive_levels": positive,
                "negative_levels": negative,
                "files": {
                    "dark_pdb": take("dark_pdb", "omit-model coordinates"),
                    "light_pdb": take("light_pdb", "light coordinates"),
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
