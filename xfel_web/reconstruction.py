"""Rebuild an omit model by placing reference atoms on omit-map peaks.

Each iteration refines the current omit model, calculates the omit map, and
inserts the closest reference atom at the highest peak. The next iteration
starts from those refined coordinates with the new atom included. The peak
search is the local-maximum and optional spline placement in ``xfel_web.peaks``.
"""

from __future__ import annotations

import math
import re
import shutil
from dataclasses import dataclass
from pathlib import Path

from xfel_web.peaks import DensityMap, choose_peak, spline_maximum_index

WARN_DISTANCE_A = 1.0
STOP_DISTANCE_A = 3.0
SECTION_HALF_A = 2.0
SECTION_STEP_A = 0.05


@dataclass(frozen=True)
class PdbAtom:
    chain: str
    resseq: int
    icode: str
    name: str
    alt: str
    resname: str
    element: str
    occ: float
    xyz: tuple[float, float, float]
    line: str

    @property
    def identity(self) -> tuple:
        """Chain, residue, insertion code, and atom name. Altloc is not included."""
        return (self.chain, self.resseq, self.icode, self.name)

    @property
    def record_key(self) -> tuple:
        return self.identity + (self.alt,)


def _element_symbol(line: str) -> str:
    if len(line) >= 78 and line[76:78].strip():
        return line[76:78].strip().upper()
    letters = "".join(ch for ch in line[12:16] if ch.isalpha())
    return letters[:1].upper()


def read_pdb_atoms(path: Path) -> list[PdbAtom]:
    atoms = []
    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.startswith(("ATOM", "HETATM")) or len(raw) < 54:
                continue
            try:
                resseq = int(raw[22:26])
                xyz = (float(raw[30:38]), float(raw[38:46]), float(raw[46:54]))
            except ValueError:
                continue
            occ_text = raw[54:60].strip() if len(raw) >= 60 else ""
            try:
                occ = float(occ_text) if occ_text else 1.0
            except ValueError:
                occ = 1.0
            atoms.append(
                PdbAtom(
                    chain=raw[21],
                    resseq=resseq,
                    icode=raw[26] if len(raw) > 26 else " ",
                    name=raw[12:16].strip(),
                    alt=raw[16],
                    resname=raw[17:20].strip(),
                    element=_element_symbol(raw),
                    occ=occ,
                    xyz=xyz,
                    line=raw.rstrip("\n"),
                )
            )
    return atoms


def describe_atom(atom: PdbAtom) -> str:
    icode = atom.icode.strip()
    residue = f"{atom.resseq}{icode}"
    alt = atom.alt.strip() or "blank"
    return (
        f"{atom.name} {atom.resname} chain {atom.chain} residue {residue} "
        f"altloc {alt} element {atom.element}"
    )


def missing_atoms(
    reference: list[PdbAtom],
    omit: list[PdbAtom],
    excluded: set[tuple],
    *,
    include_present: bool = False,
) -> list[PdbAtom]:
    """Heavy reference atoms that may be placed on a peak.

    Hydrogens and non-positive occupancies are left out. Unless
    ``include_present`` is set, an atom whose name is already in the omit
    model is left out as well. ``excluded`` holds atoms already placed in
    this run.
    """
    present = set() if include_present else {atom.identity for atom in omit}
    missing = []
    for atom in reference:
        if atom.element in {"H", "D"} or atom.occ <= 0:
            continue
        if atom.identity in present or atom.record_key in excluded:
            continue
        missing.append(atom)
    return missing


def closest_atom(
    atoms: list[PdbAtom], xyz: tuple[float, float, float]
) -> tuple[PdbAtom | None, float | None]:
    best = None
    best_distance = None
    for atom in atoms:
        distance = math.dist(atom.xyz, xyz)
        if best_distance is None or distance < best_distance:
            best = atom
            best_distance = distance
    return best, best_distance


def place_record_at(record: str, xyz: tuple[float, float, float]) -> str:
    """Copy a PDB atom record and set columns 31–54 to orthogonal coordinates."""
    text = record.rstrip("\n")
    if len(text) < 54 or not text.startswith(("ATOM", "HETATM")):
        raise ValueError("The reference atom is not an ATOM or HETATM record.")
    x, y, z = (float(value) for value in xyz)
    return f"{text[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{text[54:]}".rstrip()


def _refined_added_name(dark_id: str, atom_name: str, iteration: int) -> str:
    """File name for the refined model after one atom has been inserted."""
    label = "".join(ch if ch.isalnum() else "_" for ch in atom_name.strip()) or "atom"
    return f"{dark_id}_refined_added_{label}_step{iteration:02d}.pdb"


def insert_atom_record(pdb_text: str, record: str) -> str:
    """Insert one PDB atom, with a fresh serial, beside its residue or before END."""
    record = record.rstrip("\n")
    if len(record) < 54 or not record.startswith(("ATOM", "HETATM")):
        raise ValueError("The reference atom is not an ATOM or HETATM record.")
    lines = pdb_text.splitlines()
    serial = 0
    for line in lines:
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 11:
            try:
                serial = max(serial, int(line[6:11]))
            except ValueError:
                continue
    padded = record.ljust(80)
    updated = f"{padded[:6]}{serial + 1:5d}{padded[11:]}".rstrip()
    chain = record[21]
    resseq = int(record[22:26])
    icode = record[26] if len(record) > 26 else " "
    last_match = None
    end_at = None
    for index, line in enumerate(lines):
        if line.startswith(("ATOM", "HETATM")) and len(line) > 26:
            try:
                same = line[21] == chain and int(line[22:26]) == resseq and line[26] == icode
            except ValueError:
                same = False
            if same:
                last_match = index
        if end_at is None and line.startswith("END"):
            end_at = index
    if last_match is not None:
        lines.insert(last_match + 1, updated)
    elif end_at is not None:
        lines.insert(end_at, updated)
    else:
        lines.append(updated)
    return "\n".join(lines) + "\n"


def locate_peak(map_path: Path, neighborhood_size: int, spline_seed: bool, spline_order: int):
    """Return the highest local-maximum voxel, or None when the map has no finite peak.

    ``spline_seed`` moves the reported point to the spline maximum inside the
    neighborhood. The density is the value on the voxel ``choose_peak`` selected.
    """
    import numpy as np

    density_map = DensityMap(map_path)
    rho = density_map.original
    valid = np.isfinite(rho)
    exclusion = np.zeros(rho.shape, dtype=bool)
    found = choose_peak(rho, valid, exclusion, int(neighborhood_size))
    if found is None:
        return None
    index, density = found
    point = np.asarray(index, dtype=float)
    if spline_seed:
        point = np.asarray(
            spline_maximum_index(rho, index, int(neighborhood_size), int(spline_order)),
            dtype=float,
        )
    xyz = np.asarray(density_map.to_cart(point), dtype=float).reshape(3)
    grid_index = tuple(int(value) for value in np.asarray(index).reshape(3))
    return {
        "index": grid_index,
        "point": tuple(float(value) for value in point.reshape(3)),
        "xyz": tuple(float(value) for value in xyz),
        "density": float(density),
        "map": density_map,
    }


def _sample_density(values, index):
    """Trilinear density at a continuous array index. Outside the grid is NaN."""
    import numpy as np

    point = np.asarray(index, dtype=float).reshape(3)
    shape = np.asarray(values.shape, dtype=int)
    if np.any(point < 0) or np.any(point > shape - 1):
        return float("nan")
    base = np.floor(point).astype(int)
    upper = np.minimum(base + 1, shape - 1)
    weight = point - base
    total = 0.0
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                corner = (
                    int(upper[0] if dx else base[0]),
                    int(upper[1] if dy else base[1]),
                    int(upper[2] if dz else base[2]),
                )
                factor = (
                    (weight[0] if dx else 1.0 - weight[0])
                    * (weight[1] if dy else 1.0 - weight[1])
                    * (weight[2] if dz else 1.0 - weight[2])
                )
                total += factor * float(values[corner])
    return float(total)


def _section_directions() -> list[tuple[str, tuple[float, float, float]]]:
    """Unit vectors for the Cartesian cube through the peak.

    The three axes, six face diagonals, and four body diagonals are lines
    through the same point. Each vector has length 1 Å, so an offset along
    it is a distance in angstroms. A face diagonal is divided by √2 and a
    body diagonal by √3.
    """
    face = math.sqrt(2.0)
    body = math.sqrt(3.0)
    directions = [
        ("density_x", (1.0, 0.0, 0.0)),
        ("density_y", (0.0, 1.0, 0.0)),
        ("density_z", (0.0, 0.0, 1.0)),
        ("density_xy_p", (1.0 / face, 1.0 / face, 0.0)),
        ("density_xy_m", (1.0 / face, -1.0 / face, 0.0)),
        ("density_xz_p", (1.0 / face, 0.0, 1.0 / face)),
        ("density_xz_m", (1.0 / face, 0.0, -1.0 / face)),
        ("density_yz_p", (0.0, 1.0 / face, 1.0 / face)),
        ("density_yz_m", (0.0, 1.0 / face, -1.0 / face)),
    ]
    for signs in ((1, 1, 1), (1, 1, -1), (1, -1, 1), (1, -1, -1)):
        label = "".join("p" if sign > 0 else "m" for sign in signs)
        directions.append(
            (
                f"density_xyz_{label}",
                tuple(sign / body for sign in signs),
            )
        )
    return directions


def write_cross_sections(density_map, xyz, path: Path, *, header: list[str]) -> None:
    """Write density profiles through ``xyz``, from −2 Å to +2 Å along each line."""
    import numpy as np

    origin = np.asarray(xyz, dtype=float)
    directions = _section_directions()
    count = int(round((2 * SECTION_HALF_A) / SECTION_STEP_A))
    offsets = [-SECTION_HALF_A + index * SECTION_STEP_A for index in range(count + 1)]
    columns = []
    for _name, direction in directions:
        vector = np.asarray(direction, dtype=float)
        samples = []
        for offset in offsets:
            point = origin + float(offset) * vector
            samples.append(_sample_density(density_map.original, density_map.to_index(point)))
        columns.append(samples)
    lines = list(header)
    lines.append(
        "# offset_A is the signed distance in angstroms along each line. "
        "Face diagonals are unit vectors with components ±1/sqrt(2). "
        "Body diagonals are unit vectors with components ±1/sqrt(3). "
        "p is positive and m is negative, in x, y, z order."
    )
    lines.append("# offset_A " + " ".join(name for name, _direction in directions))

    def _fmt(value: float) -> str:
        if not math.isfinite(value):
            return "nan"
        return f"{value:.6f}"

    for row, offset in enumerate(offsets):
        values = " ".join(_fmt(column[row]) for column in columns)
        lines.append(f"{float(offset):7.3f} {values}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


class _IterationLog:
    def __init__(self, emit, path: Path):
        self._emit = emit
        self._path = path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("", encoding="utf-8")

    def __call__(self, message: str = "") -> None:
        text = str(message).rstrip("\n")
        self._emit(text)
        with self._path.open("a", encoding="utf-8") as handle:
            handle.write(text + "\n")


def _omit_config(cfg: dict, omit_pdb: Path, output_dir: Path, added: list[dict]) -> dict:
    """Omit-map settings for one iteration, plus harmonic restraints on atoms already added."""
    child = _iteration_config(cfg, omit_pdb, output_dir)
    if cfg.get("harmonic_added") and added:
        child["harmonic_atoms"] = [
            {
                "chain": item["chain"],
                "resi": int(item["resi"]),
                "name": item["name"],
                "alt": item.get("alt") or "",
                "icode": item.get("icode") or "",
            }
            for item in added
        ]
        child["harmonic_sigma"] = float(cfg.get("harmonic_sigma", 0.1))
    return child


def _iteration_config(cfg: dict, omit_pdb: Path, output_dir: Path) -> dict:
    child = dict(cfg)
    child["output_dir"] = output_dir
    child["files"] = dict(cfg["files"])
    child["files"]["dark_pdb"] = omit_pdb
    child["recompute_map"] = True
    child["recompute_fragment"] = True
    child["ray_trace"] = True
    child["show_oec_bonds"] = False
    child["geometry_bonds"] = True
    child["viewer_light_label"] = "Reference site"
    return child


def _present(
    result: dict,
    *,
    view_id: str,
    iteration: int,
    checkpoint: bool,
    caption: str,
    extras: list,
    steps: list | None = None,
) -> dict:
    shown = dict(result)
    shown["view_id"] = view_id
    shown["iteration"] = iteration
    shown["checkpoint"] = checkpoint
    shown["caption"] = caption
    files = list(result.get("files") or [])
    for role, path in extras:
        if Path(path).is_file():
            files.append({"role": role, "path": str(path)})
    shown["files"] = files
    if steps is not None:
        shown["steps"] = steps
    return shown


def _xyz_text(xyz) -> str:
    return f"{xyz[0]:.3f} {xyz[1]:.3f} {xyz[2]:.3f}"


def _atom_table_label(atom: PdbAtom) -> str:
    """Short residue label for one row of the reconstruction table."""
    icode = atom.icode.strip()
    label = f"{atom.name} {atom.resname} {atom.chain}{atom.resseq}{icode}"
    alt = atom.alt.strip()
    if alt:
        label += f" alt {alt}"
    return label


def _refmac_log(directory: Path, cfg: dict) -> Path:
    cycles = int(cfg.get("refmac_ncyc") or 0)
    return directory / f"{cfg['dark_id']}_vs_{cfg['light_id']}_refmac{cycles}.log"


_FINAL_BLOCK = re.compile(r"Final results\s*\$\$(.*?)\n\s*\$\$", re.DOTALL)
_STAT_LINE = re.compile(
    r"^\s*(R factor|R free)\s+([+-]?\d+(?:\.\d+)?)(?:\s+([+-]?\d+(?:\.\d+)?))?\s*$",
    re.MULTILINE,
)
_CYCLE_ROW = re.compile(
    r"^\s*\d+\s+([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)\b",
    re.MULTILINE,
)


def read_refmac_rfactors(log_path: Path) -> tuple[float | None, float | None]:
    """Final working and free R factors from a refmac5 log.

    The Final column of the ``Final results`` block is preferred. When that
    block is absent, the last row of the cycle table is used. A missing log
    or an unrecognised log returns ``(None, None)``.
    """
    if not log_path.is_file():
        return None, None
    text = log_path.read_text(encoding="utf-8", errors="replace")
    block = _FINAL_BLOCK.search(text)
    if block:
        found: dict[str, float] = {}
        for kind, first, second in _STAT_LINE.findall(block.group(1)):
            found[kind] = float(second if second else first)
        if "R factor" in found:
            return found.get("R factor"), found.get("R free")
    start = text.find("Rfactor analysis, stats vs cycle")
    if start < 0:
        return None, None
    end = text.find("Final results", start)
    section = text[start:] if end < 0 else text[start:end]
    rows = _CYCLE_ROW.findall(section)
    if not rows:
        return None, None
    r_factor, r_free = rows[-1]
    return float(r_factor), float(r_free)


def _step_row(
    step: int,
    r_factor: float | None,
    r_free: float | None,
    *,
    atom: str = "",
    coordinates=None,
    distance: float | None = None,
    density: float | None = None,
    note: str = "",
) -> dict:
    return {
        "step": step,
        "r_factor": r_factor,
        "r_free": r_free,
        "atom": atom,
        "coordinates": None if coordinates is None else [float(value) for value in coordinates],
        "distance": None if distance is None else float(distance),
        "density": None if density is None else float(density),
        "note": note,
    }


_STAMPED_SUFFIXES = {".pdb", ".map", ".ccp4", ".png"}
_BINARY_SUFFIXES = _STAMPED_SUFFIXES | {".ent", ".mrc", ".mtz", ".jpg", ".jpeg", ".gif", ".webp"}


def _step_tag(iteration: int) -> str:
    return f"_step{iteration:02d}"


def _rewrite_mentions(directory: Path, replacements: list[tuple[str, str]]) -> None:
    """Point logs and scripts at the renamed iteration files."""
    for path in directory.rglob("*"):
        if not path.is_file() or path.suffix.lower() in _BINARY_SUFFIXES:
            continue
        data = path.read_bytes()
        if b"\0" in data:
            continue
        try:
            text = data.decode("utf-8")
        except UnicodeDecodeError:
            continue
        updated = text
        for old, new in replacements:
            updated = updated.replace(old, new)
        if updated != text:
            path.write_text(updated, encoding="utf-8")


def _retarget_paths(value, replacements: list[tuple[str, str]]):
    if isinstance(value, str):
        for old, new in replacements:
            value = value.replace(old, new)
        return value
    if isinstance(value, dict):
        return {key: _retarget_paths(item, replacements) for key, item in value.items()}
    if isinstance(value, list):
        return [_retarget_paths(item, replacements) for item in value]
    return value


def _stamp_iteration_outputs(directory: Path, iteration: int, result: dict) -> dict:
    """Append ``_stepNN`` to pdb, map, ccp4, png, and density-section names.

    The number matches the ``r-iteration-NN`` directory. Paths stored in the
    omit-map result, and the same paths written into logs and scripts, follow
    the files.
    """
    tag = _step_tag(iteration)
    renamed: list[tuple[Path, Path]] = []
    for path in list(directory.iterdir()):
        if not path.is_file():
            continue
        if path.suffix.lower() not in _STAMPED_SUFFIXES and path.name != "density_sections.txt":
            continue
        dest = path.with_name(f"{path.stem}{tag}{path.suffix}")
        path.rename(dest)
        renamed.append((path, dest))
    if not renamed:
        return result
    replacements: list[tuple[str, str]] = []
    for old, new in renamed:
        replacements.append((str(old), str(new)))
        replacements.append((old.name, new.name))
    replacements.sort(key=lambda item: len(item[0]), reverse=True)
    _rewrite_mentions(directory, replacements)
    return _retarget_paths(result, replacements)


def run_reconstruct(cfg: dict, emit, pause=None, omit_runner=None, report_step=None) -> dict:
    """Add reference atoms onto omit-map peaks, one atom per iteration."""
    if omit_runner is None:
        from xfel_web.pipelines import run_refinement

        omit_runner = run_refinement

    output_dir = Path(cfg["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    n_atoms = int(cfg["n_atoms"])
    n_check = int(cfg["n_check_map"])
    reference_path = Path(cfg["files"]["light_pdb"])
    reference = read_pdb_atoms(reference_path)
    if not reference:
        raise ValueError(f"No atoms were found in the reference model {reference_path.name}.")

    working = output_dir / f"{cfg['dark_id']}_working.pdb"
    shutil.copy2(cfg["files"]["dark_pdb"], working)
    if not read_pdb_atoms(working):
        raise ValueError(f"No atoms were found in the omit model {Path(cfg['files']['dark_pdb']).name}.")

    include_present = cfg.get("candidate_atoms", "missing") == "all"
    missing = missing_atoms(reference, read_pdb_atoms(working), set(), include_present=include_present)
    if include_present:
        emit(
            f"Heavy atoms in {reference_path.name}: {len(missing)}. "
            "Atoms already in the omit model are eligible. Hydrogens are not candidates."
        )
    else:
        emit(
            f"Heavy atoms in {reference_path.name} that are absent from the omit model: {len(missing)}. "
            "Hydrogens are not candidates."
        )
    if not missing:
        if include_present:
            raise ValueError("The reference model has no heavy atom, so there is nothing to reconstruct.")
        raise ValueError(
            "The reference model has no heavy atom that is absent from the omit model, "
            "so there is nothing to reconstruct."
        )
    if n_atoms > len(missing):
        if include_present:
            emit(
                f"Warning: {n_atoms} atoms were requested, and the reference model has only {len(missing)} heavy atoms. "
                "The run stops when none remain."
            )
        else:
            emit(
                f"Warning: {n_atoms} atoms were requested, and only {len(missing)} heavy atoms are missing. "
                "The run stops when none remain."
            )
    if int(cfg.get("refmac_ncyc") or 0) > 0:
        emit(
            "REFMAC cycles is greater than zero. Each map is calculated after refinement. "
            "The new atom is inserted into those refined coordinates, "
            "and the next iteration starts from that model."
        )

    added: list[dict] = []
    excluded: set[tuple] = set()
    extras: list[tuple[str, Path]] = []
    steps: list[dict] = []
    iteration = 0
    result = None
    stop_reason = ""

    while len(added) < n_atoms:
        iteration += 1
        if report_step is not None:
            report_step(iteration)
        iter_dir = output_dir / f"r-iteration-{iteration:02d}"
        iter_dir.mkdir(parents=True, exist_ok=True)
        log = _IterationLog(emit, iter_dir / "reconstruction.log")
        snapshot = iter_dir / f"{cfg['dark_id']}_omit.pdb"
        shutil.copy2(working, snapshot)
        log("")
        log(f"Reconstruction iteration {iteration}: omit model {snapshot.name}, {len(added)} atoms added so far")
        result = _stamp_iteration_outputs(
            iter_dir, iteration, omit_runner(_omit_config(cfg, snapshot, iter_dir, added), log)
        )
        fragment = Path(result["fragment"])
        if not fragment.is_file():
            raise FileNotFoundError(f"The omit-map clip was not written: {fragment}")
        r_factor, r_free = read_refmac_rfactors(_refmac_log(iter_dir, cfg))

        if n_check and iteration % n_check == 0:
            if pause is None:
                log("A map check was requested, but this run cannot pause, so the search continues.")
            else:
                log("Pausing so the current residual map and omit model can be inspected.")
                caption = (
                    f"Iteration {iteration}. Inspect the residual map and the current omit model, "
                    f"then continue. {len(added)} of {n_atoms} atoms added. {result['caption']}"
                )
                pause(
                    _present(
                        result,
                        view_id=f"checkpoint-{iteration}",
                        iteration=iteration,
                        checkpoint=True,
                        caption=caption,
                        extras=extras,
                    )
                )
                log("Continuing.")

        peak = locate_peak(
            fragment,
            int(cfg["neighborhood_size"]),
            bool(cfg["spline_seed"]),
            int(cfg["spline_order"]),
        )
        if peak is None:
            stop_reason = "The clipped difference map has no finite density peak."
            log(f"Stopping: {stop_reason}")
            steps.append(_step_row(iteration, r_factor, r_free, note=stop_reason))
            break

        xyz = peak["xyz"]
        log(
            f"Peak voxel index {peak['index'][0]} {peak['index'][1]} {peak['index'][2]}, "
            f"density {peak['density']:.4f}, coordinates {_xyz_text(xyz)} Å"
        )
        if cfg["spline_seed"]:
            log(
                "The spline maximum is on, so those coordinates are the maximum of a "
                f"spline of order {int(cfg['spline_order'])} inside the "
                f"{int(cfg['neighborhood_size'])}³ neighborhood. "
                "The density is the value on the selected voxel."
            )

        candidates = missing_atoms(
            reference, read_pdb_atoms(working), excluded, include_present=include_present
        )
        atom, distance = closest_atom(candidates, xyz)
        header = [
            "# Cross sections of the omit mFo-DFc density through the peak",
            f"# voxel_index {peak['index'][0]} {peak['index'][1]} {peak['index'][2]}",
            f"# voxel_density {peak['density']:.6f}",
            f"# peak_xyz {_xyz_text(xyz)}",
        ]
        action = "stopped"
        placed_atom = ""
        placed_xyz = None
        placed_distance = None
        placed_note = ""
        if atom is None or distance is None:
            if include_present:
                stop_reason = "Every heavy atom in the reference model has already been added."
            else:
                stop_reason = "The reference model has no heavy atom that is still missing from the omit model."
            header.append("# atom none")
            log(f"Stopping: {stop_reason}")
            placed_note = stop_reason
        else:
            header.append(f"# atom {describe_atom(atom)}")
            header.append(f"# distance_A {distance:.4f}")
            kind = "reference" if include_present else "missing"
            log(f"Closest {kind} atom: {describe_atom(atom)}, {distance:.2f} Å from the peak")
            if distance > STOP_DISTANCE_A:
                stop_reason = (
                    f"{describe_atom(atom)} is {distance:.2f} Å from the peak, farther than "
                    f"{STOP_DISTANCE_A:.0f} Å. The atom was not added."
                )
                header.append("# action not_added")
                log(f"Stopping: {stop_reason}")
                placed_atom = _atom_table_label(atom)
                placed_xyz = xyz
                placed_distance = distance
                placed_note = f"Not added. Farther than {STOP_DISTANCE_A:.0f} Å from the reference atom."
            else:
                if distance > WARN_DISTANCE_A:
                    log(
                        f"Warning: {describe_atom(atom)} is {distance:.2f} Å from the peak, "
                        f"farther than {WARN_DISTANCE_A:.0f} Å."
                    )
                refined_path = iter_dir / f"{cfg['dark_id']}_refined{_step_tag(iteration)}.pdb"
                if not refined_path.is_file():
                    raise FileNotFoundError(
                        f"REFMAC did not write {cfg['dark_id']}_refined.pdb in {iter_dir.name}, "
                        "so the atom cannot be inserted into the refined model."
                    )
                placed = place_record_at(atom.line, xyz)
                added_text = insert_atom_record(
                    refined_path.read_text(encoding="utf-8", errors="replace"),
                    placed,
                )
                added_path = iter_dir / _refined_added_name(cfg["dark_id"], atom.name, iteration)
                added_path.write_text(added_text, encoding="utf-8")
                working.write_text(added_text, encoding="utf-8")
                excluded.add(atom.record_key)
                added.append(
                    {
                        "name": atom.name,
                        "resname": atom.resname,
                        "chain": atom.chain,
                        "resi": atom.resseq,
                        "alt": atom.alt,
                        "icode": atom.icode,
                        "residue": f"{atom.resseq}{atom.icode.strip()}",
                        "distance": distance,
                        "xyz": list(xyz),
                        "density": peak["density"],
                    }
                )
                action = "added"
                header.append("# action added")
                log(
                    f"Added {describe_atom(atom)} at the peak coordinates {_xyz_text(xyz)} Å "
                    f"in the refined model. Next omit model: {added_path.name}."
                )
                extras.append(
                    (f"Iteration {iteration} refined model with {atom.name.strip()}", added_path)
                )
                placed_atom = _atom_table_label(atom)
                placed_xyz = xyz
                placed_distance = distance
                notes = []
                if distance > WARN_DISTANCE_A:
                    notes.append(f"Farther than {WARN_DISTANCE_A:.0f} Å from the reference atom.")
                if cfg["spline_seed"]:
                    notes.append("Coordinates are the spline maximum. Density is the voxel value.")
                placed_note = " ".join(notes)
        sections = iter_dir / f"density_sections{_step_tag(iteration)}.txt"
        write_cross_sections(peak["map"], xyz, sections, header=header)
        log(f"Density sections: {sections}")
        extras.append((f"Iteration {iteration} density sections", sections))
        steps.append(
            _step_row(
                iteration,
                r_factor,
                r_free,
                atom=placed_atom,
                coordinates=placed_xyz,
                distance=placed_distance,
                density=peak["density"],
                note=placed_note,
            )
        )
        del peak
        if action != "added":
            break

    if result is None:
        raise RuntimeError("The reconstruction did not calculate a map.")

    if len(added) == n_atoms:
        iteration += 1
        if report_step is not None:
            report_step(iteration)
        iter_dir = output_dir / f"r-iteration-{iteration:02d}"
        iter_dir.mkdir(parents=True, exist_ok=True)
        log = _IterationLog(emit, iter_dir / "reconstruction.log")
        snapshot = iter_dir / f"{cfg['dark_id']}_omit.pdb"
        shutil.copy2(working, snapshot)
        log("")
        log(f"Final residual map of the reconstructed model ({len(added)} atoms added)")
        result = _stamp_iteration_outputs(
            iter_dir, iteration, omit_runner(_omit_config(cfg, snapshot, iter_dir, added), log)
        )
        r_factor, r_free = read_refmac_rfactors(_refmac_log(iter_dir, cfg))
        steps.append(
            _step_row(
                iteration,
                r_factor,
                r_free,
                note="Residual map of the finished model.",
            )
        )
        stop_reason = f"Added all {n_atoms} requested atoms."

    refined = iter_dir / f"{cfg['dark_id']}_refined{_step_tag(iteration)}.pdb"
    if not refined.is_file():
        raise FileNotFoundError(
            f"REFMAC did not write {cfg['dark_id']}_refined.pdb in {iter_dir.name}, "
            "so the reconstructed model cannot be saved."
        )
    reconstructed = output_dir / f"{cfg['dark_id']}_reconstructed_final.pdb"
    shutil.copy2(refined, reconstructed)
    emit(f"Reconstructed model after the final refinement: {reconstructed}")
    extras.append(("Reconstructed model", reconstructed))
    caption = (
        f"Reconstruction of {cfg['dark_id']} against {cfg['light_id']}: "
        f"{len(added)} of {n_atoms} atoms added. {stop_reason} "
        f"Each iteration is in {output_dir}. {result['caption']}"
    )
    emit("")
    emit(caption)
    return _present(
        result,
        view_id=f"final-{iteration}",
        iteration=iteration,
        checkpoint=False,
        caption=caption,
        extras=extras,
        steps=steps,
    )
