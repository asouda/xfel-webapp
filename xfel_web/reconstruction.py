"""Rebuild an omit model by placing reference atoms on omit-map peaks.

Each iteration calls the omit-map workflow. The peak search is the local-maximum
and optional spline placement in ``xfel_web.peaks``.
"""

from __future__ import annotations

import math
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
) -> list[PdbAtom]:
    """Heavy reference atoms whose name is absent from the omit model."""
    present = {atom.identity for atom in omit}
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


def _present(result: dict, *, view_id: str, iteration: int, checkpoint: bool, caption: str, extras: list) -> dict:
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
    return shown


def _xyz_text(xyz) -> str:
    return f"{xyz[0]:.3f} {xyz[1]:.3f} {xyz[2]:.3f}"


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

    missing = missing_atoms(reference, read_pdb_atoms(working), set())
    emit(
        f"Heavy atoms in {reference_path.name} that are absent from the omit model: {len(missing)}. "
        "Hydrogens are not candidates."
    )
    if not missing:
        raise ValueError(
            "The reference model has no heavy atom that is absent from the omit model, "
            "so there is nothing to reconstruct."
        )
    if n_atoms > len(missing):
        emit(
            f"Warning: {n_atoms} atoms were requested, and only {len(missing)} heavy atoms are missing. "
            "The run stops when none remain."
        )
    if int(cfg.get("refmac_ncyc") or 0) > 0:
        emit(
            "REFMAC cycles is greater than zero. Each map is calculated from the refined coordinates. "
            "The next iteration still starts from the omit model with the reference atom inserted, "
            "not from the shifted coordinates."
        )

    added: list[dict] = []
    excluded: set[tuple] = set()
    extras: list[tuple[str, Path]] = []
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
        result = omit_runner(_omit_config(cfg, snapshot, iter_dir, added), log)
        fragment = Path(result["fragment"])
        if not fragment.is_file():
            raise FileNotFoundError(f"The omit-map clip was not written: {fragment}")

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

        candidates = missing_atoms(reference, read_pdb_atoms(working), excluded)
        atom, distance = closest_atom(candidates, xyz)
        header = [
            "# Cross sections of the omit mFo-DFc density through the peak",
            f"# voxel_index {peak['index'][0]} {peak['index'][1]} {peak['index'][2]}",
            f"# voxel_density {peak['density']:.6f}",
            f"# peak_xyz {_xyz_text(xyz)}",
        ]
        action = "stopped"
        if atom is None or distance is None:
            stop_reason = "The reference model has no heavy atom that is still missing from the omit model."
            header.append("# atom none")
            log(f"Stopping: {stop_reason}")
        else:
            header.append(f"# atom {describe_atom(atom)}")
            header.append(f"# distance_A {distance:.4f}")
            log(f"Closest missing atom: {describe_atom(atom)}, {distance:.2f} Å from the peak")
            if distance > STOP_DISTANCE_A:
                stop_reason = (
                    f"{describe_atom(atom)} is {distance:.2f} Å from the peak, farther than "
                    f"{STOP_DISTANCE_A:.0f} Å. The atom was not added."
                )
                header.append("# action not_added")
                log(f"Stopping: {stop_reason}")
            else:
                if distance > WARN_DISTANCE_A:
                    log(
                        f"Warning: {describe_atom(atom)} is {distance:.2f} Å from the peak, "
                        f"farther than {WARN_DISTANCE_A:.0f} Å."
                    )
                text = working.read_text(encoding="utf-8", errors="replace")
                placed = place_record_at(atom.line, xyz)
                working.write_text(insert_atom_record(text, placed), encoding="utf-8")
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
                    f"Added {describe_atom(atom)} at the peak coordinates {_xyz_text(xyz)} Å."
                )
        sections = iter_dir / "density_sections.txt"
        write_cross_sections(peak["map"], xyz, sections, header=header)
        log(f"Density sections: {sections}")
        extras.append((f"Iteration {iteration} density sections", sections))
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
        result = omit_runner(_omit_config(cfg, snapshot, iter_dir, added), log)
        stop_reason = f"Added all {n_atoms} requested atoms."

    reconstructed = output_dir / f"{cfg['dark_id']}_reconstructed.pdb"
    shutil.copy2(working, reconstructed)
    extras.append(("Reconstructed model", reconstructed))
    if added:
        detail = "; ".join(
            f"{item['name']} {item['resname']} {item['chain']}{item['residue']} ({item['distance']:.2f} Å)"
            for item in added
        )
    else:
        detail = "no atoms were added"
    caption = (
        f"Reconstruction of {cfg['dark_id']} against {cfg['light_id']}: "
        f"{len(added)} of {n_atoms} atoms added ({detail}). {stop_reason} "
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
    )
