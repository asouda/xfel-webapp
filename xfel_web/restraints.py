"""REFMAC external restraints and optional libcheck monomer libraries."""

from __future__ import annotations

import math
import shlex
from pathlib import Path

from xfel_web.ccp4 import (
    coords_from_pdb_text,
    extract_monomer_pdb,
    residue_conformers,
    residue_name,
    scan_hetatm,
)


def pair_distance(coords: dict[str, tuple[float, float, float]], atom_a: str, atom_b: str, where: str) -> float:
    if atom_a not in coords or atom_b not in coords:
        missing = atom_a if atom_a not in coords else atom_b
        raise ValueError(
            f"Atom {missing} is not in {where}. Atoms present: {', '.join(sorted(coords))}"
        )
    ax, ay, az = coords[atom_a]
    bx, by, bz = coords[atom_b]
    return ((ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2) ** 0.5


def _normalize_distance_specs(specs, default_sigma: float) -> list[dict]:
    if not specs:
        return []
    if default_sigma <= 0:
        raise ValueError(f"Distance sigma must be positive, got {default_sigma}")
    parsed = []
    seen: set[frozenset[str]] = set()
    for spec in specs:
        if len(spec) == 2:
            atom_a, atom_b = spec
            distance, sigma = None, default_sigma
        elif len(spec) == 3:
            atom_a, atom_b, distance = spec
            sigma = default_sigma
        elif len(spec) == 4:
            atom_a, atom_b, distance, sigma = spec
        else:
            raise ValueError(
                "Each distance restraint must be (atom_1, atom_2), "
                "(atom_1, atom_2, distance), or (atom_1, atom_2, distance, sigma). "
                f"Got {spec!r}"
            )
        atom_a = str(atom_a).strip().upper()
        atom_b = str(atom_b).strip().upper()
        key = frozenset((atom_a, atom_b))
        if atom_a == atom_b or key in seen:
            raise ValueError(f"Repeated or empty atom pair in distance restraints: {spec!r}")
        seen.add(key)
        sigma = float(sigma)
        if sigma <= 0:
            raise ValueError(f"Standard deviation for {atom_a}-{atom_b} must be positive")
        parsed.append(
            {
                "atom1": atom_a,
                "atom2": atom_b,
                "distance": None if distance is None else float(distance),
                "sigma": sigma,
            }
        )
    return parsed


def external_distance_keywords(
    pdb_path: Path,
    site_copies: list[tuple[str, int]],
    specs,
    default_sigma: float,
) -> str:
    """REFMAC external distances for each residue in ``site_copies``.

    REFMAC keeps the last restraint on an atom pair, so these lines replace
    the dictionary bond for that pair. An empty spec list adds nothing.
    """
    parsed = _normalize_distance_specs(specs, default_sigma)
    if not parsed:
        return ""
    lines: list[str] = []
    for chain, resi in site_copies:
        for alt, coords in residue_conformers(pdb_path, chain, int(resi)).items():
            where = f"chain {chain} residue {int(resi)}"
            alt_txt = f" altecode {alt}" if alt else ""
            for spec in parsed:
                if spec["distance"] is None:
                    distance = pair_distance(coords, spec["atom1"], spec["atom2"], where)
                else:
                    pair_distance(coords, spec["atom1"], spec["atom2"], where)
                    distance = spec["distance"]
                lines.append(
                    f"exte dist first chain {chain} resi {int(resi)} atom {spec['atom1']}{alt_txt} "
                    f"second chain {chain} resi {int(resi)} atom {spec['atom2']}{alt_txt} "
                    f"value {distance:.3f} sigma {spec['sigma']:.3f} type 0"
                )
    if not lines:
        raise ValueError("Distance restraints produced no external restraints")
    return "\n".join(lines)


def _angle_degrees(coords, atom_a: str, vertex: str, atom_c: str, where: str) -> float:
    missing = [name for name in (atom_a, vertex, atom_c) if name not in coords]
    if missing:
        raise ValueError(
            f"Atom {missing[0]} is not in {where}. Atoms present: {', '.join(sorted(coords))}"
        )
    ax, ay, az = coords[atom_a]
    vx, vy, vz = coords[vertex]
    cx, cy, cz = coords[atom_c]
    vector_a = (ax - vx, ay - vy, az - vz)
    vector_c = (cx - vx, cy - vy, cz - vz)
    norm_a = math.sqrt(sum(value * value for value in vector_a))
    norm_c = math.sqrt(sum(value * value for value in vector_c))
    if norm_a == 0 or norm_c == 0:
        same = atom_a if norm_a == 0 else atom_c
        raise ValueError(f"{same} and {vertex} are at the same position")
    cosine = sum(a * c for a, c in zip(vector_a, vector_c)) / (norm_a * norm_c)
    cosine = max(-1.0, min(1.0, cosine))
    return math.degrees(math.acos(cosine))


def external_angle_keywords(
    pdb_path: Path,
    site_copies: list[tuple[str, int]],
    specs,
    default_sigma_deg: float,
) -> str:
    """REFMAC external angles. The middle atom is the vertex.

    Target angles and sigmas are degrees. ``refmac5`` reads both in degrees.
    """
    if not specs:
        return ""
    if default_sigma_deg <= 0:
        raise ValueError(f"Angle sigma must be positive, got {default_sigma_deg}")
    parsed = []
    seen: set[tuple[str, frozenset[str]]] = set()
    for spec in specs:
        if len(spec) == 3:
            atom_a, vertex, atom_c = spec
            angle_deg, sigma_deg = None, default_sigma_deg
        elif len(spec) == 4:
            atom_a, vertex, atom_c, angle_deg = spec
            sigma_deg = default_sigma_deg
        elif len(spec) == 5:
            atom_a, vertex, atom_c, angle_deg, sigma_deg = spec
        else:
            raise ValueError(
                "Each angle restraint must be (atom_1, vertex, atom_3), "
                "(atom_1, vertex, atom_3, angle), or "
                "(atom_1, vertex, atom_3, angle, sigma). "
                f"Got {spec!r}"
            )
        atom_a = str(atom_a).strip().upper()
        vertex = str(vertex).strip().upper()
        atom_c = str(atom_c).strip().upper()
        if len({atom_a, vertex, atom_c}) < 3:
            raise ValueError(f"Angle atoms must be three different atoms: {spec!r}")
        key = (vertex, frozenset((atom_a, atom_c)))
        if key in seen:
            raise ValueError(f"Repeated angle in angle restraints: {spec!r}")
        seen.add(key)
        sigma_deg = float(sigma_deg)
        if sigma_deg <= 0:
            raise ValueError(f"Standard deviation for {atom_a}-{vertex}-{atom_c} must be positive")
        if angle_deg is not None:
            angle_deg = float(angle_deg)
            if not 0.0 < angle_deg <= 180.0:
                raise ValueError(f"Angle for {atom_a}-{vertex}-{atom_c} must be between 0 and 180 degrees")
        parsed.append(
            {
                "atom1": atom_a,
                "vertex": vertex,
                "atom3": atom_c,
                "angle_deg": angle_deg,
                "sigma_deg": sigma_deg,
            }
        )

    lines: list[str] = []
    for chain, resi in site_copies:
        for alt, coords in residue_conformers(pdb_path, chain, int(resi)).items():
            where = f"chain {chain} residue {int(resi)}"
            if alt:
                where += f" altloc {alt}"
            alt_txt = f" altecode {alt}" if alt else ""
            for spec in parsed:
                angle_deg = spec["angle_deg"]
                if angle_deg is None:
                    angle_deg = _angle_degrees(coords, spec["atom1"], spec["vertex"], spec["atom3"], where)
                lines.append(
                    f"exte angle first chain {chain} resi {int(resi)} "
                    f"atom {spec['atom1']}{alt_txt} "
                    f"next chain {chain} resi {int(resi)} "
                    f"atom {spec['vertex']}{alt_txt} "
                    f"next chain {chain} resi {int(resi)} "
                    f"atom {spec['atom3']}{alt_txt} "
                    f"value {angle_deg:.3f} sigma {spec['sigma_deg']:.3f} type 0"
                )
    if not lines:
        raise ValueError("Angle restraints produced no external restraints")
    return "\n".join(lines)


def parse_distance_restraints(specs, default_sigma: float, coords: dict) -> list[dict]:
    parsed = _normalize_distance_specs(specs, default_sigma)
    restraints = []
    for spec in parsed:
        distance = spec["distance"]
        if distance is None:
            distance = pair_distance(coords, spec["atom1"], spec["atom2"], "the site residue")
        else:
            pair_distance(coords, spec["atom1"], spec["atom2"], "the site residue")
        restraints.append(
            {"atom1": spec["atom1"], "atom2": spec["atom2"], "distance": distance, "sigma": spec["sigma"]}
        )
    return restraints


def library_monomer_names(pdb_path: Path, site_resname: str, requested: list[str] | None) -> list[str]:
    if requested is None:
        names: list[str] = []
        seen: set[str] = set()
        for rec in scan_hetatm(pdb_path):
            if rec["resname"] not in seen:
                seen.add(rec["resname"])
                names.append(rec["resname"])
    else:
        names = [str(name).strip() for name in requested if str(name).strip()]
    if site_resname not in names:
        names.insert(0, site_resname)
    else:
        names = [site_resname] + [name for name in names if name != site_resname]
    return names


def _run_libcheck(session, work: Path, keywords: str, log_path: Path) -> None:
    script = f"cd {shlex.quote(str(work))}\nlibcheck << 'EOF'\n{keywords.strip()}\n_END\nEOF\n"
    try:
        session.run_bash(script, log_path)
    except RuntimeError as exc:
        session.emit(f"  libcheck exited with an error. See {log_path.name}")
        log_path.write_text(log_path.read_text() + "\n" + str(exc))


def run_libcheck_monomer(session, pdb_text: str, resname: str, work: Path, required: bool = True) -> str | None:
    """Monomer description from coordinates, or the CCP4 entry if that fails."""
    work.mkdir(parents=True, exist_ok=True)
    (work / f"{resname}.pdb").write_text(pdb_text)
    has_h = any(
        line.startswith(("ATOM", "HETATM")) and line[76:78].strip() == "H"
        for line in pdb_text.splitlines()
    )
    hflag = "Y" if has_h else "N"
    lib_path = work / f"{resname}.lib"
    lib_path.unlink(missing_ok=True)
    log_path = work / f"{resname}_libcheck.log"
    session.emit(f"→ libcheck {resname}  ({log_path.name})")
    _run_libcheck(
        session,
        work,
        "\n".join(
            [
                "_DOC N",
                f"_FILE_PDB {resname}.pdb",
                f"_MON {resname}",
                "_COOR Y",
                "_REF N",
                f"_HFLAG {hflag}",
                "_SRCH 0",
                f"_FILE_O {resname}",
            ]
        ),
        log_path,
    )
    if not lib_path.is_file():
        session.emit(
            f"  {resname}: no dictionary from the coordinates "
            f"(see {log_path.name}). Using the CCP4 monomer entry."
        )
        _run_libcheck(
            session,
            work,
            "\n".join(["_DOC N", f"_MON {resname}", f"_FILE_O {resname}"]),
            work / f"{resname}_libcheck_ccp4.log",
        )
    if not lib_path.is_file():
        message = (
            f"libcheck did not write {lib_path}. See {log_path}. "
            "REFMAC will use its distributed dictionary for this residue."
        )
        if required:
            raise FileNotFoundError(message)
        session.emit(f"  {resname}: skipped. {message}")
        return None
    for junk in list(work.glob(f"{resname}*.ps")) + list(work.glob(f"{resname}*.odb")):
        junk.unlink()
    return lib_path.read_text()


def _split_monomer_lib(lib_text: str) -> tuple[str, str, str]:
    lines = lib_text.splitlines()
    try:
        list_at = next(i for i, line in enumerate(lines) if line.startswith("data_comp_list"))
        block_at = next(
            i for i, line in enumerate(lines)
            if line.startswith("data_comp_") and not line.startswith("data_comp_list")
        )
    except StopIteration as exc:
        raise ValueError("libcheck library is missing data_comp_list or a monomer block") from exc
    summary = None
    header_lines = lines[: list_at + 1]
    for line in lines[list_at + 1 : block_at]:
        if line.startswith(("loop_", "_")):
            header_lines.append(line)
        elif (
            summary is None
            and line.strip()
            and not line.startswith("#")
            and set(line.strip()) != {"-"}
        ):
            summary = line
    if summary is None:
        raise ValueError("libcheck library is missing the monomer summary line")
    block = "\n".join(lines[block_at:]).strip() + "\n"
    return "\n".join(header_lines).rstrip() + "\n", summary, block


def _format_bond(comp_id: str, atom1: str, atom2: str, bond_type: str, distance: float, sigma: float) -> str:
    return f" {comp_id:<6} {atom1:<6} {atom2:<6} {bond_type:<10} {distance:8.3f} {sigma:8.3f}"


def loosen_distance_restraints(lib_text: str, comp_id: str, restraints: list[dict]) -> tuple[str, list[str]]:
    """Set the chosen site distances and their ESDs in one monomer block."""
    lines = lib_text.splitlines()
    if not any(line.strip() == "_chem_comp_bond.value_dist_esd" for line in lines):
        lines.extend(
            [
                "loop_",
                "_chem_comp_bond.comp_id",
                "_chem_comp_bond.atom_id_1",
                "_chem_comp_bond.atom_id_2",
                "_chem_comp_bond.type",
                "_chem_comp_bond.value_dist",
                "_chem_comp_bond.value_dist_esd",
            ]
        )
    atom_ids: set[str] = set()
    in_atoms = False
    for line in lines:
        if line.strip() == "_chem_comp_atom.z":
            in_atoms = True
            continue
        if in_atoms and (line.startswith("loop_") or line.startswith("_") or line.startswith("#")):
            break
        if in_atoms:
            parts = line.split()
            if len(parts) >= 2 and parts[0] == comp_id:
                atom_ids.add(parts[1])
    missing = [
        atom
        for rec in restraints
        for atom in (rec["atom1"], rec["atom2"])
        if atom not in atom_ids
    ]
    if missing:
        raise ValueError(
            f"{comp_id} library has no atom(s) {', '.join(missing)}. "
            f"Library atoms: {', '.join(sorted(atom_ids))}"
        )

    wanted = {frozenset((rec["atom1"], rec["atom2"])): rec for rec in restraints}
    found: set[frozenset[str]] = set()
    out: list[str] = []
    in_bonds = False
    report: list[str] = []

    def _emit_missing() -> None:
        for key, rec in wanted.items():
            if key in found:
                continue
            out.append(_format_bond(comp_id, rec["atom1"], rec["atom2"], "single", rec["distance"], rec["sigma"]))
            report.append(
                f"{rec['atom1']}-{rec['atom2']}  {rec['distance']:.3f} A  "
                f"sigma {rec['sigma']:.3f}  (added)"
            )
            found.add(key)

    for line in lines:
        if line.strip() == "_chem_comp_bond.value_dist_esd":
            in_bonds = True
            out.append(line)
            continue
        if in_bonds and (
            line.startswith("loop_")
            or line.startswith("_")
            or line.startswith("#")
            or line.startswith("data_")
        ):
            _emit_missing()
            in_bonds = False
            out.append(line)
            continue
        if in_bonds:
            parts = line.split()
            if len(parts) >= 6 and parts[0] == comp_id:
                try:
                    old_sigma = float(parts[5])
                except ValueError:
                    out.append(line)
                    continue
                key = frozenset((parts[1], parts[2]))
                rec = wanted.get(key)
                if rec is not None:
                    out.append(_format_bond(comp_id, parts[1], parts[2], parts[3], rec["distance"], rec["sigma"]))
                    report.append(
                        f"{parts[1]}-{parts[2]}  {rec['distance']:.3f} A  "
                        f"sigma {old_sigma:.3f} -> {rec['sigma']:.3f}"
                    )
                    found.add(key)
                    continue
        out.append(line)
    if in_bonds:
        _emit_missing()
    if len(found) != len(wanted):
        raise ValueError(f"Could not write every distance restraint for {comp_id}")
    return "\n".join(out) + "\n", report


def combine_monomer_libs(header: str, summaries: list[str], blocks: list[str]) -> str:
    parts = [header.rstrip(), *summaries, "#", "# --- DESCRIPTION OF MONOMERS ---", "#", *blocks]
    return "\n".join(part.rstrip("\n") for part in parts).rstrip() + "\n"


def build_refinement_library(
    session,
    pdb_path: Path,
    library_path: Path,
    site_chain: str,
    site_resi: int,
    distance_specs,
    distance_sigma: float,
    monomer_names: list[str] | None,
) -> Path:
    """libcheck library for the site residue and the other hetero residues."""
    site_resname = residue_name(pdb_path, site_chain, site_resi)
    site_pdb = extract_monomer_pdb(pdb_path, site_resname, site_chain, site_resi)
    site_coords = coords_from_pdb_text(site_pdb)
    restraints = (
        parse_distance_restraints(distance_specs, float(distance_sigma), site_coords)
        if distance_specs
        else []
    )
    names = library_monomer_names(pdb_path, site_resname, monomer_names)
    work = session.output_dir / "libcheck"
    header = None
    summaries: list[str] = []
    blocks: list[str] = []
    skipped: list[str] = []
    for name in names:
        pdb_text = site_pdb if name == site_resname else extract_monomer_pdb(pdb_path, name)
        lib_text = run_libcheck_monomer(session, pdb_text, name, work, required=(name == site_resname))
        if lib_text is None:
            skipped.append(name)
            continue
        if name == site_resname and restraints:
            lib_text, report = loosen_distance_restraints(lib_text, name, restraints)
            session.emit(f"Distance restraints for {name} {site_chain}/{site_resi}:")
            for row in report:
                session.emit(f"  {row}")
        this_header, summary, block = _split_monomer_lib(lib_text)
        if header is None:
            header = this_header
        summaries.append(summary)
        blocks.append(block)
    if header is None:
        raise RuntimeError("libcheck produced no monomer descriptions")
    library_path.write_text(combine_monomer_libs(header, summaries, blocks))
    session.emit(f"Monomer library ({len(summaries)} residue types): {library_path}")
    if skipped:
        session.emit("Left to the CCP4 dictionary: " + ", ".join(skipped))
    return library_path
