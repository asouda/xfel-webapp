"""REFMAC external restraints."""

from __future__ import annotations

import math
from pathlib import Path

from xfel_web.ccp4 import residue_conformers


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


# Polypeptide residues. MSE, SEC, and PYL are incorporated amino acids. The
# remaining names are the protonation states used by common force fields.
_PROTEIN_RESIDUES = frozenset(
    {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "MSE", "SEC", "PYL",
        "ASH", "GLH", "LYN", "CYM", "CYX",
        "HID", "HIE", "HIP", "HSD", "HSE", "HSP",
    }
)
_NUCLEIC_RESIDUES = frozenset(
    {
        "A", "C", "G", "U", "I",
        "DA", "DC", "DG", "DT", "DU", "DI",
        "ADE", "CYT", "GUA", "THY", "URA",
    }
)
_SOLVENT_RESIDUES = frozenset({"HOH", "WAT", "DOD", "H2O", "TIP", "SOL"})


def is_protein_atom(record: str, resname: str) -> bool:
    """True when the atom belongs to the polypeptide.

    A listed amino acid counts as ATOM or as HETATM, so selenomethionine stored
    as HETATM is still protein. Any other ATOM record is protein unless it is a
    nucleotide or a solvent. Ligands, waters, and metals recorded as HETATM are
    left out.
    """
    name = resname.strip().upper()
    if name in _PROTEIN_RESIDUES:
        return True
    if not record.startswith("ATOM"):
        return False
    return name not in _NUCLEIC_RESIDUES and name not in _SOLVENT_RESIDUES


def harmonic_atom_key(atom: dict) -> tuple:
    """Chain, residue, insertion code, atom name, and altloc of one harmonic target."""
    chain = atom.get("chain")
    chain = " " if chain is None or str(chain) == "" else str(chain)
    if chain != " ":
        chain = chain.strip() or " "
    return (
        chain[:1],
        int(atom["resi"]),
        str(atom.get("icode") or "").strip(),
        str(atom.get("name") or "").strip(),
        str(atom.get("alt") or "").strip(),
    )


def site_and_additional_residues(cfg: dict) -> list[tuple[str, int]]:
    """Residues left free: the site, its additional copies, and the neighbouring residues."""
    pairs: list[tuple[str, int]] = []
    copies = list(cfg.get("site_copies") or [])
    if not copies and cfg.get("site_chain") is not None and cfg.get("site_resi") is not None:
        copies = [(cfg["site_chain"], cfg["site_resi"])]
    for chain, resi in copies:
        pairs.append((str(chain), int(resi)))
    for item in cfg.get("extra_residues") or []:
        pairs.append((str(item["chain"]), int(item["resi"])))
    unique: list[tuple[str, int]] = []
    seen: set[tuple[str, int]] = set()
    for pair in pairs:
        if pair in seen:
            continue
        seen.add(pair)
        unique.append(pair)
    return unique


def protein_atoms_outside(pdb_path: Path, free_residues) -> list[dict]:
    """Protein atoms, in every chain, whose residue is not in ``free_residues``.

    ``free_residues`` holds ``(chain, residue number)`` pairs. An insertion code
    is not part of that pair, so every insertion at the number is left free.
    Every altloc is returned, and a repeated record is returned once.
    """
    free = {(str(chain), int(resi)) for chain, resi in free_residues}
    atoms: list[dict] = []
    seen: set[tuple] = set()
    with Path(pdb_path).open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.startswith(("ATOM", "HETATM")) or len(raw) < 27:
                continue
            if not is_protein_atom(raw[:6], raw[17:20]):
                continue
            try:
                resi = int(raw[22:26])
            except ValueError:
                continue
            chain = raw[21] if len(raw) > 21 else " "
            if (chain, resi) in free:
                continue
            name = raw[12:16].strip()
            if not name:
                continue
            atom = {
                "chain": chain,
                "resi": resi,
                "name": name,
                "alt": raw[16].strip() if len(raw) > 16 else "",
                "icode": (raw[26] if len(raw) > 26 else " ").strip(),
            }
            key = harmonic_atom_key(atom)
            if key in seen:
                continue
            seen.add(key)
            atoms.append(atom)
    return atoms


def protein_harmonic_lines(pdb_path: Path, cfg: dict, skip: list[dict] | None = None) -> str:
    """Harmonic lines for protein atoms outside the site and the additional residues.

    Atoms listed in ``skip`` already carry a harmonic line, so they are omitted.
    An empty selection returns an empty string.
    """
    atoms = protein_atoms_outside(pdb_path, site_and_additional_residues(cfg))
    if skip:
        held = {harmonic_atom_key(item) for item in skip}
        atoms = [atom for atom in atoms if harmonic_atom_key(atom) not in held]
    if not atoms:
        return ""
    sigma = float(cfg.get("harmonic_protein_sigma", 0.1))
    return external_harmonic_keywords(atoms, sigma)


def external_harmonic_keywords(atoms: list[dict], sigma: float) -> str:
    """REFMAC external harmonic lines that hold atoms at their current coordinates.

    The documented keyword is ``external harmonic chain … residue … atom … sigma``.
    REFMAC uses the coordinates in the input PDB as the target. A smaller sigma
    holds the atom more tightly. An empty atom list adds nothing.
    """
    sigma = float(sigma)
    if sigma <= 0:
        raise ValueError(f"Harmonic sigma must be positive, got {sigma}")
    lines = []
    for atom in atoms:
        chain = str(atom.get("chain") or " ").strip() or "_"
        resi = int(atom["resi"])
        name = str(atom["name"]).strip()
        if not name:
            raise ValueError("A harmonic restraint needs an atom name.")
        icode = str(atom.get("icode") or "").strip()
        alt = str(atom.get("alt") or "").strip()
        parts = ["exte harm", f"chain {chain}", f"residue {resi}"]
        if icode:
            parts.append(f"insertion {icode}")
        parts.append(f"atom {name}")
        if alt:
            parts.append(f"altcode {alt}")
        parts.append(f"sigma {sigma:.3f}")
        lines.append(" ".join(parts))
    return "\n".join(lines)


def external_distance_keywords(
    pdb_path: Path,
    site_copies: list[tuple[str, int]],
    specs,
    default_sigma: float,
) -> str:
    """REFMAC external distances for each residue in ``site_copies``.

    A spec that includes a distance or a sigma uses those numbers. A missing
    distance is taken from the coordinates, and a missing sigma is
    ``default_sigma``. REFMAC keeps the last restraint on an atom pair, so
    these lines replace the monomer-library bond, including a bond from a
    custom dictionary read earlier in the same job. An empty spec list adds
    nothing.
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
    A spec that includes an angle or a sigma uses those numbers. A missing
    angle is taken from the coordinates, and a missing sigma is
    ``default_sigma_deg``. These lines replace the monomer-library angle,
    including an angle from a custom dictionary read earlier in the same job.
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


def _sub(left, right):
    return (left[0] - right[0], left[1] - right[1], left[2] - right[2])


def _dot(left, right) -> float:
    return left[0] * right[0] + left[1] * right[1] + left[2] * right[2]


def _cross(left, right):
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def _norm(vector) -> float:
    return math.sqrt(_dot(vector, vector))


def _dihedral_degrees(coords, atom_a: str, atom_b: str, atom_c: str, atom_d: str, where: str) -> float:
    """Signed torsion of planes ABC and BCD, in degrees, from −180 to 180."""
    names = (atom_a, atom_b, atom_c, atom_d)
    missing = [name for name in names if name not in coords]
    if missing:
        raise ValueError(
            f"Atom {missing[0]} is not in {where}. Atoms present: {', '.join(sorted(coords))}"
        )
    points = [coords[name] for name in names]
    toward_a = _sub(points[0], points[1])
    toward_c = _sub(points[2], points[1])
    toward_d = _sub(points[3], points[2])
    normal_abc = _cross(toward_a, toward_c)
    normal_bcd = _cross(toward_c, toward_d)
    axis = _norm(toward_c)
    if _norm(normal_abc) == 0 or axis == 0 or _norm(normal_bcd) == 0:
        raise ValueError(
            f"Dihedral {atom_a}-{atom_b}-{atom_c}-{atom_d} is undefined because three atoms are collinear"
        )
    axis_unit = tuple(value / axis for value in toward_c)
    reference = _cross(normal_abc, axis_unit)
    return math.degrees(math.atan2(_dot(reference, normal_bcd), _dot(normal_abc, normal_bcd)))


def external_dihedral_keywords(
    pdb_path: Path,
    site_copies: list[tuple[str, int]],
    specs,
    default_sigma_deg: float,
) -> str:
    """REFMAC external torsions. The atoms are first, second, third, fourth.

    The target is the angle between the plane of the first three atoms and the
    plane of the last three, in degrees. A spec that includes an angle or a
    sigma uses those numbers. A missing angle is taken from the coordinates,
    and a missing sigma is ``default_sigma_deg``. ``period 1`` is one target.
    These lines replace the monomer-library torsion, including a torsion from
    a custom dictionary read earlier in the same job.
    """
    if not specs:
        return ""
    if default_sigma_deg <= 0:
        raise ValueError(f"Dihedral sigma must be positive, got {default_sigma_deg}")
    parsed = []
    seen: set[tuple[str, str, str, str]] = set()
    for spec in specs:
        if len(spec) == 4:
            atom_a, atom_b, atom_c, atom_d = spec
            angle_deg, sigma_deg = None, default_sigma_deg
        elif len(spec) == 5:
            atom_a, atom_b, atom_c, atom_d, angle_deg = spec
            sigma_deg = default_sigma_deg
        elif len(spec) == 6:
            atom_a, atom_b, atom_c, atom_d, angle_deg, sigma_deg = spec
        else:
            raise ValueError(
                "Each dihedral restraint must be (atom_1, atom_2, atom_3, atom_4), "
                "(atom_1, atom_2, atom_3, atom_4, angle), or "
                "(atom_1, atom_2, atom_3, atom_4, angle, sigma). "
                f"Got {spec!r}"
            )
        atom_a = str(atom_a).strip().upper()
        atom_b = str(atom_b).strip().upper()
        atom_c = str(atom_c).strip().upper()
        atom_d = str(atom_d).strip().upper()
        atoms = (atom_a, atom_b, atom_c, atom_d)
        if len(set(atoms)) < 4:
            raise ValueError(f"Dihedral atoms must be four different atoms: {spec!r}")
        reverse = (atom_d, atom_c, atom_b, atom_a)
        if atoms in seen or reverse in seen:
            raise ValueError(f"Repeated dihedral in dihedral restraints: {spec!r}")
        seen.add(atoms)
        sigma_deg = float(sigma_deg)
        if sigma_deg <= 0:
            raise ValueError(
                f"Standard deviation for {atom_a}-{atom_b}-{atom_c}-{atom_d} must be positive"
            )
        if angle_deg is not None:
            angle_deg = float(angle_deg)
            if not -180.0 <= angle_deg <= 180.0:
                raise ValueError(
                    f"Dihedral for {atom_a}-{atom_b}-{atom_c}-{atom_d} "
                    "must be between -180 and 180 degrees"
                )
        parsed.append(
            {
                "atoms": atoms,
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
                atom_a, atom_b, atom_c, atom_d = spec["atoms"]
                missing = [name for name in (atom_a, atom_b, atom_c, atom_d) if name not in coords]
                if missing:
                    raise ValueError(
                        f"Atom {missing[0]} is not in {where}. Atoms present: {', '.join(sorted(coords))}"
                    )
                angle_deg = spec["angle_deg"]
                if angle_deg is None:
                    angle_deg = _dihedral_degrees(coords, atom_a, atom_b, atom_c, atom_d, where)
                parts = []
                for role, name in (("first", atom_a), ("next", atom_b), ("next", atom_c), ("next", atom_d)):
                    parts.append(
                        f"{role} chain {chain} resi {int(resi)} atom {name}{alt_txt}"
                    )
                lines.append(
                    "exte torsion "
                    + " ".join(parts)
                    + f" value {angle_deg:.3f} sigma {spec['sigma_deg']:.3f} period 1"
                )
    if not lines:
        raise ValueError("Dihedral restraints produced no external restraints")
    return "\n".join(lines)
