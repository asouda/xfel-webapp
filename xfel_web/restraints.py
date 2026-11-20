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
