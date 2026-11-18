#!/usr/bin/env python3
"""Build a residue-clipped summary of one reconstruction run.

The summary directory receives the light reference PDB, each iteration's
refined model with the newly placed atom, the extended omit-map fragment, and
the reconstructed final PDB. PDB files keep only the residues named with
``--resid``. Iteration files are stored with ``step##`` at the start of the
name, so the viewers list them as step01, step02, step03, and so on. VMD and ChimeraX scripts in that directory load the reference, then each
added model and its map, and the reconstructed final model last. The atom
placed in an iteration is bright green. The final model is not: nothing was
added in that file. A Maestro project archive has one entry per step, with
that step's PDB and map.

Example::

    python prepare_reconstruction_summary.py \\
        runs/reconstruct/9z35_confA_omitOEXsidechains6HOH_refmac20 \\
        --resid A602,A170,A189,A332,A333,A342,A344,C354,A705,A712,A738,A776,A814,A856
"""

from __future__ import annotations

import argparse
import math
import os
import pickle
import re
import shutil
import struct
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path


ADDED_RE = re.compile(r"_refined_added_.+_step(\d+)\.pdb$")
MAP_RE = re.compile(r"_extended_.+_10A_step(\d+)\.map$")
ITERATION_RE = re.compile(r"^r-iteration-(\d+)$")
RESIDUE_RE = re.compile(r"^([A-Za-z0-9])(\d+)$")

@dataclass(frozen=True)
class Residue:
    chain: str
    resseq: int

    def __str__(self) -> str:
        return f"{self.chain}{self.resseq}"


@dataclass(frozen=True)
class AtomRecord:
    chain: str
    resseq: int
    resname: str
    name: str
    serial: int | None

    def describe(self) -> str:
        serial = "?" if self.serial is None else str(self.serial)
        return f"{self.chain}{self.resseq} {self.resname} {self.name} serial {serial}"


@dataclass
class FilteredPdb:
    path: Path
    natoms: int
    # Kept atom with the highest PDB serial. The reconstruction gives that
    # serial to the atom it places, then inserts the record beside its residue,
    # so the placed atom is not always the last line of the file.
    highlight: AtomRecord | None
    highlight_index: int | None
    newest_source: AtomRecord | None
    newest_source_kept: bool


@dataclass
class MapCopy:
    path: Path
    stdev: float
    level: float
    # Contour in units of ``stdev``: one below the highest density in the map.
    sigma: float


@dataclass
class IterationSummary:
    step: int
    added: FilteredPdb | None = None
    map_file: MapCopy | None = None


def parse_residues(text: str) -> list[Residue]:
    """Parse ``A602,A170,C354`` into chain/residue pairs, preserving order."""
    residues: list[Residue] = []
    seen: set[tuple[str, int]] = set()
    for part in text.split(","):
        token = part.strip()
        if not token:
            continue
        match = RESIDUE_RE.fullmatch(token)
        if match is None:
            raise SystemExit(
                f"Invalid residue {token!r}. Use a chain id plus a residue "
                "number, for example A602 or C354."
            )
        chain = match.group(1)
        if chain.isalpha():
            chain = chain.upper()
        residue = Residue(chain, int(match.group(2)))
        key = (residue.chain, residue.resseq)
        if key in seen:
            continue
        seen.add(key)
        residues.append(residue)
    if not residues:
        raise SystemExit("No residues were given to --resid.")
    return residues


def _atom_record(line: str) -> AtomRecord | None:
    if not line.startswith(("ATOM", "HETATM")) or len(line) < 27:
        return None
    try:
        resseq = int(line[22:26])
    except ValueError:
        return None
    serial: int | None
    try:
        serial = int(line[6:11])
    except ValueError:
        serial = None
    return AtomRecord(
        chain=line[21],
        resseq=resseq,
        resname=line[17:20].strip() or "UNK",
        name=line[12:16].strip() or "?",
        serial=serial,
    )


def _remark_block(residues: list[Residue], source_name: str) -> list[str]:
    lines = [
        "REMARK   4",
        "REMARK   4 FILTERED BY prepare_reconstruction_summary.py",
        f"REMARK   4 SOURCE: {source_name}"[:80],
    ]
    kept = ",".join(str(residue) for residue in residues)
    prefix = "REMARK   4 KEPT: "
    while kept:
        room = 80 - len(prefix)
        chunk, kept = kept[:room], kept[room:]
        lines.append(prefix + chunk)
        prefix = "REMARK   4 "
    return lines


def _is_newer(candidate: AtomRecord, current: AtomRecord) -> bool:
    """Higher PDB serial is the atom most recently inserted."""
    candidate_serial = -1 if candidate.serial is None else candidate.serial
    current_serial = -1 if current.serial is None else current.serial
    return candidate_serial >= current_serial


def filter_pdb(
    source: Path,
    dest: Path,
    residues: list[Residue],
    *,
    keep_newest: bool = False,
) -> FilteredPdb:
    """Write ``dest`` with ATOM/HETATM records for ``residues``.

    Records stay in source order. The highlight is the kept atom with the
    highest PDB serial. For a refined_added model, ``keep_newest`` also
    retains that atom when its residue is not in ``residues``: it is the
    atom this iteration inserted, and it is the one drawn in bright green.
    """
    wanted = {(residue.chain, residue.resseq) for residue in residues}
    raw_lines = source.read_text(encoding="utf-8", errors="replace").splitlines()
    newest_source: AtomRecord | None = None
    newest_index: int | None = None
    for index, line in enumerate(raw_lines):
        record = _atom_record(line)
        if record is None:
            continue
        if newest_source is None or _is_newer(record, newest_source):
            newest_source = record
            newest_index = index

    header: list[str] = []
    body: list[str] = []
    highlight: AtomRecord | None = None
    highlight_index: int | None = None
    kept_index = 0
    keep_aniso = False
    for index, line in enumerate(raw_lines):
        if line.startswith(("CRYST1", "ORIGX", "SCALE", "MTRIX")):
            header.append(line)
            keep_aniso = False
            continue
        record = _atom_record(line)
        if record is not None:
            forced = keep_newest and index == newest_index
            keep_aniso = forced or (record.chain, record.resseq) in wanted
            if keep_aniso:
                if highlight is None or _is_newer(record, highlight):
                    highlight = record
                    highlight_index = kept_index
                kept_index += 1
                body.append(line)
            continue
        if keep_aniso and line.startswith("ANISOU"):
            body.append(line)
        keep_aniso = False

    newest_outside = (
        keep_newest
        and newest_source is not None
        and (newest_source.chain, newest_source.resseq) not in wanted
    )
    dest.parent.mkdir(parents=True, exist_ok=True)
    previous_chain: str | None = None
    with dest.open("w", encoding="utf-8", newline="\n") as handle:
        for line in _remark_block(residues, source.name):
            handle.write(line + "\n")
        if newest_outside and newest_source is not None:
            handle.write(f"REMARK   4 ALSO: {newest_source.describe()}"[:80] + "\n")
        for line in header:
            handle.write(line + "\n")
        for line in body:
            record = _atom_record(line)
            if record is not None:
                if previous_chain is not None and record.chain != previous_chain:
                    handle.write("TER\n")
                previous_chain = record.chain
            handle.write(line + "\n")
        if previous_chain is not None:
            handle.write("TER\n")
        handle.write("END\n")

    newest_source_kept = newest_source is not None and (
        (newest_source.chain, newest_source.resseq) in wanted or newest_outside
    )
    return FilteredPdb(
        path=dest,
        natoms=kept_index,
        highlight=highlight,
        highlight_index=highlight_index,
        newest_source=newest_source,
        newest_source_kept=newest_source_kept,
    )


def _ccp4_voxels(raw: bytes):
    """Return the float32 voxel array of a CCP4 map."""
    import numpy as np

    for endian in ("<", ">"):
        nx, ny, nz, mode = struct.unpack_from(endian + "4i", raw, 0)
        if mode != 2 or min(nx, ny, nz) <= 0 or max(nx, ny, nz) > 20000:
            continue
        nsymbt = struct.unpack_from(endian + "i", raw, 23 * 4)[0]
        if nsymbt < 0:
            continue
        count = nx * ny * nz
        offset = 1024 + nsymbt
        if offset + count * 4 > len(raw):
            continue
        return np.frombuffer(raw, dtype=endian + "f4", count=count, offset=offset)
    raise ValueError("Not a mode-2 CCP4 map, or the grid does not fit the file")


def choose_contour(stdev: float, peak: float) -> tuple[float, float]:
    """Return the contour in sigma and in map units.

    The contour is one standard deviation below the highest density.
    ``peak`` is the maximum voxel value, so the level is ``(peak / stdev) - 1``
    in sigma units.
    """
    if stdev <= 0:
        return 0.0, peak
    used = peak / stdev - 1.0
    return used, used * stdev


def chimerax_map_origin(path: Path) -> tuple[float, float, float] | None:
    """Return a map origin that ChimeraX would otherwise discard.

    A cut-out CCP4 map keeps the original grid index in ``ncstart``,
    ``nrstart``, and ``nsstart``, while ``mx``, ``my``, and ``mz`` describe
    only the cut-out. ChimeraX treats a start index larger than ten times
    that grid as uninitialized and places the map at the coordinate origin,
    away from the model. The returned point is the position of grid index
    0, 0, 0 in angstroms. ``None`` means ChimeraX already keeps the origin.
    """
    raw = path.read_bytes()[:1024]
    for endian in ("<", ">"):
        nc, nr, ns, mode = struct.unpack_from(endian + "4i", raw, 0)
        if mode != 2 or min(nc, nr, ns) <= 0 or max(nc, nr, ns) > 20000:
            continue
        ncstart, nrstart, nsstart = struct.unpack_from(endian + "3i", raw, 16)
        mx, my, mz = struct.unpack_from(endian + "3i", raw, 28)
        xlen, ylen, zlen, alpha, beta, gamma = struct.unpack_from(
            endian + "6f", raw, 40
        )
        mapc, mapr, maps = struct.unpack_from(endian + "3i", raw, 64)
        break
    else:
        return None
    if mx <= 0 or my <= 0 or mz <= 0 or min(xlen, ylen, zlen) <= 0:
        return None
    step = (xlen / mx, ylen / my, zlen / mz)
    if {mapc, mapr, maps} == {1, 2, 3}:
        crs_to_ijk = (mapc - 1, mapr - 1, maps - 1)
        ijk_to_crs = [0, 0, 0]
        for axis, crs in enumerate(crs_to_ijk):
            ijk_to_crs[crs] = axis
    else:
        ijk_to_crs = (0, 1, 2)
    crs_size = (nc, nr, ns)
    crs_start = (ncstart, nrstart, nsstart)
    ijk_start = [crs_start[axis] for axis in ijk_to_crs]
    data_size = [crs_size[axis] for axis in ijk_to_crs]
    limit = 10 * max(mx, my, mz, max(data_size))
    if all(abs(start) <= limit for start in ijk_start):
        return None
    origin = _skew_origin(ijk_start, step, (alpha, beta, gamma))
    if not all(
        math.isfinite(value) and (spacing <= 0 or value + spacing != value)
        for value, spacing in zip(origin, step)
    ):
        return None
    return origin


def _skew_origin(
    ijk: list[int],
    step: tuple[float, float, float],
    angles: tuple[float, float, float],
) -> tuple[float, float, float]:
    """Map a grid index to Cartesian angstroms, matching ChimeraX."""
    if tuple(angles) == (90.0, 90.0, 90.0):
        cosine_gamma = cosine_beta = cosine_1 = 0.0
        sine_gamma = cosine_2 = 1.0
    else:
        alpha, beta, gamma = (math.radians(angle) for angle in angles)
        cosine_gamma = math.cos(gamma)
        sine_gamma = math.sin(gamma)
        cosine_beta = math.cos(beta)
        cosine_alpha = math.cos(alpha)
        cosine_1 = (cosine_alpha - cosine_beta * cosine_gamma) / sine_gamma
        cosine_2 = math.sqrt(max(0.0, 1 - cosine_beta ** 2 - cosine_1 ** 2))
    scaled = [spacing * index for spacing, index in zip(step, ijk)]
    return (
        scaled[0] + cosine_gamma * scaled[1] + cosine_beta * scaled[2],
        sine_gamma * scaled[1] + cosine_1 * scaled[2],
        cosine_2 * scaled[2],
    )


def map_contour(raw: bytes) -> tuple[float, float, float]:
    """Return standard deviation, contour in sigma, and absolute level."""
    voxels = _ccp4_voxels(raw)
    stdev = float(voxels.std())
    peak = float(voxels.max())
    used, level = choose_contour(stdev, peak)
    return stdev, used, level


def copy_map(source: Path, dest: Path) -> MapCopy:
    raw = source.read_bytes()
    stdev, used, level = map_contour(raw)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(raw)
    return MapCopy(path=dest, stdev=stdev, level=level, sigma=used)


def summary_filename(source_name: str) -> str:
    """Move a trailing ``_step##`` to the front of a summary file name.

    ``model_refined_added_MN2_step01.pdb`` becomes
    ``step01_model_refined_added_MN2.pdb``. Names without that suffix, such
    as the light reference, are unchanged. The step number stays zero-padded
    so step01 sorts before step10.
    """
    match = re.search(r"_step(\d+)(?=\.[^.]+$)", source_name)
    if match is None:
        return source_name
    label = match.group(1).zfill(2)
    remainder = source_name[: match.start()] + source_name[match.end() :]
    return f"step{label}_{remainder}"


def _step_number(path: Path, pattern: re.Pattern[str]) -> int | None:
    match = pattern.search(path.name)
    if match is None:
        return None
    return int(match.group(1))


def collect_iterations(run_dir: Path) -> dict[int, dict[str, Path]]:
    """Map a step number to its added PDB and extended map."""
    found: dict[int, dict[str, Path]] = {}

    def take(step: int, kind: str, path: Path) -> None:
        slot = found.setdefault(step, {})
        previous = slot.get(kind)
        if previous is not None and previous != path:
            raise SystemExit(
                f"Step {step:02d} has two {kind} files:\n  {previous}\n  {path}"
            )
        slot[kind] = path

    for folder in sorted(run_dir.iterdir()):
        if not folder.is_dir() or ITERATION_RE.fullmatch(folder.name) is None:
            continue
        folder_step = int(ITERATION_RE.fullmatch(folder.name).group(1))
        for path in sorted(folder.iterdir()):
            if not path.is_file():
                continue
            added_step = _step_number(path, ADDED_RE)
            map_step = _step_number(path, MAP_RE)
            if added_step is not None:
                take(added_step, "added", path)
                file_step = added_step
            elif map_step is not None:
                take(map_step, "map", path)
                file_step = map_step
            else:
                continue
            if file_step != folder_step:
                print(
                    f"warning: {path.name} is in {folder.name} "
                    f"but its step number is {file_step:02d}",
                    file=sys.stderr,
                )
    return found


def collect_finals(run_dir: Path) -> list[Path]:
    """Reconstructed models written beside the iteration folders."""
    return sorted(
        path
        for path in run_dir.glob("*_reconstructed_final.pdb")
        if path.is_file()
    )


def _map_caption(step: int, item: MapCopy) -> str:
    """One-line description of the contour used for a summary map."""
    maximum = item.sigma + 1.0
    return (
        f"# step {step:02d} map: {item.sigma:.4g} sigma "
        f"(1 below the maximum {maximum:.4g}), "
        f"sd={item.stdev:.6g}, level={item.level:.6g}"
    )


def _data_name(path: Path) -> str:
    """File name inside the summary directory."""
    return path.name


def _tcl_quote(name: str) -> str:
    if "{" in name or "}" in name:
        raise SystemExit(f"VMD cannot quote this path safely: {name}")
    return "{" + name + "}"


def _cx_quote(name: str) -> str:
    return '"' + name.replace("\\", "\\\\").replace('"', '\\"') + '"'


def write_vmd(
    path: Path,
    references: list[FilteredPdb],
    steps: list[IterationSummary],
    finals: list[FilteredPdb],
) -> None:
    lines = [
        "# loadVMD.cmd",
        "# Generated by prepare_reconstruction_summary.py",
        "#",
        "# From this directory:",
        "#   vmd -e loadVMD.cmd",
        "#",
        "# PDB and map paths are relative to this directory. vmd -e does not",
        "# record the script location, so start VMD from this directory. Sourcing",
        "# the script by its full path does change to that directory.",
        "#",
        "# The light reference is ball-and-stick in silver, material Ghost.",
        "# Each refined_added model is ball-and-stick colored by atom name.",
        "# The atom placed in that iteration, the one with the highest PDB",
        "# serial, is bright green.",
        "# The reconstructed final model is loaded last, colored by atom name,",
        "# with no bright-green atom.",
        "# Each map is a white wireframe at one sigma below the highest",
        "# density in that map.",
        "# No boundary box.",
        "# Only the light reference is drawn. The other molecules are loaded",
        "# with Drawn off.",
        "# Edit BALL_SCALE and BOND_RADIUS to change the ball-and-stick size.",
        "",
        "set BALL_SCALE 0.45",
        "set BOND_RADIUS 0.15",
        "set SILVER [colorinfo index silver]",
        "set WHITE [colorinfo index white]",
        "set BRIGHT [colorinfo index green2]",
        "color change rgb $BRIGHT 0.0 1.0 0.0",
        "color Display Background black",
        "display projection Orthographic",
        "display depthcue off",
        "axes location Off",
        "",
        "proc clear_reps {molid} {",
        "    set nreps [molinfo $molid get numreps]",
        "    for {set i 0} {$i < $nreps} {incr i} {",
        "        mol delrep 0 $molid",
        "    }",
        "}",
        "",
        "proc show_ballstick {molid selection colorid material} {",
        "    global BALL_SCALE BOND_RADIUS",
        "    mol representation CPK $BALL_SCALE $BOND_RADIUS 16 16",
        "    mol color ColorID $colorid",
        "    mol selection $selection",
        "    mol material $material",
        "    mol addrep $molid",
        "}",
        "",
        "proc load_reference {filepath} {",
        "    global SILVER",
        "    set molid [mol new $filepath type pdb waitfor all]",
        "    clear_reps $molid",
        "    show_ballstick $molid all $SILVER Ghost",
        "    mol rename $molid [file tail $filepath]",
        "    return $molid",
        "}",
        "",
        "proc load_added {filepath highlight} {",
        "    global BRIGHT",
        "    set molid [mol new $filepath type pdb waitfor all]",
        "    clear_reps $molid",
        "    set natoms [molinfo $molid get numatoms]",
        "    if {$natoms <= 0} {",
        "        mol rename $molid [file tail $filepath]",
        "        molinfo $molid set drawn 0",
        "        return $molid",
        "    }",
        "    if {$highlight < 0 || $highlight >= $natoms} {",
        "        set highlight [expr {$natoms - 1}]",
        "    }",
        "    if {$natoms == 1} {",
        "        show_ballstick $molid all $BRIGHT Opaque",
        "    } else {",
        "        mol representation CPK $::BALL_SCALE $::BOND_RADIUS 16 16",
        "        mol color Name",
        "        mol selection \"not index $highlight\"",
        "        mol material Opaque",
        "        mol addrep $molid",
        "        show_ballstick $molid \"index $highlight\" $BRIGHT Opaque",
        "    }",
        "    mol rename $molid [file tail $filepath]",
        "    molinfo $molid set drawn 0",
        "    return $molid",
        "}",
        "",
        "proc load_colored {filepath} {",
        "    set molid [mol new $filepath type pdb waitfor all]",
        "    clear_reps $molid",
        "    set natoms [molinfo $molid get numatoms]",
        "    if {$natoms > 0} {",
        "        mol representation CPK $::BALL_SCALE $::BOND_RADIUS 16 16",
        "        mol color Name",
        "        mol selection all",
        "        mol material Opaque",
        "        mol addrep $molid",
        "    }",
        "    mol rename $molid [file tail $filepath]",
        "    molinfo $molid set drawn 0",
        "    return $molid",
        "}",
        "",
        "# Isosurface arguments: isovalue volid show draw step size",
        "# show 0 draws the surface without the boundary box.",
        "# draw 1 is wireframe.",
        "proc load_map {filepath level} {",
        "    global WHITE",
        "    set molid [mol new $filepath type ccp4 waitfor all]",
        "    clear_reps $molid",
        "    mol representation Isosurface $level 0 0 1 1 2",
        "    mol color ColorID $WHITE",
        "    mol selection all",
        "    mol material Opaque",
        "    mol addrep $molid",
        "    mol rename $molid [file tail $filepath]",
        "    molinfo $molid set drawn 0",
        "    return $molid",
        "}",
        "",
        "set _script_path [info script]",
        "if {[file pathtype $_script_path] eq \"absolute\"} {",
        "    cd [file dirname $_script_path]",
        "}",
        "",
    ]

    def load_structure(pdb: FilteredPdb, comment: str, *, highlight: bool) -> None:
        if not highlight:
            lines.append(f"# {comment}")
            lines.append(f"load_colored {_tcl_quote(_data_name(pdb.path))}")
            return
        placed = pdb.highlight
        note = placed.describe() if placed else "none"
        index = 0 if pdb.highlight_index is None else pdb.highlight_index
        lines.append(f"# {comment}: {note}")
        lines.append(f"load_added {_tcl_quote(_data_name(pdb.path))} {index}")

    for reference in references:
        lines.append(f"load_reference {_tcl_quote(_data_name(reference.path))}")
    for step in steps:
        if step.added is not None:
            load_structure(
                step.added, f"step {step.step:02d} placed atom", highlight=True
            )
        if step.map_file is not None:
            lines.append(_map_caption(step.step, step.map_file))
            lines.append(
                f"load_map {_tcl_quote(_data_name(step.map_file.path))} "
                f"{step.map_file.level:.8f}"
            )
    for final in finals:
        load_structure(
            final,
            "reconstructed final, colored by atom name, no highlight",
            highlight=False,
        )
    lines.extend([
        "",
        "display resetview",
        "display update",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def _cx_model_lines(
    model_id: int,
    *,
    color: str,
    last: AtomRecord | None = None,
) -> list[str]:
    spec = f"#{model_id}"
    lines = [
        f"style {spec} ball",
        f"size {spec} ballScale 0.175 stickRadius 0.075",
        f"hide {spec} cartoons",
        f"hide {spec} pseudobonds",
        f"show {spec} atoms",
        f"color {spec} {color}",
    ]
    if last is not None and last.serial is not None:
        lines.append(f"color {spec}@@serial_number={last.serial} lime")
    elif last is not None:
        lines.append(
            f"color {spec}/{last.chain}:{last.resseq}@{last.name} lime"
        )
    return lines


def write_chimerax(
    path: Path,
    references: list[FilteredPdb],
    steps: list[IterationSummary],
    finals: list[FilteredPdb],
) -> None:
    """Write a ChimeraX command script.

    ChimeraX only auto-detects command scripts that end in ``.cxc``. This
    file keeps the requested name, so it is opened with an explicit format.
    """
    lines = [
        "# loadChimeraX.cmd",
        "# Generated by prepare_reconstruction_summary.py",
        "#",
        "# ChimeraX command script. From this directory:",
        "#   chimerax --cmd 'open loadChimeraX.cmd format \"cmd\"'",
        "#",
        "# Or, in the ChimeraX command line:",
        '#   open loadChimeraX.cmd format cmd',
        "#",
        "# PDB and map paths are relative to this script. ChimeraX runs the",
        "# script from the directory that contains it.",
        "#",
        "# The light reference is ball-and-stick in silver.",
        "# Each refined_added model is ball-and-stick colored by element",
        "# (ChimeraX has no PDB atom-name scheme; element is the matching one).",
        "# The atom placed in that iteration, the one with the highest PDB",
        "# serial, is lime (#00ff00).",
        "# The reconstructed final model is loaded last, colored by element,",
        "# with no lime atom.",
        "# Each map is a white mesh at one sigma below the highest",
        "# density in that map.",
        "# No outline box.",
        "# The background is black so the white mesh stays visible.",
        "",
        "set bgColor black",
        "lighting soft",
        "",
    ]
    model_id = 1
    for reference in references:
        lines.append(f"open {_cx_quote(_data_name(reference.path))}")
        lines.extend(_cx_model_lines(model_id, color="silver"))
        lines.append("")
        model_id += 1
    for step in steps:
        if step.added is not None:
            placed = step.added.highlight
            note = placed.describe() if placed else "no atoms"
            lines.append(f"# step {step.step:02d} placed atom: {note}")
            lines.append(f"open {_cx_quote(_data_name(step.added.path))}")
            lines.extend(_cx_model_lines(model_id, color="byelement", last=placed))
            lines.append("")
            model_id += 1
        if step.map_file is not None:
            lines.append(_map_caption(step.step, step.map_file))
            lines.append(f"open {_cx_quote(_data_name(step.map_file.path))}")
            origin = chimerax_map_origin(step.map_file.path)
            # meshLighting false keeps the wireframe white. Lit mesh lines
            # shade toward black and disappear on the black background.
            # The origin puts a cut-out map back on the model. ChimeraX
            # drops that grid index when it is large compared with the
            # cut-out size.
            placed = ""
            if origin is not None:
                placed = " origin " + ",".join(f"{value:.6f}" for value in origin)
            lines.append(
                f"volume #{model_id} style mesh color white step 1 "
                f"level {step.map_file.level:.8f} showOutlineBox false "
                f"meshLighting false{placed}"
            )
            lines.append("")
            model_id += 1
    for final in finals:
        lines.append("# reconstructed final, colored by element, no highlight")
        lines.append(f"open {_cx_quote(_data_name(final.path))}")
        lines.extend(_cx_model_lines(model_id, color="byelement"))
        lines.append("")
        model_id += 1
    lines.append("hide pseudobonds")
    lines.append("view")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


# VMD Name coloring, as RGB 0-255. Carbon is cyan, as in VMD.
_NAME_RGB = {
    1: (255, 255, 255),
    6: (0, 255, 255),
    7: (0, 0, 255),
    8: (255, 0, 0),
    15: (255, 128, 0),
    16: (255, 255, 0),
    20: (144, 144, 160),
    25: (178, 102, 255),
}
_SILVER = (192, 192, 192)
_BRIGHT_GREEN = (0, 255, 0)


def _schrodinger_run() -> Path | None:
    """Return the Schrödinger ``run`` launcher, if one is installed."""
    candidates: list[Path] = []
    env = os.environ.get("SCHRODINGER")
    if env:
        candidates.append(Path(env) / "run")
    opt = Path("/opt/schrodinger")
    if opt.is_dir():
        candidates.extend(sorted(opt.glob("suites*/run")))
    apps = Path("/Applications")
    if apps.is_dir():
        for folder in sorted(apps.glob("SchrodingerSuites*")):
            candidates.append(folder / "run")
            candidates.append(folder / "Contents/MacOS/run")
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def write_maestro_project(
    path: Path,
    references: list[FilteredPdb],
    steps: list[IterationSummary],
    finals: list[FilteredPdb],
) -> None:
    """Write a Maestro project archive with one entry per reconstruction step."""
    try:
        import schrodinger.infra.mmproj  # noqa: F401
    except ImportError:
        _write_maestro_project_via_schrodinger(path, references, steps, finals)
        return
    _write_maestro_project(path, references, steps, finals)


def _write_maestro_project_via_schrodinger(
    path: Path,
    references: list[FilteredPdb],
    steps: list[IterationSummary],
    finals: list[FilteredPdb],
) -> None:
    runner = _schrodinger_run()
    if runner is None:
        raise SystemExit(
            "Schrödinger Maestro is required to write the project archive, "
            "and the run launcher was not found."
        )
    with tempfile.NamedTemporaryFile(suffix=".pkl", delete=False) as handle:
        pickle.dump((str(path), references, steps, finals), handle)
        payload = handle.name
    code = (
        "import os, pickle, sys\n"
        "from pathlib import Path\n"
        f"sys.path.insert(0, {str(Path(__file__).resolve().parent)!r})\n"
        "import prepare_reconstruction_summary as tool\n"
        "archive, references, steps, finals = pickle.load("
        "open(os.environ['XFEL_MAESTRO_PICKLE'], 'rb'))\n"
        "tool._write_maestro_project("
        "Path(archive), references, steps, finals)\n"
    )
    env = os.environ.copy()
    env["XFEL_MAESTRO_PICKLE"] = payload
    try:
        subprocess.check_call([str(runner), "python3", "-c", code], env=env)
    finally:
        os.unlink(payload)


def _style_structure(
    st, pdb: FilteredPdb, *, silver: bool, highlight: bool
) -> None:
    """Ball-and-stick. Silver for the reference, VMD Name colors otherwise."""
    from schrodinger import structure
    from schrodinger.infra import mm

    st.applyStyle(structure.ATOM_BALLNSTICK, structure.BOND_BALLNSTICK)
    for atom in st.atom:
        mm.mmctg_atom_set_ribbon_style(
            st.handle, atom.index, mm.MMCT_RIBBON_STYLE_NONE
        )
        if silver:
            atom.color = _SILVER
        else:
            atom.color = _NAME_RGB.get(atom.atomic_number, (200, 200, 200))
    if silver or not highlight or pdb.highlight is None:
        return
    chosen = None
    serial = pdb.highlight.serial
    if serial is not None:
        for atom in st.atom:
            if atom.property.get("i_pdb_PDB_serial") == serial:
                chosen = atom
                break
    if chosen is None:
        index = (pdb.highlight_index or 0) + 1
        if 1 <= index <= st.atom_total:
            chosen = st.atom[index]
    if chosen is not None:
        chosen.color = _BRIGHT_GREEN


def _add_map(proj, row, item: MapCopy, label: str) -> int:
    """Attach a white mesh of ``item`` to a project row. Return vertex count."""
    from schrodinger.infra import mmproj, mmsurf

    ivol, center, _kind, _comment = mmsurf.mmvisio_read_ccp4(str(item.path))
    entry_index = mmproj.mmproj_project_entry_id_to_index(
        proj.handle, int(row.entry_id)
    )
    mmproj.mmproj_index_entry_add_volume(proj.handle, entry_index, ivol)
    level = float(item.level)
    mmsurf.mmsurf_set_gen_isolines(True)
    surface = mmsurf.mmsurf_new_isosurface(
        ivol,
        level,
        False,
        0.0,
        float(center[0]),
        float(center[1]),
        float(center[2]),
        0,
    )
    mmsurf.mmsurf_generate_molecular_isolines(ivol, surface)
    mmsurf.mmsurf_set_name(surface, label)
    mmsurf.mmsurf_set_rgb_color(surface, (255, 255, 255))
    mmsurf.mmsurf_set_style(surface, mmsurf.MMSURF_STYLE_MESH)
    mmsurf.mmsurf_set_transparency(surface, 0)
    mmsurf.mmsurf_set_visibility(surface, True)
    volume_name = mmsurf.mmvol_get_name(ivol)
    mmsurf.mmsurf_set_volume_name(surface, volume_name)
    mmsurf.mmsurf_set_scheme_volume_name(surface, volume_name)
    mmsurf.mmvol_set_visible(ivol, False)
    mmproj.mmproj_index_entry_add_surface(proj.handle, entry_index, surface, None)
    return int(mmsurf.mmsurf_get_num_vertices(surface))


def _import_entry(
    proj, pdb: FilteredPdb, title: str, *, silver: bool, highlight: bool
):
    from schrodinger import structure

    model = structure.Structure.read(str(pdb.path))
    _style_structure(model, pdb, silver=silver, highlight=highlight)
    return proj.importStructure(model, name=title)


def _write_maestro_project(
    path: Path,
    references: list[FilteredPdb],
    steps: list[IterationSummary],
    finals: list[FilteredPdb],
) -> None:
    from schrodinger import project
    from schrodinger.infra import mm, mmproj, mmsurf

    path = path.expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    mmsurf.mmvisio_initialize(mm.MMERR_DEFAULT_HANDLER)
    work = Path(tempfile.mkdtemp(prefix="xfel-maestro-"))
    project_dir = work / "loadMaestro.prj"
    try:
        handle = mmproj.mmproj_project_new(str(project_dir))
        opened = project.Project(project_handle=handle)
        reference_row = None
        for reference in references:
            reference_row = _import_entry(
                opened, reference, "reference", silver=True, highlight=False
            )
            print(f"maestro: reference ({reference.path.name})")
        final_used = False
        for step in steps:
            model = step.added
            highlight = True
            if model is None and step.map_file is not None and finals:
                model = finals[-1]
                highlight = False
                final_used = True
            if model is None:
                continue
            title = f"step{step.step:02d}"
            row = _import_entry(
                opened, model, title, silver=False, highlight=highlight
            )
            vertices = 0
            if step.map_file is not None:
                vertices = _add_map(opened, row, step.map_file, f"{title} map")
            note = "" if highlight else ", no highlight"
            print(
                f"maestro: {title} ({model.path.name}{note}, "
                f"mesh vertices {vertices})"
            )
        if not final_used:
            for final in finals:
                _import_entry(
                    opened, final, "final", silver=False, highlight=False
                )
                print(f"maestro: final ({final.path.name}, no highlight)")
        if reference_row is not None:
            reference_row.includeOnly()
        opened.close()
        if path.exists():
            path.unlink()
        written = Path(project.zip_project(str(project_dir), str(path.parent)))
        if written.resolve() != path.resolve():
            written.replace(path)
        print(f"maestro: {path}")
    finally:
        shutil.rmtree(work, ignore_errors=True)



def prepare(run_dir: Path, residues: list[Residue]) -> Path:
    run_dir = run_dir.expanduser().resolve()
    if not run_dir.is_dir():
        raise SystemExit(f"Run directory not found: {run_dir}")
    inputs = run_dir / "inputs"
    if not inputs.is_dir():
        raise SystemExit(f"No inputs directory in {run_dir}")
    light_paths = sorted(inputs.glob("light_*.pdb"))
    if not light_paths:
        raise SystemExit(f"No light_*.pdb file in {inputs}")

    iterations = collect_iterations(run_dir)
    if not iterations:
        raise SystemExit(f"No r-iteration-## directories in {run_dir}")
    final_paths = collect_finals(run_dir)
    if not final_paths:
        raise SystemExit(f"No *_reconstructed_final.pdb file in {run_dir}")

    summary = run_dir / "summary"
    summary.mkdir(parents=True, exist_ok=True)
    print(f"summary: {summary}")
    for stale in summary.glob("step*_refined.pdb"):
        stale.unlink()

    references: list[FilteredPdb] = []
    for source in light_paths:
        dest = summary / source.name
        filtered = filter_pdb(source, dest, residues)
        references.append(filtered)
        print(f"reference: {dest.name} ({filtered.natoms} atoms)")

    steps: list[IterationSummary] = []
    for step in sorted(iterations):
        slot = iterations[step]
        item = IterationSummary(step=step)
        if "added" in slot:
            source = slot["added"]
            dest = summary / summary_filename(source.name)
            item.added = filter_pdb(source, dest, residues, keep_newest=True)
            placed = (
                item.added.highlight.describe()
                if item.added.highlight
                else "no atoms"
            )
            print(
                f"step {step:02d} added:   {dest.name} "
                f"({item.added.natoms} atoms, green {placed})"
            )
            outside = (
                item.added.newest_source is not None
                and (item.added.newest_source.chain, item.added.newest_source.resseq)
                not in {(residue.chain, residue.resseq) for residue in residues}
            )
            if outside:
                print(
                    "note: the atom placed in "
                    f"{source.name} is {item.added.newest_source.describe()}, "
                    "which is outside --resid. It is kept and highlighted.",
                    file=sys.stderr,
                )
        if "map" in slot:
            source = slot["map"]
            dest = summary / summary_filename(source.name)
            print(f"step {step:02d} map:     {dest.name}")
            item.map_file = copy_map(source, dest)
            print(f"         {_map_caption(step, item.map_file)[2:]}")
        if item.added is None and item.map_file is None:
            continue
        steps.append(item)

    finals: list[FilteredPdb] = []
    for source in final_paths:
        dest = summary / summary_filename(source.name)
        filtered = filter_pdb(source, dest, residues)
        finals.append(filtered)
        print(f"final: {dest.name} ({filtered.natoms} atoms, no highlight)")
        if filtered.newest_source is not None and not filtered.newest_source_kept:
            print(
                "warning: the newest atom in "
                f"{source.name} is {filtered.newest_source.describe()}, "
                "which is outside --resid and is omitted from the summary PDB.",
                file=sys.stderr,
            )

    write_vmd(summary / "loadVMD.cmd", references, steps, finals)
    write_chimerax(summary / "loadChimeraX.cmd", references, steps, finals)
    write_maestro_project(
        summary / "loadMaestro.prjzip", references, steps, finals
    )
    print(f"wrote {summary / 'loadVMD.cmd'}")
    print(f"wrote {summary / 'loadChimeraX.cmd'}")
    print(f"wrote {summary / 'loadMaestro.prjzip'}")
    return summary


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Copy a reconstruction run into a summary folder, keep selected "
            "PDB residues, and write VMD, ChimeraX, and Maestro load scripts."
        )
    )
    parser.add_argument(
        "run_dir",
        type=Path,
        help="Reconstruction directory containing inputs/ and r-iteration-##/",
    )
    parser.add_argument(
        "--resid",
        required=True,
        help=(
            "Residues to keep, as chain id plus residue number, comma-separated. "
            "Example: A602,A170,A189,C354"
        ),
    )
    args = parser.parse_args(argv)
    residues = parse_residues(args.resid)
    prepare(args.run_dir, residues)


if __name__ == "__main__":
    main()
