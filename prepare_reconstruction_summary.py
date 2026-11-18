#!/usr/bin/env python3
"""Build a residue-clipped summary of one reconstruction run.

The summary directory receives the light reference PDB, each iteration's
refined model with the newly placed atom, the extended omit-map fragment, and
the reconstructed final PDB. PDB files keep only the residues named with
``--resid``. Iteration files are stored with ``step##`` at the start of the
name, so the viewers list them as step01, step02, step03, and so on. VMD,
ChimeraX, and Maestro scripts in that directory load the reference, then each
added model and its map, and the reconstructed final model last.

Example::

    python prepare_reconstruction_summary.py \\
        runs/reconstruct/9z35_confA_omitOEXsidechains6HOH_refmac20 \\
        --resid A602,A170,A189,A332,A333,A342,A344,C354,A705,A712,A738,A776,A814,A856
"""

from __future__ import annotations

import argparse
import re
import shutil
import struct
import sys
from dataclasses import dataclass
from pathlib import Path


ADDED_RE = re.compile(r"_refined_added_.+_step(\d+)\.pdb$")
MAP_RE = re.compile(r"_extended_.+_10A_step(\d+)\.map$")
ITERATION_RE = re.compile(r"^r-iteration-(\d+)$")
RESIDUE_RE = re.compile(r"^([A-Za-z0-9])(\d+)$")

# PyMOL figures for this workflow contour the omit map at 7 sigma.
DEFAULT_MAP_SIGMA = 7.0


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


def filter_pdb(source: Path, dest: Path, residues: list[Residue]) -> FilteredPdb:
    """Write ``dest`` with only ATOM/HETATM records for ``residues``.

    Records stay in source order. The highlight is the kept atom with the
    highest PDB serial: that is the atom this iteration placed, when its
    residue is retained.
    """
    wanted = {(residue.chain, residue.resseq) for residue in residues}
    header: list[str] = []
    body: list[str] = []
    newest_source: AtomRecord | None = None
    highlight: AtomRecord | None = None
    highlight_index: int | None = None
    kept_index = 0
    keep_aniso = False

    with source.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n\r")
            if line.startswith(("CRYST1", "ORIGX", "SCALE", "MTRIX")):
                header.append(line)
                keep_aniso = False
                continue
            record = _atom_record(line)
            if record is not None:
                if newest_source is None or _is_newer(record, newest_source):
                    newest_source = record
                keep_aniso = (record.chain, record.resseq) in wanted
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

    dest.parent.mkdir(parents=True, exist_ok=True)
    previous_chain: str | None = None
    with dest.open("w", encoding="utf-8", newline="\n") as handle:
        for line in _remark_block(residues, source.name):
            handle.write(line + "\n")
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
        newest_source.chain,
        newest_source.resseq,
    ) in wanted
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


def map_stdev(raw: bytes) -> float:
    """Population standard deviation of a CCP4 map, in map units."""
    voxels = _ccp4_voxels(raw)
    stdev = float(voxels.std())
    return stdev


def copy_map(source: Path, dest: Path, sigma: float) -> MapCopy:
    raw = source.read_bytes()
    stdev = map_stdev(raw)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(raw)
    return MapCopy(path=dest, stdev=stdev, level=sigma * stdev)


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
    sigma: float,
) -> None:
    lines = [
        "# loadVMD",
        "# Generated by prepare_reconstruction_summary.py",
        "#",
        "# From this directory:",
        "#   vmd -e loadVMD",
        "#",
        "# PDB and map paths are relative to this directory. vmd -e does not",
        "# record the script location, so start VMD from this directory. Sourcing",
        "# the script by its full path does change to that directory.",
        "#",
        "# The light reference is ball-and-stick in silver, material Ghost.",
        "# Each refined_added model is ball-and-stick colored by atom name.",
        "# The atom placed in that iteration, the one with the highest PDB",
        "# serial, is bright green.",
        "# The reconstructed final model is loaded last, in that same representation.",
        f"# Each map is a white wireframe at {sigma:g} sigma, with no boundary box.",
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

    def load_structure(pdb: FilteredPdb, comment: str) -> None:
        placed = pdb.highlight
        note = placed.describe() if placed else "none"
        index = 0 if pdb.highlight_index is None else pdb.highlight_index
        lines.append(f"# {comment}: {note}")
        lines.append(f"load_added {_tcl_quote(_data_name(pdb.path))} {index}")

    for reference in references:
        lines.append(f"load_reference {_tcl_quote(_data_name(reference.path))}")
    for step in steps:
        if step.added is not None:
            load_structure(step.added, f"step {step.step:02d} placed atom")
        if step.map_file is not None:
            lines.append(
                f"# step {step.step:02d} map: {sigma:g} sigma, "
                f"sd={step.map_file.stdev:.6g}, level={step.map_file.level:.6g}"
            )
            lines.append(
                f"load_map {_tcl_quote(_data_name(step.map_file.path))} "
                f"{step.map_file.level:.8f}"
            )
    for final in finals:
        load_structure(final, "reconstructed final, highest PDB serial")
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
        f"size {spec} ballScale 0.35 stickRadius 0.15",
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
    sigma: float,
) -> None:
    """Write a ChimeraX command script.

    ChimeraX only auto-detects command scripts that end in ``.cxc``. This
    file keeps the requested name, so it is opened with an explicit format.
    """
    lines = [
        "# loadChimeraX",
        "# Generated by prepare_reconstruction_summary.py",
        "#",
        "# ChimeraX command script. From this directory:",
        "#   chimerax --cmd 'open loadChimeraX format \"cmd\"'",
        "#",
        "# Or, in the ChimeraX command line:",
        '#   open loadChimeraX format cmd',
        "#",
        "# PDB and map paths are relative to this script. ChimeraX runs the",
        "# script from the directory that contains it.",
        "#",
        "# The light reference is ball-and-stick in silver.",
        "# Each refined_added model is ball-and-stick colored by element",
        "# (ChimeraX has no PDB atom-name scheme; element is the matching one).",
        "# The atom placed in that iteration, the one with the highest PDB",
        "# serial, is lime (#00ff00).",
        "# The reconstructed final model is loaded last, in that same representation.",
        f"# Each map is a white mesh at {sigma:g} sigma, with no outline box.",
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
            lines.append(
                f"# step {step.step:02d} map: {sigma:g} sigma, "
                f"sd={step.map_file.stdev:.6g}, level={step.map_file.level:.6g}"
            )
            lines.append(f"open {_cx_quote(_data_name(step.map_file.path))}")
            lines.append(
                f"volume #{model_id} style mesh color white step 1 "
                f"level {step.map_file.level:.8f} showOutlineBox false"
            )
            lines.append("")
            model_id += 1
    for final in finals:
        placed = final.highlight
        note = placed.describe() if placed else "no atoms"
        lines.append(f"# reconstructed final, highest PDB serial: {note}")
        lines.append(f"open {_cx_quote(_data_name(final.path))}")
        lines.extend(_cx_model_lines(model_id, color="byelement", last=placed))
        lines.append("")
        model_id += 1
    lines.append("hide pseudobonds")
    lines.append("view")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


MAESTRO_SCRIPT = """\
# loadMaestro
# Generated by prepare_reconstruction_summary.py
#
# Run inside Schrödinger Maestro from this directory, in the command area:
#
#   pythoneval exec(open("loadMaestro", encoding="utf-8").read())
#
# PDB and map paths are relative to this script's directory.
#
# The light reference is ball-and-stick in silver (RGB 192,192,192).
# Each refined_added model is ball-and-stick with the Element scheme
# (Maestro has no PDB atom-name scheme). The atom placed in that
# iteration, the one with the highest PDB serial, is bright green
# (RGB 0,255,0).
# The reconstructed final model is loaded last, in that same representation.
# Each map is a white mesh at __SIGMA__ sigma, with no volume box.
# A white mesh is difficult to see on a white workspace. Set
# Preferences > Workspace > Background to black before or after loading.

from pathlib import Path

from schrodinger import maestro
from schrodinger.structutils.color import Color
from schrodinger.surface import Style


def summary_dir():
    try:
        return Path(__file__).resolve().parent
    except NameError:
        return Path.cwd()


def data_file(name):
    # Resolve a summary-relative PDB or map name.
    path = summary_dir() / name
    if not path.is_file():
        raise FileNotFoundError("Missing summary file: %s" % path)
    return str(path)

MAP_SIGMA = __SIGMA_REPR__
REFERENCES = __REFERENCES__
STEPS = __STEPS__
FINALS = __FINALS__


def maestro_option(name):
    return maestro.get_command_option("entryimport", name)


def restore_entryimport(saved):
    maestro.command(
        "entryimport format={format} wsreplace={wsreplace} "
        "wsinclude={wsinclude} all={all} creategroups={creategroups}".format(
            **saved
        )
    )


def import_structure(filepath):
    project = maestro.project_table_get()
    before = {row.entry_id for row in project.all_rows}
    title = filepath.rsplit("/", 1)[-1]
    maestro.command(
        'entryimport entry_title="{0}" "{1}"'.format(title, filepath)
    )
    project = maestro.project_table_get()
    new_rows = [row for row in project.all_rows if row.entry_id not in before]
    if not new_rows:
        raise RuntimeError("Maestro did not import %s" % filepath)
    row = new_rows[-1]
    maestro.command("entrywsinclude entry %s" % row.entry_id)
    return row


def show_ballnstick(entry_id):
    asl = "entry.id %s" % entry_id
    maestro.command("displayatom %s" % asl)
    maestro.command("ribbon style=none %s" % asl)
    maestro.command("repatom rep=ballnstick %s" % asl)
    maestro.command("repbond rep=ballnstick")
    maestro.command("repatombonds %s" % asl)


def paint(entry_id, red, green, blue, extra_asl=""):
    asl = "entry.id %s" % entry_id
    if extra_asl:
        asl = "(%s AND %s)" % (asl, extra_asl)
    maestro.command(
        "coloratomrgb red=%d green=%d blue=%d %s" % (red, green, blue, asl)
    )


def paint_last_atom(row, chain, resseq, atom_name, fallback_index):
    structure = row.getStructure()
    chosen = None
    for atom in structure.atom:
        if (
            atom.chain.strip() == chain
            and int(atom.resnum) == int(resseq)
            and atom.pdbname.strip() == atom_name
        ):
            chosen = atom.index
    if chosen is None:
        chosen = fallback_index
    paint(row.entry_id, 0, 255, 0, "atom.entrynum %s" % chosen)


def style_map(row, map_path, level):
    before = set(row.surface)
    maestro.command('visimport entry="%s" "%s"' % (row.entry_id, map_path))
    project = maestro.project_table_get()
    row = project.getRow(row.entry_id)
    names = [name for name in row.surface if name not in before]
    if not names:
        surface = row.surface.addFromFile(
            map_path, name=map_path.rsplit("/", 1)[-1]
        )
        names = [surface.name]
        row = project.getRow(row.entry_id)
    surfaces = [row.surface[name] for name in names]
    positives = [surface for surface in surfaces if surface.isovalue >= 0]
    shown = positives or surfaces[:1]
    for surface in surfaces:
        if surface not in shown:
            surface.visible = False
    for surface in shown:
        surface.style = Style.mesh
        surface.setColoring(Color("white"))
        surface.isovalue = float(level)
        surface.visible = True
        maestro.command(
            'surfacedisplay entry=%s display=true "%s"'
            % (row.entry_id, surface.name)
        )
        try:
            maestro.command(
                'surfacesetviewasl entry=%s surface="%s" use=false'
                % (row.entry_id, surface.name)
            )
        except Exception as exc:
            print("surfacesetviewasl skipped: %s" % exc)
        volume_name = surface.volume_name
        if volume_name:
            maestro.command(
                'volumedisplay entry=%s display=false bounding_box=false "%s"'
                % (row.entry_id, volume_name)
            )
    return row


def load_summary():
    saved = {
        "format": maestro_option("format"),
        "wsreplace": maestro_option("wsreplace"),
        "wsinclude": maestro_option("wsinclude"),
        "all": maestro_option("all"),
        "creategroups": maestro_option("creategroups"),
    }
    try:
        maestro.command("entryimport wsreplace=false")
        maestro.command("entryimport wsinclude=true")
        maestro.command("entryimport all=true")
        reference_row = None
        for filepath in REFERENCES:
            reference_row = import_structure(data_file(filepath))
            show_ballnstick(reference_row.entry_id)
            paint(reference_row.entry_id, 192, 192, 192)
        for step in STEPS:
            added = step["added"]
            row = reference_row
            if added is not None:
                row = import_structure(data_file(added["path"]))
                show_ballnstick(row.entry_id)
                maestro.command(
                    'colorscheme scheme="Element" entry.id %s' % row.entry_id
                )
                if added["atom"]:
                    paint_last_atom(
                        row,
                        added["chain"],
                        added["resseq"],
                        added["atom"],
                        added["fallback_index"],
                    )
            map_info = step["map"]
            if map_info is not None and row is not None:
                style_map(row, data_file(map_info["path"]), map_info["level"])
        for final in FINALS:
            row = import_structure(data_file(final["path"]))
            show_ballnstick(row.entry_id)
            maestro.command(
                'colorscheme scheme="Element" entry.id %s' % row.entry_id
            )
            if final["atom"]:
                paint_last_atom(
                    row,
                    final["chain"],
                    final["resseq"],
                    final["atom"],
                    final["fallback_index"],
                )
        maestro.command("fit all")
    finally:
        restore_entryimport(saved)


load_summary()
"""


def _entry_payload(pdb: FilteredPdb) -> dict:
    if pdb.highlight is not None:
        return {
            "path": _data_name(pdb.path),
            "chain": pdb.highlight.chain,
            "resseq": pdb.highlight.resseq,
            "atom": pdb.highlight.name,
            "fallback_index": (pdb.highlight_index or 0) + 1,
        }
    return {
        "path": _data_name(pdb.path),
        "chain": "",
        "resseq": 0,
        "atom": "",
        "fallback_index": pdb.natoms,
    }


def write_maestro(
    path: Path,
    references: list[FilteredPdb],
    steps: list[IterationSummary],
    finals: list[FilteredPdb],
    sigma: float,
) -> None:
    """Write a Python script that must be executed inside Maestro."""
    payload_steps = []
    for step in steps:
        added = _entry_payload(step.added) if step.added is not None else None
        map_info = None
        if step.map_file is not None:
            map_info = {
                "path": _data_name(step.map_file.path),
                "level": step.map_file.level,
            }
        payload_steps.append({"step": step.step, "added": added, "map": map_info})
    script = MAESTRO_SCRIPT
    script = script.replace("__SIGMA_REPR__", repr(sigma))
    script = script.replace("__SIGMA__", format(sigma, "g"))
    script = script.replace(
        "__REFERENCES__", repr([_data_name(item.path) for item in references])
    )
    script = script.replace("__STEPS__", repr(payload_steps))
    script = script.replace("__FINALS__", repr([_entry_payload(item) for item in finals]))
    path.write_text(script, encoding="utf-8")


def prepare(run_dir: Path, residues: list[Residue], sigma: float) -> Path:
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
            item.added = filter_pdb(source, dest, residues)
            placed = (
                item.added.highlight.describe()
                if item.added.highlight
                else "no atoms"
            )
            print(
                f"step {step:02d} added:   {dest.name} "
                f"({item.added.natoms} atoms, green {placed})"
            )
            if item.added.newest_source is not None and not item.added.newest_source_kept:
                print(
                    "warning: the atom placed in "
                    f"{source.name} is {item.added.newest_source.describe()}, "
                    "which is outside --resid. The bright-green atom is the "
                    f"newest kept atom ({placed}).",
                    file=sys.stderr,
                )
        if "map" in slot:
            source = slot["map"]
            dest = summary / summary_filename(source.name)
            print(f"step {step:02d} map:     {dest.name}")
            item.map_file = copy_map(source, dest, sigma)
            print(
                f"         sd={item.map_file.stdev:.6g}  "
                f"{sigma:g} sigma -> {item.map_file.level:.6g}"
            )
        if item.added is None and item.map_file is None:
            continue
        steps.append(item)

    finals: list[FilteredPdb] = []
    for source in final_paths:
        dest = summary / summary_filename(source.name)
        filtered = filter_pdb(source, dest, residues)
        finals.append(filtered)
        placed = (
            filtered.highlight.describe() if filtered.highlight else "no atoms"
        )
        print(f"final: {dest.name} ({filtered.natoms} atoms, green {placed})")
        if filtered.newest_source is not None and not filtered.newest_source_kept:
            print(
                "warning: the newest atom in "
                f"{source.name} is {filtered.newest_source.describe()}, "
                "which is outside --resid. The bright-green atom is the "
                f"newest kept atom ({placed}).",
                file=sys.stderr,
            )

    write_vmd(summary / "loadVMD", references, steps, finals, sigma)
    write_chimerax(summary / "loadChimeraX", references, steps, finals, sigma)
    write_maestro(summary / "loadMaestro", references, steps, finals, sigma)
    print(f"wrote {summary / 'loadVMD'}")
    print(f"wrote {summary / 'loadChimeraX'}")
    print(f"wrote {summary / 'loadMaestro'}")
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
    parser.add_argument(
        "--map-sigma",
        type=float,
        default=DEFAULT_MAP_SIGMA,
        help=(
            "Map mesh contour in units of the map standard deviation "
            f"(default: {DEFAULT_MAP_SIGMA:g}, the contour used for the "
            "reconstruction figures)"
        ),
    )
    args = parser.parse_args(argv)
    if args.map_sigma == 0:
        raise SystemExit("--map-sigma must be non-zero.")
    residues = parse_residues(args.resid)
    prepare(args.run_dir, residues, args.map_sigma)


if __name__ == "__main__":
    main()
