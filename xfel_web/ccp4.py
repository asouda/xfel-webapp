"""CCP4, Coot, and PyMOL helpers shared by the three map workflows."""

from __future__ import annotations

import os
import re
import shlex
import shutil
import struct
import subprocess
import threading
import time
from pathlib import Path

CCP4_PROGRAMS = (
    "cif2mtz",
    "refmac5",
    "cad",
    "scaleit",
    "fft",
    "extends",
    "mapmask",
    "mtzdump",
    "ctruncate",
    "libcheck",
)

OEC_PAIRS_WITH_O6 = [
    ("MN1", "O1"), ("MN1", "O3"), ("MN1", "O6"),
    ("MN2", "O1"), ("MN2", "O2"), ("MN2", "O3"),
    ("MN3", "O2"), ("MN3", "O3"), ("MN3", "O4"), ("MN3", "O5"),
    ("MN4", "O5"), ("MN4", "O4"),
    ("CA1", "O1"), ("CA1", "O2"), ("CA1", "O5"), ("CA1", "O6"),
]

OEC_PAIRS = [
    ("MN1", "O1"), ("MN1", "O3"),
    ("MN2", "O1"), ("MN2", "O2"), ("MN2", "O3"),
    ("MN3", "O2"), ("MN3", "O3"), ("MN3", "O4"), ("MN3", "O5"),
    ("MN4", "O5"), ("MN4", "O4"),
    ("CA1", "O1"), ("CA1", "O2"), ("CA1", "O5"),
]

POS_COLORS = ["forest", "cyan", "blue", "density"]
NEG_COLORS = ["red", "salmon", "magenta", "hotpink"]
MAP_COLORS = ["marine", "blue", "deepblue", "density"]

LIPID_PIGMENT = {
    "LMT", "LHG", "SQD", "DGD", "BCR", "CLA", "CHL", "PHO",
    "STE", "LMG", "PL9", "FME", "BCT",
}


def find_ccp4(root: Path) -> Path:
    root = Path(root).expanduser()
    setup = root / "bin" / "ccp4.setup-sh"
    if not setup.is_file():
        raise FileNotFoundError(
            f"CCP4 setup script not found at {setup}. "
            "Install CCP4 or set the CCP4 directory to the suite prefix "
            "(for example /Applications/ccp4-9)."
        )
    return root


def find_coot(ccp4_root: Path) -> Path:
    for candidate in (ccp4_root / "bin" / "coot", Path(shutil.which("coot") or "")):
        if candidate.is_file():
            return candidate
    raise FileNotFoundError("Coot not found. CCP4 installs it at $CCP4/bin/coot.")


def find_pymol() -> Path:
    exe = shutil.which("pymol")
    if not exe:
        raise FileNotFoundError(
            "PyMOL not found on PATH. Install the conda-forge pymol package "
            "or add the PyMOL executable to PATH."
        )
    return Path(exe)


def probe_software(ccp4_root: str | Path) -> dict:
    """Report which external programs the workflows can see."""
    root = Path(ccp4_root).expanduser()
    setup = root / "bin" / "ccp4.setup-sh"
    programs: list[dict] = []
    ccp4_ok = setup.is_file()
    if ccp4_ok:
        script = "source \"$1\"\n" + "\n".join(
            f'if command -v {name} >/dev/null 2>&1; then '
            f'echo "OK {name} $(command -v {name})"; '
            f'else echo "MISSING {name}"; fi'
            for name in CCP4_PROGRAMS
        )
        proc = subprocess.run(
            ["/bin/bash", "-c", script, "probe", str(setup)],
            text=True,
            capture_output=True,
        )
        found = {}
        for line in (proc.stdout or "").splitlines():
            parts = line.split(maxsplit=2)
            if len(parts) >= 2 and parts[0] in {"OK", "MISSING"}:
                found[parts[1]] = {
                    "name": parts[1],
                    "ok": parts[0] == "OK",
                    "path": parts[2] if parts[0] == "OK" and len(parts) > 2 else "",
                }
        for name in CCP4_PROGRAMS:
            programs.append(found.get(name, {"name": name, "ok": False, "path": ""}))
    else:
        programs = [{"name": name, "ok": False, "path": ""} for name in CCP4_PROGRAMS]

    coot_path = root / "bin" / "coot"
    coot_ok = coot_path.is_file() or bool(shutil.which("coot"))
    pymol = shutil.which("pymol") or ""
    return {
        "ccp4_root": str(root),
        "ccp4_ok": ccp4_ok,
        "ccp4_setup": str(setup),
        "coot_ok": coot_ok,
        "coot_path": str(coot_path if coot_path.is_file() else shutil.which("coot") or ""),
        "pymol_ok": bool(pymol),
        "pymol_path": pymol,
        "programs": programs,
    }


def choose_resolution(
    source: str,
    explicit_low: float,
    explicit_high: float,
    ranges: list[tuple[float, float]],
) -> tuple[float, float]:
    """Pick the shell used by SCALEIT / FFT.

    ``ranges`` entries are ``(low-resolution Å, high-resolution Å)``.
    ``"sf"`` uses the intersection of those shells: the tighter low-resolution
    limit and the poorer high-resolution limit. ``"explicit"`` uses the
    numbers from the form.
    """
    if source == "sf":
        if not ranges:
            raise ValueError("No structure-factor resolution ranges were read.")
        low = min(item[0] for item in ranges)
        high = max(item[1] for item in ranges)
    elif source == "explicit":
        low = float(explicit_low)
        high = float(explicit_high)
    else:
        raise ValueError(f"Resolution source must be 'sf' or 'explicit', got {source!r}")
    if not (low > high > 0):
        raise ValueError(f"Expected low-resolution > high-resolution > 0, got {low} and {high}")
    return low, high


def res_tag(low: float, high: float) -> str:
    return f"{low:.2f}-{high:.2f}".replace(".", "p")


def fft_grid_resolution(low: float, data_high: float, requested: float | None) -> float:
    """High-resolution limit that sets the fft sampling grid.

    CCP4 fft takes the maximum Miller indices from this limit and the cell,
    then samples the map at about one third of that spacing. A limit finer
    than the reflections in the MTZ extends that reciprocal-space box. Indices
    with no observation are left out of the sum, which is a zero amplitude.
    """
    if requested is None:
        return data_high
    grid = float(requested)
    if not (low > grid > 0):
        raise ValueError(
            f"FFT grid resolution must be finer than the low-resolution cutoff "
            f"of {low:.3f} Å and greater than zero. Got {grid}."
        )
    if grid > data_high + 1e-6:
        raise ValueError(
            f"FFT grid resolution {grid:.3f} Å is coarser than the "
            f"{data_high:.3f} Å limit of the coefficients. "
            f"A finer grid needs a smaller number of ångströms."
        )
    if grid >= data_high - 1e-6:
        return data_high
    return grid


def grid_res_tag(low: float, high: float, grid_high: float) -> str:
    """Filename tag for a map, with the fft grid limit when it is finer."""
    tag = res_tag(low, high)
    if grid_high < high - 1e-4:
        return tag + "_grid" + f"{grid_high:.2f}".replace(".", "p")
    return tag


def delta_b_tag(value: float) -> str:
    if float(value).is_integer():
        return f"dB{int(value):+d}"
    return f"dB{value:+.1f}".replace(".", "p")


def require_files(files: dict[str, Path]) -> None:
    missing = [f"{key}: {path}" for key, path in files.items() if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise FileNotFoundError("Missing or empty input files:\n  " + "\n  ".join(missing))


def scan_hetatm(pdb_path: Path, skip_resnames: set[str] | None = None) -> list[dict]:
    skip_resnames = skip_resnames or {"HOH", "WAT", "DOD", "H2O"}
    residues: dict[tuple, dict] = {}
    with pdb_path.open() as handle:
        for line in handle:
            if not line.startswith("HETATM"):
                continue
            resname = line[17:20].strip()
            if resname in skip_resnames:
                continue
            chain = line[21]
            try:
                resi = int(line[22:26])
            except ValueError:
                continue
            atom = line[12:16].strip()
            occ = float(line[54:60] or 0)
            key = (chain, resi, resname)
            rec = residues.setdefault(
                key,
                {"chain": chain, "resi": resi, "resname": resname, "atoms": set(), "occ": [occ, occ]},
            )
            rec["atoms"].add(atom)
            rec["occ"][0] = min(rec["occ"][0], occ)
            rec["occ"][1] = max(rec["occ"][1], occ)
    return sorted(residues.values(), key=lambda rec: (rec["chain"], rec["resi"], rec["resname"]))


def residue_centroid(pdb_path: Path, chain: str, resi: int) -> tuple[float, float, float]:
    coords = []
    with pdb_path.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[21] != chain:
                continue
            try:
                if int(line[22:26]) != int(resi):
                    continue
            except ValueError:
                continue
            coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if not coords:
        raise ValueError(f"No atoms for chain {chain} residue {resi} in {pdb_path.name}")
    count = float(len(coords))
    return (
        sum(item[0] for item in coords) / count,
        sum(item[1] for item in coords) / count,
        sum(item[2] for item in coords) / count,
    )


def residue_atom_names(pdb_path: Path, chain: str, resi: int) -> set[str]:
    names: set[str] = set()
    with pdb_path.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[21] != chain:
                continue
            try:
                if int(line[22:26]) != int(resi):
                    continue
            except ValueError:
                continue
            names.add(line[12:16].strip())
    return names


_PML_LOAD = re.compile(r'^load\s+"([^"]+)"\s*,\s*(\S+)')
_PML_SELECT = re.compile(
    r"^select\s+([A-Za-z0-9_]+)\s*,\s*xyz and chain\s+(\S+)\s+and resi\s+(\d+)"
)
_PML_BOND = re.compile(r"^bond\s+(.+?)\s*,\s*(.+)$")
_PML_DISTANCE = re.compile(r"^distance\s+(\S+)\s*,\s*(.+?)\s*,\s*(.+)$")
_PML_ATOM = re.compile(r"\bname\s+(\S+)")
_PML_ALT = re.compile(r"\balt\s+([A-Za-z0-9])")
_PML_CARVE = re.compile(r"\bcarve=([+-]?\d+(?:\.\d+)?)")


def _pml_atom(clause: str, selections: dict[str, dict]) -> dict | None:
    """Turn one side of a PyMOL bond or distance into chain, residue, atom, and altloc."""
    words = clause.split()
    if not words or words[0] not in selections:
        return None
    found = _PML_ATOM.search(clause)
    if not found:
        return None
    spec = {
        "chain": selections[words[0]]["chain"],
        "resi": selections[words[0]]["resi"],
        "atom": found.group(1),
    }
    alt = _PML_ALT.search(clause)
    if alt:
        spec["alt"] = alt.group(1)
    return spec


def pymol_scene_links(pml_path: Path) -> dict:
    """Read the map file, explicit bonds, and distance objects written in a PyMOL script."""
    selections: dict[str, dict] = {}
    bonds: list[dict] = []
    distances: list[dict] = []
    map_path = None
    carve = None
    for raw in Path(pml_path).read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        loaded = _PML_LOAD.match(line)
        if loaded:
            path, obj = loaded.group(1), loaded.group(2)
            if obj != "xyz" and Path(path).suffix.lower() in {".map", ".ccp4", ".mrc"}:
                map_path = path
            continue
        selected = _PML_SELECT.match(line)
        if selected:
            selections[selected.group(1)] = {
                "chain": selected.group(2),
                "resi": int(selected.group(3)),
            }
            continue
        bonded = _PML_BOND.match(line)
        if bonded:
            left = _pml_atom(bonded.group(1), selections)
            right = _pml_atom(bonded.group(2), selections)
            if left and right:
                bonds.append({"a": left, "b": right})
            continue
        measured = _PML_DISTANCE.match(line)
        if measured:
            left = _pml_atom(measured.group(2), selections)
            right = _pml_atom(measured.group(3), selections)
            if left and right:
                distances.append({"name": measured.group(1), "a": left, "b": right})
        if carve is None:
            carved = _PML_CARVE.search(line)
            if carved:
                carve = float(carved.group(1))
    return {"map": map_path, "bonds": bonds, "distances": distances, "carve": carve}


def residue_name(pdb_path: Path, chain: str, resi: int) -> str:
    with pdb_path.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[21] != chain:
                continue
            try:
                if int(line[22:26]) != int(resi):
                    continue
            except ValueError:
                continue
            return line[17:20].strip()
    raise ValueError(f"No atoms for chain {chain} residue {resi} in {pdb_path.name}")


def _iter_pdb_atoms(pdb_path: Path):
    with pdb_path.open() as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                yield line if line.endswith("\n") else line + "\n"


def extract_monomer_pdb(
    pdb_path: Path,
    resname: str,
    chain: str | None = None,
    resi: int | None = None,
) -> str:
    """One conformer of one residue, with the altloc column cleared."""
    groups: dict[tuple[str, int], list[str]] = {}
    for line in _iter_pdb_atoms(pdb_path):
        if line[17:20].strip() != resname:
            continue
        try:
            key = (line[21], int(line[22:26]))
        except ValueError:
            continue
        groups.setdefault(key, []).append(line)
    if not groups:
        raise ValueError(f"Residue type {resname} was not found in {pdb_path.name}")
    if chain is not None:
        key = (chain, int(resi))
        if key not in groups:
            raise ValueError(
                f"{resname} chain {chain} residue {resi} was not found in {pdb_path.name}"
            )
    else:
        key = sorted(groups)[0]
    by_alt: dict[str, list[str]] = {}
    for line in groups[key]:
        by_alt.setdefault(line[16], []).append(line)
    chosen = by_alt[" "] if " " in by_alt else max(by_alt.values(), key=len)
    cleaned = []
    for line in chosen:
        if line[16] != " ":
            line = line[:16] + " " + line[17:]
        cleaned.append(line)
    return "".join(cleaned) + "END\n"


def coords_from_pdb_text(pdb_text: str) -> dict[str, tuple[float, float, float]]:
    coords: dict[str, tuple[float, float, float]] = {}
    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[12:16].strip()
        coords[name] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return coords


def residue_conformers(
    pdb_path: Path, chain: str, resi: int
) -> dict[str, dict[str, tuple[float, float, float]]]:
    """Altloc -> atom name -> xyz. Blank atoms are included in every conformer."""
    blank: dict[str, tuple[float, float, float]] = {}
    alternates: dict[str, dict[str, tuple[float, float, float]]] = {}
    found = False
    with pdb_path.open() as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[21] != chain:
                continue
            try:
                if int(line[22:26]) != int(resi):
                    continue
            except ValueError:
                continue
            found = True
            name = line[12:16].strip()
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            alt = line[16].strip()
            if not alt:
                blank[name] = xyz
            else:
                alternates.setdefault(alt, {})[name] = xyz
    if not found:
        raise ValueError(f"No atoms for chain {chain} residue {resi} in {pdb_path.name}")
    if not alternates:
        return {"": blank}
    return {alt: {**blank, **atoms} for alt, atoms in sorted(alternates.items())}


def copy_ccp4_header_statistics(source: Path, dest: Path) -> None:
    """Copy min, max, mean, and rms from one CCP4 map header onto another.

    ``mapmask`` recomputes these on the box it writes. Sigma contours follow
    the asymmetric-unit values already stored on the ``fft`` map, which
    ``extends`` also keeps.
    """
    src = source.read_bytes()[:1024]
    if len(src) < 1024:
        raise ValueError(f"{source} has no CCP4 map header")
    with dest.open("r+b") as handle:
        header = bytearray(handle.read(1024))
        if len(header) < 1024:
            raise ValueError(f"{dest} has no CCP4 map header")
        for word in (19, 20, 21, 54):
            start = word * 4
            header[start : start + 4] = src[start : start + 4]
        handle.seek(0)
        handle.write(header)


def ccp4_map_stats(path: Path) -> dict[str, float]:
    raw = path.read_bytes()
    nx, ny, nz, mode = struct.unpack_from("<4i", raw, 0)
    nsymbt = struct.unpack_from("<i", raw, 23 * 4)[0]
    if nsymbt < 0 or 1024 + nsymbt + nx * ny * nz * 4 > len(raw):
        nsymbt = 0
    if mode != 2:
        raise ValueError(f"Unsupported CCP4 map mode {mode} in {path}")
    offset = 1024 + nsymbt
    count = nx * ny * nz
    try:
        import numpy as np

        arr = np.frombuffer(raw, dtype="<f4", count=count, offset=offset)
        rms = float(np.sqrt(np.mean(np.square(arr))))
        return {
            "rms": rms,
            "mean": float(arr.mean()),
            "min": float(arr.min()),
            "max": float(arr.max()),
            "nx": float(nx),
            "ny": float(ny),
            "nz": float(nz),
        }
    except Exception:
        arr = struct.unpack_from("<" + str(count) + "f", raw, offset)
        mean = sum(arr) / float(count)
        rms = (sum(value * value for value in arr) / float(count)) ** 0.5
        return {
            "rms": float(rms),
            "mean": float(mean),
            "min": float(min(arr)),
            "max": float(max(arr)),
            "nx": float(nx),
            "ny": float(ny),
            "nz": float(nz),
        }


def parse_scaleit_rfactor(log_path: Path) -> str | None:
    if not log_path.is_file():
        return None
    text = log_path.read_text(errors="replace")
    match = re.search(r"THE TOTALS\s+[0-9.]+\s+[0-9.]+\s+[0-9.]+\s+[0-9.]+\s+([0-9.]+)", text)
    if match:
        return match.group(1)
    return None


class Ccp4Session:
    """One calculation. Shell snippets run after ``ccp4.setup-sh``."""

    def __init__(self, output_dir: Path, ccp4_root: Path, emit):
        self.output_dir = Path(output_dir)
        self.ccp4_root = find_ccp4(ccp4_root)
        self.coot = find_coot(self.ccp4_root)
        self.pymol = find_pymol()
        self.emit = emit
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def announce(self, message: str) -> None:
        self.emit("")
        self.emit(message)

    def run_bash(self, script: str, log_path: Path | None = None) -> str:
        setup = self.ccp4_root / "bin" / "ccp4.setup-sh"
        wrapped = (
            "set -eo pipefail\n"
            'export SSL_CERT_FILE="${SSL_CERT_FILE:-}"\n'
            f"source {shlex.quote(str(setup))}\n"
            "set -eo pipefail\n"
            f"cd {shlex.quote(str(self.output_dir))}\n"
            f"{script}\n"
        )
        lines = self._stream(["/bin/bash", "-c", wrapped], cwd=self.output_dir, env=None)
        text = "\n".join(lines)
        if log_path is not None:
            log_path.write_text(text + ("\n" if text else ""))
        return text

    def run_command(self, args: list[str], log_path: Path | None = None, env: dict | None = None) -> str:
        lines = self._stream(args, cwd=self.output_dir, env=env)
        text = "\n".join(lines)
        if log_path is not None:
            log_path.write_text(text + ("\n" if text else ""))
        return text

    def _stream(self, args: list[str], cwd: Path, env: dict | None) -> list[str]:
        proc = subprocess.Popen(
            args,
            cwd=str(cwd),
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        lines: list[str] = []

        def reader() -> None:
            assert proc.stdout is not None
            for line in proc.stdout:
                cleaned = line.rstrip("\n")
                lines.append(cleaned)
                self.emit(cleaned)

        thread = threading.Thread(target=reader, daemon=True)
        thread.start()
        started = time.time()
        next_heartbeat = started + 20
        while proc.poll() is None:
            thread.join(timeout=2)
            now = time.time()
            if proc.poll() is None and now >= next_heartbeat:
                self.emit(f"… still running ({int(now - started)} s)")
                next_heartbeat = now + 20
        thread.join()
        if proc.returncode != 0:
            tail = "\n".join(lines[-40:]) if lines else "(no output)"
            raise RuntimeError(f"Command failed (exit {proc.returncode}):\n{tail}")
        return lines

    def extend_around_model(
        self,
        map_in: Path,
        map_out: Path,
        xyz_in: Path,
        border: float,
        log_path: Path,
    ) -> None:
        """Pad ``map_in`` so it covers ``xyz_in`` plus ``border`` angstroms.

        ``extends`` stores one map section in a fixed array of 200,000 points.
        A fine grid can make the model box larger than that. ``mapmask`` has
        no such limit and applies the same border; the asymmetric-unit sigma
        is then copied back onto the padded map.
        """
        keywords = f"BORDER {border:.1f}\nEND"
        try:
            self.run_ccp4(
                "extends",
                [("mapin", map_in), ("mapout", map_out), ("xyzin", xyz_in)],
                keywords,
                log_path,
            )
            return
        except RuntimeError as exc:
            if "arrays too small" not in str(exc):
                raise
        self.emit(
            "extends cannot store a section of this map. "
            "mapmask will apply the same border and keep the asymmetric-unit sigma."
        )
        map_out.unlink(missing_ok=True)
        self.run_ccp4(
            "mapmask",
            [("mapin", map_in), ("mapout", map_out), ("xyzin", xyz_in)],
            f"BORDER {border:.1f}\nEXTEND XTAL\nEND",
            log_path.with_name(log_path.stem + "_mapmask.log"),
        )
        copy_ccp4_header_statistics(map_in, map_out)

    def run_ccp4(self, program: str, kv: list[tuple[str, Path | str]], keywords: str, log_path: Path) -> str:
        parts = [program]
        for key, value in kv:
            parts.append(str(key))
            parts.append(shlex.quote(str(value)))
        body = " ".join(parts) + " << 'EOF'\n" + keywords.strip() + "\nEOF"
        self.emit(f"→ {program}  ({log_path.name})")
        return self.run_bash(body, log_path)

    def mtz_labels(self, mtz_path: Path) -> list[str]:
        text = self.run_ccp4("mtzdump", [("hklin", mtz_path)], "END", self.output_dir / f"{mtz_path.stem}_mtzdump.log")
        labels: list[str] = []
        capture = False
        for line in text.splitlines():
            if "Column Labels" in line:
                capture = True
                continue
            if capture:
                stripped = line.strip()
                if stripped.startswith("*") and labels:
                    break
                if stripped.startswith("*"):
                    continue
                labels.extend(stripped.split())
        return labels

    def mtz_resolution_angstrom(self, mtz_path: Path) -> tuple[float, float]:
        log_path = self.output_dir / f"{mtz_path.stem}_mtzdump.log"
        text = log_path.read_text(errors="replace")
        match = re.search(
            r"Resolution Range\s*:.*?([0-9.]+)\s*-\s*([0-9.]+)\s*A",
            text,
            flags=re.S,
        )
        if not match:
            raise RuntimeError(f"Could not read resolution limits from {log_path}")
        first, second = float(match.group(1)), float(match.group(2))
        return max(first, second), min(first, second)

    @staticmethod
    def pick_label(labels: list[str], candidates: list[str], what: str) -> str:
        upper = {label.upper(): label for label in labels}
        for name in candidates:
            if name.upper() in upper:
                return upper[name.upper()]
        raise KeyError(f"No {what} column among {candidates}. MTZ labels were: {labels}")

    def export_fragment(
        self,
        extended_map: Path,
        pdb_path: Path,
        site_copies: list[tuple[str, int]],
        site_label: str,
        radius: float,
        sigma_level: float,
        difference: bool,
        recompute: bool,
    ) -> Path:
        chain, resi = site_copies[0]
        fragment = self.output_dir / (
            f"{extended_map.stem}_{site_label}_{chain}{resi}_{int(radius)}A.map"
        )
        if fragment.is_file() and not recompute:
            self.emit(f"Reusing existing fragment: {fragment}")
            return fragment

        copies_py = ",\n    ".join(f'("{item_chain}", {item_resi})' for item_chain, item_resi in site_copies)
        diff_flag = 1 if difference else 0
        script = self.output_dir / "coot_site_export.py"
        script.write_text(
            f'''#!/usr/bin/env python
import os
import sys

MAP_FILE = "{extended_map}"
PDB_FILE = "{pdb_path}"
OUT_DIR = "{self.output_dir}"
RADIUS_A = {radius:.1f}
SITE_COPIES = [
    {copies_py},
]

imol_map = read_ccp4_map(MAP_FILE, {diff_flag})
if imol_map is None or imol_map < 0:
    sys.exit("Failed to read map: %s" % MAP_FILE)
set_map_is_difference_map(imol_map, {diff_flag})
try:
    set_contour_level_in_sigma(imol_map, {float(sigma_level):.1f})
    rmsd = map_sigma(imol_map)
    print("Map RMSD (1 sigma) = %.4f" % rmsd)
except Exception:
    rmsd = None

imol_model = read_pdb(PDB_FILE)
if imol_model is None or imol_model < 0:
    sys.exit("Failed to read coordinates: %s" % PDB_FILE)

base_name = os.path.splitext(os.path.basename(MAP_FILE))[0]
for chain_id, resno in SITE_COPIES:
    atoms = residue_info_py(imol_model, chain_id, resno, "")
    if not atoms:
        print("WARNING: residue %s/%d not found - skipping" % (chain_id, resno))
        continue
    coords = [a[2] for a in atoms]
    n = len(coords)
    x = sum(c[0] for c in coords) / float(n)
    y = sum(c[1] for c in coords) / float(n)
    z = sum(c[2] for c in coords) / float(n)
    set_rotation_centre(x, y, z)
    out_name = "%s_{site_label}_%s%d_%dA.map" % (base_name, chain_id, resno, int(RADIUS_A))
    out_path = os.path.join(OUT_DIR, out_name)
    status = export_map_fragment(imol_map, x, y, z, RADIUS_A, out_path)
    print("Exported %s (centre %.2f %.2f %.2f, radius %.1f A) -> status %s" % (
        out_path, x, y, z, RADIUS_A, status))

print("Done.")
try:
    coot_real_exit(0)
except NameError:
    pass
'''
        )
        self.emit("→ coot  (coot_site_export.log)")
        env = os.environ.copy()
        env["COOT_VERBOSITY_LEVEL"] = "0"
        log_path = self.output_dir / "coot_site_export.log"
        try:
            self.run_command(
                [str(self.coot), "--no-graphics", "--script", str(script)],
                log_path=log_path,
                env=env,
            )
        except RuntimeError as exc:
            raise RuntimeError(f"Coot failed. See {log_path}\n{exc}") from exc
        if not fragment.is_file():
            raise FileNotFoundError(f"Expected Coot fragment was not written: {fragment}")
        self.emit(f"Map fragment of radius {radius:.1f} Å centered on {site_label}: {fragment}")
        return fragment

    def render_figure(
        self,
        *,
        fragment: Path,
        model_pdb: Path,
        map_name: str,
        png_path: Path,
        pml_path: Path,
        site_chain: str,
        site_resi: int,
        extra_residues: list[dict],
        positive_levels: list[float],
        negative_levels: list[float] | None,
        positive_colors: list[str],
        negative_colors: list[str],
        contour_mode: str,
        surface: str,
        carve: float,
        show_oec_bonds: bool,
        draw_distances: bool,
        oec_pairs: list[tuple[str, str]],
        png_width: int,
        png_height: int,
        ray_trace: bool,
    ) -> None:
        stats = ccp4_map_stats(fragment)
        self.emit(
            f"Fragment file stats (raw voxels): RMS={stats['rms']:.5f}  "
            f"mean={stats['mean']:.5f}  min={stats['min']:.5f}  max={stats['max']:.5f}  "
            f"grid={int(stats['nx'])}x{int(stats['ny'])}x{int(stats['nz'])}"
        )
        if contour_mode == "sigma":
            normalize_flag = 1
            self.emit("Contour mode: PyMOL-normalized sigma (level n = n sigma)")
        elif contour_mode == "absolute":
            normalize_flag = 0
            self.emit("Contour mode: absolute map units (normalize_ccp4_maps off)")
        else:
            raise ValueError(f"Contour mode must be 'sigma' or 'absolute', got {contour_mode!r}")

        if surface == "mesh":
            contour_cmd = "isomesh"
            surface_transparency = 0.25
        elif surface == "solid":
            contour_cmd = "isosurface"
            surface_transparency = 0.45
        else:
            raise ValueError(f"Surface representation must be 'mesh' or 'solid', got {surface!r}")
        self.emit(f"Surface representation: {surface} ({contour_cmd})")

        site_sel = f"xyz and chain {site_chain} and resi {site_resi}"
        lines = [
            f"set normalize_ccp4_maps, {normalize_flag}",
            f'load "{fragment}", {map_name}',
            f'load "{model_pdb}", xyz',
            "hide everything",
            "util.cbay xyz",
            f"select site, {site_sel}",
            "show spheres, site",
            "show sticks, site",
            "set sphere_scale, 0.15, xyz",
            "set stick_radius, 0.12, xyz",
        ]
        for extra in extra_residues:
            label = extra.get("label") or f"{extra['chain']}{extra['resi']}"
            lines += [
                f"select {label}, xyz and chain {extra['chain']} and resi {extra['resi']}",
                f"show sticks, {label}",
                f"show spheres, {label}",
                f"set sphere_scale, 0.15, {label}",
                f"set stick_radius, 0.12, {label}",
            ]
        present = residue_atom_names(model_pdb, site_chain, site_resi)
        if show_oec_bonds:
            for atom_a, atom_b in oec_pairs:
                if atom_a in present and atom_b in present:
                    lines.append(f"bond site and name {atom_a}, site and name {atom_b}")
        if draw_distances:
            if "MN1" in present and "O5" in present:
                lines.append("distance d_mn1_o5, site and name MN1, site and name O5")
            for extra in extra_residues:
                extra_atoms = residue_atom_names(model_pdb, extra["chain"], extra["resi"])
                if "MN1" in present and "OE2" in extra_atoms:
                    label = extra.get("label") or f"{extra['chain']}{extra['resi']}"
                    lines.append(
                        f"distance d_mn1_{label}_A, site and name MN1, {label} and alt A and name OE2"
                    )
                    lines.append(
                        f"distance d_mn1_{label}_B, site and name MN1, {label} and alt B and name OE2"
                    )
        lines += [
            "viewport %d,%d" % (png_width, png_height),
            "bg_color white",
            "space cmyk",
            "set ray_opaque_background, 1",
            "set ray_trace_mode, 0",
            "set ray_trace_fog, 1",
            "set orthoscopic, on",
            "set antialias, 4",
            "set mesh_width, 0.25",
            "set dash_gap, 0.15",
            "set dash_radius, 0.05",
            "hide labels",
            "set auto_zoom, 0",
        ]
        carve_target = "site"
        for extra in extra_residues:
            carve_target += " or " + (extra.get("label") or f"{extra['chain']}{extra['resi']}")
        lines.append(f"select carve_sel, {carve_target}")
        lines += [
                  "orient carve_sel",
                  "zoom carve_sel, 2",
        ]

        if negative_levels is None:
            pos_i = 0
            neg_i = 0
            for level in positive_levels:
                token = f"{level:g}".replace("-", "m").replace(".", "p")
                name = f"map_{token}"
                if level >= 0:
                    color = positive_colors[min(pos_i, len(positive_colors) - 1)]
                    pos_i += 1
                else:
                    color = negative_colors[min(neg_i, len(negative_colors) - 1)]
                    neg_i += 1
                lines += [
                    f"{contour_cmd} {name}, {map_name}, {level:.6f}, carve_sel, carve={carve:.2f}",
                    f"color {color}, {name}",
                ]
                self.emit(f"  {level:g} sigma -> {name} at {level:.5f} ({color})")
        else:
            for index, level in enumerate(positive_levels):
                name = "mappos_" + str(level).replace(".", "p")
                color = positive_colors[min(index, len(positive_colors) - 1)]
                lines += [
                    f"{contour_cmd} {name}, {map_name}, {level:.6f}, carve_sel, carve={carve:.2f}",
                    f"color {color}, {name}",
                ]
                self.emit(f"  +{level:g} sigma -> {name} at {level:.5f} ({color})")
            for index, level in enumerate(negative_levels):
                name = "mapneg_" + str(level).replace(".", "p")
                color = negative_colors[min(index, len(negative_colors) - 1)]
                lines += [
                    f"{contour_cmd} {name}, {map_name}, {-level:.6f}, carve_sel, carve={carve:.2f}",
                    f"color {color}, {name}",
                ]
                self.emit(f"  -{level:g} sigma -> {name} at {-level:.5f} ({color})")

        if surface == "mesh":
            lines.append("set mesh_width, 0.75, map*")
        lines.append(f"set transparency, {surface_transparency:.2f}, map*")

        header = (
            "# Isodensity scene written by the XFEL map web application.\n"
            "# Open in PyMOL with: pymol "
            + pml_path.name
            + "\n"
        )
        pml_path.write_text(header + "\n".join(lines) + "\n")
        self.emit(f"Wrote {pml_path}")

        ray_flag = 1 if ray_trace else 0
        render_path = pml_path.with_name(pml_path.stem + "_render.pml")
        png_cmd = f'png {png_path}, {png_width}, {png_height}, dpi=200, ray={ray_flag}'
        render_path.write_text(f"@{pml_path.name}\n{png_cmd}\nquit\n")
        self.emit("→ pymol  (pymol_render.log)")
        log_path = self.output_dir / "pymol_render.log"

        def _render(script: Path) -> bool:
            self.run_command([str(self.pymol), "-cq", str(script)], log_path=log_path)
            return png_path.is_file()

        rendered = False
        try:
            rendered = _render(render_path)
        except RuntimeError:
            rendered = False
        if not rendered and ray_trace:
            self.emit("Ray-traced PNG failed; retrying with ray=0 …")
            render_path.write_text(
                f"@{pml_path.name}\n"
                f'png {png_path}, {png_width}, {png_height}, dpi=200, ray=0\n'
                "quit\n"
            )
            rendered = _render(render_path)
        if not rendered:
            raise RuntimeError(f"PyMOL did not write {png_path}. See {log_path}")
        self.emit(f"Saved image: {png_path}")
