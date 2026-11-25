"""Fo−Fo, omit-model mFo−DFc, sharpened 2mFo−DFc, and omit-map reconstruction."""

from __future__ import annotations

import math
import shutil
from pathlib import Path

from xfel_web.ccp4 import (
    MAP_COLORS,
    NEG_COLORS,
    OEC_PAIRS,
    OEC_PAIRS_WITH_O6,
    POS_COLORS,
    CommandFailed,
    Ccp4Session,
    choose_resolution,
    delta_b_tag,
    fft_grid_resolution,
    parse_scaleit_rfactor,
    require_files,
    reasonable_bond_specs,
    residue_atom_names,
    residue_centroid,
)
from xfel_web.direct_map import direct_site_cube, parse_fft_sampling, whole_cell_exceeds_fft
from xfel_web.restraints import (
    external_angle_keywords,
    external_dihedral_keywords,
    external_distance_keywords,
    external_harmonic_keywords,
    protein_harmonic_lines,
)

# Protocol sent on every refmac5 job. A trailing "-" continues the keyword.
# A custom dictionary is a libin file assignment on the command line. NCYC is last.
REFMAC_FIXED_KEYWORDS = """\
MAKE -
    HYDROGEN NO -
    HOUT NO -
    PEPTIDE NO -
    CISPEPTIDE YES -
    SSBRIDGE YES -
    SYMMETRY YES -
    SUGAR YES -
    CONNECTIVITY NO -
    LINK NO
REFI -
    TYPE REST -
    RESI MLKF -
    METH CGMAT -
    BREF ANIS
SCAL -
    TYPE SIMP -
    LSSC -
    ANISO -
    EXPE
SOLVENT YES
WEIGHT -
    AUTO"""


def exclude_restraint_keywords(specs) -> str:
    """``RESTRAINT EXCLUDE RESIDUE`` lines for residues left out of restraints.

    Each spec is ``(chain, residue, atoms)``. An atom list of ``*`` means every
    atom in that residue.
    """
    lines = []
    for chain, resi, atoms in specs or []:
        atom_text = " ".join(str(atom) for atom in atoms)
        lines.append(f"RESTRAINT EXCLUDE RESIDUE {int(resi)} {chain} ATOMS {atom_text}")
    return "\n".join(lines)


def refmac_keywords(
    ncyc: int,
    leading: str = "",
    bref: str = "ANIS",
    refi_type: str = "REST",
    exclude=None,
) -> str:
    """REFMAC standard input.

    External restraints, when present, come first, so they replace matching
    distances and angles from the CCP4 monomer library and from a ``libin``
    dictionary assigned on the command line. The fixed protocol follows.
    ``RESTRAINT EXCLUDE`` lines, when present, come after that protocol.
    ``NCYC`` comes last.

    ``bref`` and ``refi_type`` are the ``REFI BREF`` and ``REFI TYPE`` tokens.
    The defaults keep the Fo−Fo zero-cycle job on anisotropic restrained
    refinement. Omit maps and reconstruction pass the form values, whose
    defaults are isotropic and restrained.

    ``exclude`` is a list of ``(chain, residue, atoms)`` specs. An atom list of
    ``*`` means every atom in that residue.
    """
    bref = str(bref).upper()
    refi_type = str(refi_type).upper()
    if bref not in {"ISOT", "ANIS", "MIXE", "OVER"}:
        raise ValueError(f"REFI BREF must be ISOT, ANIS, MIXE, or OVER, got {bref!r}")
    if refi_type not in {"REST", "UNRE"}:
        raise ValueError(f"REFI TYPE must be REST or UNRE, got {refi_type!r}")
    protocol = REFMAC_FIXED_KEYWORDS.replace("    TYPE REST -\n", f"    TYPE {refi_type} -\n", 1)
    protocol = protocol.replace("    BREF ANIS\n", f"    BREF {bref}\n", 1)
    parts: list[str] = []
    if leading and leading.strip():
        parts.append(leading.strip())
    parts.append(protocol.strip())
    excluded = exclude_restraint_keywords(exclude)
    if excluded:
        parts.append(excluded)
    return "\n".join(parts) + f"\nNCYC {int(ncyc)}\nEND"


def _sampling_grid(session: Ccp4Session, low: float, high: float, requested: float | None) -> tuple[float, str]:
    """High-resolution limit of the ``fft`` grid, and a filename suffix when it is finer."""
    grid_high = fft_grid_resolution(low, high, requested)
    suffix = ""
    if grid_high < high - 1e-4:
        suffix = "_grid" + f"{grid_high:.2f}".replace(".", "p")
        session.emit(
            f"FFT sampling grid is the {grid_high:.3f} Å grid. "
            f"The coefficients stop at {high:.3f} Å, so reciprocal-space points "
            f"from {high:.3f} to {grid_high:.3f} Å contribute zero amplitude."
        )
    return grid_high, suffix


def _cartesian_center(cfg: dict) -> tuple[float, float, float] | None:
    if cfg.get("map_center") != "cartesian":
        return None
    point = cfg.get("map_center_xyz")
    if not point:
        return None
    return (float(point[0]), float(point[1]), float(point[2]))


def _additional_center(site: dict, pdb_path: Path) -> tuple[float, float, float]:
    """Center of the additional-site map: an explicit point, or the residue centroid."""
    if site.get("map_center") == "cartesian" and site.get("map_center_xyz"):
        x, y, z = site["map_center_xyz"]
        return (float(x), float(y), float(z))
    try:
        return residue_centroid(pdb_path, str(site["chain"]), int(site["resi"]))
    except ValueError as exc:
        raise ValueError(f"Additional site: {exc}") from exc


def _additional_cube(session, cfg: dict, pdb_path: Path, extended_map: Path, direct_cube: dict) -> Path:
    """A fine-grid cube around the additional site when the whole cell was not written."""
    site = cfg["additional_site"]
    center = _additional_center(site, pdb_path)
    cube = extended_map.with_name(
        f"{extended_map.stem}_{site['label']}_{site['chain']}{int(site['resi'])}_cube.ccp4"
    )
    session.emit(
        "The whole-cell transform was replaced by a cube around the main site. "
        f"The additional site is written as its own cube: {cube.name}"
    )
    text = direct_site_cube(
        direct_cube["mtz"],
        cube,
        direct_cube["f"],
        direct_cube["phi"],
        direct_cube["low"],
        direct_cube["high"],
        direct_cube["sampling"],
        center,
        cfg["export_radius"],
        None,
        emit=session.emit,
    )
    cube.with_suffix(".log").write_text(text)
    return cube


def _export_additional_fragment(
    session,
    cfg: dict,
    source_map: Path,
    pdb_path: Path,
    sigma: float,
    difference: bool,
) -> Path | None:
    """Clip and save the additional-site map. The peak search keeps the main clip."""
    site = cfg.get("additional_site")
    if not site:
        return None
    session.emit(
        f"Additional site: chain {site['chain']}  residue {site['resi']}  label {site['label']}"
    )
    center = None
    if site.get("map_center") == "cartesian" and site.get("map_center_xyz"):
        center = _additional_center(site, pdb_path)
        session.emit(
            "Additional map center is the Cartesian point: "
            f"{center[0]:.2f} {center[1]:.2f} {center[2]:.2f} Å"
        )
    else:
        centroid = _additional_center(site, pdb_path)
        session.emit(
            "Additional map center is the site centroid: "
            f"{centroid[0]:.2f} {centroid[1]:.2f} {centroid[2]:.2f} Å"
        )
    session.announce("Export the additional site")
    return session.export_fragment(
        source_map,
        pdb_path,
        [(str(site["chain"]), int(site["resi"]))],
        str(site["label"]),
        cfg["export_radius"],
        sigma,
        difference=difference,
        recompute=cfg["recompute_fragment"],
        center=center,
        script_stem="coot_additional_export",
    )


def _check_site(pdb_path: Path, chain: str, resi: int, label: str):
    atoms = residue_atom_names(pdb_path, chain, resi)
    if not atoms:
        raise ValueError(
            f"Residue {chain}/{resi} was not found in {pdb_path.name}. "
            "Use “List hetero residues” and set the chain and residue."
        )
    cx, cy, cz = residue_centroid(pdb_path, chain, resi)
    print_site = f"Configured site: chain {chain}  residue {resi}  label {label}"
    return atoms, (cx, cy, cz), print_site


def _prepare_resolution(session: Ccp4Session, source: str, lowres: float, dmin: float, mtzs: list[tuple[str, Path]]):
    ranges = []
    for _label, mtz in mtzs:
        columns = session.mtz_labels(mtz)
        low, high = session.mtz_resolution_angstrom(mtz)
        session.emit(f"Structure factors {mtz.name}: {low:.2f}–{high:.2f} Å  columns={columns}")
        ranges.append((low, high))
    if len(ranges) > 1:
        session.emit(
            "Common shell: "
            f"{min(item[0] for item in ranges):.2f}–{max(item[1] for item in ranges):.2f} Å"
        )
    low, high = choose_resolution(source, lowres, dmin, ranges)
    if source == "sf":
        session.emit(f"Using MTZ resolution limits: low {low:.3f} Å, high {high:.3f} Å")
    else:
        session.emit(f"Using explicit resolution: low {low:.3f} Å, high {high:.3f} Å")
    return low, high


def _cif2mtz(session: Ccp4Session, pdb_id: str, cif_path: Path, mtz_path: Path, kind: str) -> None:
    session.run_ccp4(
        "cif2mtz",
        [("hklin", cif_path), ("hklout", mtz_path)],
        "END",
        session.output_dir / f"{pdb_id}_cif2mtz_{kind}.log",
    )


def is_reflection_mtz(path: Path) -> bool:
    """True when the observed-reflection input is already an MTZ."""
    return path.name.lower().endswith(".mtz")


def install_reflection_mtz(source: Path, dest: Path) -> None:
    """Copy a supplied reflection MTZ into the job directory."""
    if source.resolve() == dest.resolve():
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, dest)


_MEAN_INTENSITIES = ("I", "SIGI")
_FRIEDEL_INTENSITIES = ("I(+)", "SIGI(+)", "I(-)", "SIGI(-)")
_AMPLITUDE_PAIRS = (
    ("F", "SIGF"),
    ("FP", "SIGFP"),
    ("FOBS", "SIGFOBS"),
    ("FMEAN", "SIGFMEAN"),
)
_FREE_LABELS = ("FREE", "FREER", "FREERFLAG", "FREER_FLAG")


def _labels_by_name(labels: list[str]) -> dict[str, str]:
    return {label.upper(): label for label in labels}


def _column_group(found: dict[str, str], names: tuple[str, ...]) -> str:
    inner = ",".join(found[name] for name in names)
    return f"'/*/*/[{inner}]'"


def intensity_truncate_args(labels: list[str]) -> str:
    """ctruncate column selection for mean intensities or a Friedel pair.

    Depositions such as 6w1v store ``I`` and ``SIGI``. Depositions such as
    9z35 store only ``I(+)``, ``SIGI(+)``, ``I(-)``, and ``SIGI(-)``.
    """
    found = _labels_by_name(labels)
    if all(name in found for name in _MEAN_INTENSITIES):
        return "-colin " + _column_group(found, _MEAN_INTENSITIES)
    if all(name in found for name in _FRIEDEL_INTENSITIES):
        return "-colano " + _column_group(found, _FRIEDEL_INTENSITIES)
    raise ValueError(
        "Observed structure factors have no mean intensities (I, SIGI) "
        "and no complete Friedel pair (I(+), SIGI(+), I(-), SIGI(-)). "
        f"MTZ labels were: {labels}"
    )


def free_column(labels: list[str]) -> str | None:
    """FREE, or a FreeR_flag-style column, in the CCP4 convention REFMAC reads."""
    found = _labels_by_name(labels)
    for name in _FREE_LABELS:
        if name in found:
            return found[name]
    return None


def has_observed_intensities(labels: list[str]) -> bool:
    """Mean intensities, or a complete Friedel pair, as deposited structure factors store them."""
    found = _labels_by_name(labels)
    if all(name in found for name in _MEAN_INTENSITIES):
        return True
    return all(name in found for name in _FRIEDEL_INTENSITIES)


def amplitude_pair(labels: list[str]) -> tuple[str, str] | None:
    """First recognized mean-amplitude pair. ``F`` is preferred to ``FP`` when both exist."""
    found = _labels_by_name(labels)
    for amplitude, sigma in _AMPLITUDE_PAIRS:
        if amplitude in found and sigma in found:
            return found[amplitude], found[sigma]
    return None


def refmac_amplitude_cad_keywords(amplitude: str, sigma: str, free: str | None) -> str | None:
    """Rename amplitudes to ``F`` and ``SIGF``, keeping a free-set column from the same file.

    ``None`` means the columns are already ``F`` and ``SIGF``, so REFMAC can read the file unchanged.
    """
    if amplitude.upper() == "F" and sigma.upper() == "SIGF":
        return None
    labin = f"LABIN FILE 1 E1={amplitude} E2={sigma}"
    labout = "LABOUT FILE 1 E1=F E2=SIGF"
    if free:
        labin += f" E3={free}"
        labout += " E3=FREE"
    return f"{labin}\n{labout}\nEND"


def amplitude_cad_keywords(labels: list[str], free: str = "FREE") -> tuple[str, str, str]:
    """Mean amplitudes from ctruncate, and the cad script that adds the free-set column.

    ``-colin`` writes ``F`` and ``SIGF``. ``-colano`` writes ``FMEAN`` and
    ``SIGFMEAN`` alongside Friedel amplitudes. REFMAC reads ``F`` and ``SIGF``,
    so a Friedel run renames the mean pair and leaves the anomalous columns
    out of the file REFMAC sees. ``free`` is the column name in the intensity MTZ.
    """
    found = _labels_by_name(labels)
    if "F" in found and "SIGF" in found:
        return found["F"], found["SIGF"], f"LABIN FILE 1 ALL\nLABIN FILE 2 E1={free}\nEND"
    if "FMEAN" in found and "SIGFMEAN" in found:
        amplitude, sigma = found["FMEAN"], found["SIGFMEAN"]
        keywords = (
            f"LABIN FILE 1 E1={amplitude} E2={sigma}\n"
            "LABOUT FILE 1 E1=F E2=SIGF\n"
            f"LABIN FILE 2 E1={free}\n"
            "END"
        )
        return amplitude, sigma, keywords
    raise ValueError(
        "ctruncate did not write mean amplitudes F and SIGF, or FMEAN and SIGFMEAN. "
        f"MTZ labels were: {labels}"
    )


def _adopt_refmac_amplitudes(
    session: Ccp4Session,
    source: Path,
    amplitude: str,
    sigma: str,
    free: str | None,
    dest: Path,
    log_name: str,
) -> Path:
    """Return an MTZ whose mean amplitudes are ``F`` and ``SIGF``.

    A file that already uses those names is returned unchanged, including its ``FREE``
    column when that column is present. Any other recognized pair is renamed with ``cad``.
    """
    keywords = refmac_amplitude_cad_keywords(amplitude, sigma, free)
    if keywords is None:
        return source
    session.run_ccp4(
        "cad",
        [("hklin1", source), ("hklout", dest)],
        keywords,
        session.output_dir / log_name,
    )
    return dest


def run_fofo(cfg: dict, emit) -> dict:
    session = Ccp4Session(cfg["output_dir"], cfg["ccp4_root"], emit)
    files = cfg["files"]
    session.announce("1. Check software and input files")
    require_files(files)
    session.emit(f"CCP4 root : {session.ccp4_root}")
    session.emit(f"Coot      : {session.coot}")
    session.emit(f"PyMOL     : {session.pymol}")
    for key, path in files.items():
        session.emit(f"  {key:18s} {path}  ({path.stat().st_size / 1e6:.1f} MB)")
    session.run_bash('echo "CCP4=$CCP4"; which cif2mtz refmac5 cad scaleit fft mtzdump')

    session.announce("2. Inspect the light coordinates and the chosen site")
    atoms, centroid, site_line = _check_site(files["light_pdb"], cfg["site_chain"], cfg["site_resi"], cfg["site_label"])
    session.emit(site_line)
    session.emit("Atom names at the site: " + ", ".join(sorted(atoms)))
    session.emit(f"Centroid: {centroid[0]:.2f} {centroid[1]:.2f} {centroid[2]:.2f} Å")

    session.announce("3. Build the Fo−Fo difference map")
    dark_id = cfg["dark_id"]
    light_id = cfg["light_id"]
    out = session.output_dir
    dark_sf_mtz = out / f"{dark_id}-sf.mtz"
    light_sf_mtz = out / f"{light_id}-sf.mtz"
    if (not dark_sf_mtz.is_file()) or cfg["recompute_map"]:
        _cif2mtz(session, dark_id, files["dark_sf"], dark_sf_mtz, "sf")
    if (not light_sf_mtz.is_file()) or cfg["recompute_map"]:
        _cif2mtz(session, light_id, files["light_sf"], light_sf_mtz, "sf")
    low, high = _prepare_resolution(
        session, cfg["resolution_source"], cfg["lowres"], cfg["dmin"],
        [(dark_id, dark_sf_mtz), (light_id, light_sf_mtz)],
    )
    grid_high, grid_suffix = _sampling_grid(session, low, high, cfg.get("fft_dmin"))
    diff_map = out / f"{light_id}_minus_{dark_id}{grid_suffix}.ccp4"
    extended_map = out / f"{light_id}_minus_{dark_id}{grid_suffix}_extended.ccp4"

    if extended_map.is_file() and not cfg["recompute_map"]:
        session.emit(f"Reusing existing map: {extended_map}")
    else:
        dark_val_mtz = out / f"{dark_id}_validation_2fo-fc_map_coef.mtz"
        light_val_mtz = out / f"{light_id}_validation_2fo-fc_map_coef.mtz"
        dark_calc_mtz = out / f"{dark_id}_calc.mtz"
        light_calc_mtz = out / f"{light_id}_calc.mtz"
        _cif2mtz(session, dark_id, files["dark_validation"], dark_val_mtz, "val")
        _cif2mtz(session, light_id, files["light_validation"], light_val_mtz, "val")
        for pdb_id, sf_mtz, pdb_file, calc_mtz in (
            (dark_id, dark_sf_mtz, files["dark_pdb"], dark_calc_mtz),
            (light_id, light_sf_mtz, files["light_pdb"], light_calc_mtz),
        ):
            session.run_ccp4(
                "refmac5",
                [
                    ("hklin", sf_mtz),
                    ("hklout", calc_mtz),
                    ("xyzin", pdb_file),
                    ("xyzout", out / f"{pdb_id}_checked.pdb"),
                ],
                refmac_keywords(0),
                out / f"{pdb_id}_refmac0.log",
            )
        light_fp = session.pick_label(session.mtz_labels(light_calc_mtz), ["FP", "F", "FOBS", "FOSC"], "LIGHT amplitude")
        light_sig = session.pick_label(session.mtz_labels(light_calc_mtz), ["SIGFP", "SIGF", "SIGFOBS"], "LIGHT sigma")
        dark_fp = session.pick_label(session.mtz_labels(dark_calc_mtz), ["FP", "F", "FOBS", "FOSC"], "DARK amplitude")
        dark_sig = session.pick_label(session.mtz_labels(dark_calc_mtz), ["SIGFP", "SIGF", "SIGFOBS"], "DARK sigma")
        phase_map = {
            "light_validation": (light_val_mtz, ["PHWT", "PHIC", "PHI"]),
            "dark_validation": (dark_val_mtz, ["PHWT", "PHIC", "PHI"]),
            "light_refmac": (light_calc_mtz, ["PHWT", "PHIC", "PHI"]),
            "dark_refmac": (dark_calc_mtz, ["PHWT", "PHIC", "PHI"]),
        }
        phase_mtz, phase_candidates = phase_map[cfg["phase_source"]]
        phi_lab = session.pick_label(session.mtz_labels(phase_mtz), phase_candidates, "phase")
        session.emit(
            f"CAD labels: LIGHT {light_fp}/{light_sig}, DARK {dark_fp}/{dark_sig}, "
            f"PHI {phi_lab} from {phase_mtz.name}"
        )
        cad_mtz = out / f"cad_{light_id}_{dark_id}.mtz"
        session.run_ccp4(
            "cad",
            [("hklin1", light_calc_mtz), ("hklin2", dark_calc_mtz), ("hklin3", phase_mtz), ("hklout", cad_mtz)],
            (
                f"LABIN FILE 1 E1={light_fp} E2={light_sig}\n"
                "LABOUT FILE 1 E1=F_LIGHT E2=SIGF_LIGHT\n"
                "CTYPE FILE 1 E1=F E2=Q\n"
                f"LABIN FILE 2 E1={dark_fp} E2={dark_sig}\n"
                "LABOUT FILE 2 E1=F_DARK E2=SIGF_DARK\n"
                "CTYPE FILE 2 E1=F E2=Q\n"
                f"LABIN FILE 3 E1={phi_lab}\n"
                "LABOUT FILE 3 E1=PHI_REF\n"
                "CTYPE FILE 3 E1=P\n"
                "END"
            ),
            out / f"cad_{light_id}_{dark_id}.log",
        )
        scale_mtz = out / f"scaleit_{light_id}_{dark_id}.mtz"
        scale_log = out / f"scaleit_{light_id}_{dark_id}.log"
        session.run_ccp4(
            "scaleit",
            [("hklin", cad_mtz), ("hklout", scale_mtz)],
            (
                f"RESOLUTION {low:.2f} {high:.2f}\n"
                "REFINE ISOTROPIC\n"
                "LABIN FP=F_DARK SIGFP=SIGF_DARK FPH1=F_LIGHT SIGFPH1=SIGF_LIGHT\n"
                "END"
            ),
            scale_log,
        )
        rfac = parse_scaleit_rfactor(scale_log)
        if rfac:
            session.emit(
                f"SCALEIT overall R-factor = {rfac}  "
                "(≲ 0.10–0.15 is typical for a usable isomorphous pair)"
            )
        session.run_ccp4(
            "fft",
            [("hklin", scale_mtz), ("mapout", diff_map)],
            (
                f"RESOLUTION {low:.2f} {grid_high:.2f}\n"
                "LABIN F1=F_LIGHT SIG1=SIGF_LIGHT F2=F_DARK SIG2=SIGF_DARK PHI=PHI_REF\n"
                "END"
            ),
            out / f"fft_{light_id}_{dark_id}.log",
        )
        session.extend_around_model(
            diff_map,
            extended_map,
            files["light_pdb"],
            cfg["map_border"],
            out / f"pad_{light_id}_{dark_id}.log",
        )
        session.emit(f"Created: {extended_map}")
    session.emit(f"Interpretation: positive density is stronger in {light_id}; negative is stronger in {dark_id}.")

    session.announce("4. Export a site fragment with Coot")
    fragment = session.export_fragment(
        extended_map,
        files["light_pdb"],
        cfg["site_copies"],
        cfg["site_label"],
        cfg["export_radius"],
        cfg["positive_levels"][0],
        difference=True,
        recompute=cfg["recompute_fragment"],
    )
    additional_fragment = _export_additional_fragment(
        session,
        cfg,
        extended_map,
        files["light_pdb"],
        cfg["positive_levels"][0],
        True,
    )

    session.announce("5. Contour isodensity surfaces in PyMOL")
    png_path = out / f"{light_id}_minus_{dark_id}_{cfg['site_label']}_{cfg['site_chain']}{cfg['site_resi']}_diffmap.png"
    pml_path = out / f"pymol_{light_id}_minus_{dark_id}_{cfg['site_label']}.pml"
    session.render_figure(
        fragment=fragment,
        model_pdb=files["light_pdb"],
        map_name="DIFF",
        png_path=png_path,
        pml_path=pml_path,
        site_chain=cfg["site_chain"],
        site_resi=cfg["site_resi"],
        extra_residues=cfg["extra_residues"],
        positive_levels=cfg["positive_levels"],
        negative_levels=cfg["negative_levels"],
        positive_colors=POS_COLORS,
        negative_colors=NEG_COLORS,
        contour_mode=cfg["contour_mode"],
        surface=cfg["surface"],
        carve=cfg["carve"],
        show_oec_bonds=cfg["show_oec_bonds"],
        draw_distances=cfg["draw_distances"],
        oec_pairs=OEC_PAIRS_WITH_O6,
        png_width=cfg["png_width"],
        png_height=cfg["png_height"],
        ray_trace=cfg["ray_trace"],
    )
    sigma_txt = ", ".join(f"±{level:g}" for level in cfg["positive_levels"])
    caption = (
        f"{light_id} − {dark_id} difference density around {cfg['site_label']} "
        f"{cfg['site_chain']}/{cfg['site_resi']}, contoured at {sigma_txt} σ. "
        "Positive density (green, cyan, blue) is stronger in the light dataset."
    )
    return _result(
        caption,
        png_path,
        extended_map,
        fragment,
        pml_path,
        out,
        additional_fragment=additional_fragment,
        viewer=_viewer_scene(
            model_pdb=files["light_pdb"],
            fragment=fragment,
            cfg=cfg,
            positive_colors=POS_COLORS,
            negative_colors=NEG_COLORS,
            additional_fragment=additional_fragment,
        ),
    )


def run_refinement(cfg: dict, emit) -> dict:
    session = Ccp4Session(cfg["output_dir"], cfg["ccp4_root"], emit)
    files = cfg["files"]
    session.announce("1. Check software and input files")
    require_files(files)
    session.emit(f"CCP4 root : {session.ccp4_root}")
    session.emit(f"Coot      : {session.coot}")
    session.emit(f"PyMOL     : {session.pymol}")
    for key, path in files.items():
        session.emit(f"  {key:18s} {path}  ({path.stat().st_size / 1e6:.1f} MB)")
    tools = "ctruncate cad refmac5 fft mtzdump"
    if not is_reflection_mtz(files["light_sf"]):
        tools = "cif2mtz " + tools
    session.run_bash(f'echo "CCP4=$CCP4"; which {tools}')

    session.announce("2. Inspect the omit-model coordinates and the chosen site")
    explicit_center = _cartesian_center(cfg)
    atoms = residue_atom_names(files["dark_pdb"], cfg["site_chain"], cfg["site_resi"])
    session.emit(
        f"Configured site: chain {cfg['site_chain']}  residue {cfg['site_resi']}  label {cfg['site_label']}"
    )
    if atoms:
        session.emit("Atom names at the site: " + ", ".join(sorted(atoms)))
    elif explicit_center is None:
        raise ValueError(
            f"Residue {cfg['site_chain']}/{cfg['site_resi']} was not found in {files['dark_pdb'].name}. "
            "Use “List hetero residues” and set the chain and residue, or set the map center to a Cartesian point."
        )
    else:
        session.emit(
            f"Residue {cfg['site_chain']}/{cfg['site_resi']} has no atoms in {files['dark_pdb'].name}."
        )
    if explicit_center is None:
        centroid = residue_centroid(files["dark_pdb"], cfg["site_chain"], cfg["site_resi"])
        session.emit(f"Map center is the site centroid: {centroid[0]:.2f} {centroid[1]:.2f} {centroid[2]:.2f} Å")
    else:
        session.emit(
            "Map center is the Cartesian point: "
            f"{explicit_center[0]:.2f} {explicit_center[1]:.2f} {explicit_center[2]:.2f} Å"
        )
    session.emit(f"REFMAC cycles: {cfg['refmac_ncyc']}")

    session.announce("3. Build the mFo−DFc map")
    dark_id = cfg["dark_id"]
    light_id = cfg["light_id"]
    out = session.output_dir
    light_sf_mtz = out / f"{light_id}-sf.mtz"
    if (not light_sf_mtz.is_file()) or cfg["recompute_map"]:
        if is_reflection_mtz(files["light_sf"]):
            session.emit(f"Using the supplied reflection MTZ: {files['light_sf']}")
            install_reflection_mtz(files["light_sf"], light_sf_mtz)
        else:
            _cif2mtz(session, light_id, files["light_sf"], light_sf_mtz, "sf")
    low, high = _prepare_resolution(
        session, cfg["resolution_source"], cfg["lowres"], cfg["dmin"], [(light_id, light_sf_mtz)]
    )
    grid_high, grid_suffix = _sampling_grid(session, low, high, cfg.get("fft_dmin"))
    diff_map = out / f"{light_id}_DELFWT_{dark_id}{grid_suffix}.ccp4"
    extended_map = out / f"{light_id}_DELFWT_{dark_id}{grid_suffix}_extended.ccp4"
    refined_pdb = out / f"{dark_id}_refined.pdb"
    direct_cube = None

    if extended_map.is_file() and not cfg["recompute_map"]:
        session.emit(f"Reusing existing map: {extended_map}")
        if not refined_pdb.is_file():
            raise FileNotFoundError(
                f"Refined coordinates {refined_pdb.name} are missing. Turn map recomputation back on."
            )
        session.emit(f"Reusing refined model: {refined_pdb}")
    else:
        light_sf_truncate = out / f"{light_id}-sf-truncate.mtz"
        light_sf_truncate_free = out / f"{light_id}-sf-truncate-free.mtz"
        observed_labels = session.mtz_labels(light_sf_mtz)
        free = free_column(observed_labels)
        if has_observed_intensities(observed_labels):
            truncate_log = out / f"{light_id}_ctruncate_sf.log"
            column_args = intensity_truncate_args(observed_labels)
            session.emit(f"→ ctruncate  ({truncate_log.name})")
            session.emit(f"ctruncate {column_args}")
            session.run_bash(
                (
                    f"ctruncate -hklin { _q(light_sf_mtz) } "
                    f"-hklout { _q(light_sf_truncate) } "
                    f"{column_args}"
                ),
                truncate_log,
            )
            amplitude, sigma, cad_keywords = amplitude_cad_keywords(
                session.mtz_labels(light_sf_truncate), free or "FREE"
            )
            session.emit(f"Amplitude columns for REFMAC: {amplitude}, {sigma}")
            if free:
                session.run_ccp4(
                    "cad",
                    [
                        ("hklin1", light_sf_truncate),
                        ("hklin2", light_sf_mtz),
                        ("hklout", light_sf_truncate_free),
                    ],
                    cad_keywords,
                    out / f"{light_id}_cad_free.log",
                )
                refmac_hklin = light_sf_truncate_free
            else:
                session.emit("No FREE column was found, so every reflection is used.")
                refmac_hklin = _adopt_refmac_amplitudes(
                    session, light_sf_truncate, amplitude, sigma, None, light_sf_truncate_free, f"{light_id}_cad_free.log"
                )
        else:
            pair = amplitude_pair(observed_labels)
            if pair is None:
                raise ValueError(
                    "Observed reflections have neither intensities (I, SIGI or a complete Friedel pair) "
                    "nor amplitudes (F/SIGF, FP/SIGFP, FOBS/SIGFOBS, or FMEAN/SIGFMEAN). "
                    f"MTZ labels were: {observed_labels}"
                )
            amplitude, sigma = pair
            session.emit(
                "The reflection MTZ already contains amplitudes, so ctruncate is not run."
            )
            session.emit(f"Amplitude columns for REFMAC: {amplitude}, {sigma}")
            if free:
                session.emit(f"Free-set column: {free}")
            else:
                session.emit("No FREE column was found, so every reflection is used.")
            refmac_hklin = _adopt_refmac_amplitudes(
                session, light_sf_mtz, amplitude, sigma, free, light_sf_truncate_free, f"{light_id}_cad_free.log"
            )
        omit_mtz = out / f"{dark_id}_vs_{light_id}.mtz"
        refmac_files = [
            ("hklin", refmac_hklin),
            ("hklout", omit_mtz),
            ("xyzin", files["dark_pdb"]),
            ("xyzout", refined_pdb),
        ]
        custom_dictionary = cfg.get("custom_dictionary")
        if custom_dictionary and int(cfg["refmac_ncyc"]) > 0:
            refmac_files.append(("libin", custom_dictionary))
            session.emit(
                f"Custom dictionary: {custom_dictionary}. "
                "External restraints are read after this file and replace matching distances, angles, and dihedrals."
            )
        elif custom_dictionary:
            session.emit(
                "Custom dictionary was not used. REFMAC cycles is 0, "
                "so the CCP4 monomer library is used."
            )
        site_missing = explicit_center is not None and not atoms

        def _restraints(label: str, build) -> str:
            try:
                return build()
            except ValueError as exc:
                if not site_missing:
                    raise
                session.emit(f"{label} were not applied. {exc}")
                return ""

        distance_text = _restraints(
            "Distance restraints",
            lambda: external_distance_keywords(
                files["dark_pdb"], cfg["distance_restraints"], cfg["distance_sigma"]
            ),
        )
        angle_text = _restraints(
            "Angle restraints",
            lambda: external_angle_keywords(
                files["dark_pdb"], cfg["angle_restraints"], cfg["angle_sigma"]
            ),
        )
        dihedral_text = _restraints(
            "Dihedral restraints",
            lambda: external_dihedral_keywords(
                files["dark_pdb"],
                cfg.get("dihedral_restraints") or [],
                cfg.get("dihedral_sigma", 10.0),
            ),
        )
        geometry_text = "\n".join(part for part in (distance_text, angle_text, dihedral_text) if part)
        harmonic_atoms = list(cfg.get("harmonic_atoms") or [])
        harmonic_text = ""
        if harmonic_atoms and int(cfg["refmac_ncyc"]) > 0:
            harmonic_sigma = float(cfg.get("harmonic_sigma", 0.1))
            harmonic_text = external_harmonic_keywords(harmonic_atoms, harmonic_sigma)
            session.emit(
                f"Harmonic restraints on {len(harmonic_atoms)} atoms, sigma {harmonic_sigma:.3f} Å. "
                "Each atom is held at the coordinates in the input model."
            )
        elif harmonic_atoms:
            session.emit(
                "Harmonic restraints were not applied. REFMAC cycles is 0."
            )
        protein_text = ""
        if cfg.get("harmonic_protein") and int(cfg["refmac_ncyc"]) > 0:
            protein_sigma = float(cfg.get("harmonic_protein_sigma", 0.1))
            protein_text = protein_harmonic_lines(files["dark_pdb"], cfg, harmonic_atoms)
            if protein_text:
                protein_count = protein_text.count("\n") + 1
                session.emit(
                    f"Harmonic restraints on {protein_count} protein atoms outside the site "
                    f"and the additional residues, sigma {protein_sigma:.3f} Å. "
                    "Each atom is held at the coordinates in the input model."
                )
            else:
                session.emit(
                    "Protein harmonic restraints were requested, but the model has no protein atom "
                    "outside the site and the additional residues."
                )
        elif cfg.get("harmonic_protein"):
            session.emit(
                "Protein harmonic restraints were not applied. REFMAC cycles is 0."
            )
        restraint_text = "\n".join(
            part for part in (geometry_text, harmonic_text, protein_text) if part
        )
        exclude_specs = list(cfg.get("exclude_restraints") or [])
        exclude_text = exclude_restraint_keywords(exclude_specs)
        if exclude_text:
            session.emit("REFMAC will leave these atoms out of restraints:\n" + exclude_text)
        if restraint_text:
            restraint_path = out / f"{dark_id}_external_restraints.txt"
            restraint_path.write_text(restraint_text + "\n")
            session.emit(f"External restraints: {restraint_path}")
            # The protein block can be one line per atom of a photosystem. The
            # file keeps every line. The log keeps the short restraints and the count above.
            logged = "\n".join(part for part in (geometry_text, harmonic_text) if part)
            if logged:
                session.emit(logged)
        session.run_ccp4(
            "refmac5",
            refmac_files,
            refmac_keywords(
                int(cfg["refmac_ncyc"]),
                restraint_text,
                bref=cfg.get("refmac_bref", "ISOT"),
                refi_type=cfg.get("refmac_type", "REST"),
                exclude=exclude_specs,
            ),
            out / f"{dark_id}_vs_{light_id}_refmac{int(cfg['refmac_ncyc'])}.log",
        )
        omit_labels = session.mtz_labels(omit_mtz)
        fo_lab = session.pick_label(omit_labels, ["FP", "F", "FOBS", "FOSC"], "Fo")
        fc_lab = session.pick_label(omit_labels, ["FC_ALL", "FC"], "Fc")
        delfwt_lab = session.pick_label(omit_labels, ["DELFWT"], "DELFWT")
        phdelwt_lab = session.pick_label(omit_labels, ["PHDELWT"], "PHDELWT")
        session.emit(f"REFMAC {omit_mtz.name}: Fo={fo_lab}, Fc={fc_lab}, map={delfwt_lab}/{phdelwt_lab}")
        fft_log = out / f"fft_{light_id}_DELFWT_{dark_id}.log"
        try:
            session.run_ccp4(
                "fft",
                [("hklin", omit_mtz), ("mapout", diff_map)],
                (
                    # The high limit sets the sampling grid as well as the coefficient
                    # cutoff. A finer limit than the MTZ extends the reciprocal box;
                    # missing reflections stay out of the sum and contribute zero.
                    f"RESOLUTION {low:.2f} {grid_high:.2f}\n"
                    f"LABIN F1={delfwt_lab} PHI={phdelwt_lab}\n"
                    "END"
                ),
                fft_log,
            )
        except CommandFailed as exc:
            text = "\n".join(exc.lines)
            sampling = parse_fft_sampling(text)
            # A product past 2^31-1 either stops in CCPALC or, once the wrapped
            # length is positive, segfaults. Both mean the whole cell was never written.
            if sampling is None or not whole_cell_exceeds_fft(sampling):
                raise
            if diff_map.is_file():
                diff_map.unlink()
            center = explicit_center or residue_centroid(refined_pdb, cfg["site_chain"], cfg["site_resi"])
            direct_cube = {
                "mtz": omit_mtz,
                "f": delfwt_lab,
                "phi": phdelwt_lab,
                "low": low,
                "high": grid_high,
                "sampling": sampling,
            }
            direct_log = out / f"direct_{light_id}_DELFWT_{dark_id}.log"
            text = direct_site_cube(
                omit_mtz,
                extended_map,
                delfwt_lab,
                phdelwt_lab,
                low,
                grid_high,
                sampling,
                center,
                cfg["export_radius"],
                None,
                emit=session.emit,
            )
            direct_log.write_text(text)
        else:
            session.extend_around_model(
                diff_map,
                extended_map,
                refined_pdb,
                cfg["map_border"],
                out / f"pad_{light_id}_DELFWT_{dark_id}.log",
            )
        session.emit(f"Created: {extended_map}")
        session.emit(f"Refined model ({int(cfg['refmac_ncyc'])} cycles): {refined_pdb}")
    session.emit(
        f"Interpretation: positive density is present in the {light_id} observations "
        f"and absent from the {dark_id} model."
    )

    session.announce("4. Export a site fragment with Coot")
    fragment = session.export_fragment(
        extended_map,
        refined_pdb,
        cfg["site_copies"],
        cfg["site_label"],
        cfg["export_radius"],
        cfg["positive_levels"][0],
        difference=True,
        recompute=cfg["recompute_fragment"],
        center=explicit_center,
    )
    additional_source = extended_map
    if direct_cube and cfg.get("additional_site"):
        additional_source = _additional_cube(session, cfg, refined_pdb, extended_map, direct_cube)
    additional_fragment = _export_additional_fragment(
        session,
        cfg,
        additional_source,
        refined_pdb,
        cfg["positive_levels"][0],
        True,
    )

    session.announce("5. Contour isodensity surfaces in PyMOL")
    png_path = out / (
        f"{light_id}_DELFWT_{dark_id}_{cfg['site_label']}_{cfg['site_chain']}{cfg['site_resi']}_diffmap.png"
    )
    pml_path = out / f"pymol_{light_id}_DELFWT_{dark_id}_{cfg['site_label']}.pml"
    session.render_figure(
        fragment=fragment,
        model_pdb=refined_pdb,
        map_name="DIFF",
        png_path=png_path,
        pml_path=pml_path,
        site_chain=cfg["site_chain"],
        site_resi=cfg["site_resi"],
        extra_residues=cfg["extra_residues"],
        positive_levels=cfg["positive_levels"],
        negative_levels=cfg["negative_levels"],
        positive_colors=POS_COLORS,
        negative_colors=NEG_COLORS,
        contour_mode=cfg["contour_mode"],
        surface=cfg["surface"],
        carve=cfg["carve"],
        show_oec_bonds=cfg["show_oec_bonds"],
        draw_distances=cfg["draw_distances"],
        oec_pairs=OEC_PAIRS,
        png_width=cfg["png_width"],
        png_height=cfg["png_height"],
        ray_trace=cfg["ray_trace"],
        geometry_bonds=bool(cfg.get("geometry_bonds")),
        render_png=bool(cfg.get("render_png", True)),
    )
    sigma_txt = ", ".join(f"±{level:g}" for level in cfg["positive_levels"])
    grid_sentence = ""
    if grid_high < high - 1e-4:
        grid_sentence = (
            f" Sampled on the {grid_high:.2f} Å grid, with zero amplitude "
            f"from {high:.2f} Å to that limit."
        )
    if explicit_center is None:
        place = f"{cfg['site_chain']}/{cfg['site_resi']}"
    else:
        place = "at {:.2f}, {:.2f}, {:.2f} Å".format(*explicit_center)
    caption = (
        f"mFo−DFc density of {dark_id} against {light_id} around {cfg['site_label']} "
        f"{place}, {cfg['refmac_ncyc']} REFMAC cycles, "
        f"contoured at {sigma_txt} σ. Positive density is in the data and missing from the omit model."
        f"{grid_sentence}"
    )
    return _result(
        caption,
        png_path,
        extended_map,
        fragment,
        pml_path,
        out,
        extra=[refined_pdb],
        additional_fragment=additional_fragment,
        viewer=_viewer_scene(
            model_pdb=refined_pdb,
            fragment=fragment,
            cfg=cfg,
            positive_colors=POS_COLORS,
            negative_colors=NEG_COLORS,
            light_pdb=files["light_pdb"],
            additional_fragment=additional_fragment,
        ),
    )


def run_twofofc(cfg: dict, emit) -> dict:
    session = Ccp4Session(cfg["output_dir"], cfg["ccp4_root"], emit)
    files = cfg["files"]
    session.announce("1. Check software and input files")
    require_files(files)
    session.emit(f"CCP4 root : {session.ccp4_root}")
    session.emit(f"Coot      : {session.coot}")
    session.emit(f"PyMOL     : {session.pymol}")
    for key, path in files.items():
        session.emit(f"  {key:18s} {path}  ({path.stat().st_size / 1e6:.1f} MB)")
    session.run_bash('echo "CCP4=$CCP4"; which cif2mtz fft mtzdump')

    session.announce("2. Inspect the coordinates and the chosen site")
    atoms, centroid, site_line = _check_site(files["pdb"], cfg["site_chain"], cfg["site_resi"], cfg["site_label"])
    session.emit(site_line)
    session.emit("Atom names at the site: " + ", ".join(sorted(atoms)))
    session.emit(f"Centroid: {centroid[0]:.2f} {centroid[1]:.2f} {centroid[2]:.2f} Å")

    session.announce("3. Build the sharpened 2mFo−DFc map")
    pdb_id = cfg["pdb_id"]
    out = session.output_dir
    sf_mtz = out / f"{pdb_id}-sf.mtz"
    val_mtz = out / f"{pdb_id}_validation_2fo-fc_map_coef.mtz"
    b_tag = delta_b_tag(cfg["delta_b"])
    if (not sf_mtz.is_file()) or cfg["recompute_map"]:
        _cif2mtz(session, pdb_id, files["sf"], sf_mtz, "sf")
    low, high = _prepare_resolution(
        session, cfg["resolution_source"], cfg["lowres"], cfg["dmin"], [(pdb_id, sf_mtz)]
    )
    grid_high, grid_suffix = _sampling_grid(session, low, high, cfg.get("fft_dmin"))
    density_map = out / f"{pdb_id}_2fofc_{b_tag}{grid_suffix}.ccp4"
    extended_map = out / f"{pdb_id}_2fofc_{b_tag}{grid_suffix}_extended.ccp4"

    if extended_map.is_file() and not cfg["recompute_map"]:
        session.emit(f"Reusing existing map: {extended_map}")
    else:
        _cif2mtz(session, pdb_id, files["validation"], val_mtz, "val")
        val_labels = session.mtz_labels(val_mtz)
        fwt_lab = session.pick_label(val_labels, ["FWT"], "FWT")
        phwt_lab = session.pick_label(val_labels, ["PHWT"], "PHWT")
        val_low, val_high = session.mtz_resolution_angstrom(val_mtz)
        session.emit(
            f"Validation coefficients {val_mtz.name}: {val_low:.2f}–{val_high:.2f} Å  "
            f"map={fwt_lab}/{phwt_lab}  columns={val_labels}"
        )
        if high + 1e-3 < val_high:
            session.emit(
                f"WARNING: high-resolution cutoff {high:.2f} Å is finer than the {fwt_lab} "
                f"limit {val_high:.2f} Å. Those extra shells are not in the coefficient file."
            )
        stol2 = 1.0 / (4.0 * high * high)
        factor = math.exp(-cfg["delta_b"] * stol2)
        session.emit(
            f"Sharpening: {fwt_lab} *= exp(-ΔB (sinθ/λ)^2) with ΔB = {cfg['delta_b']:.2f} Å^2"
        )
        session.emit(f"  amplitude multiplier at {high:.2f} Å = {factor:.3f}")
        session.run_ccp4(
            "fft",
            [("hklin", val_mtz), ("mapout", density_map)],
            (
                f"TITLE {pdb_id} weighted 2mFo-DFc sharpened by dB={cfg['delta_b']:.2f}\n"
                f"RESOLUTION {low:.2f} {grid_high:.2f}\n"
                f"LABIN F1={fwt_lab} PHI={phwt_lab}\n"
                f"SCALE F1 1.0 {cfg['delta_b']:.2f}\n"
                "END"
            ),
            out / f"fft_{pdb_id}_2fofc_{b_tag}.log",
        )
        session.extend_around_model(
            density_map,
            extended_map,
            files["pdb"],
            cfg["map_border"],
            out / f"pad_{pdb_id}_2fofc_{b_tag}.log",
        )
        session.emit(f"Created: {extended_map}")
    session.emit(
        f"Interpretation: blue mesh is sharpened 2mFo−DFc density for {pdb_id} "
        f"(ΔB = {cfg['delta_b']:.1f} Å²), not a difference map."
    )

    session.announce("4. Export a site fragment with Coot")
    fragment = session.export_fragment(
        extended_map,
        files["pdb"],
        cfg["site_copies"],
        cfg["site_label"],
        cfg["export_radius"],
        cfg["sigma_levels"][0],
        difference=False,
        recompute=cfg["recompute_fragment"],
    )
    additional_fragment = _export_additional_fragment(
        session,
        cfg,
        extended_map,
        files["pdb"],
        cfg["sigma_levels"][0],
        False,
    )

    session.announce("5. Contour isodensity surfaces in PyMOL")
    png_path = out / f"{pdb_id}_2fofc_{b_tag}_{cfg['site_label']}_{cfg['site_chain']}{cfg['site_resi']}.png"
    pml_path = out / f"pymol_{pdb_id}_2fofc_{b_tag}_{cfg['site_label']}.pml"
    session.render_figure(
        fragment=fragment,
        model_pdb=files["pdb"],
        map_name="TWOFOFC",
        png_path=png_path,
        pml_path=pml_path,
        site_chain=cfg["site_chain"],
        site_resi=cfg["site_resi"],
        extra_residues=cfg["extra_residues"],
        positive_levels=cfg["sigma_levels"],
        negative_levels=None,
        positive_colors=MAP_COLORS,
        negative_colors=NEG_COLORS,
        contour_mode=cfg["contour_mode"],
        surface=cfg["surface"],
        carve=cfg["carve"],
        show_oec_bonds=cfg["show_oec_bonds"],
        draw_distances=cfg["draw_distances"],
        oec_pairs=OEC_PAIRS,
        png_width=cfg["png_width"],
        png_height=cfg["png_height"],
        ray_trace=cfg["ray_trace"],
    )
    sigma_txt = ", ".join(f"{level:g}" for level in cfg["sigma_levels"])
    caption = (
        f"{pdb_id} sharpened 2mFo−DFc density (ΔB = {cfg['delta_b']:.0f} Å²) around "
        f"{cfg['site_label']} {cfg['site_chain']}/{cfg['site_resi']}, contoured at {sigma_txt} σ."
    )
    return _result(
        caption,
        png_path,
        extended_map,
        fragment,
        pml_path,
        out,
        additional_fragment=additional_fragment,
        viewer=_viewer_scene(
            model_pdb=files["pdb"],
            fragment=fragment,
            cfg=cfg,
            positive_colors=MAP_COLORS,
            negative_colors=[],
            additional_fragment=additional_fragment,
        ),
    )


def initial_viewer_sigma(
    minimum: float,
    maximum: float,
    levels: list[float],
    negative: list[float] | None = None,
) -> float:
    """Pick the starting contour: the level closest to zero that sits inside the slider range."""
    candidates = [float(level) for level in levels]
    if negative:
        candidates.extend(-abs(float(level)) for level in negative)
    inside = [level for level in candidates if minimum - 1e-9 <= level <= maximum + 1e-9]
    if not inside:
        return (minimum + maximum) / 2.0
    return min(inside, key=lambda level: (abs(level), level < 0, level))


def _q(path: Path) -> str:
    import shlex
    return shlex.quote(str(path))


def _viewer_scene(
    *,
    model_pdb: Path,
    fragment: Path,
    cfg: dict,
    positive_colors: list[str],
    negative_colors: list[str],
    light_pdb: Path | None = None,
    additional_fragment: Path | None = None,
) -> dict:
    """Return the data needed to recreate the local PyMOL scene in WebGL.

    The PNG still uses every sigma level. The viewer draws one surface and
    slides its contour between ``viewer_sigma_min`` and ``viewer_sigma_max``.
    ``light_pdb`` is the omit workflow's observed-dataset coordinates. The
    viewer can draw that site in ball and stick; the map is unchanged.
    """
    levels = list(cfg["positive_levels"] if "positive_levels" in cfg else cfg["sigma_levels"])
    negative = cfg.get("negative_levels")
    sigma_min = float(cfg["viewer_sigma_min"])
    sigma_max = float(cfg["viewer_sigma_max"])
    scene = {
        "model": str(model_pdb),
        "model_ext": model_pdb.suffix.lstrip(".").lower() or "pdb",
        "map": str(fragment),
        "site_label": cfg.get("site_label") or "",
        "site": {"chain": cfg["site_chain"], "resi": cfg["site_resi"]},
        "extra_residues": cfg["extra_residues"],
        "positive_levels": levels,
        "negative_levels": negative,
        "positive_colors": positive_colors,
        "negative_colors": negative_colors,
        "contour_mode": cfg["contour_mode"],
        "surface": cfg["surface"],
        "sigma_min": sigma_min,
        "sigma_max": sigma_max,
        "sigma": initial_viewer_sigma(sigma_min, sigma_max, levels, negative),
    }
    if cfg.get("viewer_light_label"):
        scene["light_label"] = str(cfg["viewer_light_label"])
    if light_pdb is not None:
        scene["light_model"] = str(light_pdb)
        scene["light_model_ext"] = light_pdb.suffix.lstrip(".").lower() or "pdb"
        # The reference site has atoms the omit model lacks. Its sticks are
        # calculated from these coordinates, so a missing omit atom still
        # shows the bonds around it.
        if cfg.get("geometry_bonds"):
            scene["light_bonds"] = reasonable_bond_specs(
                light_pdb,
                cfg["site_chain"],
                int(cfg["site_resi"]),
                cfg.get("extra_residues"),
            )
    site = cfg.get("additional_site")
    if additional_fragment is not None and site:
        additional = {
            "label": site["label"],
            "map": str(additional_fragment),
            "site": {"chain": site["chain"], "resi": int(site["resi"])},
            "extra_residues": site.get("extra_residues") or [],
        }
        if light_pdb is not None and cfg.get("geometry_bonds"):
            additional["light_bonds"] = reasonable_bond_specs(
                light_pdb,
                str(site["chain"]),
                int(site["resi"]),
                site.get("extra_residues"),
            )
        scene["additional"] = additional
    return scene


def _result(
    caption,
    png_path,
    extended_map,
    fragment,
    pml_path,
    output_dir,
    extra=None,
    viewer=None,
    additional_fragment: Path | None = None,
) -> dict:
    named = [
        ("Figure", png_path),
        ("Extended map", extended_map),
        ("Site fragment", fragment),
        ("Additional site map", additional_fragment),
        ("PyMOL script", pml_path),
        ("PyMOL render script", pml_path.with_name(pml_path.stem + "_render.pml")),
    ]
    if extra:
        named.extend(("Refined coordinates" if path.name.endswith("_refined.pdb") else "Output", path) for path in extra)
    files = [{"role": role, "path": str(path)} for role, path in named if path is not None and path.is_file()]
    manifest = {
        "caption": caption,
        "output_dir": str(output_dir),
        "files": files,
    }
    (output_dir / "webapp_manifest.json").write_text(
        __import__("json").dumps(manifest, indent=2) + "\n"
    )
    return {
        "caption": caption,
        "figure": str(png_path),
        "final_map": str(extended_map),
        "fragment": str(fragment),
        "pymol_script": str(pml_path),
        "output_dir": str(output_dir),
        "files": files,
        "viewer": viewer,
        "additional_fragment": str(additional_fragment) if additional_fragment else None,
    }


def run_reconstruct(cfg: dict, emit, pause=None, report_step=None, report_row=None) -> dict:
    from xfel_web.reconstruction import run_reconstruct as _run

    return _run(
        cfg,
        emit,
        pause=pause,
        omit_runner=run_refinement,
        report_step=report_step,
        report_row=report_row,
    )


RUNNERS = {
    "fofo": run_fofo,
    "refinement": run_refinement,
    "twofofc": run_twofofc,
    "reconstruct": run_reconstruct,
}
