"""The three notebook calculations: Fo−Fo, omit-model mFo−DFc, and 2Fo−Fc."""

from __future__ import annotations

import math
from pathlib import Path

from xfel_web.ccp4 import (
    MAP_COLORS,
    NEG_COLORS,
    OEC_PAIRS,
    OEC_PAIRS_WITH_O6,
    POS_COLORS,
    Ccp4Session,
    choose_resolution,
    delta_b_tag,
    parse_scaleit_rfactor,
    require_files,
    res_tag,
    residue_atom_names,
    residue_centroid,
)
from xfel_web.restraints import (
    build_refinement_library,
    external_angle_keywords,
    external_distance_keywords,
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
    for label, mtz in mtzs:
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


_MEAN_INTENSITIES = ("I", "SIGI")
_FRIEDEL_INTENSITIES = ("I(+)", "SIGI(+)", "I(-)", "SIGI(-)")


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


def amplitude_cad_keywords(labels: list[str]) -> tuple[str, str, str]:
    """Mean amplitudes from ctruncate, and the cad script that adds FREE.

    ``-colin`` writes ``F`` and ``SIGF``. ``-colano`` writes ``FMEAN`` and
    ``SIGFMEAN`` alongside Friedel amplitudes. REFMAC reads ``F`` and ``SIGF``,
    so a Friedel run renames the mean pair and leaves the anomalous columns
    out of the file REFMAC sees.
    """
    found = _labels_by_name(labels)
    if "F" in found and "SIGF" in found:
        return found["F"], found["SIGF"], "LABIN FILE 1 ALL\nLABIN FILE 2 E1=FREE\nEND"
    if "FMEAN" in found and "SIGFMEAN" in found:
        amplitude, sigma = found["FMEAN"], found["SIGFMEAN"]
        keywords = (
            f"LABIN FILE 1 E1={amplitude} E2={sigma}\n"
            "LABOUT FILE 1 E1=F E2=SIGF\n"
            "LABIN FILE 2 E1=FREE\n"
            "END"
        )
        return amplitude, sigma, keywords
    raise ValueError(
        "ctruncate did not write mean amplitudes F and SIGF, or FMEAN and SIGFMEAN. "
        f"MTZ labels were: {labels}"
    )


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
    session.run_bash('echo "CCP4=$CCP4"; which cif2mtz refmac5 cad scaleit fft extends mtzdump')

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
    tag = res_tag(low, high)
    diff_map = out / f"{light_id}_minus_{dark_id}_{tag}.ccp4"
    extended_map = out / f"{light_id}_minus_{dark_id}_{tag}_extended.ccp4"

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
                "NCYC 0\nEND",
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
                f"RESOLUTION {low:.2f} {high:.2f}\n"
                "LABIN F1=F_LIGHT SIG1=SIGF_LIGHT F2=F_DARK SIG2=SIGF_DARK PHI=PHI_REF\n"
                "END"
            ),
            out / f"fft_{light_id}_{dark_id}.log",
        )
        session.run_ccp4(
            "extends",
            [("mapin", diff_map), ("mapout", extended_map), ("xyzin", files["light_pdb"])],
            f"BORDER {cfg['map_border']:.1f}\nEND",
            out / f"extends_{light_id}_{dark_id}.log",
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
        viewer=_viewer_scene(
            model_pdb=files["light_pdb"],
            fragment=fragment,
            cfg=cfg,
            positive_colors=POS_COLORS,
            negative_colors=NEG_COLORS,
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
    tools = "cif2mtz ctruncate cad refmac5 fft extends mtzdump"
    if cfg["use_libcheck"]:
        tools += " libcheck"
    session.run_bash(f'echo "CCP4=$CCP4"; which {tools}')

    session.announce("2. Inspect the omit-model coordinates and the chosen site")
    atoms, centroid, site_line = _check_site(files["dark_pdb"], cfg["site_chain"], cfg["site_resi"], cfg["site_label"])
    session.emit(site_line)
    session.emit("Atom names at the site: " + ", ".join(sorted(atoms)))
    session.emit(f"Centroid: {centroid[0]:.2f} {centroid[1]:.2f} {centroid[2]:.2f} Å")
    session.emit(f"REFMAC cycles: {cfg['refmac_ncyc']}")

    session.announce("3. Build the mFo−DFc map")
    dark_id = cfg["dark_id"]
    light_id = cfg["light_id"]
    out = session.output_dir
    light_sf_mtz = out / f"{light_id}-sf.mtz"
    if (not light_sf_mtz.is_file()) or cfg["recompute_map"]:
        _cif2mtz(session, light_id, files["light_sf"], light_sf_mtz, "sf")
    low, high = _prepare_resolution(
        session, cfg["resolution_source"], cfg["lowres"], cfg["dmin"], [(light_id, light_sf_mtz)]
    )
    tag = res_tag(low, high)
    diff_map = out / f"{light_id}_DELFWT_{dark_id}_{tag}.ccp4"
    extended_map = out / f"{light_id}_DELFWT_{dark_id}_{tag}_extended.ccp4"
    refined_pdb = out / f"{dark_id}_refined.pdb"

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
        truncate_log = out / f"{light_id}_ctruncate_sf.log"
        column_args = intensity_truncate_args(session.mtz_labels(light_sf_mtz))
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
        amplitude, sigma, cad_keywords = amplitude_cad_keywords(session.mtz_labels(light_sf_truncate))
        session.emit(f"Amplitude columns for REFMAC: {amplitude}, {sigma}")
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
        omit_mtz = out / f"{dark_id}_vs_{light_id}.mtz"
        refmac_files = [
            ("hklin", light_sf_truncate_free),
            ("hklout", omit_mtz),
            ("xyzin", files["dark_pdb"]),
            ("xyzout", refined_pdb),
        ]
        angle_text = external_angle_keywords(
            files["dark_pdb"], cfg["site_copies"], cfg["angle_restraints"], cfg["angle_sigma"]
        )
        if cfg["use_libcheck"]:
            library_cif = out / f"{dark_id}_lib.cif"
            build_refinement_library(
                session,
                files["dark_pdb"],
                library_cif,
                cfg["site_chain"],
                cfg["site_resi"],
                cfg["distance_restraints"],
                cfg["distance_sigma"],
                cfg["lib_monomers"],
            )
            refmac_files.append(("libin", library_cif))
            restraint_text = angle_text
        else:
            distance_text = external_distance_keywords(
                files["dark_pdb"], cfg["site_copies"], cfg["distance_restraints"], cfg["distance_sigma"]
            )
            restraint_text = "\n".join(part for part in (distance_text, angle_text) if part)
        if restraint_text:
            restraint_path = out / f"{dark_id}_external_restraints.txt"
            restraint_path.write_text(restraint_text + "\n")
            session.emit(f"External restraints: {restraint_path}")
            session.emit(restraint_text)
        session.run_ccp4(
            "refmac5",
            refmac_files,
            (
                (restraint_text + "\n" if restraint_text else "")
                + f"NCYC {int(cfg['refmac_ncyc'])}\n"
                "MAKE HYDR Y\n"
                "MAKE HOUT Y\n"
                "END"
            ),
            out / f"{dark_id}_vs_{light_id}_refmac{int(cfg['refmac_ncyc'])}.log",
        )
        omit_labels = session.mtz_labels(omit_mtz)
        fo_lab = session.pick_label(omit_labels, ["FP", "F", "FOBS", "FOSC"], "Fo")
        fc_lab = session.pick_label(omit_labels, ["FC_ALL", "FC"], "Fc")
        delfwt_lab = session.pick_label(omit_labels, ["DELFWT"], "DELFWT")
        phdelwt_lab = session.pick_label(omit_labels, ["PHDELWT"], "PHDELWT")
        session.emit(f"REFMAC {omit_mtz.name}: Fo={fo_lab}, Fc={fc_lab}, map={delfwt_lab}/{phdelwt_lab}")
        session.run_ccp4(
            "fft",
            [("hklin", omit_mtz), ("mapout", diff_map)],
            (
                f"RESOLUTION {low:.2f} {high:.2f}\n"
                f"LABIN F1={delfwt_lab} PHI={phdelwt_lab}\n"
                "END"
            ),
            out / f"fft_{light_id}_DELFWT_{dark_id}.log",
        )
        session.run_ccp4(
            "extends",
            [("mapin", diff_map), ("mapout", extended_map), ("xyzin", refined_pdb)],
            f"BORDER {cfg['map_border']:.1f}\nEND",
            out / f"extends_{light_id}_DELFWT_{dark_id}.log",
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
    )
    sigma_txt = ", ".join(f"±{level:g}" for level in cfg["positive_levels"])
    caption = (
        f"mFo−DFc density of {dark_id} against {light_id} around {cfg['site_label']} "
        f"{cfg['site_chain']}/{cfg['site_resi']}, {cfg['refmac_ncyc']} REFMAC cycles, "
        f"contoured at {sigma_txt} σ. Positive density is in the data and missing from the omit model."
    )
    return _result(
        caption,
        png_path,
        extended_map,
        fragment,
        pml_path,
        out,
        extra=[refined_pdb],
        viewer=_viewer_scene(
            model_pdb=refined_pdb,
            fragment=fragment,
            cfg=cfg,
            positive_colors=POS_COLORS,
            negative_colors=NEG_COLORS,
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
    session.run_bash('echo "CCP4=$CCP4"; which cif2mtz fft extends mtzdump')

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
    tag = res_tag(low, high)
    density_map = out / f"{pdb_id}_2fofc_{b_tag}_{tag}.ccp4"
    extended_map = out / f"{pdb_id}_2fofc_{b_tag}_{tag}_extended.ccp4"

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
                f"RESOLUTION {low:.2f} {high:.2f}\n"
                f"LABIN F1={fwt_lab} PHI={phwt_lab}\n"
                f"SCALE F1 1.0 {cfg['delta_b']:.2f}\n"
                "END"
            ),
            out / f"fft_{pdb_id}_2fofc_{b_tag}.log",
        )
        session.run_ccp4(
            "extends",
            [("mapin", density_map), ("mapout", extended_map), ("xyzin", files["pdb"])],
            f"BORDER {cfg['map_border']:.1f}\nEND",
            out / f"extends_{pdb_id}_2fofc_{b_tag}.log",
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
        viewer=_viewer_scene(
            model_pdb=files["pdb"],
            fragment=fragment,
            cfg=cfg,
            positive_colors=MAP_COLORS,
            negative_colors=[],
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


def _viewer_scene(*, model_pdb: Path, fragment: Path, cfg: dict, positive_colors: list[str], negative_colors: list[str]) -> dict:
    """Return the data needed to recreate the local PyMOL scene in WebGL.

    The PNG still uses every sigma level. The viewer draws one surface and
    slides its contour between ``viewer_sigma_min`` and ``viewer_sigma_max``.
    """
    levels = list(cfg["positive_levels"] if "positive_levels" in cfg else cfg["sigma_levels"])
    negative = cfg.get("negative_levels")
    sigma_min = float(cfg["viewer_sigma_min"])
    sigma_max = float(cfg["viewer_sigma_max"])
    return {
        "model": str(model_pdb),
        "model_ext": model_pdb.suffix.lstrip(".").lower() or "pdb",
        "map": str(fragment),
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


def _result(caption, png_path, extended_map, fragment, pml_path, output_dir, extra=None, viewer=None) -> dict:
    named = [
        ("Figure", png_path),
        ("Extended map", extended_map),
        ("Site fragment", fragment),
        ("PyMOL script", pml_path),
        ("PyMOL render script", pml_path.with_name(pml_path.stem + "_render.pml")),
    ]
    if extra:
        named.extend(("Refined coordinates" if path.name.endswith("_refined.pdb") else "Output", path) for path in extra)
    files = [{"role": role, "path": str(path)} for role, path in named if path.is_file()]
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
    }


RUNNERS = {
    "fofo": run_fofo,
    "refinement": run_refinement,
    "twofofc": run_twofofc,
}
