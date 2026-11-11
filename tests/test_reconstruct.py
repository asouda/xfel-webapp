"""Reconstruction workflow: peak placement, atom insertion, and the form."""

import math
import shutil
import tempfile
import unittest
from pathlib import Path

import gemmi
import numpy as np

from xfel_web.ccp4 import reasonable_bond_commands, reasonable_bond_specs
from xfel_web.pipelines import _viewer_scene
from xfel_web.restraints import external_harmonic_keywords
from xfel_web.forms import FormError, parse_job
from xfel_web.peaks import DensityMap
from xfel_web.reconstruction import (
    _section_directions,
    closest_atom,
    describe_atom,
    insert_atom_record,
    place_record_at,
    locate_peak,
    missing_atoms,
    read_pdb_atoms,
    run_reconstruct,
    write_cross_sections,
)


def _pdb_line(serial, name, x, y, z, element, resname="OEX", resi=416, alt=" "):
    line = [" "] * 80
    line[0:6] = list("HETATM")
    line[6:11] = list(f"{serial:5d}")
    line[12:16] = list(f"{name:>4}")
    line[16] = alt
    line[17:20] = list(f"{resname:>3}")
    line[21] = "A"
    line[22:26] = list(f"{resi:4d}")
    line[30:38] = list(f"{x:8.3f}")
    line[38:46] = list(f"{y:8.3f}")
    line[46:54] = list(f"{z:8.3f}")
    line[54:60] = list(f"{1.0:6.2f}")
    line[60:66] = list(f"{20.0:6.2f}")
    line[76:78] = list(f"{element:>2}")
    return "".join(line)


def _write_map(path: Path, grid: gemmi.FloatGrid) -> None:
    ccp4 = gemmi.Ccp4Map()
    ccp4.grid = grid
    ccp4.update_ccp4_header(2)
    ccp4.write_ccp4_map(str(path))


def _spike_map(path: Path) -> None:
    grid = gemmi.FloatGrid(16, 16, 16)
    grid.set_unit_cell(gemmi.UnitCell(32, 32, 32, 90, 90, 90))
    grid.spacegroup = gemmi.SpaceGroup("P 1")
    grid.set_value(8, 8, 8, 5.0)
    _write_map(path, grid)


def _x_ramp_map(path: Path) -> None:
    """Orthogonal map whose value at each node is that node's Cartesian x, in Å."""
    grid = gemmi.FloatGrid(16, 16, 16)
    grid.set_unit_cell(gemmi.UnitCell(32, 32, 32, 90, 90, 90))
    grid.spacegroup = gemmi.SpaceGroup("P 1")
    for i in range(16):
        for j in range(16):
            for k in range(16):
                grid.set_value(i, j, k, 2.0 * i)
    _write_map(path, grid)


class HarmonicTests(unittest.TestCase):
    def test_a_harmonic_line_holds_the_added_atom_with_the_given_sigma(self):
        text = external_harmonic_keywords(
            [{"chain": "A", "resi": 416, "name": "O5", "alt": "B", "icode": ""}],
            0.1,
        )
        self.assertEqual(text, "exte harm chain A residue 416 atom O5 altcode B sigma 0.100")
        self.assertEqual(external_harmonic_keywords([], 0.1), "")


class BondTests(unittest.TestCase):
    def test_nearby_heavy_atoms_are_bonded_and_distant_or_hydrogen_atoms_are_not(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "site.pdb"
            path.write_text(
                "\n".join(
                    [
                        _pdb_line(1, "MN1", 0, 0, 0, "MN"),
                        _pdb_line(2, "O5", 1.9, 0, 0, "O"),
                        _pdb_line(3, "MN4", 2.8, 0, 0, "MN"),
                        _pdb_line(4, "H1", 1.0, 0, 0, "H"),
                        _pdb_line(5, "O1", 0, 1.8, 0, "O", alt="A"),
                        _pdb_line(6, "O1", 0, 1.8, 0, "O", alt="B"),
                        "END",
                    ]
                )
                + "\n"
            )
            commands = reasonable_bond_commands(path, "A", 416, [])
            specs = reasonable_bond_specs(path, "A", 416, [])
        self.assertIn("bond site and name MN1, site and name O5", commands)
        self.assertFalse(any("MN4" in command and "MN1" in command for command in commands))
        self.assertFalse(any("name H1" in command for command in commands))
        self.assertNotIn("bond site and alt A and name O1, site and alt B and name O1", commands)
        self.assertIn({"a": {"chain": "A", "resi": 416, "atom": "MN1"}, "b": {"chain": "A", "resi": 416, "atom": "O5"}}, specs)
        self.assertFalse(any(spec["a"]["atom"] == "MN1" and spec["b"]["atom"] == "MN4" for spec in specs))

    def test_the_reference_scene_bonds_the_reference_coordinates(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "reference.pdb"
            path.write_text(
                "\n".join([_pdb_line(1, "MN1", 0, 0, 0, "MN"), _pdb_line(2, "O5", 1.9, 0, 0, "O"), "END"]) + "\n"
            )
            cfg = {
                "positive_levels": [3.5],
                "negative_levels": None,
                "viewer_sigma_min": -5,
                "viewer_sigma_max": 5,
                "contour_mode": "sigma",
                "surface": "mesh",
                "site_chain": "A",
                "site_resi": 416,
                "extra_residues": [],
                "geometry_bonds": True,
                "viewer_light_label": "Reference site",
            }
            scene = _viewer_scene(
                model_pdb=path,
                fragment=path,
                cfg=cfg,
                positive_colors=["blue"],
                negative_colors=["red"],
                light_pdb=path,
            )
            self.assertEqual(scene["light_label"], "Reference site")
            self.assertEqual(len(scene["light_bonds"]), 1)
            plain = dict(cfg)
            plain["geometry_bonds"] = False
            other = _viewer_scene(
                model_pdb=path,
                fragment=path,
                cfg=plain,
                positive_colors=["blue"],
                negative_colors=["red"],
                light_pdb=path,
            )
            self.assertNotIn("light_bonds", other)


class AtomTests(unittest.TestCase):
    def test_missing_atoms_skip_hydrogens_and_atoms_already_present(self):
        reference = read_pdb_atoms_from(
            [
                _pdb_line(1, "MN1", 0, 0, 0, "MN"),
                _pdb_line(2, "O5", 1, 0, 0, "O"),
                _pdb_line(3, "H1", 1.2, 0, 0, "H"),
            ]
        )
        omit = read_pdb_atoms_from([_pdb_line(1, "MN1", 0, 0, 0, "MN")])
        missing = missing_atoms(reference, omit, set())
        self.assertEqual([atom.name for atom in missing], ["O5"])
        atom, distance = closest_atom(missing, (1.2, 0, 0))
        self.assertEqual(atom.name, "O5")
        self.assertAlmostEqual(distance, 0.2, places=3)
        self.assertIn("O5 OEX chain A residue 416", describe_atom(atom))

    def test_inserted_record_keeps_its_coordinates_and_gets_a_new_serial(self):
        original = _pdb_line(1, "MN1", 0, 0, 0, "MN") + "\nEND\n"
        added = insert_atom_record(original, _pdb_line(7, "O5", 1.5, 2.5, 3.5, "O"))
        lines = [line for line in added.splitlines() if line.startswith("HETATM")]
        self.assertEqual(len(lines), 2)
        self.assertEqual(lines[1][6:11], "    2")
        self.assertEqual(lines[1][12:16], "  O5")
        self.assertIn("   1.500", lines[1])
        self.assertTrue(added.strip().endswith("END"))

    def test_a_placed_record_uses_the_peak_coordinates(self):
        source = _pdb_line(7, "O5", 1.5, 2.5, 3.5, "O")
        placed = place_record_at(source, (8.25, -0.5, 12))
        self.assertEqual(placed[12:30], source[12:30])
        self.assertEqual(placed[30:54], "   8.250  -0.500  12.000")
        self.assertEqual(placed[54:].rstrip(), source[54:].rstrip())


def read_pdb_atoms_from(lines: list[str]):
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "model.pdb"
        path.write_text("\n".join(lines) + "\nEND\n")
        return read_pdb_atoms(path)


class PeakTests(unittest.TestCase):
    def test_spike_is_the_selected_voxel_and_the_section_passes_through_it(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "spike.ccp4"
            _spike_map(path)
            peak = locate_peak(path, 3, False, 3)
            self.assertIsNotNone(peak)
            self.assertAlmostEqual(peak["density"], 5.0, places=4)
            sections = Path(tmp) / "density_sections.txt"
            write_cross_sections(peak["map"], peak["xyz"], sections, header=["# test"])
            text = sections.read_text()
            self.assertIn("density_xy_p", text)
            self.assertIn("density_xyz_ppp", text)
            rows = [line.split() for line in text.splitlines() if line[:1].isdigit() or line.startswith(" ")]
            center = [row for row in rows if abs(float(row[0])) < 1e-6][0]
            self.assertEqual(len(center), 14)
            for value in center[1:]:
                self.assertAlmostEqual(float(value), 5.0, places=4)
            refined = locate_peak(path, 3, True, 3)
            self.assertAlmostEqual(refined["density"], 5.0, places=4)
            self.assertTrue(np.all(np.isfinite(refined["xyz"])))

    def test_diagonal_offsets_are_distances_in_angstroms(self):
        directions = _section_directions()
        self.assertEqual(
            [name for name, _vector in directions],
            [
                "density_x",
                "density_y",
                "density_z",
                "density_xy_p",
                "density_xy_m",
                "density_xz_p",
                "density_xz_m",
                "density_yz_p",
                "density_yz_m",
                "density_xyz_ppp",
                "density_xyz_ppm",
                "density_xyz_pmp",
                "density_xyz_pmm",
            ],
        )
        for _name, vector in directions:
            self.assertAlmostEqual(float(np.linalg.norm(vector)), 1.0, places=8)

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "ramp.ccp4"
            _x_ramp_map(path)
            density_map = DensityMap(path)
            origin = density_map.to_cart((8, 8, 8))
            sections = Path(tmp) / "density_sections.txt"
            write_cross_sections(density_map, origin, sections, header=["# ramp"])
            rows = [
                line.split()
                for line in sections.read_text().splitlines()
                if line[:1].isdigit() or line.startswith(" ") or line.startswith("-")
            ]
            row = [item for item in rows if abs(float(item[0]) - 1.0) < 1e-6][0]
            by_name = {name: float(row[index + 1]) for index, (name, _vector) in enumerate(directions)}
            self.assertAlmostEqual(by_name["density_x"], 17.0, places=5)
            self.assertAlmostEqual(by_name["density_y"], 16.0, places=5)
            self.assertAlmostEqual(by_name["density_xy_p"], 16.0 + 1.0 / math.sqrt(2.0), places=5)
            self.assertAlmostEqual(by_name["density_xy_m"], 16.0 + 1.0 / math.sqrt(2.0), places=5)
            self.assertAlmostEqual(by_name["density_yz_p"], 16.0, places=5)
            self.assertAlmostEqual(by_name["density_xyz_ppp"], 16.0 + 1.0 / math.sqrt(3.0), places=5)
            self.assertAlmostEqual(by_name["density_xyz_pmm"], 16.0 + 1.0 / math.sqrt(3.0), places=5)


class LoopTests(unittest.TestCase):
    def _run(self, root: Path, xyz, extra_atoms, n_atoms=1, n_check=0, pause=None, report_step=None, harmonic=False):
        source = root / "spike.ccp4"
        _spike_map(source)
        omit = root / "omit.pdb"
        reference = root / "reference.pdb"
        keeper = _pdb_line(1, "MN1", 0, 0, 0, "MN")
        omit.write_text(keeper + "\nEND\n")
        reference_lines = [keeper, _pdb_line(2, "O5", *xyz, "O"), *extra_atoms, "END"]
        reference.write_text("\n".join(reference_lines) + "\n")
        calls = []

        def omit_runner(cfg, emit):
            emit("omit map")
            dest = Path(cfg["output_dir"]) / "fragment.ccp4"
            shutil.copy2(source, dest)
            calls.append(Path(cfg["files"]["dark_pdb"]).read_text())
            self.last_harmonics.append(cfg.get("harmonic_atoms"))
            return {
                "caption": "omit caption",
                "figure": str(dest),
                "final_map": str(dest),
                "fragment": str(dest),
                "pymol_script": str(dest),
                "output_dir": str(cfg["output_dir"]),
                "files": [],
                "viewer": {"model": str(cfg["files"]["dark_pdb"]), "map": str(dest)},
            }

        logs = []
        self.last_harmonics = []
        result = run_reconstruct(
            {
                "output_dir": root / "out",
                "dark_id": "omit",
                "light_id": "ref",
                "n_atoms": n_atoms,
                "n_check_map": n_check,
                "neighborhood_size": 3,
                "spline_seed": False,
                "spline_order": 3,
                "harmonic_added": harmonic,
                "harmonic_sigma": 0.1,
                "refmac_ncyc": 0,
                "files": {"dark_pdb": omit, "light_pdb": reference, "light_sf": reference},
            },
            logs.append,
            pause=pause,
            omit_runner=omit_runner,
            report_step=report_step,
        )
        return result, calls, logs

    def test_a_peak_adds_the_reference_atom_and_a_final_map_is_calculated(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "spike.ccp4"
            _spike_map(source)
            peak = locate_peak(source, 3, False, 3)
            steps = []
            result, calls, logs = self._run(root, peak["xyz"], [], report_step=steps.append)
            self.assertEqual(steps, [1, 2])
            self.assertEqual(len(calls), 2)
            self.assertIn("O5", calls[1])
            self.assertNotIn("O5", calls[0])
            self.assertEqual(result["view_id"], "final-2")
            self.assertFalse(result["checkpoint"])
            self.assertIn("O5", (root / "out" / "omit_reconstructed.pdb").read_text())
            self.assertTrue((root / "out" / "r-iteration-01" / "density_sections.txt").is_file())
            self.assertTrue(any("Added" in line and "O5" in line for line in logs))
            self.assertIn("omit caption", result["caption"])

    def test_an_atom_farther_than_3_angstroms_stops_the_search(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "spike.ccp4"
            _spike_map(source)
            peak = locate_peak(source, 3, False, 3)
            far = (peak["xyz"][0] + 4.0, peak["xyz"][1], peak["xyz"][2])
            result, calls, logs = self._run(root, far, [])
            self.assertEqual(len(calls), 1)
            self.assertNotIn("O5", (root / "out" / "omit_reconstructed.pdb").read_text())
            self.assertTrue(any("farther than 3" in line for line in logs))
            self.assertIn("not added", result["caption"])

    def test_an_atom_farther_than_1_angstrom_is_added_with_a_warning(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "spike.ccp4"
            _spike_map(source)
            peak = locate_peak(source, 3, False, 3)
            near = (peak["xyz"][0] + 1.5, peak["xyz"][1], peak["xyz"][2])
            _result, calls, logs = self._run(root, near, [])
            self.assertEqual(len(calls), 2)
            self.assertTrue(any(line.startswith("Warning:") and "1.50" in line for line in logs))
            rebuilt = read_pdb_atoms(root / "out" / "omit_reconstructed.pdb")
            added = [atom for atom in rebuilt if atom.name == "O5"]
            self.assertEqual(len(added), 1)
            self.assertAlmostEqual(added[0].xyz[0], peak["xyz"][0], places=3)
            self.assertAlmostEqual(added[0].xyz[1], peak["xyz"][1], places=3)
            self.assertAlmostEqual(added[0].xyz[2], peak["xyz"][2], places=3)
            self.assertGreater(abs(added[0].xyz[0] - near[0]), 1.0)

    def test_a_map_check_pauses_before_the_atom_is_added(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "spike.ccp4"
            _spike_map(source)
            peak = locate_peak(source, 3, False, 3)
            pauses = []

            def pause(partial):
                pauses.append(partial)
                self.assertNotIn("O5", (root / "out" / "omit_working.pdb").read_text())

            _result, calls, _logs = self._run(root, peak["xyz"], [], n_check=1, pause=pause)
            self.assertEqual([item["view_id"] for item in pauses], ["checkpoint-1"])
            self.assertTrue(pauses[0]["checkpoint"])
            self.assertEqual(len(calls), 2)
            self.assertIn("O5", (root / "out" / "omit_reconstructed.pdb").read_text())

    def test_atoms_added_earlier_are_sent_as_harmonic_restraints(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "spike.ccp4"
            _spike_map(source)
            peak = locate_peak(source, 3, False, 3)
            self._run(root, peak["xyz"], [], harmonic=True)
            self.assertIsNone(self.last_harmonics[0])
            self.assertEqual(self.last_harmonics[1][0]["name"], "O5")
            self.assertEqual(self.last_harmonics[1][0]["resi"], 416)
            self.assertEqual(self.last_harmonics[1][0]["chain"], "A")


class ReconstructFormTests(unittest.TestCase):
    def _form(self, root: Path, **overrides):
        pdb = root / "model.pdb"
        reference = root / "reference.pdb"
        sf = root / "obs.cif"
        pdb.write_text("END\n")
        reference.write_text("END\n")
        sf.write_text("data_test\n")
        form = {
            "output_dir": str(root / "job-out"),
            "ccp4_root": "/Applications/ccp4-9",
            "resolution_source": "explicit",
            "lowres": "30",
            "dmin": "1.89",
            "map_border": "5",
            "site_chain": "A",
            "site_resi": "416",
            "site_label": "OEC",
            "extra_residues": "",
            "extra_copies": "",
            "export_radius": "10",
            "carve": "3",
            "sigma_levels": "3.5",
            "contour_mode": "sigma",
            "surface": "mesh",
            "viewer_sigma_min": "-5",
            "viewer_sigma_max": "5",
            "png_width": "800",
            "png_height": "600",
            "dark_id": "omit",
            "light_id": "obs",
            "refmac_ncyc": "0",
            "distance_sigma": "0.5",
            "angle_sigma": "25",
            "distance_restraints": "",
            "angle_restraints": "",
            "n_atoms": "2",
            "n_check_map": "0",
            "neighborhood_size": "3",
            "spline_order": "3",
            "path_dark_pdb": str(pdb),
            "path_light_pdb": str(reference),
            "path_light_sf": str(sf),
        }
        form.update(overrides)
        return form

    def test_defaults_match_the_reconstruction_form(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cfg = parse_job("reconstruct", self._form(root), {}, root)
        self.assertEqual(cfg["n_atoms"], 2)
        self.assertEqual(cfg["n_check_map"], 0)
        self.assertEqual(cfg["neighborhood_size"], 3)
        self.assertFalse(cfg["spline_seed"])
        self.assertEqual(cfg["spline_order"], 3)
        self.assertEqual(cfg["distance_restraints"], [])
        self.assertEqual(cfg["angle_restraints"], [])
        self.assertIsNone(cfg["custom_dictionary"])
        self.assertTrue(cfg["geometry_bonds"])
        self.assertFalse(cfg["show_oec_bonds"])
        self.assertTrue(cfg["ray_trace"])
        self.assertEqual(cfg["viewer_light_label"], "Reference site")
        self.assertFalse(cfg["harmonic_added"])
        self.assertEqual(cfg["harmonic_sigma"], 0.1)

    def test_a_non_positive_harmonic_sigma_is_rejected_when_the_option_is_on(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaises(FormError):
                parse_job("reconstruct", self._form(root, harmonic_added="on", harmonic_sigma="0"), {}, root)
            cfg = parse_job("reconstruct", self._form(root, harmonic_added="on", harmonic_sigma="0.1"), {}, root)
        self.assertTrue(cfg["harmonic_added"])
        self.assertEqual(cfg["harmonic_sigma"], 0.1)

    def test_even_neighborhood_and_a_non_pdb_model_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaises(FormError):
                parse_job("reconstruct", self._form(root, neighborhood_size="4"), {}, root)
            cif = root / "model.cif"
            cif.write_text("data_x\n")
            with self.assertRaises(FormError):
                parse_job("reconstruct", self._form(root, path_dark_pdb=str(cif)), {}, root)

    def test_the_omit_map_form_does_not_gain_reconstruction_fields(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            form = self._form(root)
            form["path_light_pdb"] = form["path_dark_pdb"]
            cfg = parse_job("refinement", form, {}, root)
        self.assertNotIn("n_atoms", cfg)
        self.assertNotIn("geometry_bonds", cfg)
        self.assertNotIn("n_check_map", cfg)
        self.assertNotIn("harmonic_added", cfg)
