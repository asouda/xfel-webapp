import gzip
import struct
import tempfile
import unittest
from io import BytesIO
from pathlib import Path

from werkzeug.datastructures import FileStorage

from xfel_web.ccp4 import (
    ccp4_map_center,
    choose_resolution,
    copy_ccp4_header_statistics,
    delta_b_tag,
    fft_grid_resolution,
    grid_res_tag,
    pad_map_around_model,
    pymol_scene_links,
    residue_conformers,
    selection_coordinates,
    view_frame_half_widths,
)
from xfel_web.forms import (
    FormError,
    decompress_gzip,
    parse_angle_restraints,
    parse_copies,
    parse_distance_restraints,
    parse_extra_residues,
    parse_job,
    parse_viewer_sigma,
    resolve_input_file,
)
from xfel_web.pipelines import (
    _viewer_scene,
    amplitude_cad_keywords,
    initial_viewer_sigma,
    intensity_truncate_args,
)
from xfel_web.restraints import external_angle_keywords, external_distance_keywords


def _pdb_line(serial, name, x, y, z, element):
    line = [" "] * 80
    line[0:6] = list("HETATM")
    line[6:11] = list(f"{serial:5d}")
    line[12:16] = list(f"{name:>4}")
    line[17:20] = list("OEY")
    line[21] = "A"
    line[22:26] = list(f"{416:4d}")
    line[30:38] = list(f"{x:8.3f}")
    line[38:46] = list(f"{y:8.3f}")
    line[46:54] = list(f"{z:8.3f}")
    line[54:60] = list(f"{1.0:6.2f}")
    line[60:66] = list(f"{20.0:6.2f}")
    line[76:78] = list(f"{element:>2}")
    return "".join(line) + "\n"


class ParseTests(unittest.TestCase):
    def test_resolution_intersection(self):
        low, high = choose_resolution("sf", 99, 1, [(30.0, 2.0), (20.0, 2.1)])
        self.assertEqual((low, high), (20.0, 2.1))

    def test_explicit_resolution_rejects_inverted_shell(self):
        with self.assertRaises(ValueError):
            choose_resolution("explicit", 1.5, 2.0, [])

    def test_fft_grid_extends_only_to_a_finer_limit(self):
        self.assertEqual(fft_grid_resolution(30, 1.89, None), 1.89)
        self.assertEqual(fft_grid_resolution(30, 1.89, 1.2), 1.2)
        self.assertEqual(fft_grid_resolution(30, 1.89, 1.89), 1.89)
        self.assertEqual(grid_res_tag(30, 1.89, 1.89), "30p00-1p89")
        self.assertEqual(grid_res_tag(30, 1.89, 1.2), "30p00-1p89_grid1p20")
        with self.assertRaises(ValueError):
            fft_grid_resolution(30, 1.89, 2.5)
        with self.assertRaises(ValueError):
            fft_grid_resolution(30, 1.89, 40)

    def test_mapmask_keeps_asymmetric_unit_sigma(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "asu.ccp4"
            dest = root / "padded.ccp4"
            header = bytearray(1024)
            header[19 * 4 : 22 * 4] = b"aminamaxmean"
            header[54 * 4 : 55 * 4] = b"rms!"
            source.write_bytes(header)
            dest.write_bytes(b"\x11" * 1024)
            copy_ccp4_header_statistics(source, dest)
            written = dest.read_bytes()
        self.assertEqual(written[19 * 4 : 22 * 4], b"aminamaxmean")
        self.assertEqual(written[54 * 4 : 55 * 4], b"rms!")
        self.assertEqual(written[0:4], b"\x11" * 4)

    def test_map_center_matches_the_viewer_grid_box(self):
        def write_map(path, *, ncstart, nrstart, nsstart, mapc, mapr, maps, origin):
            raw = bytearray(1024)
            ints = {
                0: 10, 1: 10, 2: 10, 3: 2,
                4: ncstart, 5: nrstart, 6: nsstart,
                7: 100, 8: 100, 9: 100,
                16: mapc, 17: mapr, 18: maps,
            }
            for index, value in ints.items():
                struct.pack_into("<i", raw, index * 4, value)
            for index, value in enumerate((100.0, 200.0, 300.0, 90.0, 90.0, 90.0), start=10):
                struct.pack_into("<f", raw, index * 4, value)
            for index, value in enumerate(origin, start=49):
                struct.pack_into("<f", raw, index * 4, value)
            path.write_bytes(raw)

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            straight = root / "straight.ccp4"
            write_map(straight, ncstart=0, nrstart=20, nsstart=0, mapc=1, mapr=2, maps=3, origin=(0.0, 0.0, 0.0))
            straight_center = ccp4_map_center(straight)
            self.assertAlmostEqual(straight_center[0], 4.5)
            self.assertAlmostEqual(straight_center[1], 49.0)
            self.assertAlmostEqual(straight_center[2], 13.5)
            swapped = root / "swapped.ccp4"
            write_map(
                swapped,
                ncstart=0, nrstart=20, nsstart=0,
                mapc=2, mapr=3, maps=1,
                origin=(1.5, -2.0, 0.25),
            )
            center = ccp4_map_center(swapped)
        self.assertAlmostEqual(center[0], 4.5 + 1.5)
        self.assertAlmostEqual(center[1], 4.5 * 2 - 2.0)
        self.assertAlmostEqual(center[2], 24.5 * 3 + 0.25)

    def test_figure_frame_stays_centered_on_the_map(self):
        self.assertEqual(view_frame_half_widths([], (1.0, 2.0, 3.0), 3.0), (3.0, 3.0, 3.0))
        self.assertEqual(
            view_frame_half_widths([(0.0, 5.0, 1.0), (4.0, 0.0, 1.0)], (1.0, 2.0, 3.0), 3.0),
            (3.0, 3.0, 2.0),
        )
        with tempfile.TemporaryDirectory() as tmp:
            pdb = Path(tmp) / "site.pdb"
            pdb.write_text(
                "\n".join([
                    "HETATM    1  C   OEC A 416       0.000   0.000   0.000  1.00 20.00           C",
                    "HETATM    2  N   OEC A 416       2.000   0.000   0.000  1.00 20.00           N",
                    "HETATM    3  O   GLU A 189       9.000   1.000   0.000  1.00 20.00           O",
                    "HETATM    4  C   OEC B 416       8.000   8.000   8.000  1.00 20.00           C",
                    "END",
                    "",
                ])
            )
            coords = selection_coordinates(pdb, "A", 416, [{"chain": "A", "resi": 189, "label": "GLU"}])
        self.assertEqual(coords, [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (9.0, 1.0, 0.0)])

    def test_delta_b_tag(self):
        self.assertEqual(delta_b_tag(-40), "dB-40")
        self.assertEqual(delta_b_tag(-12.5), "dB-12p5")

    def test_pymol_script_links(self):
        path = Path(__file__).with_name("_scene.pml")
        path.write_text(
            "\n".join([
                'load "/tmp/out/extended_OEC_A416_10A.map", DIFF',
                'load "/tmp/out/model.pdb", xyz',
                "select site, xyz and chain A and resi 416",
                "select GLU189, xyz and chain A and resi 189",
                "bond site and name MN1, site and name O1",
                "distance d_mn1_o5, site and name MN1, site and name O5",
                "distance d_mn1_GLU189_B, site and name MN1, GLU189 and alt B and name OE2",
                "isosurface map_2, TWOFOFC, 2.000000, carve_sel, carve=3.00",
            ]) + "\n"
        )
        try:
            links = pymol_scene_links(path)
        finally:
            path.unlink(missing_ok=True)
        self.assertEqual(links["map"], "/tmp/out/extended_OEC_A416_10A.map")
        self.assertEqual(links["bonds"], [{
            "a": {"chain": "A", "resi": 416, "atom": "MN1"},
            "b": {"chain": "A", "resi": 416, "atom": "O1"},
        }])
        self.assertEqual(links["distances"][1]["b"]["alt"], "B")
        self.assertEqual(links["distances"][1]["b"]["atom"], "OE2")
        self.assertEqual(links["carve"], 3.0)

    def test_residue_lines(self):
        extras = parse_extra_residues("A 189 GLU189\n# comment\n")
        self.assertEqual(extras, [{"chain": "A", "resi": 189, "label": "GLU189"}])
        self.assertEqual(parse_copies("A", 416, "D 416\nA 416"), [("A", 416), ("D", 416)])
        self.assertEqual(parse_distance_restraints("MN1 O5\nMN4 O5 2.4 0.5"), [("MN1", "O5"), ("MN4", "O5", 2.4, 0.5)])
        self.assertEqual(parse_angle_restraints("MN3 O5 MN4 109.5"), [("MN3", "O5", "MN4", 109.5)])
        with self.assertRaises(FormError):
            parse_extra_residues("AB 189")

    def test_viewer_sigma_range(self):
        self.assertEqual(parse_viewer_sigma({"viewer_sigma_min": "-5", "viewer_sigma_max": "5"}), (-5.0, 5.0))
        self.assertEqual(initial_viewer_sigma(-5, 5, [3, 4, 5], [3, 4, 5]), 3)
        self.assertEqual(initial_viewer_sigma(1, 5, [2, 3.5]), 2)
        self.assertEqual(initial_viewer_sigma(4.5, 8, [3, 4, 5]), 5)
        self.assertEqual(initial_viewer_sigma(-4, -1, [3, 4, 5]), -2.5)
        with self.assertRaises(FormError):
            parse_viewer_sigma({"viewer_sigma_min": "5", "viewer_sigma_max": "-5"})
        with self.assertRaises(FormError):
            parse_viewer_sigma({"viewer_sigma_min": "2", "viewer_sigma_max": "2"})

    def test_viewer_scene_uses_sigma_range(self):
        scene = _viewer_scene(
            model_pdb=Path("model.pdb"),
            fragment=Path("fragment.map"),
            cfg={
                "site_chain": "A",
                "site_resi": 416,
                "extra_residues": [],
                "positive_levels": [3, 4, 5],
                "negative_levels": [3, 4, 5],
                "contour_mode": "sigma",
                "surface": "mesh",
                "viewer_sigma_min": -5,
                "viewer_sigma_max": 5,
            },
            positive_colors=["blue"],
            negative_colors=["red"],
        )
        self.assertEqual(scene["sigma_min"], -5)
        self.assertEqual(scene["sigma_max"], 5)
        self.assertEqual(scene["sigma"], 3)
        self.assertEqual(scene["positive_levels"], [3, 4, 5])
        self.assertNotIn("light_model", scene)
        light = _viewer_scene(
            model_pdb=Path("model.pdb"),
            fragment=Path("fragment.map"),
            cfg={
                "site_chain": "A",
                "site_resi": 416,
                "extra_residues": [],
                "positive_levels": [3.5],
                "negative_levels": [3.5],
                "contour_mode": "sigma",
                "surface": "mesh",
                "viewer_sigma_min": -5,
                "viewer_sigma_max": 5,
            },
            positive_colors=["blue"],
            negative_colors=["red"],
            light_pdb=Path("/data/6w1v.cif"),
        )
        self.assertEqual(light["light_model"], "/data/6w1v.cif")
        self.assertEqual(light["light_model_ext"], "cif")

    def test_gzip_inputs_are_decompressed_with_gunzip(self):
        plain = b"data_test\n_entry.id 9z35\n"
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "9z35-sf.cif.gz"
            with gzip.open(source, "wb") as handle:
                handle.write(plain)
            plain_pdb = root / "model.pdb"
            plain_pdb.write_bytes(b"END\n")
            dest = decompress_gzip(source, root / "out", "structure factors", dest_name=source.name)
            self.assertEqual(dest.name, "9z35-sf.cif")
            self.assertEqual(dest.read_bytes(), plain)
            self.assertTrue(source.is_file())
            owned = root / "out" / "owned.cif.gz"
            owned.parent.mkdir(exist_ok=True)
            with gzip.open(owned, "wb") as handle:
                handle.write(plain)
            expanded = decompress_gzip(owned, owned.parent, "structure factors")
            self.assertEqual(expanded.read_bytes(), plain)
            self.assertFalse(owned.exists())
            kept = decompress_gzip(plain_pdb, root / "out", "coordinates")
            self.assertEqual(kept, plain_pdb)
            output = root / "job"
            output.mkdir()
            plain_cif = root / "9z35-sf.cif"
            plain_cif.write_bytes(plain)
            self.assertEqual(
                resolve_input_file({"path_sf": str(plain_cif)}, {}, "sf", "structure factors", output, root),
                plain_cif.resolve(),
            )
            self.assertEqual(
                resolve_input_file({"path_pdb": str(plain_pdb)}, {}, "pdb", "coordinates", output, root),
                plain_pdb.resolve(),
            )
            gzip_named_pdb = root / "still.pdb"
            with gzip.open(gzip_named_pdb, "wb") as handle:
                handle.write(b"END\n")
            compressed_bytes = gzip_named_pdb.read_bytes()
            self.assertEqual(
                resolve_input_file({"path_pdb": str(gzip_named_pdb)}, {}, "pdb", "coordinates", output, root),
                gzip_named_pdb.resolve(),
            )
            self.assertEqual(gzip_named_pdb.read_bytes(), compressed_bytes)
            upload = FileStorage(stream=BytesIO(b"END\n"), filename="model.pdb")
            saved = resolve_input_file({}, {"file_pdb": upload}, "pdb", "coordinates", output, root)
            self.assertEqual(saved.read_bytes(), b"END\n")
            self.assertTrue(saved.name.endswith(".pdb"))
            self.assertFalse(saved.name.lower().endswith(".gz"))
            resolved = resolve_input_file(
                {"path_sf": str(source)},
                {},
                "sf",
                "structure factors",
                output,
                root,
            )
            self.assertEqual(resolved.read_bytes(), plain)
            self.assertEqual(resolved.parent, output / "inputs")
            broken = root / "broken.pdb.gz"
            broken.write_bytes(b"not gzip")
            with self.assertRaises(FormError):
                decompress_gzip(broken, root / "out", "coordinates")

    def test_external_restraints_from_coordinates(self):
        path = Path(__file__).with_name("_tiny.pdb")
        path.write_text(
            _pdb_line(1, "MN1", 0, 0, 0, "MN")
            + _pdb_line(2, "O5", 1, 0, 0, "O")
            + _pdb_line(3, "MN4", 1, 1, 0, "MN")
            + "END\n"
        )
        try:
            text = external_distance_keywords(path, [("A", 416)], [("MN1", "O5")], 0.5)
            self.assertIn("exte dist first chain A resi 416 atom MN1", text)
            self.assertIn("value 1.000 sigma 0.500", text)
            angle = external_angle_keywords(path, [("A", 416)], [("MN1", "O5", "MN4")], 25)
            self.assertIn("exte angle first chain A resi 416 atom MN1", angle)
            self.assertIn("value 90.000 sigma 25.000", angle)
            conformers = residue_conformers(path, "A", 416)
            self.assertIn("MN1", conformers[""])
        finally:
            path.unlink(missing_ok=True)


class RefinementFormTests(unittest.TestCase):
    def _form(self, root: Path, **overrides):
        pdb = root / "model.pdb"
        sf = root / "obs.cif"
        pdb.write_text("END\n")
        sf.write_text("data_test\n")
        output = root / "job-out"
        form = {
            "output_dir": str(output),
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
            "path_dark_pdb": str(pdb),
            "path_light_pdb": str(pdb),
            "path_light_sf": str(sf),
        }
        form.update(overrides)
        return form

    def test_blank_grid_resolution_leaves_the_data_cutoff(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cfg = parse_job("refinement", self._form(root), {}, root)
        self.assertIsNone(cfg["fft_dmin"])
        self.assertEqual(cfg["map_center"], "residue")
        self.assertIsNone(cfg["map_center_xyz"])

    def test_cartesian_map_center_is_kept(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cfg = parse_job(
                "refinement",
                self._form(root, map_center="cartesian", center_x="1.5", center_y="-2", center_z="3.25"),
                {},
                root,
            )
        self.assertEqual(cfg["map_center"], "cartesian")
        self.assertEqual(cfg["map_center_xyz"], (1.5, -2.0, 3.25))

    def test_cartesian_map_center_requires_coordinates(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaises(FormError):
                parse_job("refinement", self._form(root, map_center="cartesian"), {}, root)

    def test_finer_grid_resolution_is_kept(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cfg = parse_job("refinement", self._form(root, fft_dmin="1.20"), {}, root)
        self.assertEqual(cfg["fft_dmin"], 1.2)

    def test_custom_dictionary_is_optional(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cfg = parse_job("refinement", self._form(root), {}, root)
        self.assertIsNone(cfg["custom_dictionary"])

    def test_custom_dictionary_path_is_kept_for_any_cycle_count(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            library = root / "ligand.cif"
            library.write_text("data_comp_list\n")
            cfg = parse_job(
                "refinement",
                self._form(root, refmac_ncyc="0", path_custom_dictionary=str(library)),
                {},
                root,
            )
        self.assertEqual(cfg["custom_dictionary"], library.resolve())
        self.assertEqual(cfg["refmac_ncyc"], 0)

    def test_missing_custom_dictionary_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaises(FormError) as caught:
                parse_job(
                    "refinement",
                    self._form(root, path_custom_dictionary=str(root / "missing.cif")),
                    {},
                    root,
                )
        self.assertIn("custom dictionary", str(caught.exception))

    def test_coarser_grid_resolution_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaises(FormError) as caught:
                parse_job("refinement", self._form(root, fft_dmin="2.50"), {}, root)
        self.assertIn("coarser", str(caught.exception))


class TruncateColumnTests(unittest.TestCase):
    def test_mean_intensities_use_colin(self):
        args = intensity_truncate_args(["H", "K", "L", "FREE", "I", "SIGI"])
        self.assertEqual(args, "-colin '/*/*/[I,SIGI]'")

    def test_friedel_intensities_use_colano(self):
        args = intensity_truncate_args(
            ["H", "K", "L", "FREE", "I(+)", "SIGI(+)", "I(-)", "SIGI(-)"]
        )
        self.assertEqual(args, "-colano '/*/*/[I(+),SIGI(+),I(-),SIGI(-)]'")

    def test_mean_intensities_win_when_both_layouts_are_present(self):
        args = intensity_truncate_args(
            ["I", "SIGI", "I(+)", "SIGI(+)", "I(-)", "SIGI(-)"]
        )
        self.assertEqual(args, "-colin '/*/*/[I,SIGI]'")

    def test_incomplete_columns_are_rejected(self):
        with self.assertRaises(ValueError) as caught:
            intensity_truncate_args(["H", "K", "L", "F(+)", "SIGF(+)"])
        self.assertIn("F(+)", str(caught.exception))

    def test_mean_amplitudes_keep_every_truncate_column(self):
        amplitude, sigma, keywords = amplitude_cad_keywords(["F", "SIGF", "I", "SIGI"])
        self.assertEqual((amplitude, sigma), ("F", "SIGF"))
        self.assertEqual(keywords, "LABIN FILE 1 ALL\nLABIN FILE 2 E1=FREE\nEND")

    def test_friedel_amplitudes_are_renamed_to_mean_columns(self):
        amplitude, sigma, keywords = amplitude_cad_keywords(
            ["FMEAN", "SIGFMEAN", "F(+)", "SIGF(+)", "F(-)", "SIGF(-)"]
        )
        self.assertEqual((amplitude, sigma), ("FMEAN", "SIGFMEAN"))
        self.assertEqual(
            keywords,
            "LABIN FILE 1 E1=FMEAN E2=SIGFMEAN\n"
            "LABOUT FILE 1 E1=F E2=SIGF\n"
            "LABIN FILE 2 E1=FREE\n"
            "END",
        )


def _write_map(path, values, *, nxyz, starts, axes, cell, symops, arms):
    """values is indexed [section, row, column] and matches axes (mapc, mapr, maps)."""
    import numpy as np

    ns, nr, nc = values.shape
    header = bytearray(1024)
    struct.pack_into("<10i", header, 0, nc, nr, ns, 2, starts[0], starts[1], starts[2], *nxyz)
    struct.pack_into("<6f", header, 10 * 4, *cell)
    struct.pack_into("<3i", header, 16 * 4, *axes)
    sym = b"".join(op.encode("ascii").ljust(80) for op in symops)
    struct.pack_into("<i", header, 23 * 4, len(sym))
    header[52 * 4 : 52 * 4 + 4] = b"MAP "
    header[53 * 4 : 53 * 4 + 4] = b"DA\x00\x00"
    struct.pack_into("<f", header, 54 * 4, arms)
    path.write_bytes(bytes(header) + sym + np.asarray(values, dtype="<f4").tobytes())


class PadMapTests(unittest.TestCase):
    def test_full_cell_wraps_and_keeps_fft_sigma(self):
        import numpy as np

        # Fast Y, medium X, slow Z, the axis order fft uses for these maps.
        values = np.zeros((4, 4, 4), dtype=np.float32)
        for z in range(4):
            for x in range(4):
                for y in range(4):
                    values[z, x, y] = x + 10 * y + 100 * z
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "cell.ccp4"
            dest = root / "padded.ccp4"
            pdb = root / "model.pdb"
            _write_map(
                source,
                values,
                nxyz=(4, 4, 4),
                starts=(0, 0, 0),
                axes=(2, 1, 3),
                cell=(10, 10, 10, 90, 90, 90),
                symops=("X, Y, Z",),
                arms=0.07091,
            )
            pdb.write_text(_pdb_line(1, "CA", -2.5, 2.5, 22.5, "C"))
            pad_map_around_model(source, dest, pdb, 0.0)
            raw = dest.read_bytes()
            nc, nr, ns = struct.unpack_from("<3i", raw, 0)
            ncstart, nrstart, nsstart = struct.unpack_from("<3i", raw, 4 * 4)
            arms = struct.unpack_from("<f", raw, 54 * 4)[0]
            self.assertEqual((nc, nr, ns), (1, 1, 1))
            self.assertEqual((ncstart, nrstart, nsstart), (1, -1, 9))
            self.assertAlmostEqual(arms, 0.07091, places=5)
            nsymbt = struct.unpack_from("<i", raw, 23 * 4)[0]
            value = struct.unpack_from("<f", raw, 1024 + nsymbt)[0]
            self.assertEqual(value, 113.0)

    def test_partial_cell_uses_symmetry(self):
        import numpy as np

        values = np.zeros((4, 4, 2), dtype=np.float32)
        for z in range(4):
            for y in range(4):
                for x in range(2):
                    values[z, y, x] = x + 10 * y + 100 * z
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "asu.ccp4"
            dest = root / "padded.ccp4"
            pdb = root / "model.pdb"
            _write_map(
                source,
                values,
                nxyz=(4, 4, 4),
                starts=(0, 0, 0),
                axes=(1, 2, 3),
                cell=(8, 8, 8, 90, 90, 90),
                symops=("X, Y, Z", "-X, Y, Z"),
                arms=1.25,
            )
            pdb.write_text(_pdb_line(1, "CA", 6.0, 2.0, 2.0, "C"))
            pad_map_around_model(source, dest, pdb, 0.0)
            raw = dest.read_bytes()
            nsymbt = struct.unpack_from("<i", raw, 23 * 4)[0]
            value = struct.unpack_from("<f", raw, 1024 + nsymbt)[0]
            self.assertEqual(value, 111.0)
            self.assertAlmostEqual(struct.unpack_from("<f", raw, 54 * 4)[0], 1.25, places=5)


if __name__ == "__main__":
    unittest.main()
