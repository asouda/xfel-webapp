import gzip
import tempfile
import unittest
from io import BytesIO
from pathlib import Path

from werkzeug.datastructures import FileStorage

from xfel_web.ccp4 import choose_resolution, delta_b_tag, pymol_scene_links, residue_conformers
from xfel_web.forms import (
    FormError,
    decompress_gzip,
    parse_angle_restraints,
    parse_copies,
    parse_distance_restraints,
    parse_extra_residues,
    parse_viewer_sigma,
    resolve_input_file,
)
from xfel_web.pipelines import _viewer_scene, initial_viewer_sigma
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


if __name__ == "__main__":
    unittest.main()
