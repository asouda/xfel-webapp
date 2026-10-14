import unittest
from pathlib import Path

from xfel_web.ccp4 import choose_resolution, delta_b_tag, pymol_scene_links, residue_conformers
from xfel_web.forms import (
    FormError,
    parse_angle_restraints,
    parse_copies,
    parse_distance_restraints,
    parse_extra_residues,
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
