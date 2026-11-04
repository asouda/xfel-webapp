"""The site-cube Fourier sum matches CCP4 fft where both can run."""

import struct
import subprocess
import tempfile
import unittest
from pathlib import Path

import gemmi
import numpy as np

from xfel_web.direct_map import (
    _write_map,
    direct_site_cube,
    expand_coefficients,
    parse_fft_sampling,
    synthesize_cube,
    whole_cell_exceeds_fft,
)


def _write_mtz(path: Path, spacegroup: int, cell, rows):
    mtz = gemmi.Mtz()
    mtz.spacegroup = gemmi.find_spacegroup_by_number(spacegroup)
    mtz.cell = gemmi.UnitCell(*cell)
    mtz.add_dataset("data")
    for label, typ in (("H", "H"), ("K", "H"), ("L", "H"), ("F", "F"), ("PHI", "P")):
        mtz.add_column(label, typ)
    mtz.set_data(np.asarray(rows, dtype=np.float32))
    mtz.write_to_file(str(path))
    return mtz


def _map_array(path: Path) -> np.ndarray:
    raw = path.read_bytes()
    nc, nr, ns = struct.unpack_from("<3i", raw, 0)
    nsymbt = struct.unpack_from("<i", raw, 23 * 4)[0]
    count = nc * nr * ns
    data = np.frombuffer(raw, dtype="<f4", count=count, offset=1024 + nsymbt)
    return data.reshape(ns, nr, nc)


def _brute(hkl, coeff, nxyz, bounds, volume):
    rho = np.zeros(tuple(stop - start + 1 for start, stop in bounds), dtype=np.float64)
    ix = np.arange(bounds[0][0], bounds[0][1] + 1)
    iy = np.arange(bounds[1][0], bounds[1][1] + 1)
    iz = np.arange(bounds[2][0], bounds[2][1] + 1)
    for (h, k, l), value in zip(hkl, coeff):
        phase = (
            h * ix[:, None, None] / nxyz[0]
            + k * iy[None, :, None] / nxyz[1]
            + l * iz[None, None, :] / nxyz[2]
        )
        rho += (value * np.exp(-2j * np.pi * phase)).real
    return rho / volume


class DirectMapTests(unittest.TestCase):
    def test_sampling_line_is_the_crystallographic_grid(self):
        text = """
        Before permutation :
           Sampling intervals on xyz ............. 1176 2200 3080
        After permutation :
           Sampling intervals on xyz ............. 1176 3080 2200
        """
        self.assertEqual(parse_fft_sampling(text), (1176, 2200, 3080))

    def test_allocator_limit_follows_the_true_grid_product(self):
        # 0.5 Å on the photosystem cell fits. 0.3 Å wraps negative (CCPALC).
        # 0.2 Å wraps positive and fft then segfaults.
        self.assertFalse(whole_cell_exceeds_fft((704, 1320, 1848)))
        self.assertTrue(whole_cell_exceeds_fft((1176, 2200, 3080)))
        self.assertTrue(whole_cell_exceeds_fft((1764, 3344, 4624)))

    def test_cube_matches_a_direct_sum_and_fft(self):
        rows = [
            [1, 2, 3, 50.0, 30.0],
            [2, 0, 1, 20.0, -40.0],
            [0, 1, 0, 10.0, 0.0],
            [4, -1, 2, 15.0, 80.0],
        ]
        cell = (30, 40, 50, 90, 90, 90)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mtz_path = root / "small.mtz"
            mtz = _write_mtz(mtz_path, 1, cell, rows)
            asu = mtz.get_f_phi("F", "PHI")
            miller = np.array(asu.miller_array, dtype=np.int32)
            values = np.array(asu.value_array, dtype=np.complex128)
            spacing = np.asarray(mtz.cell.calculate_d_array(miller), dtype=np.float64)
            # RESOLUTION 20 2 drops the 40 Å reflection, as fft does.
            keep = (spacing <= 20 + 1e-4) & (spacing + 1e-4 >= 2)
            hkl, coeff = expand_coefficients(miller[keep], values[keep], mtz.spacegroup)
            nxyz = (44, 60, 76)
            bounds = ((5, 12), (8, 20), (20, 30))
            got = synthesize_cube(hkl, coeff, nxyz, bounds) / mtz.cell.volume
            brute = _brute(hkl, coeff, nxyz, bounds, mtz.cell.volume)
            self.assertLess(np.max(np.abs(got - brute)), 1e-9)

            setup = Path("/Applications/ccp4-9/bin/ccp4.setup-sh")
            if not setup.is_file():
                return
            map_path = root / "small.ccp4"
            script = f"""
            set -eo pipefail
            source {setup}
            fft hklin {mtz_path} mapout {map_path} << 'EOF'
            RESOLUTION 20 2
            LABIN F1=F PHI=PHI
            END
            EOF
            """
            subprocess.run(["/bin/bash", "-c", script], check=True, cwd=root)
            # P1 fft: sections Y, rows X, columns Z.
            data = _map_array(map_path)
            ix0, ix1 = bounds[0]
            iy0, iy1 = bounds[1]
            iz0, iz1 = bounds[2]
            sub = data[iy0:iy1 + 1, ix0:ix1 + 1, iz0:iz1 + 1]
            pred = np.transpose(got, (1, 0, 2))
            self.assertEqual(sub.shape, pred.shape)
            self.assertLess(float(np.max(np.abs(sub - pred))), 2e-6)

    def test_p212121_cube_matches_fft(self):
        setup = Path("/Applications/ccp4-9/bin/ccp4.setup-sh")
        if not setup.is_file():
            self.skipTest("CCP4 fft is not installed")
        cell = (40, 50, 60, 90, 90, 90)
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mtz_path = root / "sg19.mtz"
            _write_mtz(mtz_path, 19, cell, [[3, 4, 5, 80.0, 25.0]])
            map_path = root / "sg19.ccp4"
            script = f"""
            set -eo pipefail
            source {setup}
            fft hklin {mtz_path} mapout {map_path} << 'EOF'
            RESOLUTION 30 3
            LABIN F1=F PHI=PHI
            END
            EOF
            """
            subprocess.run(["/bin/bash", "-c", script], check=True, cwd=root)
            header = map_path.read_bytes()[:1024]
            nxyz = struct.unpack_from("<3i", header, 7 * 4)
            cube_path = root / "cube.ccp4"
            direct_site_cube(
                mtz_path,
                cube_path,
                "F",
                "PHI",
                30,
                3,
                nxyz,
                (8.0, 10.0, 12.0),
                4.0,
            )
            whole = _map_array(map_path)
            cube = _map_array(cube_path)
            cube_header = cube_path.read_bytes()[:1024]
            nc, nr, ns = struct.unpack_from("<3i", cube_header, 0)
            cstart, rstart, sstart = struct.unpack_from("<3i", cube_header, 4 * 4)
            mapc, mapr, maps = struct.unpack_from("<3i", cube_header, 16 * 4)
            # sg 19 fft stores sections along Z, rows along X, columns along Y.
            # The cube keeps that axis order. Its header sampling is the cube,
            # in P1, with the cell scaled so a grid index keeps its coordinate.
            self.assertEqual((mapc, mapr, maps), (2, 1, 3))
            self.assertEqual(struct.unpack_from("<i", cube_header, 22 * 4)[0], 1)
            header_nxyz = struct.unpack_from("<3i", cube_header, 7 * 4)
            header_cell = struct.unpack_from("<6f", cube_header, 10 * 4)
            counts = [None, None, None]
            for length, axis in zip((nc, nr, ns), (mapc, mapr, maps)):
                counts[axis - 1] = length
            self.assertEqual(header_nxyz, tuple(counts))
            for axis in range(3):
                self.assertAlmostEqual(
                    header_cell[axis] / header_nxyz[axis],
                    cell[axis] / nxyz[axis],
                    places=5,
                )
            ix = int(round(8.0 / cell[0] * nxyz[0]))
            iy = int(round(10.0 / cell[1] * nxyz[1]))
            iz = int(round(12.0 / cell[2] * nxyz[2]))
            file_index = {
                mapc: (iy - cstart, nc),
                mapr: (ix - rstart, nr),
                maps: (iz - sstart, ns),
            }
            # whole and cube are both [section, row, column] = [iz, ix, iy]
            self.assertAlmostEqual(
                float(whole[iz, ix, iy]),
                float(cube[file_index[maps][0], file_index[mapr][0], file_index[mapc][0]]),
                places=5,
            )
            arms = struct.unpack_from("<f", cube_header, 54 * 4)[0]
            self.assertGreater(arms, 0)
            self.assertLess(nc * nr * ns, struct.unpack_from("<3i", header, 0)[0]
                            * struct.unpack_from("<3i", header, 0)[1]
                            * struct.unpack_from("<3i", header, 0)[2])

    def test_coot_places_the_cube_on_the_model_coordinate(self):
        coot = Path("/Applications/ccp4-9/bin/coot")
        if not coot.is_file():
            self.skipTest("Coot is not installed")
        # A negative crystallographic index, stored on a cube-sized P1 grid.
        cell = gemmi.UnitCell(100, 200, 300, 90, 90, 90)
        nxyz = (100, 200, 300)
        bounds = ((-20, -16), (100, 104), (200, 204))
        rho = np.zeros((5, 5, 5), dtype=np.float32)
        rho[2, 2, 2] = 7.0
        index = (-18, 102, 202)
        point = (
            index[0] / nxyz[0] * cell.a,
            index[1] / nxyz[1] * cell.b,
            index[2] / nxyz[2] * cell.c,
        )
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            map_path = root / "cube.ccp4"
            _write_map(
                map_path,
                rho,
                nxyz=nxyz,
                bounds=bounds,
                cell=cell,
                spacegroup=gemmi.find_spacegroup_by_number(19),
                rms=0.5,
            )
            header = map_path.read_bytes()[:1024]
            self.assertEqual(struct.unpack_from("<3i", header, 7 * 4), (5, 5, 5))
            script = root / "probe.py"
            script.write_text(
                "\n".join([
                    "imol = read_ccp4_map(%r, 0)" % str(map_path),
                    "val = density_at_point(imol, %.10f, %.10f, %.10f)" % point,
                    "open(%r, 'w').write('%%.4f\\n' %% val)" % str(root / "value.txt"),
                    "try:",
                    "    coot_real_exit(0)",
                    "except NameError:",
                    "    pass",
                    "",
                ])
            )
            subprocess.run(
                [str(coot), "--no-graphics", "--script", str(script)],
                check=True,
                cwd=root,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
            )
            self.assertEqual((root / "value.txt").read_text().strip(), "7.0000")


if __name__ == "__main__":
    unittest.main()
