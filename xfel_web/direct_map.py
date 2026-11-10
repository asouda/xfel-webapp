"""Direct Fourier sum on a cube around the site.

CCP4 ``fft`` allocates one sample for every grid point of the unit cell.
A fine grid (for example 0.3 Å on a photosystem cell) overflows the 32-bit
length of that allocation. The values inside a small cube are the same sum,
taken only at the fine-grid points that fall in the cube. Coefficients that
the extended reciprocal lattice would add beyond the measured data have
amplitude zero, so the sum uses the reflections in the MTZ.
"""

from __future__ import annotations

import math
import struct
from pathlib import Path

import gemmi
import numpy as np

# A 0.1 Å cube about 40 Å on a side. Larger than the default 20 Å clip.
_MAX_CUBE_POINTS = 80_000_000
# Dense coefficient box for a ~2 Å photosystem is about 10 million cells.
_MAX_HKL_BOX = 80_000_000


# ``fft`` multiplies the three sampling intervals in a signed 32-bit integer.
# A product past this limit wraps. A negative wrap stops in CCPALC; a positive
# wrap is handed to the allocator and the process segfaults.
_FFT_GRID_LIMIT = 2**31 - 1


def parse_fft_sampling(text: str) -> tuple[int, int, int] | None:
    """Crystallographic (nx, ny, nz) printed by ``fft`` before it permutes axes."""
    import re

    match = re.search(
        r"Before permutation\s*:.*?Sampling intervals on xyz\s*\.*\s*(\d+)\s+(\d+)\s+(\d+)",
        text,
        flags=re.S,
    )
    if not match:
        return None
    return int(match.group(1)), int(match.group(2)), int(match.group(3))


def whole_cell_exceeds_fft(nxyz: tuple[int, int, int]) -> bool:
    """True when the whole-cell grid cannot be allocated by CCP4 ``fft``."""
    nx, ny, nz = nxyz
    return nx * ny * nz > _FFT_GRID_LIMIT


def _unique_first(hkl: np.ndarray, coeff: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Drop duplicate Miller indices, keeping the earliest row."""
    offset = 1 << 20
    key = (
        ((hkl[:, 0].astype(np.int64) + offset) << 42)
        | ((hkl[:, 1].astype(np.int64) + offset) << 21)
        | (hkl[:, 2].astype(np.int64) + offset)
    )
    order = np.argsort(key, kind="stable")
    key = key[order]
    hkl = hkl[order]
    coeff = coeff[order]
    keep = np.empty(len(key), dtype=bool)
    keep[0] = True
    keep[1:] = key[1:] != key[:-1]
    return hkl[keep], coeff[keep]


def expand_coefficients(miller: np.ndarray, values: np.ndarray, spacegroup: gemmi.SpaceGroup):
    """Symmetry-expand complex F, then add any missing Friedel mate as F*."""
    miller = np.asarray(miller, dtype=np.float64)
    values = np.asarray(values, dtype=np.complex128)
    blocks_h = []
    blocks_c = []
    for op in spacegroup.operations():
        matrix = np.array(op.float_seitz(), dtype=np.float64)
        rotation, translation = matrix[:3, :3], matrix[:3, 3]
        images = np.rint(miller @ rotation).astype(np.int32)
        phased = values * np.exp(-2j * np.pi * (miller @ translation))
        blocks_h.append(images)
        blocks_c.append(phased)
    hkl = np.concatenate(blocks_h)
    coeff = np.concatenate(blocks_c)
    hkl = np.concatenate([hkl, -hkl])
    coeff = np.concatenate([coeff, np.conjugate(coeff)])
    hkl, coeff = _unique_first(hkl, coeff)
    keep = np.any(hkl != 0, axis=1)
    return hkl[keep], coeff[keep]


def _bluestein(samples: np.ndarray, period: int, count: int) -> np.ndarray:
    """Y[q] = sum_p A[p] exp(-2πi p q / period), q = 0 .. count-1."""
    width = samples.shape[-1]
    index = np.arange(width)
    query = np.arange(count)
    phased = samples * np.exp(-1j * np.pi * index**2 / period)
    lags = np.arange(-(width - 1), count)
    kernel = np.exp(1j * np.pi * lags**2 / period)
    length = 1 << (int(lags.size) - 1).bit_length()
    transformed = np.fft.fft(phased, n=length, axis=-1)
    kernel_fft = np.fft.fft(kernel, n=length)
    convolved = np.fft.ifft(transformed * kernel_fft, axis=-1)
    picked = convolved[..., width - 1 : width - 1 + count]
    return picked * np.exp(-1j * np.pi * query**2 / period)


def _partial_dft(volume: np.ndarray, hmin: int, period: int, i0: int, i1: int) -> np.ndarray:
    """Replace the last axis (Miller hmin, hmin+1, ...) with grid indices i0..i1."""
    offset = np.arange(volume.shape[-1])
    query = np.arange(i1 - i0 + 1)
    phased = volume * np.exp(-2j * np.pi * offset * i0 / period)
    inner = _bluestein(phased, period, query.size)
    return inner * np.exp(-2j * np.pi * hmin * (i0 + query) / period)


def synthesize_cube(
    hkl: np.ndarray,
    coeff: np.ndarray,
    nxyz: tuple[int, int, int],
    bounds: tuple[tuple[int, int], tuple[int, int], tuple[int, int]],
) -> np.ndarray:
    """Real density on grid indices [i0, i1] inclusive, divided by the cell volume later.

    The returned array is indexed ``[ix, iy, iz]`` and is still multiplied by the
    cell volume. Callers divide by the volume.
    """
    if len(hkl) == 0:
        raise ValueError("No reflections left for the Fourier sum")
    lowers = hkl.min(axis=0)
    uppers = hkl.max(axis=0)
    spans = uppers - lowers + 1
    box = int(spans[0]) * int(spans[1]) * int(spans[2])
    if box > _MAX_HKL_BOX:
        raise RuntimeError(
            "The reflections span "
            f"{int(spans[0])} × {int(spans[1])} × {int(spans[2])} Miller indices. "
            "That box is too large to sum on the site cube."
        )
    volume = np.zeros(tuple(int(n) for n in spans), dtype=np.complex128)
    volume[hkl[:, 0] - lowers[0], hkl[:, 1] - lowers[1], hkl[:, 2] - lowers[2]] += coeff
    # Last axis is l, then k, then h. The result is indexed (iz, iy, ix).
    volume = _partial_dft(volume, int(lowers[2]), nxyz[2], bounds[2][0], bounds[2][1])
    volume = np.moveaxis(volume, 1, -1)
    volume = _partial_dft(volume, int(lowers[1]), nxyz[1], bounds[1][0], bounds[1][1])
    volume = np.moveaxis(volume, 0, -1)
    volume = _partial_dft(volume, int(lowers[0]), nxyz[0], bounds[0][0], bounds[0][1])
    return np.moveaxis(volume.real, (2, 1, 0), (0, 1, 2))


def _cube_half_side(center: np.ndarray, radius: float, extra_centers: list[np.ndarray]) -> float:
    half = float(radius)
    for other in extra_centers:
        offset = np.max(np.abs(other - center))
        half = max(half, float(offset) + float(radius))
    return half


def _grid_bounds(cell: gemmi.UnitCell, nxyz: tuple[int, int, int], center: np.ndarray, half: float):
    corners = []
    for sign_x in (-1.0, 1.0):
        for sign_y in (-1.0, 1.0):
            for sign_z in (-1.0, 1.0):
                position = gemmi.Position(
                    float(center[0] + sign_x * half),
                    float(center[1] + sign_y * half),
                    float(center[2] + sign_z * half),
                )
                fractional = cell.fractionalize(position)
                corners.append((fractional.x, fractional.y, fractional.z))
    fractional = np.array(corners, dtype=np.float64)
    bounds = []
    for axis, period in enumerate(nxyz):
        low = float(fractional[:, axis].min())
        high = float(fractional[:, axis].max())
        # One grid step outside the cube so the clip sphere can be interpolated.
        i0 = math.floor(low * period + 1e-8) - 1
        i1 = math.ceil(high * period - 1e-8) + 1
        if i1 < i0:
            i1 = i0
        bounds.append((i0, i1))
    return tuple(bounds)


def _axis_order(spacegroup_number: int) -> tuple[int, int, int]:
    """CCP4 fast, medium, slow axes. Matches fft for P1 and for P212121 and above."""
    if spacegroup_number >= 19:
        return (2, 1, 3)
    return (3, 1, 2)


def _write_map(
    path: Path,
    rho_xyz: np.ndarray,
    *,
    nxyz: tuple[int, int, int],
    bounds: tuple[tuple[int, int], tuple[int, int], tuple[int, int]],
    cell: gemmi.UnitCell,
    spacegroup: gemmi.SpaceGroup,
    rms: float,
) -> None:
    mapc, mapr, maps = _axis_order(spacegroup.number)
    file_rho = np.transpose(rho_xyz, (maps - 1, mapr - 1, mapc - 1)).astype("<f4", copy=False)
    starts = [bounds[0][0], bounds[1][0], bounds[2][0]]
    # Clipper allocates the header sampling (NX, NY, NZ) before it looks at the
    # stored window. Keep the crystallographic grid indices, and scale each
    # cell edge so index/sampling*edge is unchanged. The sampling is then the
    # cube itself. Symmetry stays P1: these fractional coordinates belong to
    # the cube, not to the crystal cell.
    counts = tuple(stop - start + 1 for start, stop in bounds)
    header_cell = (
        counts[0] * cell.a / nxyz[0],
        counts[1] * cell.b / nxyz[1],
        counts[2] * cell.c / nxyz[2],
        cell.alpha,
        cell.beta,
        cell.gamma,
    )
    header = bytearray(1024)
    struct.pack_into(
        "<3i",
        header,
        0,
        file_rho.shape[2],
        file_rho.shape[1],
        file_rho.shape[0],
    )
    struct.pack_into("<i", header, 3 * 4, 2)
    struct.pack_into(
        "<3i",
        header,
        4 * 4,
        starts[mapc - 1],
        starts[mapr - 1],
        starts[maps - 1],
    )
    struct.pack_into("<3i", header, 7 * 4, counts[0], counts[1], counts[2])
    struct.pack_into("<6f", header, 10 * 4, *header_cell)
    struct.pack_into("<3i", header, 16 * 4, mapc, mapr, maps)
    struct.pack_into("<f", header, 19 * 4, float(file_rho.min()))
    struct.pack_into("<f", header, 20 * 4, float(file_rho.max()))
    struct.pack_into("<f", header, 21 * 4, 0.0)
    struct.pack_into("<i", header, 22 * 4, 1)
    symmetry = b"X,Y,Z".ljust(80)
    struct.pack_into("<i", header, 23 * 4, len(symmetry))
    header[52 * 4 : 52 * 4 + 4] = b"MAP "
    header[53 * 4 : 53 * 4 + 4] = b"DA\x00\x00"
    struct.pack_into("<f", header, 54 * 4, float(rms))
    struct.pack_into("<i", header, 55 * 4, 1)
    title = b"site cube, full-cell sigma".ljust(80)
    header[56 * 4 : 56 * 4 + 80] = title
    path.write_bytes(bytes(header) + symmetry + file_rho.tobytes(order="C"))


def direct_site_cube(
    mtz_path: Path,
    map_out: Path,
    f_label: str,
    phi_label: str,
    d_low: float,
    d_high: float,
    nxyz: tuple[int, int, int],
    center: tuple[float, float, float],
    radius: float,
    extra_centers: list[tuple[float, float, float]] | None = None,
    emit=None,
) -> str:
    """Write the fine-grid density in a cube centered on ``center``.

    ``radius`` is the clip radius in ångströms. The cube's half-side is at
    least that, and grows so every extra site center keeps its own clip
    inside the same cube. The written cell is the cube, in P1, with the
    crystallographic grid indices kept so model coordinates still land on
    the same density. The header rms is the full-cell root-mean-square
    implied by the coefficients, which is the sigma ``fft`` would store.
    """
    mtz = gemmi.read_mtz_file(str(mtz_path))
    asu = mtz.get_f_phi(f_label, phi_label)
    miller = np.array(asu.miller_array, dtype=np.int32)
    values = np.array(asu.value_array, dtype=np.complex128)
    finite = np.isfinite(values.real) & np.isfinite(values.imag)
    miller = miller[finite]
    values = values[finite]
    spacing = np.asarray(mtz.cell.calculate_d_array(miller), dtype=np.float64)
    keep = (spacing <= d_low + 1e-4) & (spacing + 1e-4 >= d_high)
    miller = miller[keep]
    values = values[keep]
    if len(miller) == 0:
        raise ValueError(
            f"No {f_label} reflections between {d_low:.2f} and {d_high:.2f} Å in {mtz_path.name}"
        )
    hkl, coeff = expand_coefficients(miller, values, mtz.spacegroup)
    origin = np.asarray(center, dtype=np.float64)
    others = [np.asarray(point, dtype=np.float64) for point in (extra_centers or [])]
    half = _cube_half_side(origin, radius, others)
    bounds = _grid_bounds(mtz.cell, nxyz, origin, half)
    counts = [stop - start + 1 for start, stop in bounds]
    points = counts[0] * counts[1] * counts[2]
    if points > _MAX_CUBE_POINTS:
        raise RuntimeError(
            f"The site cube is {counts[0]} × {counts[1]} × {counts[2]} grid points "
            f"({points / 1e6:.1f} million) on the {d_high:.2f} Å sampling. "
            "Reduce the clip radius or use a coarser FFT grid resolution for a final map."
        )
    lines = [
        "CCP4 fft cannot allocate this whole-cell grid, so the density is summed "
        "directly on the site cube.",
        f"Sampling intervals X Y Z: {nxyz[0]} {nxyz[1]} {nxyz[2]}",
        (
            f"Cube center {origin[0]:.2f} {origin[1]:.2f} {origin[2]:.2f} Å, "
            f"half-side {half:.2f} Å (clip radius {radius:.2f} Å)."
        ),
        "Grid limits X Y Z: " + " ".join(f"{start}:{stop}" for start, stop in bounds),
        f"Expanded reflections: {len(hkl)}",
    ]

    def say(message: str) -> None:
        lines.append(message)
        if emit is not None:
            emit(message)

    for line in lines:
        if emit is not None:
            emit(line)
    say("Summing the site cube.")
    density = np.array(synthesize_cube(hkl, coeff, nxyz, bounds), dtype=np.float64, copy=True)
    volume = float(mtz.cell.volume)
    density /= volume
    rms = math.sqrt(float(np.sum(np.abs(coeff) ** 2))) / volume
    partial = map_out.with_name(map_out.name + ".partial")
    if partial.is_file():
        partial.unlink()
    try:
        _write_map(
            partial,
            density.astype(np.float32, copy=False),
            nxyz=nxyz,
            bounds=bounds,
            cell=mtz.cell,
            spacegroup=mtz.spacegroup,
            rms=rms,
        )
        partial.replace(map_out)
    except Exception:
        partial.unlink(missing_ok=True)
        raise
    say(
        f"Full-cell sigma {rms:.6f}. Wrote {map_out.name} "
        f"({map_out.stat().st_size / 1e6:.1f} MB)."
    )
    return "\n".join(lines) + "\n"
