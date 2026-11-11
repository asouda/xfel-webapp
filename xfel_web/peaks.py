"""Local density maxima on a CCP4 map.

The map reader, neighborhood maximum, and spline placement are the voxel
search from ``xfel-gaussian`` ``sequential_peaks_gaussian``. They live here so
the reconstruction does not import that package.
"""

import gemmi
import numpy as np
from scipy.ndimage import map_coordinates, maximum_filter, spline_filter
from scipy.optimize import minimize


class DensityMap:
    """CCP4 grid with Cartesian coordinates for each stored index."""

    def __init__(self, path):
        self.ccp4 = gemmi.read_ccp4_map(str(path), setup=False)
        self.original = np.array(self.ccp4.grid, dtype=np.float64, copy=True)
        self.axis = np.asarray(self.ccp4.axis_positions(), dtype=int)
        self.start = np.asarray([self.ccp4.header_i32(i) for i in (5, 6, 7)], dtype=float)
        self.sampling = np.asarray([self.ccp4.header_i32(i) for i in (8, 9, 10)], dtype=float)
        if np.any(self.sampling <= 0):
            raise ValueError("Invalid CCP4 sampling intervals")
        orth = np.array(self.ccp4.grid.unit_cell.orth.mat, dtype=float)
        # The 3 array-index directions are not necessarily x,y,z crystallographic directions.
        self.origin = self.to_cart(np.zeros(3), orth)
        self.step = np.stack(
            [self.to_cart(np.eye(3)[i], orth) - self.origin for i in range(3)],
            axis=1,
        )
        self.invstep = np.linalg.inv(self.step)
        self.voxel_A3 = abs(np.linalg.det(self.step))
        self.residual = self.original.copy()

    def to_cart(self, idx, orth=None):
        if orth is None:
            orth = np.array(self.ccp4.grid.unit_cell.orth.mat, dtype=float)
        fractional = (np.asarray(idx, dtype=float)[..., self.axis] + self.start[self.axis]) / self.sampling
        return fractional @ orth.T

    def to_index(self, xyz):
        return (np.asarray(xyz) - self.origin) @ self.invstep.T


def choose_peak(rho, original_mask, exclusion, neighborhood_size=3):
    """Strongest voxel that is a local maximum in a cubic neighborhood.

    A candidate equals the strongest valid value in a cube of side
    ``neighborhood_size``. Outside the array the filter uses −∞, so the border
    needs no separate margin. Returns ``(index, density)``, or None when no
    finite candidate remains.
    """
    valid = np.isfinite(rho) & original_mask & ~exclusion
    candidates = valid & (
        rho
        == maximum_filter(
            np.where(valid, rho, -np.inf),
            size=neighborhood_size,
            mode="constant",
            cval=-np.inf,
        )
    )
    if not np.any(candidates):
        return None
    ind = np.argmax(np.where(candidates, rho, -np.inf))
    idx = np.unravel_index(ind, rho.shape)
    if not np.isfinite(rho[idx]):
        return None
    return np.array(idx), float(rho[idx])


def _density_on_spline(coeffs, order, point):
    sampled = map_coordinates(
        coeffs,
        np.asarray(point, dtype=float).reshape(3, 1),
        order=order,
        prefilter=False,
        mode="nearest",
    )
    return float(sampled[0])


def _maximize_spline(patch, order, anchor):
    """Index of the spline maximum inside a neighborhood patch.

    ``anchor`` is the detected voxel in patch coordinates. It is returned when
    the interpolant is no higher anywhere else in the patch.
    """
    coeffs = spline_filter(np.asarray(patch, dtype=np.float64), order=order)
    shape = np.asarray(patch.shape, dtype=int)
    anchor = np.clip(np.asarray(anchor, dtype=float), 0.0, shape - 1.0)
    axes = []
    for length in shape:
        if length <= 1:
            axes.append(np.zeros(1))
            continue
        count = min((int(length) - 1) * 8 + 1, 49)
        axes.append(np.linspace(0.0, float(length - 1), count))
    grids = np.meshgrid(*axes, indexing="ij")
    coords = np.vstack([grid.ravel() for grid in grids])
    sampled = map_coordinates(coeffs, coords, order=order, prefilter=False, mode="nearest")
    if not np.isfinite(sampled).any():
        raise ValueError("Spline evaluation is not finite")
    start = coords[:, int(np.argmax(np.where(np.isfinite(sampled), sampled, -np.inf)))]
    free = [axis for axis, length in enumerate(shape) if length > 1]
    polished = start.copy()
    if free:
        def negative(reduced):
            point = start.copy()
            point[free] = reduced
            return -_density_on_spline(coeffs, order, point)

        bounds = [(0.0, float(shape[axis] - 1)) for axis in free]
        result = minimize(negative, start[free], method="L-BFGS-B", bounds=bounds)
        if result.x is not None and np.all(np.isfinite(result.x)):
            polished = start.copy()
            polished[free] = np.clip(np.asarray(result.x, dtype=float), 0.0, shape[free] - 1.0)
    anchor_value = _density_on_spline(coeffs, order, anchor)
    best = polished if _density_on_spline(coeffs, order, polished) >= _density_on_spline(coeffs, order, start) else start
    best_value = _density_on_spline(coeffs, order, best)
    tolerance = 1e-8 * max(1.0, abs(anchor_value))
    if best_value <= anchor_value + tolerance:
        return anchor
    return best


def spline_maximum_index(rho, idx, neighborhood_size, spline_order):
    """Continuous stored index where a spline of the neighborhood is highest.

    The window has side ``neighborhood_size`` and is centered on ``idx``.
    It is clipped to the array. The interpolant is the ``scipy.ndimage``
    spline of ``spline_order`` on stored-index coordinates: 3 is cubic, and
    orders run from 0 through 5. Orders below 2 reach their maximum at a
    sample. A window that is not entirely finite, or a spline that cannot
    be built, returns ``idx``.
    """
    idx = np.asarray(idx, dtype=int)
    if int(spline_order) < 2:
        return idx.astype(float)
    half = int(neighborhood_size) // 2
    lo = np.maximum(idx - half, 0)
    hi = np.minimum(idx + half + 1, rho.shape)
    patch = np.asarray(rho[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]], dtype=np.float64)
    if patch.size == 0 or not np.all(np.isfinite(patch)):
        return idx.astype(float)
    try:
        anchor = idx.astype(float) - lo.astype(float)
        point = _maximize_spline(patch, int(spline_order), anchor)
    except (RuntimeError, ValueError):
        return idx.astype(float)
    return lo.astype(float) + np.asarray(point, dtype=float)
