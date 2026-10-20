"""Synthetic meshes and surface volumes with planted index offsets / shifts / depth offsets (self-tests only)."""
from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates

from .conventions import exporter_grid_map
from .exporter import F32, apply_affine_f32, resample_linear_preserving_invalids, transform_tifxyz
from .slicecentre import layer_depth_um, slice_centre


def curved_sheet(rows=120, cols=150, spacing=20.0, radius=1800.0, holes=6, seed=0):
    """A wavy cylinder-section sheet in coarse-volume voxels, grid spacing ~20 vx, with invalid holes and a ragged
    border. Returns float32 (rows, cols, 3) with -1 for invalid."""
    rng = np.random.default_rng(seed)
    r, c = np.mgrid[:rows, :cols].astype(np.float64)
    theta = c * spacing / radius
    wob = 25.0 * np.sin(r / 17.0) * np.cos(c / 23.0)
    R = radius + wob
    x = 3000 + R * np.cos(theta)
    y = 2500 + R * np.sin(theta)
    z = 12000 + r * spacing + 6.0 * np.sin(c / 9.0)
    P = np.stack([x, y, z], -1).astype(F32)
    bad = np.zeros((rows, cols), bool)
    for _ in range(holes):
        cy, cx, rad = rng.integers(10, rows - 10), rng.integers(10, cols - 10), rng.integers(2, 6)
        bad |= (r - cy) ** 2 + (c - cx) ** 2 <= rad ** 2
    edge = rng.integers(0, 4, size=rows)
    bad |= c < edge[:, None]
    bad[:2] = True
    P[bad] = F32(-1)
    return P


def demo_affine(scale=3.8976, rot_deg=0.4, t=(340.0, 215.0, -49800.0)):
    a = np.deg2rad(rot_deg)
    Rz = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
    L = scale * Rz @ np.diag([1.0, 1.0002, 0.9998])
    return np.c_[L, np.asarray(t)]


def planted_resample(points: np.ndarray, new_h: int, new_w: int, shift_cells=(0.0, 0.0)):
    """Like the exporter's resampler but with an extra planted source-index shift (source cells)."""
    rows, cols = points.shape[:2]
    out = np.full((new_h, new_w, 3), F32(-1), dtype=F32)
    # half-pixel index rule (destination grid -> source grid) from the clean conventions module
    sx = exporter_grid_map(np.arange(new_w), new_w, cols) + shift_cells[1]
    sy = exporter_grid_map(np.arange(new_h), new_h, rows) + shift_cells[0]
    ok_x = (sx >= 0) & (sx <= cols - 1)
    ok_y = (sy >= 0) & (sy <= rows - 1)
    sx = np.clip(sx, 0, cols - 1); sy = np.clip(sy, 0, rows - 1)
    x0 = np.minimum(np.floor(sx).astype(int), cols - 2); y0 = np.minimum(np.floor(sy).astype(int), rows - 2)
    fx = (sx - x0)[None, :, None]; fy = (sy - y0)[:, None, None]
    P = points.astype(np.float64)
    v = (points[..., 0] != -1)
    ok = v[np.ix_(y0, x0)] & v[np.ix_(y0, x0 + 1)] & v[np.ix_(y0 + 1, x0)] & v[np.ix_(y0 + 1, x0 + 1)]
    ok &= ok_y[:, None] & ok_x[None, :]
    val = (P[np.ix_(y0, x0)] * (1 - fx) + P[np.ix_(y0, x0 + 1)] * fx) * (1 - fy) + \
          (P[np.ix_(y0 + 1, x0)] * (1 - fx) + P[np.ix_(y0 + 1, x0 + 1)] * fx) * fy
    out[ok] = val[ok].astype(F32)
    return out


def texture_world(shape_wuv, sigma=(2.0, 1.2, 1.6), seed=0):
    rng = np.random.default_rng(seed)
    t = gaussian_filter(rng.standard_normal(shape_wuv).astype(np.float32), sigma)
    t -= t.mean(); t /= t.std()
    return t


def synth_sv_pair(true_map, px_f=2.403, px_c=9.366, n_f=49, n_c=12, Hc=150, Wc=150, depth_true_um=0.0,
                  lattice_um=2.403, seed=0, shift_field=None):
    """Returns (fine level-0 array (n_f, Hf, Wf), coarse array (n_c, Hc, Wc)).
    Physical sheet coords: fine pixel (k, i, j) <-> (w, u, v) = ((k - cF) px_f, i px_f, j px_f).
    true_map = (ay, by, ax, bx): coarse canvas = a + b * fine canvas (+ shift_field(p_f_y, p_f_x) -> (dy, dx) coarse px).
    Coarse layer k at depth d shows the fine content at d + depth_true_um (sign convention of ct.py)."""
    ay, by, ax, bx = true_map
    # fine extent needed to cover the coarse canvas
    Hf = int(np.ceil((Hc + 4 - ay) / by)) + 8
    Wf = int(np.ceil((Wc + 4 - ax) / bx)) + 8
    cF, cC = slice_centre(n_f), slice_centre(n_c)
    wmax = max(cF * px_f, cC * px_c + abs(depth_true_um) + px_c)
    nw = int(np.ceil(2 * wmax / lattice_um)) + 6
    w0 = -wmax - 3 * lattice_um
    nu = int(np.ceil(Hf * px_f / lattice_um)) + 6
    nv = int(np.ceil(Wf * px_f / lattice_um)) + 6
    T = texture_world((nw, nu, nv), seed=seed)
    # fine: sample at voxel centres (texture already smooth on the 2.4 um lattice)
    k = np.arange(n_f); i = np.arange(Hf); j = np.arange(Wf)
    W_ = (layer_depth_um(k, n_f, px_f) - w0) / lattice_um
    U_ = (i * px_f) / lattice_um
    V_ = (j * px_f) / lattice_um
    G = np.meshgrid(W_, U_, V_, indexing="ij")
    fine = map_coordinates(T, G, order=1, mode="nearest").astype(np.float32)
    # coarse: texture blurred to the coarse voxel (box ~ Gaussian 0.42 px_c), sampled at mapped centres
    sb = 0.42 * px_c / lattice_um
    Tb = gaussian_filter(T, (sb, sb, sb))
    kc = np.arange(n_c); qy = np.arange(Hc, dtype=np.float64); qx = np.arange(Wc, dtype=np.float64)
    QY, QX = np.meshgrid(qy, qx, indexing="ij")
    # invert the true map (iterate if a shift field is present)
    PY = (QY - ay) / by; PX = (QX - ax) / bx
    if shift_field is not None:
        for _ in range(5):
            dy, dx = shift_field(PY, PX)
            PY = (QY - dy - ay) / by; PX = (QX - dx - ax) / bx
    coarse = np.zeros((n_c, Hc, Wc), np.float32)
    for kk in kc:
        w = layer_depth_um(kk, n_c, px_c, depth_true_um)
        coarse[kk] = map_coordinates(Tb, [np.full_like(PY, (w - w0) / lattice_um), PY * px_f / lattice_um,
                                          PX * px_f / lattice_um], order=1, mode="nearest")
    # 8-bit-like positive range, no exact zeros (zeros mean 'no data' in ct.py validity)
    fine = np.clip(128 + 30 * fine, 1, 255).astype(np.float32)
    coarse = np.clip(128 + 30 * coarse / coarse.std(), 1, 255).astype(np.float32)
    return fine, coarse


def block_mean(a: np.ndarray, f: int):
    """In-plane 2^l block mean (pyramid level convention assumed in ct.py)."""
    n, H, W = a.shape
    H2, W2 = H // f * f, W // f * f
    return a[:, :H2, :W2].reshape(n, H2 // f, f, W2 // f, f).mean((2, 4)).astype(np.float32)
