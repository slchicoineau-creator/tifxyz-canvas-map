"""Reverse route (TIFXYZ-REVERSE, prereg/2026-10-13_partner_tifxyz_reverse_route.md): coarse vertex -> fine surface.

For every (sampled) coarse vertex: map it into the fine volume with the catalog affine M (coarse -> fine), start at the
exporter formula position on the fine grid, project onto the fine bilinear surface (Gauss-Newton), then measure
 (a) the normal distance to the fine surface (coarse voxels),
 (b) the in-plane residual of the found position vs the formula (coarse canvas px).
The forward route (fit.inverse_map) is untouched.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import binary_erosion

from . import conventions as cv
from .conventions import MeshGeom
from .fit import linear_fit_axis

ON_SURFACE_COARSE_VX = 0.05
INPLANE_P95_PX = 0.003
FRAC_ON_MIN = 0.95
TANG_CONV_VX = 1e-3
IN_FILES_MIN_PX = 3.0


def _bilinear_gather(points, valid, gy, gx):
    """Like fit.bilinear but gathers the 4 corners first (no float64 copy of the whole mesh)."""
    rows, cols = points.shape[:2]
    y0 = np.floor(gy).astype(np.int64)
    x0 = np.floor(gx).astype(np.int64)
    inside = (y0 >= 0) & (x0 >= 0) & (y0 + 1 <= rows - 1) & (x0 + 1 <= cols - 1)
    y0c = np.clip(y0, 0, rows - 2)
    x0c = np.clip(x0, 0, cols - 2)
    fy = (gy - y0c)[:, None]
    fx = (gx - x0c)[:, None]
    p00, p01 = points[y0c, x0c].astype(np.float64), points[y0c, x0c + 1].astype(np.float64)
    p10, p11 = points[y0c + 1, x0c].astype(np.float64), points[y0c + 1, x0c + 1].astype(np.float64)
    ok = inside & valid[y0c, x0c] & valid[y0c, x0c + 1] & valid[y0c + 1, x0c] & valid[y0c + 1, x0c + 1]
    val = (p00 * (1 - fx) + p01 * fx) * (1 - fy) + (p10 * (1 - fx) + p11 * fx) * fy
    du = (p01 - p00) * (1 - fy) + (p11 - p10) * fy
    dv = (p10 - p00) * (1 - fx) + (p11 - p01) * fx
    return val, du, dv, ok


def _project(points, valid, Y, gy, gx, iters=12):
    for _ in range(iters):
        val, du, dv, ok = _bilinear_gather(points, valid, gy, gx)
        r = Y - val
        a11 = (du * du).sum(1); a12 = (du * dv).sum(1); a22 = (dv * dv).sum(1)
        b1 = (du * r).sum(1); b2 = (dv * r).sum(1)
        det = a11 * a22 - a12 * a12
        det = np.where(np.abs(det) < 1e-12, np.nan, det)
        gx = gx + np.clip(np.nan_to_num((a22 * b1 - a12 * b2) / det), -2, 2)
        gy = gy + np.clip(np.nan_to_num((a11 * b2 - a12 * b1) / det), -2, 2)
    val, du, dv, ok = _bilinear_gather(points, valid, gy, gx)
    r = Y - val
    nrm = np.cross(du, dv)
    nrm /= np.maximum(np.linalg.norm(nrm, axis=1, keepdims=True), 1e-12)
    nd = (r * nrm).sum(1)
    tang = np.linalg.norm(r - nd[:, None] * nrm, axis=1)
    return gy, gx, nd, tang, ok


def reverse_route(fine_pts, fine_valid, coarse_pts, coarse_valid, gF: MeshGeom, gC: MeshGeom, M,
                  px_fine: float, px_coarse: float, n_max: int = 100000, seed: int = 0, erode: int = 2, iters: int = 12):
    """Returns a dict of per-pair summary numbers (see the prereg, section 2)."""
    v = binary_erosion(coarse_valid, iterations=erode) if erode else coarse_valid
    rr, cc = np.nonzero(v)
    n_valid_eroded = int(rr.size)
    if rr.size > n_max:
        pick = np.random.default_rng(seed).choice(rr.size, n_max, replace=False)
        rr, cc = rr[pick], cc[pick]
    pCy = gC.grid_to_canvas(rr.astype(np.float64), 0)
    pCx = gC.grid_to_canvas(cc.astype(np.float64), 1)
    M = np.asarray(M, np.float64)
    X = coarse_pts[rr, cc].astype(np.float64)
    Y = X @ M[:, :3].T + M[:, 3]
    pFy = cv.exact_canvas_map_inv(pCy, gF, gC, 0)
    pFx = cv.exact_canvas_map_inv(pCx, gF, gC, 1)
    gy0 = gF.canvas_to_grid(pFy, 0)
    gx0 = gF.canvas_to_grid(pFx, 1)
    gy, gx, nd_f, tg, ok = _project(fine_pts, fine_valid, Y, gy0, gx0, iters)
    nd = nd_f * (px_fine / px_coarse)
    conv = ok & (tg < TANG_CONV_VX)
    on = conv & (np.abs(nd) <= ON_SURFACE_COARSE_VX)
    qFy = gF.grid_to_canvas(gy, 0)
    qFx = gF.grid_to_canvas(gx, 1)
    ey = cv.exact_canvas_map(qFy, gF, gC, 0) - pCy
    ex = cv.exact_canvas_map(qFx, gF, gC, 1) - pCx
    inpl = np.hypot(ey, ex)
    out = dict(n_valid_eroded=n_valid_eroded, n_sampled=int(rr.size), n_converged=int(conv.sum()), n_on_surface=int(on.sum()),
               frac_on=float(on.mean()) if rr.size else float("nan"))
    if conv.any():
        out.update(median_normal_coarse_vx=float(np.median(np.abs(nd[conv]))),
                   p95_normal_coarse_vx=float(np.percentile(np.abs(nd[conv]), 95)),
                   inplane_p95_px=float(np.percentile(inpl[conv], 95)), inplane_median_px=float(np.median(inpl[conv])),
                   inplane_max_px=float(inpl[conv].max()))
    else:
        out.update(median_normal_coarse_vx=float("nan"), p95_normal_coarse_vx=float("nan"),
                   inplane_p95_px=float("nan"), inplane_median_px=float("nan"), inplane_max_px=float("nan"))
    out["pair_pass"] = bool(out["frac_on"] >= FRAC_ON_MIN and np.isfinite(out["inplane_p95_px"])
                            and out["inplane_p95_px"] <= INPLANE_P95_PX)
    # offset vs naive (readable-pair verdict line)
    sel = on if on.sum() >= 100 else conv
    if sel.sum() >= 100:
        k = px_fine / px_coarse
        rrF, ccF = np.nonzero(fine_valid)   # fine valid box (full mask)
        y0, y1 = gF.grid_to_canvas(rrF.min(), 0), gF.grid_to_canvas(rrF.max(), 0)
        x0, x1 = gF.grid_to_canvas(ccF.min(), 1), gF.grid_to_canvas(ccF.max(), 1)
        ay, by, _ = linear_fit_axis(qFy[sel], pCy[sel])
        ax, bx, _ = linear_fit_axis(qFx[sel], pCx[sel])
        yc, xc = 0.5 * (y0 + y1), 0.5 * (x0 + x1)
        oy, ox = (ay + by * yc) - k * yc, (ax + bx * xc) - k * xc
        out.update(offset_vs_naive_centre_px=[float(oy), float(ox)], offset_norm_px=float(np.hypot(oy, ox)),
                   in_files=bool(np.hypot(oy, ox) >= IN_FILES_MIN_PX))
    else:
        out.update(offset_vs_naive_centre_px=None, offset_norm_px=float("nan"), in_files=False)
    return out
