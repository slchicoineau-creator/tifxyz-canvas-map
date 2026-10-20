"""Model-free index correspondence between a fine and a coarse mesh of one segment (mesh route; no CT).

Route I  (inverse map, uses the published affine M: coarse-volume coords -> fine-volume coords):
    for sampled fine vertices X_F, Y = M^-1 X_F, then find the coarse grid position g_C whose bilinear surface point is
    closest to Y (vectorised Gauss-Newton on the bilinear patch). Gives a per-vertex correspondence p_C(p_F) in canvas
    px, without assuming any resampling convention, plus the normal (depth) distance |C(g_C) - Y|.
Route F  (free affine; does not use the published affine):
    per-axis linear canvas map p_C = a + b p_F and a free 3-D affine fitted jointly by variable projection.
    Identifiable only through the sheet's curvature; reported as a cross-check of route I and of the catalog affine.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import binary_erosion
from scipy.optimize import least_squares

from .conventions import MeshGeom


def bilinear(points: np.ndarray, valid: np.ndarray, gy: np.ndarray, gx: np.ndarray):
    """Bilinear sample of a (rows, cols, 3) grid at fractional (gy, gx); ok = all 4 corners valid and inside."""
    rows, cols = points.shape[:2]
    y0 = np.floor(gy).astype(np.int64)
    x0 = np.floor(gx).astype(np.int64)
    inside = (y0 >= 0) & (x0 >= 0) & (y0 + 1 <= rows - 1) & (x0 + 1 <= cols - 1)
    y0c = np.clip(y0, 0, rows - 2)
    x0c = np.clip(x0, 0, cols - 2)
    fy = (gy - y0c)[:, None]
    fx = (gx - x0c)[:, None]
    P = points.astype(np.float64)
    p00, p01 = P[y0c, x0c], P[y0c, x0c + 1]
    p10, p11 = P[y0c + 1, x0c], P[y0c + 1, x0c + 1]
    ok = inside & valid[y0c, x0c] & valid[y0c, x0c + 1] & valid[y0c + 1, x0c] & valid[y0c + 1, x0c + 1]
    val = (p00 * (1 - fx) + p01 * fx) * (1 - fy) + (p10 * (1 - fx) + p11 * fx) * fy
    du = (p01 - p00) * (1 - fy) + (p11 - p10) * fy        # d/dgx
    dv = (p10 - p00) * (1 - fx) + (p11 - p01) * fx        # d/dgy
    return val, du, dv, ok


def sample_fine_vertices(fine_pts: np.ndarray, fine_valid: np.ndarray, n: int, seed: int = 0, erode: int = 2):
    v = binary_erosion(fine_valid, iterations=erode) if erode else fine_valid
    rr, cc = np.nonzero(v)
    rng = np.random.default_rng(seed)
    if rr.size > n:
        pick = rng.choice(rr.size, n, replace=False)
        rr, cc = rr[pick], cc[pick]
    return rr, cc


def _gn_project(coarse_pts, coarse_valid, Y, gy, gx, iters):
    for _ in range(iters):
        val, du, dv, ok = bilinear(coarse_pts, coarse_valid, gy, gx)
        r = Y - val
        a11 = (du * du).sum(1); a12 = (du * dv).sum(1); a22 = (dv * dv).sum(1)
        b1 = (du * r).sum(1); b2 = (dv * r).sum(1)
        det = a11 * a22 - a12 * a12
        det = np.where(np.abs(det) < 1e-12, np.nan, det)
        dgx = np.clip(np.nan_to_num((a22 * b1 - a12 * b2) / det), -2, 2)
        dgy = np.clip(np.nan_to_num((a11 * b2 - a12 * b1) / det), -2, 2)
        gx = gx + dgx
        gy = gy + dgy
    val, du, dv, ok = bilinear(coarse_pts, coarse_valid, gy, gx)
    r = Y - val
    nrm = np.cross(du, dv)
    nrm /= np.maximum(np.linalg.norm(nrm, axis=1, keepdims=True), 1e-12)
    normal_d = (r * nrm).sum(1)
    tang = np.linalg.norm(r - normal_d[:, None] * nrm, axis=1)
    return gy, gx, normal_d, tang, ok


def inverse_map(fine_pts, fine_valid, coarse_pts, coarse_valid, fine: MeshGeom, coarse: MeshGeom, M: np.ndarray,
                init_canvas_map, n: int = 60000, seed: int = 0, iters: int = 12, on_surface_vx: float = 0.05,
                restart_cells=(-1, 0, 1)):
    """Route I. M: 3x4 affine coarse-volume -> fine-volume. init_canvas_map(py, px) -> (cy, cx) start guess.
    A sample is 'ok' when the projection converged tangentially (< 1e-3 coarse vx) AND lies on the coarse surface
    (|normal distance| <= on_surface_vx). Samples failing from the start guess are retried from start guesses shifted by
    whole coarse cells (restart_cells, model-free), keeping the smallest residual. Returns per-sample arrays."""
    rr, cc = sample_fine_vertices(fine_pts, fine_valid, n, seed)
    py = fine.grid_to_canvas(rr.astype(np.float64), 0)
    px = fine.grid_to_canvas(cc.astype(np.float64), 1)
    X = fine_pts[rr, cc].astype(np.float64)
    L = np.asarray(M, np.float64)[:, :3]
    t = np.asarray(M, np.float64)[:, 3]
    Y = np.linalg.solve(L, (X - t).T).T
    cy0, cx0 = init_canvas_map(py, px)
    gy0 = coarse.canvas_to_grid(cy0, 0)
    gx0 = coarse.canvas_to_grid(cx0, 1)
    gy, gx, nd, tg, ok = _gn_project(coarse_pts, coarse_valid, Y, gy0, gx0, iters)
    good = ok & (tg < 1e-3) & (np.abs(nd) <= on_surface_vx)
    n_first = int(good.sum())
    best = np.where(ok, np.hypot(nd, tg), np.inf)
    for dy in restart_cells:
        for dx in restart_cells:
            if dy == 0 and dx == 0:
                continue
            idx = np.nonzero(~good)[0]
            if idx.size == 0:
                break
            gy2, gx2, nd2, tg2, ok2 = _gn_project(coarse_pts, coarse_valid, Y[idx], gy0[idx] + dy, gx0[idx] + dx, iters)
            sc = np.where(ok2, np.hypot(nd2, tg2), np.inf)
            better = sc < best[idx]
            j = idx[better]
            gy[j], gx[j], nd[j], tg[j], ok[j], best[j] = gy2[better], gx2[better], nd2[better], tg2[better], ok2[better], sc[better]
            good = ok & (tg < 1e-3) & (np.abs(nd) <= on_surface_vx)
    return dict(py=py, px=px, cy=coarse.grid_to_canvas(gy, 0), cx=coarse.grid_to_canvas(gx, 1),
                normal_d=nd, tang_resid=tg, ok=good, n_ok_first_start=n_first,
                tangential_converged=ok & (tg < 1e-3))


def linear_fit_axis(p_f: np.ndarray, p_c: np.ndarray):
    A = np.c_[np.ones_like(p_f), p_f]
    coef, *_ = np.linalg.lstsq(A, p_c, rcond=None)
    res = p_c - A @ coef
    return float(coef[0]), float(coef[1]), res


def split_half_noise(py, px, cy, cx):
    """/2 rule: fit (a, b) per axis on left/right and top/bottom halves; noise of the whole-mesh value at the canvas
    centre = |half_a - half_b| / 2 (each half carries half the samples)."""
    out = {}
    yc, xc = np.median(py), np.median(px)
    for name, m in (("lr", px < xc), ("tb", py < yc)):
        vals = []
        for sel in (m, ~m):
            ay, by, _ = linear_fit_axis(py[sel], cy[sel])
            ax, bx, _ = linear_fit_axis(px[sel], cx[sel])
            vals.append((ay + by * yc, ax + bx * xc))
        d = np.subtract(vals[0], vals[1])
        out[name] = float(np.linalg.norm(d) / 2)
    return out


def free_affine_fit(fine_pts, fine_valid, coarse_pts, coarse_valid, fine: MeshGeom, coarse: MeshGeom, init_ab,
                    n: int = 30000, seed: int = 1):
    """Route F. init_ab = (ay, by, ax, bx). Returns fitted (ay, by, ax, bx), the 3x4 affine (coarse -> fine) and
    residual stats (fine vx). Samples kept fixed: those whose start position lies >= 3 cells inside valid coarse."""
    rr, cc = sample_fine_vertices(fine_pts, fine_valid, n, seed)
    py = fine.grid_to_canvas(rr.astype(np.float64), 0)
    px = fine.grid_to_canvas(cc.astype(np.float64), 1)
    X = fine_pts[rr, cc].astype(np.float64)
    cv_er = binary_erosion(coarse_valid, iterations=3)
    ay, by, ax, bx = init_ab
    gy0 = coarse.canvas_to_grid(ay + by * py, 0)
    gx0 = coarse.canvas_to_grid(ax + bx * px, 1)
    iy = np.clip(np.rint(gy0).astype(int), 0, coarse_pts.shape[0] - 1)
    ix = np.clip(np.rint(gx0).astype(int), 0, coarse_pts.shape[1] - 1)
    keep = cv_er[iy, ix]
    py, px, X = py[keep], px[keep], X[keep]
    sc = np.array([1.0, 1e-3, 1.0, 1e-3])      # parameter scaling (b is ~0.25 +- 1e-3)

    def resid(theta):
        ay, by, ax, bx = theta * sc
        gy = coarse.canvas_to_grid(ay + by * py, 0)
        gx = coarse.canvas_to_grid(ax + bx * px, 1)
        Q, _, _, ok = bilinear(coarse_pts, coarse_valid, gy, gx)
        Q1 = np.c_[Q, np.ones(len(Q))]
        W, *_ = np.linalg.lstsq(Q1[ok], X[ok], rcond=None)
        r = np.zeros_like(X)
        r[ok] = X[ok] - Q1[ok] @ W
        return r.ravel()

    th0 = np.asarray(init_ab, np.float64) / sc
    sol = least_squares(resid, th0, x_scale="jac", diff_step=1e-6, xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=200)
    ay, by, ax, bx = sol.x * sc
    gy = coarse.canvas_to_grid(ay + by * py, 0)
    gx = coarse.canvas_to_grid(ax + bx * px, 1)
    Q, _, _, ok = bilinear(coarse_pts, coarse_valid, gy, gx)
    Q1 = np.c_[Q, np.ones(len(Q))]
    W, *_ = np.linalg.lstsq(Q1[ok], X[ok], rcond=None)
    e = np.linalg.norm(X[ok] - Q1[ok] @ W, axis=1)
    return dict(ay=ay, by=by, ax=ax, bx=bx, affine=W.T, n=int(ok.sum()), resid_median_vx=float(np.median(e)),
                resid_p95_vx=float(np.percentile(e, 95)), status=int(sol.status), nfev=int(sol.nfev))


def affine_point_deviation(M1: np.ndarray, M2: np.ndarray, pts: np.ndarray) -> dict:
    """Distance between two 3x4 affines applied to the same points (target-volume voxels)."""
    a = pts @ M1[:, :3].T + M1[:, 3]
    b = pts @ M2[:, :3].T + M2[:, 3]
    d = np.linalg.norm(a - b, axis=1)
    return dict(median=float(np.median(d)), p95=float(np.percentile(d, 95)), max=float(d.max()))
