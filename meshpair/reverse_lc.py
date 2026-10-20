"""Lattice-clean route (TIFXYZ-REVERSE amendment 1, prereg/2026-10-13_partner_tifxyz_reverse_amend1.md).

Why (reverse_route/DIAGNOSIS.md): both published meshes are exporter resamples of one source mesh S (bilinear per source
cell). Inside a reference-mesh cell that lies wholly within ONE source cell (on both axes) the reference bilinear patch IS
the source surface, so a query vertex (itself a bilinear sample of S) lies on it to float rounding. Inside a reference cell
that straddles a source grid line the source surface has a kink and the reference chord departs from it (the S1 off-surface
vertices: 99.9-100 % of them sit in straddling cells). So the original bars (0.05 coarse vx, 0.95, 0.003 px) are applied
only to query vertices whose FORMULA position falls in a lattice-clean reference cell. Bars are unchanged.

Direction: reverse (coarse vertex -> fine surface) unless the coarse mesh is a plain copy of the source; then every
reverse query sits on a source node, i.e. inside a straddling fine cell, and the forward direction (fine vertex -> coarse
surface, formula start) is the lattice-clean one (the coarse cells ARE source cells).
The frozen reverse.py is not changed.
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import binary_erosion

from . import conventions as cv
from .conventions import MeshGeom
from .reverse import _project, ON_SURFACE_COARSE_VX, INPLANE_P95_PX, FRAC_ON_MIN, TANG_CONV_VX, IN_FILES_MIN_PX
from .fit import linear_fit_axis

N_CLEAN_MIN = 100
EPS = 1e-6
F32 = np.float32


def ulp32(x) -> float:
    return float(np.spacing(F32(abs(float(x)))))


def float_floor_px(F_pts, F_valid, C_pts, C_valid, px_f, px_c, n_src, gC: MeshGeom, fine_resampled: bool, coarse_resampled: bool):
    """float32 resolution of the published files, in COARSE canvas px (~ coarse voxels along the sheet):
    u_coord = one ulp of the largest stored fine coordinate (converted to coarse vx) + one ulp of the largest coarse one;
    u_index = for each mesh made by the float32 resampler, one ulp of the largest source index (src = (dst+0.5)*scale-0.5
              is computed in float32), converted from source cells to coarse canvas px: x (n_C / n_S) / s_C."""
    mf = float(np.abs(F_pts[F_valid]).max()) if F_valid.any() else 0.0
    mc = float(np.abs(C_pts[C_valid]).max()) if C_valid.any() else 0.0
    u_coord = ulp32(mf) * (px_f / px_c) + ulp32(mc)
    R = gC.render_scale
    per = max(ulp32(n_src[0]) * (gC.rows / n_src[0]) * R / gC.scale_y, ulp32(n_src[1]) * (gC.cols / n_src[1]) * R / gC.scale_x)
    u_index = per * (int(fine_resampled) + int(coarse_resampled))
    return u_coord, u_index


def source_index(g, n_ref: int, n_src: int):
    """Exporter resample relation, reference grid -> source grid (one axis): s = (g + 0.5) n_src / n_ref - 0.5
    (the clean conventions module's exporter_grid_map)."""
    return cv.exporter_grid_map(g, n_ref, n_src)


def clean_axis(g, n_ref: int, n_src: int):
    """True where the reference cell [floor(g), floor(g)+1] contains no source grid line strictly inside."""
    j = np.floor(np.asarray(g, np.float64))
    s0 = source_index(j, n_ref, n_src)
    s1 = source_index(j + 1, n_ref, n_src)
    m = np.floor(s0 + EPS) + 1.0          # first integer strictly above s0
    return ~(m < s1 - EPS)


def lattice_clean_route(ref_pts, ref_valid, q_pts, q_valid, g_ref: MeshGeom, g_q: MeshGeom, M_q_to_ref,
                        n_src_ref_axes, px_ref: float, px_q: float, px_coarse: float, q_is_coarse: bool, inplane_bar: float = INPLANE_P95_PX,
                        n_max: int = 100000, seed: int = 0, erode: int = 2, iters: int = 12):
    """Query vertices of mesh q, mapped into the reference volume by M, projected (GN, formula start) onto the
    reference bilinear surface. Normal distance and in-plane residual reported in COARSE voxels / COARSE canvas px.
    n_src_ref_axes = (source rows, source cols) of the common source mesh (shape only)."""
    v = binary_erosion(q_valid, iterations=erode) if erode else q_valid
    rr, cc = np.nonzero(v)
    n_valid_eroded = int(rr.size)
    if rr.size > n_max:
        pick = np.random.default_rng(seed).choice(rr.size, n_max, replace=False)
        rr, cc = rr[pick], cc[pick]
    pQy = g_q.grid_to_canvas(rr.astype(np.float64), 0)
    pQx = g_q.grid_to_canvas(cc.astype(np.float64), 1)
    M = np.asarray(M_q_to_ref, np.float64)
    X = q_pts[rr, cc].astype(np.float64)
    Y = X @ M[:, :3].T + M[:, 3]
    # formula position on the reference canvas / grid
    if q_is_coarse:      # reference = fine
        pRy = cv.exact_canvas_map_inv(pQy, g_ref, g_q, 0); pRx = cv.exact_canvas_map_inv(pQx, g_ref, g_q, 1)
    else:                # reference = coarse
        pRy = cv.exact_canvas_map(pQy, g_q, g_ref, 0); pRx = cv.exact_canvas_map(pQx, g_q, g_ref, 1)
    gy0 = g_ref.canvas_to_grid(pRy, 0)
    gx0 = g_ref.canvas_to_grid(pRx, 1)
    clean = clean_axis(gy0, g_ref.rows, n_src_ref_axes[0]) & clean_axis(gx0, g_ref.cols, n_src_ref_axes[1])
    gy, gx, nd_r, tg, ok = _project(ref_pts, ref_valid, Y, gy0.copy(), gx0.copy(), iters)
    nd = nd_r * (px_ref / px_coarse)
    conv = ok & (tg < TANG_CONV_VX)     # same rule as reverse.py: tangential residual < 1e-3 reference-volume voxels
    on = conv & (np.abs(nd) <= ON_SURFACE_COARSE_VX)
    # in-plane residual in COARSE canvas px
    if q_is_coarse:
        ey = cv.exact_canvas_map(g_ref.grid_to_canvas(gy, 0), g_ref, g_q, 0) - pQy
        ex = cv.exact_canvas_map(g_ref.grid_to_canvas(gx, 1), g_ref, g_q, 1) - pQx
    else:
        ey = g_ref.grid_to_canvas(gy, 0) - pRy
        ex = g_ref.grid_to_canvas(gx, 1) - pRx
    inpl = np.hypot(ey, ex)
    n_clean = int(clean.sum())
    cc_ = conv & clean
    out = dict(direction="reverse" if q_is_coarse else "forward", inplane_bar_px=float(inplane_bar), n_valid_eroded=n_valid_eroded, n_sampled=int(rr.size),
               n_clean=n_clean, frac_clean=float(clean.mean()) if rr.size else float("nan"),
               n_clean_converged=int(cc_.sum()), n_clean_on=int((on & clean).sum()),
               frac_on_clean=float((on & clean).sum() / n_clean) if n_clean else float("nan"),
               frac_on_all=float(on.mean()) if rr.size else float("nan"))
    if cc_.any():
        out.update(median_normal_clean=float(np.median(np.abs(nd[cc_]))), p95_normal_clean=float(np.percentile(np.abs(nd[cc_]), 95)),
                   max_normal_clean=float(np.abs(nd[cc_]).max()),
                   inplane_p95_clean=float(np.percentile(inpl[cc_], 95)), inplane_median_clean=float(np.median(inpl[cc_])),
                   inplane_max_clean=float(inpl[cc_].max()))
    else:
        out.update(median_normal_clean=float("nan"), p95_normal_clean=float("nan"), max_normal_clean=float("nan"),
                   inplane_p95_clean=float("nan"), inplane_median_clean=float("nan"), inplane_max_clean=float("nan"))
    nst = conv & ~clean
    if nst.any():   # descriptive: the straddling (lattice-affected) vertices
        out.update(p95_normal_straddle=float(np.percentile(np.abs(nd[nst]), 95)),
                   inplane_p95_straddle=float(np.percentile(inpl[nst], 95)))
    # offset vs naive (same line as reverse.py), from the clean on-surface vertices: fine canvas <-> coarse canvas
    sel = on & clean
    if sel.sum() >= N_CLEAN_MIN:
        if q_is_coarse:
            pf_y, pf_x = g_ref.grid_to_canvas(gy, 0)[sel], g_ref.grid_to_canvas(gx, 1)[sel]; pc_y, pc_x = pQy[sel], pQx[sel]
            gFi, fvalid = g_ref, ref_valid
        else:
            pf_y, pf_x = pQy[sel], pQx[sel]; pc_y, pc_x = g_ref.grid_to_canvas(gy, 0)[sel], g_ref.grid_to_canvas(gx, 1)[sel]
            gFi, fvalid = g_q, q_valid
        k = (px_ref / px_q) if q_is_coarse else (px_q / px_ref)     # px_fine / px_coarse
        rF, cF = np.nonzero(fvalid)
        yc = 0.5 * (gFi.grid_to_canvas(rF.min(), 0) + gFi.grid_to_canvas(rF.max(), 0))
        xc = 0.5 * (gFi.grid_to_canvas(cF.min(), 1) + gFi.grid_to_canvas(cF.max(), 1))
        ay, by, _ = linear_fit_axis(pf_y, pc_y); ax, bx, _ = linear_fit_axis(pf_x, pc_x)
        oy, ox = (ay + by * yc) - k * yc, (ax + bx * xc) - k * xc
        out.update(offset_vs_naive_centre_px=[float(oy), float(ox)], offset_norm_px=float(np.hypot(oy, ox)),
                   in_files=bool(np.hypot(oy, ox) >= IN_FILES_MIN_PX))
    else:
        out.update(offset_vs_naive_centre_px=None, offset_norm_px=float("nan"), in_files=False)
    if n_clean < N_CLEAN_MIN:
        out["verdict"] = "NOT_EVALUABLE"
    else:
        ok_pair = (out["frac_on_clean"] >= FRAC_ON_MIN and np.isfinite(out["inplane_p95_clean"])
                   and out["inplane_p95_clean"] <= inplane_bar)
        out["verdict"] = "PASS" if ok_pair else "FAIL"
    return out


def run_pair(F_pts, F_valid, C_pts, C_valid, gF: MeshGeom, gC: MeshGeom, M_c2f, n_src, px_f, px_c, coarse_is_copy: bool,
             fine_is_copy: bool = False, **kw):
    """Choose the lattice-clean direction and run it. n_src = (rows, cols) of the common source mesh."""
    if fine_is_copy:
        n_src = (gF.rows, gF.cols) if n_src is None else n_src
    if coarse_is_copy:
        n_src = (gC.rows, gC.cols) if n_src is None else n_src
    u_coord, u_index = float_floor_px(F_pts, F_valid, C_pts, C_valid, px_f, px_c, n_src, gC, not fine_is_copy, not coarse_is_copy)
    bar = max(INPLANE_P95_PX, u_coord + u_index)
    if coarse_is_copy:
        Minv = np.linalg.inv(np.vstack([np.asarray(M_c2f, np.float64), [0, 0, 0, 1]]))[:3]
        out = lattice_clean_route(C_pts, C_valid, F_pts, F_valid, gC, gF, Minv, n_src, px_c, px_f, px_c, False, bar, **kw)
    else:
        out = lattice_clean_route(F_pts, F_valid, C_pts, C_valid, gF, gC, M_c2f, n_src, px_f, px_c, px_c, True, bar, **kw)
    out.update(u_coord_px=u_coord, u_index_px=u_index, n_src=[int(n_src[0]), int(n_src[1])])
    return out
