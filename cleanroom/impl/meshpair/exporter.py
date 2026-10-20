# SPDX-License-Identifier: MIT
# Authors: Chris Thomas; clean-room implementation by Claude Code (Opus 5.5), implementer session 2026-10-18.
"""Python model of the tifxyz mesh transfer step (float32 affine, grid-spacing measure, size rule, half-pixel
bilinear resample that keeps invalids) and a mesh comparison.

clean-room implementation from cleanroom/SPEC.md v1 (2026-10-18), Part B.
"""
from __future__ import annotations

import math

import numpy as np

F32 = np.float32


def _minus1(arr):
    dt = arr.dtype
    return dt.type(-1) if np.issubdtype(dt, np.floating) else F32(-1)


def _point_valid(points):
    """Point validity (SPEC B.1): x != -1 and all three components finite."""
    m1 = _minus1(points)
    return (points[..., 0] != m1) & np.isfinite(points).all(axis=-1)


def _sample_valid(points):
    """Sample validity (SPEC B.1): no component == -1 and all three components finite."""
    m1 = _minus1(points)
    return (points != m1).all(axis=-1) & np.isfinite(points).all(axis=-1)


def _median_cpp(values):
    v = np.sort(np.asarray(values, dtype=np.float64))
    n = v.size
    if n == 0:
        return float("nan")
    m = v[n // 2]
    if n % 2 == 0:
        m = 0.5 * (m + v[n // 2 - 1])
    return float(m)


def _dist(a, b):
    return np.linalg.norm(a.astype(np.float64) - b.astype(np.float64), axis=-1)


def measure_grid_axis_spacing(points):
    """(med_h, med_v, n_h, n_v): median neighbour distance along cols (h) and rows (v)."""
    points = np.asarray(points)
    rows, cols = points.shape[:2]
    if rows < 2 or cols < 2:
        return (float("nan"), float("nan"), 0, 0)
    if rows < 20 or cols < 20:
        rs = np.arange(1, rows)
        cs = np.arange(1, cols)
    else:
        r0 = max(1, int(rows * 0.1))
        r1 = int(rows * 0.9)
        c0 = max(1, int(cols * 0.1))
        c1 = int(cols * 0.9)
        r1 = max(r0 + 1, r1)
        c1 = max(c0 + 1, c1)
        rs = np.arange(r0, r1, 4)
        cs = np.arange(c0, c1, 4)
    R, C = np.meshgrid(rs, cs, indexing="ij")
    R = R.ravel()
    C = C.ravel()
    valid = _point_valid(points)
    out = []
    for (dr, dc) in ((0, 1), (1, 0)):
        ok = valid[R, C] & valid[R - dr, C - dc]
        d = _dist(points[R[ok], C[ok]], points[R[ok] - dr, C[ok] - dc])
        d = d[np.isfinite(d) & (d > 0)]
        out.append(d)
    dh, dv = out
    return (_median_cpp(dh), _median_cpp(dv), int(dh.size), int(dv.size))


def apply_affine_f32(points, A, scale_before=1.0, scale_after=1.0):
    """Apply a 3x4 affine [L | t] to the valid vertices in float32/float64 as SPEC B.4; invalid vertices are copied."""
    out = np.array(points, dtype=np.float32, copy=True)
    A = np.asarray(A, dtype=np.float64)
    L, t = A[:, :3], A[:, 3]
    V = _point_valid(out)
    p = out[V] * np.float32(scale_before)
    q = p.astype(np.float64) @ L.T + t
    q = q.astype(np.float32) * np.float32(scale_after)
    bad = ~np.isfinite(q).all(axis=-1) | (q[:, 0] == np.float32(-1))
    q[bad] = np.float32(-1)
    out[V] = q
    return out


def _axis_taps(n_src, n_dst):
    scale = np.float32(n_src) / np.float32(n_dst)
    d = np.arange(n_dst, dtype=np.float32)
    s = (d + np.float32(0.5)) * scale - np.float32(0.5)
    s = np.clip(s, np.float32(0), np.float32(n_src - 1))
    i0 = np.floor(s).astype(np.int64)
    i1 = np.minimum(i0 + 1, n_src - 1)
    f = s - i0.astype(np.float32)
    return i0, i1, f


def resample_linear_preserving_invalids(points, new_h, new_w):
    """Half-pixel bilinear resample (SPEC B.5); a target vertex is written only if all 4 corners are valid samples."""
    points = np.asarray(points)
    rows, cols = points.shape[:2]
    out = np.full((new_h, new_w, 3), np.float32(-1), dtype=np.float32)
    y0, y1, fy = _axis_taps(rows, new_h)
    x0, x1, fx = _axis_taps(cols, new_w)
    one = np.float32(1)
    p00 = points[y0[:, None], x0[None, :]]
    p01 = points[y0[:, None], x1[None, :]]
    p10 = points[y1[:, None], x0[None, :]]
    p11 = points[y1[:, None], x1[None, :]]
    wx1 = fx[None, :, None]
    wx0 = (one - fx)[None, :, None]
    wy1 = fy[:, None, None]
    wy0 = (one - fy)[:, None, None]
    top = (p00 * wx0) + (p01 * wx1)
    bot = (p10 * wx0) + (p11 * wx1)
    res = ((top * wy0) + (bot * wy1)).astype(np.float32)
    ok = _sample_valid(p00) & _sample_valid(p01) & _sample_valid(p10) & _sample_valid(p11)
    out[ok] = res[ok]
    return out


def transform_tifxyz(points, A, scale_before=1.0, scale_after=1.0):
    """Affine + spacing-ratio size rule + resample (SPEC B.6). Returns (out float32, info dict)."""
    base = measure_grid_axis_spacing(points)
    A_use = np.c_[np.eye(3), np.zeros(3)] if A is None else A
    out = apply_affine_f32(points, A_use, scale_before, scale_after)
    after = measure_grid_axis_spacing(out)

    def _ratio(i, ni):
        if base[ni] > 0 and after[ni] > 0 and base[i] > 0:
            return after[i] / base[i]
        return float("nan")

    ratio_x = _ratio(0, 2)
    ratio_y = _ratio(1, 3)
    fb = scale_before * scale_after
    if A is not None:
        sv = np.linalg.svd(np.asarray(A, dtype=np.float64)[:, :3], compute_uv=False)
        mu = np.mean(sv)
        if np.max(np.abs(sv - mu)) / mu <= 1e-4:
            fb = fb * mu
    if not (math.isfinite(ratio_x) and ratio_x > 0):
        ratio_x = fb
    if not (math.isfinite(ratio_y) and ratio_y > 0):
        ratio_y = fb
    info = {"spacing_before": (base[0], base[1]), "spacing_after": (after[0], after[1]),
            "ratio_x": ratio_x, "ratio_y": ratio_y, "resampled": False}
    if abs(ratio_x - 1.0) > 1e-4 or abs(ratio_y - 1.0) > 1e-4:
        fx, fy = np.float32(ratio_x), np.float32(ratio_y)
        if abs(fx - np.float32(1)) < np.float32(0.001) and abs(fy - np.float32(1)) < np.float32(0.001):
            return out, info
        rows, cols = out.shape[:2]
        new_w = max(1, int(np.round(np.float32(cols) * fx)))
        new_h = max(1, int(np.round(np.float32(rows) * fy)))
        out = resample_linear_preserving_invalids(out, new_h, new_w)
        info["resampled"] = True
        info["new_shape"] = (new_h, new_w)
    return out, info


def compare_meshes(pred, pub) -> dict:
    pred = np.asarray(pred)
    pub = np.asarray(pub)
    r = {"pred_shape": list(pred.shape[:2]), "pub_shape": list(pub.shape[:2]),
         "same_shape": pred.shape == pub.shape}
    if not r["same_shape"]:
        return r
    vp, vq = _point_valid(pred), _point_valid(pub)
    both = vp & vq
    e = _dist(pred[both], pub[both])
    mism = int((vp ^ vq).sum())
    r["n_valid_pub"] = int(vq.sum())
    r["n_valid_pred"] = int(vp.sum())
    r["n_both"] = int(both.sum())
    r["mask_mismatch"] = mism
    r["mask_agreement"] = float(1 - mism / max(1, int((vp | vq).sum())))
    if e.size:
        r["err_median_vx"] = float(np.median(e))
        r["err_p99_vx"] = float(np.percentile(e, 99))
        r["err_max_vx"] = float(e.max())
    else:
        r["err_median_vx"] = r["err_p99_vx"] = r["err_max_vx"] = float("nan")
    return r
