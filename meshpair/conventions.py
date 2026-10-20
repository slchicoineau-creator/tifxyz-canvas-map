# SPDX-License-Identifier: MIT
# Authors: Chris Thomas; clean-room implementation by Claude Code (Opus 5.5), implementer session 2026-10-18.
"""Renderer canvas <-> grid conventions and the exporter's half-pixel index map.

clean-room implementation from cleanroom/SPEC.md v1 (2026-10-18), Part A.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


def _lround(x: float) -> int:
    """Round half away from zero, evaluated as floor(x + 0.5) / -floor(-x + 0.5) in float64 (SPEC §0)."""
    x = float(x)
    if x >= 0:
        return int(math.floor(x + 0.5))
    return -int(math.floor(-x + 0.5))


def _f64(v):
    return np.asarray(v, dtype=np.float64)


@dataclass(frozen=True)
class MeshGeom:
    rows: int
    cols: int
    scale_x: float = 0.05
    scale_y: float = 0.05
    render_scale: float = 1.0
    px_um: float = float("nan")

    @property
    def canvas_w(self) -> int:
        return max(1, _lround((self.cols * self.render_scale) / self.scale_x))

    @property
    def canvas_h(self) -> int:
        return max(1, _lround((self.rows * self.render_scale) / self.scale_y))

    def _params(self, axis: int):
        if axis == 0:
            return self.rows, self.scale_y, self.canvas_h
        if axis == 1:
            return self.cols, self.scale_x, self.canvas_w
        raise ValueError(f"axis must be 0 or 1, got {axis!r}")

    def canvas_to_grid(self, p, axis: int):
        n, s, W = self._params(axis)
        p = _f64(p)
        R = self.render_scale
        return ((s * ((p - 0.5 * W) + 0.5)) / R) + (n / 2.0)

    def grid_to_canvas(self, g, axis: int):
        n, s, W = self._params(axis)
        g = _f64(g)
        R = self.render_scale
        return ((((g - n / 2.0) * R) / s) + 0.5 * W) - 0.5


def exporter_grid_map(g_f, n_f, n_c):
    q = n_c / n_f
    return (_f64(g_f) + 0.5) * q - 0.5


def exporter_grid_map_inv(g_c, n_f, n_c):
    q = n_f / n_c
    return (_f64(g_c) + 0.5) * q - 0.5


def _n(geom: MeshGeom, axis: int) -> int:
    return geom.rows if axis == 0 else geom.cols


def exact_canvas_map(p_f, fine: MeshGeom, coarse: MeshGeom, axis: int):
    g_f = fine.canvas_to_grid(p_f, axis)
    g_c = exporter_grid_map(g_f, _n(fine, axis), _n(coarse, axis))
    return coarse.grid_to_canvas(g_c, axis)


def exact_canvas_map_inv(p_c, fine: MeshGeom, coarse: MeshGeom, axis: int):
    g_c = coarse.canvas_to_grid(p_c, axis)
    g_f = exporter_grid_map_inv(g_c, _n(fine, axis), _n(coarse, axis))
    return fine.grid_to_canvas(g_f, axis)


def naive_canvas_map(p_f, k):
    return _f64(p_f) * k


def exact_linear_coeffs(fine: MeshGeom, coarse: MeshGeom, axis: int):
    """(a, b) of p_C = a + b p_F, evaluated from the map at p_F = 0 and 1 (not a closed form)."""
    p0 = float(exact_canvas_map(0.0, fine, coarse, axis))
    p1 = float(exact_canvas_map(1.0, fine, coarse, axis))
    return p0, p1 - p0


def offset_summary(fine: MeshGeom, coarse: MeshGeom, k=None, valid_box=None) -> dict:
    if k is None:
        k = fine.px_um / coarse.px_um
    if valid_box is None:
        valid_box = (0, fine.canvas_h - 1, 0, fine.canvas_w - 1)
    y0, y1, x0, x1 = valid_box
    out = {"k_naive": k}
    for axis, pre, lo, hi in ((0, "y", y0, y1), (1, "x", x0, x1)):
        a, b = exact_linear_coeffs(fine, coarse, axis)
        ends = np.array([lo, hi], dtype=np.float64)
        off = (a + b * ends) - k * ends
        mid = 0.5 * (lo + hi)
        out[f"{pre}_a"] = a
        out[f"{pre}_b"] = b
        out[f"{pre}_b_minus_k"] = float(b - k)
        out[f"{pre}_offset_at_origin_px"] = a
        out[f"{pre}_offset_box_min_px"] = float(off.min())
        out[f"{pre}_offset_box_max_px"] = float(off.max())
        out[f"{pre}_offset_box_centre_px"] = float((a + b * mid) - k * mid)
    if math.isfinite(coarse.px_um):
        for name in [n for n in out if n.endswith("_px")]:
            out[name[:-3] + "_um"] = float(out[name] * coarse.px_um)
    return out
