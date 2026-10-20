# SPDX-License-Identifier: MIT
# Authors: Chris Thomas; clean-room implementation by Claude Code (Opus 5.5), implementer session 2026-10-18.
"""Slice centre of a surface volume and the default coarse-layer choice (helpers for ct.py).

clean-room implementation from cleanroom/SPEC.md v1 (2026-10-18), Part D.
"""
from __future__ import annotations


def slice_centre(n_layers: int) -> float:
    """Fractional layer index at the mesh surface: c(N) = (N - 1) / 2."""
    return (n_layers - 1) / 2


def layer_depth_um(k, n_layers: int, px_um: float, delta_um: float = 0.0):
    """Signed depth of layer k along the normal: ((k - c(N)) * px) + delta."""
    return ((k - slice_centre(n_layers)) * px_um) + delta_um


def depth_to_layer(d_um, px_um: float, n_layers0: int):
    """Fractional level-0 layer of a volume for depth d: (d / px) + c(N0)."""
    return (d_um / px_um) + slice_centre(n_layers0)


def coarse_layers(n_c: int, count=None) -> list:
    """The `count` layers nearest the surface (default 8 for even n_c, 7 for odd), ties -> smaller k, sorted."""
    if count is None:
        count = 8 if n_c % 2 == 0 else 7
    c = slice_centre(n_c)
    ranked = sorted(range(n_c), key=lambda k: (abs(k - c), k))
    return sorted(ranked[:count])
