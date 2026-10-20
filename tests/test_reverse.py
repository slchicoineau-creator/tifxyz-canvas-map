"""S0 self-tests of the reverse route (synthetic only). prereg/2026-10-13_partner_tifxyz_reverse_route.md"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from meshpair.exporter import F32, transform_tifxyz                        # noqa: E402
from meshpair.reverse import reverse_route                                  # noqa: E402
from meshpair.synth import curved_sheet, planted_resample, demo_affine      # noqa: E402
from meshpair.tifxyz import Mesh                                            # noqa: E402

PXF, PXC = 2.403, 9.366
META = {"scale": [0.05, 0.05]}


def A_coarse():
    """~1/3.9 affine (source = fine volume -> coarse volume), mild rotation/anisotropy, shift."""
    A = demo_affine(scale=1 / 3.8976, rot_deg=0.4, t=(1210.0, 540.0, -3100.0))
    return A


def fine_is_source(seed, shift_px=(0.0, 0.0)):
    S = curved_sheet(rows=300, cols=360, spacing=20.0, radius=9000.0, holes=10, seed=seed)
    A = A_coarse()
    Cpts, info = transform_tifxyz(S, A)
    if shift_px != (0.0, 0.0):
        from meshpair.exporter import apply_affine_f32
        T = apply_affine_f32(S, A)
        nh, nw = info["new_shape"]
        # shift in OUTPUT cells = px * 0.05 (grid = s * canvas px); in source cells x (n_old / n_new)
        sh = (shift_px[0] * 0.05 * (T.shape[0] / nh), shift_px[1] * 0.05 * (T.shape[1] / nw))
        Cpts = planted_resample(T, nh, nw, shift_cells=sh)
    M = np.linalg.inv(np.vstack([A, [0, 0, 0, 1]]))[:3]          # coarse -> fine volume
    return Mesh(S, META), Mesh(Cpts, META), M


def run(F, C, M):
    return reverse_route(F.points, F.valid, C.points, C.valid, F.geom(PXF), C.geom(PXC), M, PXF, PXC)


def test_fine_is_source_normal_distance_and_pass():
    F, C, M = fine_is_source(0)
    r = run(F, C, M)
    assert r["frac_on"] > 0.99
    assert r["p95_normal_coarse_vx"] <= 0.01
    assert r["pair_pass"], r


@pytest.mark.parametrize("seed", range(0, 100, 7))
def test_normal_p95_below_gate(seed):
    F, C, M = fine_is_source(seed)
    assert run(F, C, M)["p95_normal_coarse_vx"] <= 0.01


@pytest.mark.parametrize("shift", [0.05, 0.25, 1.0])
def test_planted_inplane_shifts_fail_pair_rule(shift):
    fails = 0
    n = 20
    for seed in range(n):
        d = [(shift, 0.0), (0.0, shift), (shift / np.sqrt(2), shift / np.sqrt(2))][seed % 3]
        F, C, M = fine_is_source(seed, shift_px=d)
        r = run(F, C, M)
        fails += int(not r["pair_pass"])
    assert fails / n >= 0.95


def test_readable_orientation_coarse_is_source_has_offset():
    """coarse = source (copy), fine = 3.9x resample: the readable configuration; offset >= 3 px must be reported."""
    S = curved_sheet(rows=110, cols=130, seed=1)
    A = demo_affine()
    Fpts, info = transform_tifxyz(S, A)
    r = run(Mesh(Fpts, META), Mesh(S, META), A)
    assert r["offset_norm_px"] >= 3.0 and r["in_files"]
