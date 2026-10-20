"""Realistic synthetic model of the diagnosed cause (amendment 1): one jittered source mesh S, both published meshes made
from it by the exporter replica (float32 affine + float32 bilinear resample), in the four catalog configurations."""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from meshpair.exporter import F32, transform_tifxyz, apply_affine_f32       # noqa: E402
from meshpair.synth import curved_sheet, demo_affine, planted_resample      # noqa: E402

PXF, PXC = 2.399, 9.362
META = {"scale": [0.05, 0.05]}


def jittered_source(seed, rows=200, cols=240, spacing=20.0, jitter=1.0, z_shift=0.0):
    S = curved_sheet(rows=rows, cols=cols, spacing=spacing, radius=9000.0, holes=8, seed=seed)
    v = S[..., 0] != -1
    rng = np.random.default_rng(10_000 + seed)
    J = rng.standard_normal(S.shape) * jitter
    S2 = S.astype(np.float64) + J
    S2[..., 2] += z_shift
    S2 = S2.astype(F32); S2[~v] = F32(-1)
    return S2


def h44(A):
    return np.vstack([np.asarray(A, np.float64), [0, 0, 0, 1]])


CONFIGS = {
    # name: (A_F or None for copy, A_C or None for copy)
    "both_ratio4": (np.c_[4.0 * np.eye(3), np.zeros(3)], demo_affine(scale=1.0285, rot_deg=-1.3, t=(6677.0, 6689.0, 624.0))),
    "both_ratio1.95": (demo_affine(scale=1.949, rot_deg=0.4, t=(-18.0, 40.0, 25.0)), demo_affine(scale=0.4625, rot_deg=-1.1, t=(1200.0, 500.0, -300.0))),
    "fine_copy": (None, demo_affine(scale=1 / 3.88, rot_deg=0.4, t=(1210.0, 540.0, -3100.0))),
    "coarse_copy": (demo_affine(scale=3.8976, rot_deg=0.4, t=(340.0, 215.0, -45000.0)), None),
}


def make_pair(S, AF, AC, shift_px=(0.0, 0.0)):
    """Returns (F_pts, C_pts, M_c2f, fine_is_copy, coarse_is_copy). shift_px: planted in-plane shift in COARSE canvas px,
    planted in the mesh that is the QUERY of the lattice-clean route (coarse unless the coarse is the copy)."""
    I = np.c_[np.eye(3), np.zeros(3)]
    fine_copy, coarse_copy = AF is None, AC is None
    AFm = I if fine_copy else AF
    ACm = I if coarse_copy else AC
    F = S.copy() if fine_copy else transform_tifxyz(S, AFm)[0]
    if coarse_copy:
        C = S.copy()
    else:
        C, info = transform_tifxyz(S, ACm)
    if shift_px != (0.0, 0.0):
        if not coarse_copy:
            T = apply_affine_f32(S, ACm); nh, nw = C.shape[:2]
            C = planted_resample(T, nh, nw, shift_cells=(shift_px[0] * 0.05 * T.shape[0] / nh, shift_px[1] * 0.05 * T.shape[1] / nw))
        else:   # coarse = source: 1 coarse px = 0.05 source cells; plant in the fine (query of the forward route)
            T = apply_affine_f32(S, AFm); nh, nw = F.shape[:2]
            F = planted_resample(T, nh, nw, shift_cells=(-shift_px[0] * 0.05, -shift_px[1] * 0.05))
    M = (h44(AFm) @ np.linalg.inv(h44(ACm)))[:3]
    return F, C, M, fine_copy, coarse_copy
