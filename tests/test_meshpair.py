"""Self-tests on synthetic data (no network, no real scroll data). Run: .venv/bin/python -m pytest -q tests"""
import os
import sys

import numpy as np
import pytest
from scipy.ndimage import map_coordinates

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from meshpair import conventions as cv                                    # noqa: E402
from meshpair.catalog import pairs, coarse_to_fine_affine                 # noqa: E402
from meshpair.ct import CanvasMap, TileSpec, _peak, array_sv, measure_tile, fit_field  # noqa: E402
from meshpair.exporter import transform_tifxyz, resample_linear_preserving_invalids, compare_meshes  # noqa: E402
from meshpair.fit import free_affine_fit, inverse_map                     # noqa: E402
from meshpair.meshcheck import check_pair                                 # noqa: E402
from meshpair.synth import block_mean, curved_sheet, demo_affine, planted_resample, synth_sv_pair  # noqa: E402
from meshpair.tifxyz import Mesh                                          # noqa: E402

PXF, PXC = 2.403, 9.366
K = PXF / PXC


# ---------------------------------------------------------------- conventions
def test_canvas_grid_roundtrip_and_renderer_rule():
    g = cv.MeshGeom(rows=211, cols=238, px_um=PXC)
    p = np.array([0.0, 17.3, 4759.0])
    for ax in (0, 1):
        assert np.allclose(g.grid_to_canvas(g.canvas_to_grid(p, ax), ax), p)
    # renderer rule with R = 1, s = 0.05: grid = s (p + 0.5)
    assert np.allclose(g.canvas_to_grid(p, 1), 0.05 * (p + 0.5))


def test_exact_map_inverse_and_w00_numbers():
    F = cv.MeshGeom(rows=823, cols=928, px_um=PXF)
    C = cv.MeshGeom(rows=211, cols=238, px_um=PXC)
    p = np.array([0.0, 5000.0, 16000.0])
    for ax in (0, 1):
        q = cv.exact_canvas_map(p, F, C, ax)
        assert np.allclose(cv.exact_canvas_map_inv(q, F, C, ax), p)
    a, b = cv.exact_linear_coeffs(F, C, 0)
    # p_C = (p_F + 10.5) rho - 10.5 with rho = 211/823
    assert abs(b - 211 / 823) < 1e-12 and abs(a - (10.5 * 211 / 823 - 10.5)) < 1e-9


# ---------------------------------------------------------------- exporter replica
def test_resampler_matches_independent_float64_halfpixel_bilinear():
    S = curved_sheet(rows=60, cols=70, holes=0, seed=3)
    S[:2] = S[2]          # fill the invalid top rows (no invalids for this comparison)
    S[:, :4] = S[:, 4:5]
    out = resample_linear_preserving_invalids(S, 233, 271)
    gy = np.clip((np.arange(233) + 0.5) * (60 / 233) - 0.5, 0, 59)
    gx = np.clip((np.arange(271) + 0.5) * (70 / 271) - 0.5, 0, 69)
    GY, GX = np.meshgrid(gy, gx, indexing="ij")
    ref = np.stack([map_coordinates(S[..., i].astype(np.float64), [GY, GX], order=1) for i in range(3)], -1)
    assert np.abs(out - ref).max() < 2e-3          # float32 vs float64 on coordinates ~1e4


def test_transform_size_rule_and_invalid_propagation():
    S = curved_sheet(rows=80, cols=90, seed=1)
    A = demo_affine()
    out, info = transform_tifxyz(S, A)
    assert info["resampled"]
    h, w = info["new_shape"]
    assert h == int(np.round(np.float32(80) * np.float32(info["ratio_y"])))
    assert w == int(np.round(np.float32(90) * np.float32(info["ratio_x"])))
    # identity affine: no resample, points unchanged
    out2, info2 = transform_tifxyz(S, np.c_[np.eye(3), np.zeros(3)])
    assert not info2["resampled"] and np.array_equal(out2, S)
    # an invalid source vertex invalidates every target vertex whose 4 corners touch it
    assert (out[..., 0] == -1).sum() > (S[..., 0] == -1).sum() * (h * w) / (80 * 90) * 0.9


# ---------------------------------------------------------------- mesh routes
def _pair(shift_cells=(0.0, 0.0), seed=0):
    S = curved_sheet(rows=110, cols=130, seed=seed)
    A = demo_affine()
    Fpts, info = transform_tifxyz(S, A)
    if shift_cells != (0.0, 0.0):
        from meshpair.exporter import apply_affine_f32
        Fpts = planted_resample(apply_affine_f32(S, A), *info["new_shape"], shift_cells=shift_cells)
    meta = {"scale": [0.05, 0.05]}
    return Mesh(Fpts, meta), Mesh(S, meta), A


def test_route_I_recovers_exporter_formula_without_knowing_it():
    Fm, Cm, A = _pair()
    rec = check_pair(Fm, Cm, PXF, PXC, A, n_inv=20000, do_free=False)
    ri = rec["route_I"]
    assert ri["frac_ok"] > 0.99
    assert ri["vs_exporter_model_px"]["p95"] < 0.01
    em = rec["exporter_model"]
    assert abs(ri["offset_vs_naive_centre_px"][0] - em["y_offset_box_centre_px"]) < 0.02
    assert abs(ri["offset_vs_naive_centre_px"][1] - em["x_offset_box_centre_px"]) < 0.02
    assert rec["lines"]["M1_in_files"] == "IN_FILES"          # the half-pixel rule alone gives ~7.8 px


@pytest.mark.parametrize("plant", [(0.3, -0.6), (-0.45, 0.2)])
def test_route_I_and_F_recover_planted_index_shift(plant):
    Fm, Cm, A = _pair(shift_cells=plant)
    rec = check_pair(Fm, Cm, PXF, PXC, A, n_inv=20000, n_free=15000, do_free=True)
    gF, gC = Fm.geom(PXF), Cm.geom(PXC)
    # a planted source-index shift of d cells moves the coarse canvas position by d / s = 20 d px
    want = np.array(plant) * 20.0
    ri = rec["route_I"]
    yc = 0.5 * (rec["fine_canvas_valid_box"][0] + rec["fine_canvas_valid_box"][1])
    xc = 0.5 * (rec["fine_canvas_valid_box"][2] + rec["fine_canvas_valid_box"][3])
    got_I = np.array([ri["ay"] + ri["by"] * yc - cv.exact_canvas_map(yc, gF, gC, 0),
                      ri["ax"] + ri["bx"] * xc - cv.exact_canvas_map(xc, gF, gC, 1)])
    assert np.abs(got_I - want).max() < 0.02
    rf = rec["route_F"]
    got_F = np.array([rf["ay"] + rf["by"] * yc - cv.exact_canvas_map(yc, gF, gC, 0),
                      rf["ax"] + rf["bx"] * xc - cv.exact_canvas_map(xc, gF, gC, 1)])
    assert np.abs(got_F - want).max() < 0.05
    assert not rec["lines"]["M2b_model_matches_route_I"]       # the planted file must NOT look like the exporter


def test_reproduction_flags_a_planted_file():
    Fm, Cm, A = _pair(shift_cells=(0.25, 0.0))
    rec = check_pair(Fm, Cm, PXF, PXC, A, source=Cm, A_src_to_fine=A, coarse_is_copy=True, n_inv=5000, do_free=False)
    assert rec["lines"]["M2a_reproduced"] is False
    Fm, Cm, A = _pair()
    rec = check_pair(Fm, Cm, PXF, PXC, A, source=Cm, A_src_to_fine=A, coarse_is_copy=True, n_inv=5000, do_free=False)
    assert rec["lines"]["M2a_reproduced"] is True


# ---------------------------------------------------------------- CT stage
W00_LIKE = (-7.808, 211 / 823, -7.807, 238 / 928)


@pytest.mark.parametrize("d,dep", [((0.0, 0.0), 0.0), ((1.3, -0.7), 4.8), ((-2.5, 3.2), -7.0)])
def test_ct_recovers_planted_shift_and_depth(d, dep):
    m = W00_LIKE
    tm = (m[0] + d[0], m[1], m[2] + d[1], m[3])
    f, c = synth_sv_pair(tm, Hc=140, Wc=140, n_c=12, depth_true_um=dep, seed=2)
    spec = TileSpec(T=96, S=10, depth_range_um=19.0)
    r = measure_tile(array_sv(f, PXF), array_sv(c, PXC), CanvasMap(*m), (70, 70), spec)
    assert r["valid"]
    assert np.abs(np.subtract(r["map_correction_px"], d)).max() < 0.2
    assert abs(r["depth_offset_um"] - dep) < 0.2 * PXC          # 0.2 coarse layer
    # pyramid level 2 (block mean) gives the same answer
    r2 = measure_tile(array_sv(block_mean(f, 4), PXF, level=2, n_layers0=f.shape[0]), array_sv(c, PXC),
                      CanvasMap(*m), (70, 70), spec)
    assert np.abs(np.subtract(r2["map_correction_px"], d)).max() < 0.2


def test_ct_differential_plant_and_naive_map():
    m = W00_LIKE
    f, c = synth_sv_pair(m, Hc=140, Wc=140, n_c=12, seed=4)
    F, C = array_sv(f, PXF), array_sv(c, PXC)
    spec = TileSpec(T=96, S=12, depth_range_um=10.0)
    base = measure_tile(F, C, CanvasMap(*m), (70, 70), spec)
    for pl in ((4.0, 0.0), (0.0, -8.0)):
        r = measure_tile(F, C, CanvasMap(*m).planted(*pl), (70, 70), spec)
        delta = np.subtract(r["map_correction_px"], base["map_correction_px"])
        assert np.abs(delta - np.array(pl)).max() < 0.5          # planted roll recovered (item stage D bar)
    # the naive map is off by the exporter offset; the CT stage must report that offset
    naive = measure_tile(F, C, CanvasMap(0.0, K, 0.0, K), (70, 70), spec)
    fy, fx = CanvasMap(*m).inv(70.0, 70.0)
    want = (m[0] + m[1] * fy - K * fy, m[2] + m[3] * fx - K * fx)
    assert np.abs(np.subtract(naive["map_correction_px"], want)).max() < 0.3


def test_peak_runner_up_rules():
    yy, xx = np.mgrid[-10:11, -10:11]
    one = np.exp(-((yy - 1) ** 2 + xx ** 2) / 18.0)              # broad single peak: skirt but no 2nd peak
    pk = _peak(one, 10)
    assert pk["runner_up"] < 0 and pk["runner_up_any"] > 0.5
    two = one + 0.8 * np.exp(-((yy + 6) ** 2 + (xx - 5) ** 2) / 2.0)
    pk2 = _peak(two, 10)
    assert pk2["runner_up"] > 0.7


def test_fit_field_affine():
    rng = np.random.default_rng(0)
    tiles = []
    for _ in range(30):
        y, x = rng.uniform(0, 4000, 2)
        tiles.append(dict(centre_c=[y, x], map_correction_px=[0.5 + 1e-4 * y, -0.3 + 2e-4 * x]))
    ff = fit_field(tiles)
    assert ff["loo_p95_px"] < 1e-6


# ---------------------------------------------------------------- catalog
def test_catalog_pairs_and_affine_composition():
    A = demo_affine().tolist()
    cat = {"samples": {"PHercX": {
        "sample": {"properties": {}},
        "volumes": {"v24": {"long_id": "v24", "properties": {"pixel_size_um": 2.4, "shape": [10, 10, 10]}},
                    "v94": {"long_id": "v94", "properties": {"pixel_size_um": 9.4, "shape": [5, 5, 5]}}},
        "segments": {"s1": {"long_id": "s1-w", "original_volume_id": "v94", "data": [
            {"type": "tifxyz-transformed", "origins": [{"path": "PHercX/segments/s1/mesh/a.tifxyz/"}],
             "parameters": {"target_volume": "v94"},
             "creation_info": {"process": "copy-normalized-for-original",
                               "derived_from": {"source_path": "s3://vesuvius-challenge/PHercX/segments/s1/mesh/intermediate/tifxyz_normalized/"}}},
            {"type": "tifxyz-transformed", "origins": [{"path": "PHercX/segments/s1/mesh/b.tifxyz/"}],
             "parameters": {"target_volume": "v24"},
             "creation_info": {"process": "transform-tifxyz-normalized",
                               "derived_from": {"source_path": "s3://vesuvius-challenge/PHercX/segments/s1/mesh/intermediate/tifxyz_normalized/"},
                               "provenance": {"parameters": {"commit": "abc", "affine-matrix": str(A)}}}}]}}}}}
    P = pairs(cat)
    assert len(P) == 1 and P[0]["fine"]["px"] == 2.4 and P[0]["coarse"]["is_copy"]
    assert P[0]["source_key"] == "PHercX/segments/s1/mesh/intermediate/tifxyz_normalized"
    M = coarse_to_fine_affine(P[0])
    assert np.allclose(M, np.asarray(A))
