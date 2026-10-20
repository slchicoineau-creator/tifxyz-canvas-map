"""Tests added for the public release: the canvasmap helper, the pair-map table, and the CLI changes (CHANGES.md)."""
import json
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import check                                                             # noqa: E402
from meshpair import canvasmap as cm                                     # noqa: E402
from meshpair import conventions as cv                                   # noqa: E402
from meshpair.tifxyz import write_tifxyz                                 # noqa: E402

F32_005 = float(np.float32(0.05))


def _mesh(tmp, name, rows, cols, scale=0.05):
    pts = np.zeros((rows, cols, 3), np.float32)
    write_tifxyz(tmp / name, pts, {"scale": [scale, scale], "format": "tifxyz"})
    return tmp / name


def test_helper_matches_published_formula(tmp_path):
    """With s = 0.05 and render scale 1 the map is p_C = (p_F + 10.5) * n_C/n_F - 10.5 on each axis."""
    F = _mesh(tmp_path, "fine", 823, 928)
    C = _mesh(tmp_path, "coarse", 211, 238)
    p = np.array([0.0, 17.25, 9000.0, 18559.0])
    yc, xc = cm.fine_to_coarse_px(p, p, F, C)
    np.testing.assert_allclose(yc, (p + 10.5) * 211 / 823 - 10.5, atol=1e-9)
    np.testing.assert_allclose(xc, (p + 10.5) * 238 / 928 - 10.5, atol=1e-9)


def test_helper_round_trip_and_levels(tmp_path):
    F = _mesh(tmp_path, "fine", 300, 410, F32_005)
    C = _mesh(tmp_path, "coarse", 77, 105)
    rng = np.random.default_rng(0)
    y, x = rng.uniform(0, 6000, 50), rng.uniform(0, 8200, 50)
    yc, xc = cm.fine_to_coarse_px(y, x, F, C)
    yb, xb = cm.coarse_to_fine_px(yc, xc, F, C)
    np.testing.assert_allclose(yb, y, atol=1e-8)
    np.testing.assert_allclose(xb, x, atol=1e-8)
    # a fine level-2 pixel maps to the same place as its level-0 centre 4 i + 1.5
    i = np.array([3.0, 250.0])
    y2, x2 = cm.fine_to_coarse_px(i, i, F, C, fine_level=2)
    y0, x0 = cm.fine_to_coarse_px(4 * i + 1.5, 4 * i + 1.5, F, C)
    np.testing.assert_allclose([y2, x2], [y0, x0], atol=1e-9)
    # coarse level 1 is half the level-0 index after the centre shift
    y1, _ = cm.fine_to_coarse_px(y, x, F, C, coarse_level=1)
    np.testing.assert_allclose(y1, (yc - 0.5) / 2, atol=1e-9)


def test_helper_reads_float32_scale(tmp_path):
    """Published meshes store 0.05 as a double or as float32; the helper uses the stored value."""
    F = _mesh(tmp_path, "fine", 926, 738, F32_005)
    C = _mesh(tmp_path, "coarse", 257, 205, 0.05)
    g = cm.mesh_geom(F)
    assert g.scale_x == F32_005 and g.rows == 926 and g.cols == 738
    c = cm.linear_coeffs(F, C)
    a, b = cv.exact_linear_coeffs(cv.MeshGeom(926, 738, F32_005, F32_005), cv.MeshGeom(257, 205), 0)
    assert (c["y_a"], c["y_b"]) == (a, b)


def test_warp_moves_content_to_mapped_place(tmp_path):
    """A blob drawn on the fine canvas lands where fine_to_coarse_px says, not where index scaling says."""
    F = _mesh(tmp_path, "fine", 80, 90)
    C = _mesh(tmp_path, "coarse", 21, 23)
    gF = cm.mesh_geom(F)
    img = np.zeros((gF.canvas_h, gF.canvas_w), np.float32)
    yy, xx = np.mgrid[0:gF.canvas_h, 0:gF.canvas_w]
    cy, cx = 900.0, 1000.0
    img[:] = np.exp(-((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * 40.0 ** 2))
    out = cm.warp_canvas_image(img, F, C, "fine_to_coarse")
    gC = cm.mesh_geom(C)
    assert out.shape == (gC.canvas_h, gC.canvas_w)
    py, px = np.unravel_index(np.argmax(out), out.shape)
    ey, ex = cm.fine_to_coarse_px(cy, cx, F, C)
    assert abs(py - ey) <= 1 and abs(px - ex) <= 1
    back = cm.warp_canvas_image(out, F, C, "coarse_to_fine")
    by, bx = np.unravel_index(np.argmax(back), back.shape)
    assert abs(by - cy) <= 4 and abs(bx - cx) <= 4


def test_pair_maps_table_reproduces_records():
    """Every stored per-pair map equals the conventions code run on that pair's grid shapes (with each mesh's stored
    scale, 0.05 or float32 0.05): the table, the records and the code agree."""
    recs = [json.loads(l) for f in ("out/mesh_records.jsonl", "out/dev_w00_mesh.jsonl")
            for l in (ROOT / f).read_text().splitlines() if l.strip()]
    assert len(recs) == 207
    for r in recs:
        em = r["exporter_model"]
        best = min(
            max(abs(a - em[k]) for a, k in zip(
                [*cv.exact_linear_coeffs(cv.MeshGeom(*r["fine_shape"][:2], sf, sf), cv.MeshGeom(*r["coarse_shape"][:2], sc, sc), 0),
                 *cv.exact_linear_coeffs(cv.MeshGeom(*r["fine_shape"][:2], sf, sf), cv.MeshGeom(*r["coarse_shape"][:2], sc, sc), 1)],
                ["y_a", "y_b", "x_a", "x_b"]))
            for sf in (0.05, F32_005) for sc in (0.05, F32_005))
        assert best == 0.0, (r["sample"], r["segment"], best)
    rows = (ROOT / "results" / "pair_maps.tsv").read_text().splitlines()
    assert len(rows) == 208
    head = rows[0].split("\t")
    by = {(r["sample"], r["segment"]): r["exporter_model"] for r in recs}
    for line in rows[1:]:
        row = dict(zip(head, line.split("\t")))
        em = by[(row["sample"], row["segment"])]
        for k in ("y_a", "y_b", "x_a", "x_b"):
            assert abs(float(row[k]) - em[k]) < 1e-8, (row["segment"], k)


def _args(**kw):
    d = dict(samples=None, segments=None, exclude_dev=False, sample_n=None, seed=0)
    d.update(kw)
    return SimpleNamespace(**d)


def test_segments_accepts_list_or_file(tmp_path):
    P = [dict(sample="A", segment=s, long_id=s + "-x") for s in ("1", "2", "3")]
    assert [p["segment"] for p in check._select(P, _args(segments="1,3"))] == ["1", "3"]
    f = tmp_path / "segs.txt"
    f.write_text("2\n3\n")
    assert [p["segment"] for p in check._select(P, _args(segments=str(f)))] == ["2", "3"]
    f.write_text("1,2")
    assert [p["segment"] for p in check._select(P, _args(segments=str(f)))] == ["1", "2"]


def test_prereg_guard_is_opt_in(monkeypatch):
    sel = [dict(sample="PHerc0139", segment="x")]
    monkeypatch.setattr(check, "REQUIRE_PUSHED_PREREG", False)
    check._guard(sel, None)                       # off by default: no exit
    monkeypatch.setattr(check, "REQUIRE_PUSHED_PREREG", True)
    with pytest.raises(SystemExit):
        check._guard(sel, None)
    check._guard([dict(sample="PHerc0841", segment="20260220213127")], None)   # the development pair is always allowed


def test_canvasmap_cli(tmp_path):
    F = _mesh(tmp_path, "fine", 823, 928)
    C = _mesh(tmp_path, "coarse", 211, 238)
    out = subprocess.run([sys.executable, "-m", "meshpair.canvasmap", str(F), str(C), "--fine-px", "0", "0"],
                         cwd=ROOT, capture_output=True, text=True, check=True).stdout
    assert '"y_b"' in out and "coarse px:" in out
