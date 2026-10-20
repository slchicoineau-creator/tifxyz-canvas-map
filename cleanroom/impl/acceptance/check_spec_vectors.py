"""Acceptance vectors of cleanroom/SPEC.md section 7 (run: python cr02_accept.py <package root>)."""
import sys, json, hashlib, tempfile
from pathlib import Path
import numpy as np
ROOT = Path(sys.argv[1]); sys.path.insert(0, str(ROOT))
from meshpair import conventions as cv
from meshpair.exporter import (measure_grid_axis_spacing, transform_tifxyz, apply_affine_f32,
                               resample_linear_preserving_invalids, compare_meshes, _median_cpp)
from meshpair.tifxyz import read_tifxyz, write_tifxyz
from meshpair.slicecentre import coarse_layers
from meshpair.synth import curved_sheet, demo_affine
F32_005 = float(np.float32(0.05))
fails = []


def chk(name, got, want):
    ok = got == want
    if not ok:
        fails.append(name)
    print(("OK  " if ok else "FAIL"), name, "" if ok else f"got {got!r} want {want!r}")


def sha(a):
    return hashlib.sha256(np.ascontiguousarray(a, dtype="<f4").tobytes()).hexdigest()


# 7.1
tab = {(0.05, 0.05): [-7.808019441069064, 0.2563791008501539, -7.80711206896558, 0.25646551724184974],
       (0.05, F32_005): [-7.807987890721506, 0.2563790970298214, -7.807076495318142, 0.256465513420153],
       (F32_005, 0.05): [-7.808050880609244, 0.25637910467094116, -7.807147531818373, 0.25646552106309173],
       (F32_005, F32_005): [-7.8080193302612315, 0.25637910085060867, -7.807111958170481, 0.256465517241395]}
for (sf, sc), w in tab.items():
    F = cv.MeshGeom(823, 928, sf, sf); C = cv.MeshGeom(211, 238, sc, sc)
    chk(f"coeffs {sf},{sc}", [*cv.exact_linear_coeffs(F, C, 0), *cv.exact_linear_coeffs(F, C, 1)], w)
    chk(f"canvas {sf}", (F.canvas_h, F.canvas_w), (16460, 18560))
g = cv.MeshGeom(211, 238); p = [0.0, 17.3, 4759.0]
chk("c2g ax1", g.canvas_to_grid(p, 1).tolist(), [0.024999999999991473, 0.8900000000000006, 237.97500000000002])
chk("c2g ax0", g.canvas_to_grid(p, 0).tolist(), [0.024999999999991473, 0.8900000000000006, 237.975])
chk("ecm", float(cv.exact_canvas_map(5000.0, cv.MeshGeom(823, 928), g, 0)), 1274.0874848116646)
chk("ecm_inv", float(cv.exact_canvas_map_inv(1270.0, cv.MeshGeom(823, 928), g, 1)), 4982.373949579833)
chk("cw R.5", cv.MeshGeom(10, 10, 0.05, 0.05, render_scale=0.5).canvas_w, 100)
chk("cw R0", cv.MeshGeom(1, 1, 0.05, 0.05, render_scale=0.0).canvas_w, 1)
chk("cw 7x5", cv.MeshGeom(7, 5, 0.04, 0.04).canvas_w, 125)
chk("lround", [cv._lround(2.5), cv._lround(-2.5), cv._lround(3.5)], [3, -3, 4])
os_ = cv.offset_summary(cv.MeshGeom(823, 928, px_um=2.403), cv.MeshGeom(211, 238, px_um=9.366))
want = dict(k_naive=0.25656630365150546, y_a=-7.808019441069064, y_b=0.2563791008501539, y_b_minus_k=-0.0001872028013515381,
            y_offset_at_origin_px=-7.808019441069064, y_offset_box_min_px=-10.889190348513694, y_offset_box_max_px=-7.808019441069064,
            y_offset_box_centre_px=-9.348604894791606, x_a=-7.80711206896558, x_b=0.25646551724184974, x_b_minus_k=-0.00010078640965571761,
            x_offset_at_origin_px=-7.80711206896558, x_offset_box_min_px=-9.677607045766308, x_offset_box_max_px=-7.80711206896558,
            x_offset_box_centre_px=-8.742359557365944, y_offset_at_origin_um=-73.12991008505284, y_offset_box_min_um=-101.98815680417925,
            y_offset_box_max_um=-73.12991008505284, y_offset_box_centre_um=-87.55903344461818, x_offset_at_origin_um=-73.12141163793162,
            x_offset_box_min_um=-90.64046759064723, x_offset_box_max_um=-73.12141163793162, x_offset_box_centre_um=-81.88093961428943)
chk("offset_summary keys", list(os_), list(want)); chk("offset_summary vals", os_, want)
o2 = cv.offset_summary(cv.MeshGeom(823, 928), cv.MeshGeom(211, 238), k=0.25)
chk("os2", (len(o2), o2["y_b_minus_k"], o2["y_offset_box_max_px"], o2["y_offset_box_centre_px"], o2["x_offset_box_max_px"],
            o2["x_offset_box_centre_px"]),
    # SPEC A.3/section 7.1 say 14 keys; the key list itself gives 15 (erratum, see IMPLEMENTER_LOG.md gap 1)
    (15, 0.006379100850153918, 97.18560145161427, 44.68879100527238, 112.18642242252372, 52.18965517677907))
recs = [json.loads(l) for f in ("out/mesh_records.jsonl", "out/dev_w00_mesh.jsonl")
        for l in (ROOT / f).read_text().splitlines() if l.strip()]
n = 0; split = {}
for r in recs:
    hit = None
    for sf, sn in ((0.05, "0.05"), (F32_005, "f32")):
        for sc, cn in ((0.05, "0.05"), (F32_005, "f32")):
            gF = cv.MeshGeom(*r["fine_shape"][:2], sf, sf, px_um=r["px_fine"])
            gC = cv.MeshGeom(*r["coarse_shape"][:2], sc, sc, px_um=r["px_coarse"])
            m = cv.offset_summary(gF, gC, r["px_fine"] / r["px_coarse"], valid_box=tuple(r["fine_canvas_valid_box"]))
            if list(m) == list(r["exporter_model"]) and m == r["exporter_model"]:
                hit = (sn, cn); break
        if hit:
            break
    if hit:
        n += 1; split[hit] = split.get(hit, 0) + 1
chk("records", len(recs), 207)
chk("full-record parity", n, 207); print("split", split)
# 7.2
S80 = curved_sheet(rows=80, cols=90, seed=1)
chk("S80 sha", sha(S80), "d640295bbb43a0cf86688a5f2274031c76c5744b3b7fd96d20891941eacac4c9")
chk("SP1", measure_grid_axis_spacing(S80), (20.00825044360449, 20.00935154363681, 271, 270))
o, i = transform_tifxyz(S80, demo_affine())
chk("T1", (o.shape, sha(o), int((o[..., 0] == -1).sum()), i["spacing_after"], i["ratio_x"], i["ratio_y"], i["resampled"], i["new_shape"]),
    ((312, 351, 3), "1df4347843d55fe734a2e47f48c02f63dd6fd77bfca21b4033ad2ea2ef806b08", 11165,
     (77.99673724418912, 77.97294373053788), 3.898228756383859, 3.89682512001915, True, (312, 351)))
I = np.c_[np.eye(3), np.zeros(3)]
o, i = transform_tifxyz(S80, I)
chk("T-id", (sha(o), i["ratio_x"], i["ratio_y"], i["resampled"], "new_shape" in i), (sha(S80), 1.0, 1.0, False, False))
o, i = transform_tifxyz(S80, np.c_[1.0005 * np.eye(3), np.zeros(3)])
chk("T-win", (i["ratio_x"], i["ratio_y"], i["resampled"], o.shape[:2], sha(o)),
    (1.000502824758029, 1.000527017154908, False, (80, 90), "aa8d590a495b330fd85464bbb7ef9be4395d9dcf2e31e1ce9aa806b8d7088ec1"))
o, i = transform_tifxyz(S80, np.c_[1.002 * np.eye(3), np.zeros(3)]); X = o
chk("T-win2", (i["ratio_x"], i["ratio_y"], i["resampled"], i["new_shape"], int((o[..., 0] == -1).sum()), sha(o)),
    (1.0019980692007522, 1.0020020491101451, True, (80, 90), 641, "eff60163cd8b709b19f2d14791e7269837bd47c13db7353b56e93d414e3ae1e8"))
S60 = curved_sheet(rows=60, cols=70, holes=0, seed=3); S60[:2] = S60[2]; S60[:, :4] = S60[:, 4:5]
chk("S60 sha", sha(S60), "1c9f80ffcacaf50a7081557d737fb04635f332a0465b7366aa107c05c57b7335")
chk("R1", sha(resample_linear_preserving_invalids(S60, 233, 271)), "b0c4a8d3f3022d31c51b364e7339964ba5bf2f817091f048b6cf12c98d801109")
r2 = resample_linear_preserving_invalids(S80, 233, 271)
chk("R2", (sha(r2), int((r2[..., 0] == -1).sum())), ("fe29e885c128671f071f47c00b161040b8ef1b8297f93c78b8efb3bd453616c3", 6376))
r, c = np.mgrid[0:6, 0:7].astype(np.float64)
V1 = np.stack([c * 20 + 0.1 * r * r, r * 20 + 0.05 * c, 1000 + 0.3 * r * c], -1).astype(np.float32)
V1[3, 4] = -1; V1[1, 2, 1] = -1
chk("V1 sha", sha(V1), "5b8b07e9a68c09a31681ebd3e964b1bee7548791cb627ad79ec17770c4d0b14e")
A = np.array([[3.9, 0, 0, 10.25], [0, 3.9, 0, -3.5], [0, 0, 3.9, 7.0]])
chk("SP-V1", measure_grid_axis_spacing(V1), (20.02816521491387, 20.037216479929768, 28, 28))
af = apply_affine_f32(V1, A, 2.0, 0.5)
chk("AF1", (sha(af), af[0, 0].tolist(), af[5, 6].tolist(), af[1, 2].tolist()),
    ("a60dc5dac8cea5e87d02503454750df836d50c5e75fa152573224be1d432cb1c", [5.125, -1.75, 3903.5],
     [482.875, 389.4200134277344, 3938.60009765625], [161.51499938964844, -5.650000095367432, 3905.83984375]))
o, i = transform_tifxyz(V1, A)
chk("TV1", (o.shape[:2], sha(o), int((o[..., 0] == -1).sum()), i["ratio_x"], i["ratio_y"], i["new_shape"], o[0, 0].tolist(), o[22, 26].tolist()),
    ((23, 27), "7492b7162611290571e302dc63c92a104a129914f034b37401709ce411e8cdf2", 56, 3.8999999112942025, 3.899999914065268,
     (23, 27), [10.25, -3.5, 3907.0], [488.0, 387.6700134277344, 3942.10009765625]))
L = np.array([[[0, 0, 0], [20, 0, 0], [40, 0, 0]]], np.float32)
o, i = transform_tifxyz(L, np.c_[1.5 * np.eye(3), np.zeros(3)])
chk("TL", (i["new_shape"], o[..., 0].tolist(), o[..., 1].tolist(), o[..., 2].tolist()),
    ((2, 4), [[0.0, 18.75, 41.25, 60.0]] * 2, [[0.0] * 4] * 2, [[0.0] * 4] * 2))
o, i = transform_tifxyz(L, np.c_[np.diag([1.5, 1.5, 1.6]), np.zeros(3)])
chk("TL-nu", (i["ratio_x"], i["ratio_y"], i["resampled"], o.shape[:2]), (1.0, 1.0, False, (1, 3)))
e = apply_affine_f32(np.array([[(1, 2, 3), (np.nan, 2, 3), (5, -1, 7), (-1, 0, 0)]], np.float32), np.c_[2 * np.eye(3), np.zeros(3)])
chk("AF-edge", repr(e.tolist()), repr([[[2.0, 4.0, 6.0], [float("nan"), 2.0, 3.0], [10.0, -2.0, 14.0], [-1.0, 0.0, 0.0]]]))
chk("AF-m1", apply_affine_f32(np.array([[(2, 0, 0)]], np.float32), np.c_[np.eye(3), [-3, 0, 0]]).tolist(), [[[-1.0, -1.0, -1.0]]])
chk("CM", compare_meshes(S80, S80[:-1]), {"pred_shape": [80, 90], "pub_shape": [79, 90], "same_shape": False})
cm = compare_meshes(X, X)
chk("CM2", (cm["n_valid_pub"], cm["n_valid_pred"], cm["n_both"], cm["mask_mismatch"], cm["mask_agreement"],
            cm["err_median_vx"], cm["err_p99_vx"], cm["err_max_vx"]), (6559, 6559, 6559, 0, 1.0, 0.0, 0.0, 0.0))
chk("median", (_median_cpp([3, 1, 2, 4]), _median_cpp([3, 1, 2])), (2.5, 2.0))
# 7.3
with tempfile.TemporaryDirectory() as t:
    P = np.ones((2, 2, 3), np.float32); P[0, 0, 2] = 0; P[0, 1, 2] = -1; P[1, 0, 2] = 5; P[1, 1, 2] = np.nan
    write_tifxyz(Path(t) / "m", P, {"scale": [0.04, 0.05]})
    m = read_tifxyz(Path(t) / "m", vc_load=True)
    chk("vc_load", repr(m.points.tolist()),
        repr([[[-1.0, -1.0, -1.0], [-1.0, -1.0, -1.0]], [[1.0, 1.0, 5.0], [1.0, 1.0, float("nan")]]]))
    g = m.geom(); chk("geom scale", (g.scale_x, g.scale_y, g.rows, g.cols), (0.04, 0.05, 2, 2))
    chk("no vc_load", float(read_tifxyz(Path(t) / "m").points[0, 1, 2]), -1.0)
# D
chk("coarse_layers", [coarse_layers(12), coarse_layers(7), coarse_layers(13), coarse_layers(12, 3), coarse_layers(5, 9),
                      coarse_layers(1), coarse_layers(2)],
    [list(range(2, 10)), list(range(7)), list(range(3, 10)), [4, 5, 6], list(range(5)), [0], [0, 1]])
print("FAILS:", fails)
sys.exit(1 if fails else 0)
