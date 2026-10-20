"""Post-verdict, descriptive (NOT a gate): are the in-plane residuals of the 2 S2 pairs above 0.003 px the float32 source-index
rounding of the exporter's resampler? Prediction per vertex: (src32 - src64) of the resampled mesh's own index, in coarse px."""
import json, sys
from pathlib import Path
import numpy as np
from scipy.ndimage import binary_erosion
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from meshpair.tifxyz import read_tifxyz
from meshpair.catalog import coarse_to_fine_affine
from meshpair import conventions as cv
from meshpair.reverse import _project
from meshpair.reverse_lc import clean_axis
F32 = np.float32
CACHE = HERE.parents[1] / "cache" / "meshes"


def src32(n_new, n_old):
    sc = F32(n_old) / F32(n_new)
    return ((np.arange(n_new, dtype=F32) + F32(0.5)) * sc - F32(0.5)).astype(np.float64)


def src64(n_new, n_old):
    return (np.arange(n_new) + 0.5) * (n_old / n_new) - 0.5


P = {(p["sample"], p["segment"]): p for p in json.load(open(HERE.parents[1] / "out" / "pairs.json"))["pairs"]}
out = []
for key in [("PHerc1667", "20260612121456"), ("PHercParis4", "20260623153216")]:
    p = P[key]; f, c = p["fine"], p["coarse"]
    F = read_tifxyz(CACHE / f["key"]); C = read_tifxyz(CACHE / c["key"])
    import tifffile
    with tifffile.TiffFile(CACHE / p["source_key"] / "x.tif") as tf:
        nS = tf.pages[0].shape[:2]
    gF, gC = F.geom(f["px"]), C.geom(c["px"]); M = coarse_to_fine_affine(p)
    nF, nC = F.points.shape[:2], C.points.shape[:2]
    if c["is_copy"]:   # forward: fine vertex -> coarse (= source) surface
        q, qg, rg, Mq = F, gF, gC, np.linalg.inv(np.vstack([M, [0, 0, 0, 1]]))[:3]
    else:
        q, qg, rg, Mq = C, gC, gF, M
    v = binary_erosion(q.valid, iterations=2); rr, cc = np.nonzero(v)
    pick = np.random.default_rng(0).choice(rr.size, min(100000, rr.size), replace=False); rr, cc = rr[pick], cc[pick]
    pQy, pQx = qg.grid_to_canvas(rr.astype(float), 0), qg.grid_to_canvas(cc.astype(float), 1)
    Y = q.points[rr, cc].astype(float) @ Mq[:, :3].T + Mq[:, 3]
    if c["is_copy"]:
        pRy, pRx = cv.exact_canvas_map(pQy, gF, gC, 0), cv.exact_canvas_map(pQx, gF, gC, 1)
    else:
        pRy, pRx = cv.exact_canvas_map_inv(pQy, gF, gC, 0), cv.exact_canvas_map_inv(pQx, gF, gC, 1)
    g0y, g0x = rg.canvas_to_grid(pRy, 0), rg.canvas_to_grid(pRx, 1)
    nref = (nC if c["is_copy"] else nF)
    clean = clean_axis(g0y, nref[0], nS[0]) & clean_axis(g0x, nref[1], nS[1])
    gy, gx, nd, tg, ok = _project((C if c["is_copy"] else F).points, (C if c["is_copy"] else F).valid, Y,
                                  rg.canvas_to_grid(pRy, 0), rg.canvas_to_grid(pRx, 1))
    # measured in-plane residual in coarse px (signed, per axis)
    if c["is_copy"]:
        ey, ex = rg.grid_to_canvas(gy, 0) - pRy, rg.grid_to_canvas(gx, 1) - pRx
    else:
        ey = cv.exact_canvas_map(gF.grid_to_canvas(gy, 0), gF, gC, 0) - pQy
        ex = cv.exact_canvas_map(gF.grid_to_canvas(gx, 1), gF, gC, 1) - pQx
    # predicted: float32 minus float64 source index of each resampled mesh, in source cells -> coarse px (x n_C/n_S / s_C)
    to_px_y, to_px_x = (nC[0] / nS[0]) / gC.scale_y, (nC[1] / nS[1]) / gC.scale_x
    py = np.zeros(rr.size); px = np.zeros(rr.size)
    if not f["is_copy"]:   # fine resampled: its vertex index (query in forward; found cell in reverse)
        if c["is_copy"]:
            py += (src32(nF[0], nS[0]) - src64(nF[0], nS[0]))[rr] * to_px_y
            px += (src32(nF[1], nS[1]) - src64(nF[1], nS[1]))[cc] * to_px_x
        else:
            jy, jx = np.clip(np.round(gy).astype(int), 0, nF[0] - 1), np.clip(np.round(gx).astype(int), 0, nF[1] - 1)
            py -= (src32(nF[0], nS[0]) - src64(nF[0], nS[0]))[jy] * to_px_y
            px -= (src32(nF[1], nS[1]) - src64(nF[1], nS[1]))[jx] * to_px_x
    if not c["is_copy"]:   # coarse resampled: its own vertex index (query in reverse)
        py += (src32(nC[0], nS[0]) - src64(nC[0], nS[0]))[rr] * to_px_y
        px += (src32(nC[1], nS[1]) - src64(nC[1], nS[1]))[cc] * to_px_x
    good = ok & (tg < 1e-3) & clean   # lattice-clean vertices only, as the rule
    res = np.hypot(ey - py, ex - px)[good]
    rec = dict(sample=key[0], segment=key[1], n=int(good.sum()),
               inplane_p95_measured=float(np.percentile(np.hypot(ey, ex)[good], 95)),
               inplane_p95_after_removing_float32_index=float(np.percentile(res, 95)),
               corr_x=float(np.corrcoef(ex[good], px[good])[0, 1]), corr_y=float(np.corrcoef(ey[good], py[good])[0, 1]),
               pred_p95=float(np.percentile(np.hypot(py, px)[good], 95)))
    print(rec, flush=True); out.append(rec)
json.dump(out, open(HERE / "post_float_check.json", "w"), indent=1)
