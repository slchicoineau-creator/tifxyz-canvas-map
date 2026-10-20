"""Step-1 diagnosis of the TIFXYZ-REVERSE off-surface vertices (post-hoc, S1 pairs + synthetic only; NOT a gate).

Competing explanations, each with a distinguishing prediction:
 H-LAT  lattice: the fine bilinear surface != the source surface inside fine cells that straddle a SOURCE grid line
        (the bilinear source has a kink there). Predicts: (i) coarse vertices lie on the source surface A_F(S) to float noise;
        (ii) off-surface vertices sit only in straddling fine cells, ~0 elsewhere; (iii) the exact chain
        e = A_F(bilerp_S(s_C)) - bilerp_F(g_F formula) predicts the measured signed normal distance per vertex.
 H-SMO  fine (or coarse) re-meshed / smoothed after the transform: predicts replica(S) != published mesh, and off vertices
        not tied to the straddling phase.
 H-AFF  catalog affine M imprecise: predicts a smooth field |Y - A_F(bilerp_S(s_C))| >> float, no phase dependence.
 H-MSK  masked-cell interpolation: predicts off vertices next to invalid fine corners.
 H-FLT  float32 rounding of stored coordinates: predicts |nd| <= a few ulp (~1e-3 coarse vx), no phase dependence.
Source meshes (tifxyz_normalized, public) are fetched into cache/meshes.
"""
import json, sys, time
from pathlib import Path
import numpy as np
from scipy.ndimage import binary_erosion, distance_transform_edt
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from meshpair.tifxyz import read_tifxyz, fetch_tifxyz
from meshpair.catalog import coarse_to_fine_affine
from meshpair.exporter import transform_tifxyz, compare_meshes
from meshpair import reverse as R, conventions as cv

CACHE = HERE.parent / "cache" / "meshes"


def bilerp(P, gy, gx):
    rows, cols = P.shape[:2]
    y0 = np.clip(np.floor(gy).astype(np.int64), 0, rows - 2); x0 = np.clip(np.floor(gx).astype(np.int64), 0, cols - 2)
    fy = (gy - y0)[:, None]; fx = (gx - x0)[:, None]
    q = lambda a, b: P[a, b].astype(np.float64)
    return (q(y0, x0) * (1 - fx) + q(y0, x0 + 1) * fx) * (1 - fy) + (q(y0 + 1, x0) * (1 - fx) + q(y0 + 1, x0 + 1) * fx) * fy


def kink_pred(P, gy, gx):
    """1-D kink model from the FINE mesh alone: along an axis, within one source cell the fine samples are collinear,
    so the second difference D_j is non-zero only next to a source line; for a kink at fraction a of cell [j, j+1]:
    D_j = Delta (1-a), D_{j+1} = Delta a and the chord error at t is -t D_j (t<=a) or -(1-t) D_{j+1} (t>a)
    = -(the smaller of t|D_j|, (1-t)|D_{j+1}|). Returns the predicted error vector (true - interp), x-part + y-part."""
    rows, cols = P.shape[:2]
    P = P.astype(np.float64)
    i = np.clip(np.floor(gy).astype(np.int64), 1, rows - 3); j = np.clip(np.floor(gx).astype(np.int64), 1, cols - 3)
    ty = (gy - i)[:, None]; tx = (gx - j)[:, None]
    def d2x(r, c): return P[r, c + 1] - 2 * P[r, c] + P[r, c - 1]
    def d2y(r, c): return P[r + 1, c] - 2 * P[r, c] + P[r - 1, c]
    def e1(Dl, Dr, t):
        a = t * Dl; b = (1 - t) * Dr
        pick = (np.linalg.norm(a, axis=1) <= np.linalg.norm(b, axis=1))[:, None]
        return -np.where(pick, a, b)
    ex = (1 - ty) * e1(d2x(i, j), d2x(i, j + 1), tx) + ty * e1(d2x(i + 1, j), d2x(i + 1, j + 1), tx)
    ey = (1 - tx) * e1(d2y(i, j), d2y(i + 1, j), ty) + tx * e1(d2y(i, j + 1), d2y(i + 1, j + 1), ty)
    return ex + ey


def straddle(g, nF, nS):
    j = np.floor(g)
    s0 = (j + 0.5) * nS / nF - 0.5; s1 = (j + 1.5) * nS / nF - 0.5
    return np.floor(s0) != np.floor(s1)


def one(p):
    f, c = p["fine"], p["coarse"]
    for k in (f["key"], c["key"], p["source_key"]):
        fetch_tifxyz(k, CACHE)
    F = read_tifxyz(CACHE / f["key"]); C = read_tifxyz(CACHE / c["key"]); S = read_tifxyz(CACHE / p["source_key"])
    gF, gC = F.geom(f["px"]), C.geom(c["px"]); M = coarse_to_fine_affine(p)
    AF = np.asarray(f["affine_src_to_this"], float); AC = np.asarray(c["affine_src_to_this"], float)
    # H-SMO: replica vs published
    repF, iF = transform_tifxyz(S.points, AF); repC, iC = transform_tifxyz(S.points, AC)
    cmpF = compare_meshes(repF, F.points); cmpC = compare_meshes(repC, C.points)
    nS = S.points.shape[:2]; nF = F.points.shape[:2]; nC = C.points.shape[:2]
    v = binary_erosion(C.valid, iterations=2); rr, cc = np.nonzero(v)
    if rr.size > 100000:
        pick = np.random.default_rng(0).choice(rr.size, 100000, replace=False); rr, cc = rr[pick], cc[pick]
    pCy = gC.grid_to_canvas(rr.astype(float), 0); pCx = gC.grid_to_canvas(cc.astype(float), 1)
    X = C.points[rr, cc].astype(float); Y = X @ M[:, :3].T + M[:, 3]
    gy0 = gF.canvas_to_grid(cv.exact_canvas_map_inv(pCy, gF, gC, 0), 0); gx0 = gF.canvas_to_grid(cv.exact_canvas_map_inv(pCx, gF, gC, 1), 1)
    gy, gx, nd, tg, ok = R._project(F.points, F.valid, Y, gy0, gx0)
    k = f["px"] / c["px"]; nd = nd * k; conv = ok & (tg < 1e-3)
    off = conv & (np.abs(nd) > 0.05)
    qy = cv.exact_canvas_map(gF.grid_to_canvas(gy, 0), gF, gC, 0) - pCy
    qx = cv.exact_canvas_map(gF.grid_to_canvas(gx, 1), gF, gC, 1) - pCx
    inpl = np.hypot(qy, qx)
    # H-AFF / (i): source-surface chain. s_C from the coarse resample formula; true point in fine volume.
    sy = (rr + 0.5) * nS[0] / nC[0] - 0.5; sx = (cc + 0.5) * nS[1] / nC[1] - 0.5
    Sv = S.valid
    T_src = bilerp(S.points, sy, sx)
    T = T_src @ AF[:, :3].T + AF[:, 3]
    dYT = np.linalg.norm(Y - T, axis=1) * k     # coarse vx
    # lattice-free control: project Y onto the SOURCE surface carried into the fine volume (A_F(S), no resample)
    SF = (np.where(Sv[..., None], S.points.astype(np.float64) @ AF[:, :3].T + AF[:, 3], -1.0)).astype(np.float64)
    sgy, sgx, snd, stg, sok = R._project(SF, Sv, Y, sy.copy(), sx.copy())
    snd = snd * k; sconv = sok & (stg < 1e-3)
    # (ii) straddle classes of the found fine cell
    sty = straddle(gy, nF[0], nS[0]); stx = straddle(gx, nF[1], nS[1]); st = sty | stx
    # (iii) exact chain prediction at the formula position and the found-point normal
    I0 = bilerp(F.points, gy0, gx0)
    e = T - I0
    _, du, dv, _ = R._bilinear_gather(F.points, F.valid, gy, gx)
    n = np.cross(du, dv); n /= np.linalg.norm(n, axis=1, keepdims=True)
    pred_exact = (e * n).sum(1) * k
    ek = kink_pred(F.points, gy0, gx0)
    pred_kink = (ek * n).sum(1) * k
    # H-MSK: distance (fine cells) to the nearest invalid fine vertex
    dF = distance_transform_edt(F.valid)[np.clip(np.round(gy).astype(int), 0, nF[0] - 1), np.clip(np.round(gx).astype(int), 0, nF[1] - 1)]
    # H-FLT: float32 ulp of the fine coordinates, in coarse vx
    ulp = float(np.spacing(np.float32(np.abs(F.points[F.valid]).max())) * k)
    q = lambda a, qs=(.5, .9, .99): [float(np.quantile(a, x)) for x in qs] if a.size else None
    c_ = conv
    r = dict(sample=p["sample"], segment=p["segment"], nS=list(nS), nF=list(nF), nC=list(nC),
             ratio_F=[nF[0] / nS[0], nF[1] / nS[1]], ratio_C=[nC[0] / nS[0], nC[1] / nS[1]],
             replica_F=dict((kk, cmpF.get(kk)) for kk in ("same_shape", "err_median_vx", "err_p99_vx", "err_max_vx", "mask_agreement")),
             replica_C=dict((kk, cmpC.get(kk)) for kk in ("same_shape", "err_median_vx", "err_p99_vx", "err_max_vx", "mask_agreement")),
             n=int(rr.size), frac_conv=float(conv.mean()), frac_off=float(off.mean()),
             Y_minus_T_coarse_vx_q50_90_99=q(dYT[np.isfinite(dYT)]),
             source_surface_absnd_q50_90_99=q(np.abs(snd[sconv])), source_surface_frac_on=float((sconv & (np.abs(snd) <= 0.05)).mean()),
             source_surface_frac_conv=float(sconv.mean()),
             frac_found_cell_straddles=float(st[c_].mean()),
             frac_off_among_straddle=float((np.abs(nd[c_ & st]) > 0.05).mean()),
             frac_off_among_nonstraddle=float((np.abs(nd[c_ & ~st]) > 0.05).mean()),
             frac_of_off_in_straddle=float(st[off].mean()) if off.any() else None,
             absnd_nonstraddle_q50_90_99=q(np.abs(nd[c_ & ~st])), absnd_straddle_q50_90_99=q(np.abs(nd[c_ & st])),
             inpl_nonstraddle_q50_95_99=q(inpl[c_ & ~st], (.5, .95, .99)), inpl_straddle_q50_95_99=q(inpl[c_ & st], (.5, .95, .99)),
             exact_pred_resid_abs_q50_90_99=q(np.abs(nd - pred_exact)[c_]),
             exact_pred_corr=float(np.corrcoef(nd[c_], pred_exact[c_])[0, 1]),
             kink_pred_resid_abs_q50_90_99=q(np.abs(nd - pred_kink)[c_]),
             kink_pred_corr=float(np.corrcoef(nd[c_], pred_kink[c_])[0, 1]),
             off_dist_to_invalid_fine_q10_50=q(dF[off], (.1, .5)) if off.any() else None,
             all_dist_to_invalid_fine_q10_50=q(dF[c_], (.1, .5)),
             ulp_fine_coord_in_coarse_vx=ulp)
    return r


if __name__ == "__main__":
    P = json.load(open(HERE.parent / "out" / "pairs.json"))["pairs"]
    want = sys.argv[1:]
    sel = [p for p in P if p["segment"] in want and p["sample"] in ("PHerc0139", "PHerc0500P2")]
    out = HERE / "diag" / "diagnose_s1.jsonl"
    with open(out, "a") as fo:
        for p in sel:
            t = time.time(); r = one(p); r["seconds"] = round(time.time() - t, 1)
            fo.write(json.dumps(r) + "\n"); fo.flush(); print(json.dumps(r), flush=True)
