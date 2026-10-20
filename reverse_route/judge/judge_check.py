"""Independent review checks of TIFXYZ-REVERSE amendment 1 (a separate AI review session that did not do the work).

CPU only, cached public meshes only (no downloads: a pair whose files are not all cached is skipped). Not a gate.
  A  re-run the frozen rule (meshpair.reverse_lc.run_pair) on spot-check pairs; compare with s2/s2.jsonl.
  B  independent re-implementation (own code, no Gauss-Newton): evaluate the reference bilinear surface AT the exporter
     formula position of every sampled query vertex; residual r = bilerp_ref(g_formula) - Y. Normal part and in-plane
     part (least squares on the two grid tangents) on clean vs straddling vertices, own clean classification.
  C  planted errors on REAL meshes, measured with the frozen rule:
     C1 global in-plane shift of the query mesh (first-order along its own grid tangents) of 0.01/0.05/0.25/1 coarse px;
     C2 local normal bump of 0.2 coarse vx on a contiguous block holding ~3 % / ~10 % of the eroded query vertices;
     C3 normal bump of 0.2 coarse vx on EVERY query vertex whose formula cell straddles a source line (blind-spot probe).
  D  (application preview, labelled, not a test) the frozen rule on 2 S1 pairs and 2 readable pairs, cached files only.
Usage (from the repository root):  python reverse_route/judge/judge_check.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import tifffile
from scipy.ndimage import binary_erosion

HERE = Path(__file__).resolve().parent
W = HERE.parents[1]
sys.path.insert(0, str(W))
from meshpair.tifxyz import read_tifxyz                      # noqa: E402
from meshpair.catalog import coarse_to_fine_affine            # noqa: E402
from meshpair.reverse_lc import run_pair                      # noqa: E402

CACHE = W / "cache" / "meshes"
PAIRS = {(p["sample"], p["segment"]): p for p in json.load(open(W / "out" / "pairs.json"))["pairs"]}
S2 = {(r["sample"], r["segment"]): r for r in map(json.loads, open(W / "reverse_route" / "s2" / "s2.jsonl"))}

SPOT = [("PHerc0814", "20260226000000"),        # fine is copy (largest 0814)
        ("PHercParis4", "20260604223808"),      # fine is copy
        ("PHercParis4", "20260623153216"),      # both transformed; the float-floor-decided pair
        ("PHercParis4", "20260603185441"),      # both transformed, small
        ("PHercParis4", "20260701183126"),      # both transformed, large
        ("PHerc1667", "20260612121456")]        # coarse is copy (forward), float-floor-decided
PLANT = [("PHercParis4", "20260603185441"), ("PHerc0814", "20260226000000")]
PREVIEW = [("PHerc0139", "20260311000000"), ("PHerc0500P2", "20250611171318")]


def cached(p):
    return all((CACHE / k / "x.tif").exists() for k in (p["fine"]["key"], p["coarse"]["key"]))


def load(p):
    F = read_tifxyz(CACHE / p["fine"]["key"]); C = read_tifxyz(CACHE / p["coarse"]["key"])
    sx = CACHE / p["source_key"] / "x.tif"
    if sx.exists():
        with tifffile.TiffFile(sx) as tf:
            n_src = tuple(int(x) for x in tf.pages[0].shape[:2])
    elif p["coarse"]["is_copy"]:
        n_src = tuple(C.points.shape[:2])
    elif p["fine"]["is_copy"]:
        n_src = tuple(F.points.shape[:2])
    else:
        n_src = None
    return F, C, n_src, coarse_to_fine_affine(p)


def official(p, F, C, n_src, M, Fpts=None, Cpts=None):
    return run_pair(F.points if Fpts is None else Fpts, F.valid, C.points if Cpts is None else Cpts, C.valid,
                    F.geom(p["fine"]["px"]), C.geom(p["coarse"]["px"]), M, n_src, p["fine"]["px"], p["coarse"]["px"],
                    coarse_is_copy=p["coarse"]["is_copy"], fine_is_copy=p["fine"]["is_copy"], n_max=100000, seed=0)


# ---------- independent implementation (B) ----------
# The half-pixel relation below is written out as ((g + 0.5) * n_src) / n_ref on purpose (the judge's own order).
# meshpair.conventions.exporter_grid_map computes (g + 0.5) * (n_src / n_ref): 1-ulp differences on 27 % of the indices
# and changed clean flags on 4 real shape pairs (checked 2026-10-18), so swapping it in would change judge_check.json.
def own_clean(g, n_ref, n_src):
    j = np.floor(g)
    s0 = (j + 0.5) * n_src / n_ref - 0.5
    s1 = (j + 1.5) * n_src / n_ref - 0.5
    k = np.floor(s0) + 1.0                      # smallest integer > s0 (s0 integer -> s0+1)
    return ~(k < s1 - 1e-9)


def own_bilerp(P, V, gy, gx):
    rows, cols = P.shape[:2]
    y0 = np.floor(gy).astype(np.int64); x0 = np.floor(gx).astype(np.int64)
    ins = (y0 >= 0) & (x0 >= 0) & (y0 <= rows - 2) & (x0 <= cols - 2)
    y0 = np.clip(y0, 0, rows - 2); x0 = np.clip(x0, 0, cols - 2)
    ty = (gy - y0)[:, None]; tx = (gx - x0)[:, None]
    a = P[y0, x0].astype(np.float64); b = P[y0, x0 + 1].astype(np.float64)
    c = P[y0 + 1, x0].astype(np.float64); d = P[y0 + 1, x0 + 1].astype(np.float64)
    ok = ins & V[y0, x0] & V[y0, x0 + 1] & V[y0 + 1, x0] & V[y0 + 1, x0 + 1]
    val = (a * (1 - tx) + b * tx) * (1 - ty) + (c * (1 - tx) + d * tx) * ty
    tx_ = (b - a) * (1 - ty) + (d - c) * ty     # d/d col
    ty_ = (c - a) * (1 - tx) + (d - b) * tx     # d/d row
    return val, tx_, ty_, ok


def independent(p, F, C, n_src, M, seed=0):
    px_f, px_c = p["fine"]["px"], p["coarse"]["px"]
    if p["coarse"]["is_copy"]:      # forward: query fine -> reference coarse
        Q, R = F, C; Mq = np.linalg.inv(np.vstack([M, [0, 0, 0, 1]]))[:3]; px_r = px_c
    else:
        Q, R = C, F; Mq = M; px_r = px_f
    nq = Q.points.shape[:2]; nr = R.points.shape[:2]
    v = binary_erosion(Q.valid, iterations=2)
    rr, cc = np.nonzero(v)
    if rr.size > 100000:
        pick = np.random.default_rng(seed).choice(rr.size, 100000, replace=False); rr, cc = rr[pick], cc[pick]
    Y = Q.points[rr, cc].astype(np.float64) @ Mq[:, :3].T + Mq[:, 3]
    gy = (rr + 0.5) * nr[0] / nq[0] - 0.5        # exporter relation between two resamples of one source
    gx = (cc + 0.5) * nr[1] / nq[1] - 0.5
    clean = own_clean(gy, nr[0], n_src[0]) & own_clean(gx, nr[1], n_src[1])
    val, tx_, ty_, ok = own_bilerp(R.points, R.valid, gy, gx)
    r = val - Y
    n = np.cross(tx_, ty_); n /= np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-12)
    dn = (r * n).sum(1) * px_r / px_c                                   # coarse vx
    # in-plane: r_t = a*tx_ + b*ty_  (a,b in reference grid cells)
    A11 = (tx_ * tx_).sum(1); A12 = (tx_ * ty_).sum(1); A22 = (ty_ * ty_).sum(1)
    b1 = (tx_ * r).sum(1); b2 = (ty_ * r).sum(1); det = A11 * A22 - A12 * A12
    da = (A22 * b1 - A12 * b2) / det; db = (A11 * b2 - A12 * b1) / det   # grid cells (col, row) of the reference
    gC = C.geom(px_c)
    if R is F:   # reference grid cells -> coarse grid cells -> coarse canvas px
        ex = da * (C.points.shape[1] / nr[1]) / gC.scale_x; ey = db * (C.points.shape[0] / nr[0]) / gC.scale_y
    else:
        ex = da / gC.scale_x; ey = db / gC.scale_y
    inpl = np.hypot(ex, ey)
    out = dict(n=int(rr.size), frac_clean=float(clean.mean()))
    fin = np.isfinite(dn) & np.isfinite(inpl)
    out["n_nonfinite_dropped"] = int((ok & ~fin).sum())
    for name, m in (("clean", clean & ok & fin), ("straddle", ~clean & ok & fin)):
        if m.sum() == 0:
            continue
        out[name] = dict(n=int(m.sum()), frac_on=float((np.abs(dn[m]) <= 0.05).mean()),
                         p95_abs_normal=float(np.percentile(np.abs(dn[m]), 95)), max_abs_normal=float(np.abs(dn[m]).max()),
                         inplane_p95=float(np.percentile(inpl[m], 95)), inplane_median=float(np.median(inpl[m])))
    return out


# ---------- plants (C) ----------
def query_side(p):
    return "F" if p["coarse"]["is_copy"] else "C"


def grid_tangents(P):
    P = P.astype(np.float64)
    tc = np.zeros_like(P); tr = np.zeros_like(P)
    tc[:, 1:-1] = 0.5 * (P[:, 2:] - P[:, :-2]); tr[1:-1] = 0.5 * (P[2:] - P[:-2])
    return tc, tr


def normals(P):
    tc, tr = grid_tangents(P)
    n = np.cross(tc, tr); n /= np.maximum(np.linalg.norm(n, axis=-1, keepdims=True), 1e-12)
    return n


def cells_per_coarse_px(p, Q, C, axis):
    """query grid cells per coarse canvas px along an axis."""
    gC = C.geom(p["coarse"]["px"]); s = gC.scale_y if axis == 0 else gC.scale_x
    return s * (Q.points.shape[axis] / C.points.shape[axis])


def planted(p, F, C, n_src, M):
    side = query_side(p); Q = F if side == "F" else C
    px_q = p["fine"]["px"] if side == "F" else p["coarse"]["px"]; px_c = p["coarse"]["px"]
    res = []
    tc, tr = grid_tangents(Q.points)
    for shift in (0.01, 0.05, 0.25, 1.0):
        for ax in ("x", "y", "xy"):
            dc = shift * cells_per_coarse_px(p, Q, C, 1) * (1 if "x" in ax else 0) / (np.sqrt(2) if ax == "xy" else 1)
            dr = shift * cells_per_coarse_px(p, Q, C, 0) * (1 if "y" in ax else 0) / (np.sqrt(2) if ax == "xy" else 1)
            Qp = (Q.points.astype(np.float64) + dc * tc + dr * tr).astype(np.float32)
            Qp[~Q.valid] = -1
            o = official(p, F, C, n_src, M, **({"Fpts": Qp} if side == "F" else {"Cpts": Qp}))
            res.append(dict(kind="C1_shift", shift_px=shift, axis=ax, verdict=o["verdict"], frac_on_clean=o["frac_on_clean"],
                            inplane_median_clean=o["inplane_median_clean"], inplane_p95_clean=o["inplane_p95_clean"],
                            bar=o["inplane_bar_px"], p95_normal_clean=o["p95_normal_clean"]))
    del tc, tr
    nrm = normals(Q.points)
    h = 0.2 * px_c / px_q                                                # 0.2 coarse vx in query-volume voxels
    v = binary_erosion(Q.valid, iterations=2); rr, cc = np.nonzero(v)
    r0, c0 = int(np.median(rr)), int(np.median(cc))
    for target in (0.03, 0.10):
        # square block centred on the median vertex grown until it holds the target share of eroded vertices
        lo, hi = 1, max(Q.points.shape)
        while lo < hi:
            mid = (lo + hi) // 2
            share = ((np.abs(rr - r0) <= mid) & (np.abs(cc - c0) <= mid)).mean()
            lo, hi = (mid + 1, hi) if share < target else (lo, mid)
        blk = np.zeros(Q.valid.shape, bool); blk[max(0, r0 - lo):r0 + lo + 1, max(0, c0 - lo):c0 + lo + 1] = True
        blk &= Q.valid
        Qp = Q.points.astype(np.float64); Qp[blk] += h * nrm[blk]; Qp = Qp.astype(np.float32)
        o = official(p, F, C, n_src, M, **({"Fpts": Qp} if side == "F" else {"Cpts": Qp}))
        res.append(dict(kind="C2_local_bump", target_share=target, share_eroded=float(((np.abs(rr - r0) <= lo) & (np.abs(cc - c0) <= lo)).mean()),
                        bump_coarse_vx=0.2, verdict=o["verdict"], frac_on_clean=o["frac_on_clean"],
                        p95_normal_clean=o["p95_normal_clean"]))
    # C3: bump only the vertices whose formula cell straddles (both-transformed reverse pairs only)
    if side == "C" and not p["fine"]["is_copy"]:
        nr = F.points.shape[:2]; nq = C.points.shape[:2]
        R_, Cc = np.meshgrid(np.arange(nq[0]), np.arange(nq[1]), indexing="ij")
        gy = (R_ + 0.5) * nr[0] / nq[0] - 0.5; gx = (Cc + 0.5) * nr[1] / nq[1] - 0.5
        st = ~(own_clean(gy, nr[0], n_src[0]) & own_clean(gx, nr[1], n_src[1])) & Q.valid
        Qp = Q.points.astype(np.float64); Qp[st] += h * nrm[st]; Qp = Qp.astype(np.float32)
        o = official(p, F, C, n_src, M, Cpts=Qp)
        res.append(dict(kind="C3_straddle_only_bump", share_bumped=float(st[v].mean()), bump_coarse_vx=0.2, verdict=o["verdict"],
                        frac_on_clean=o["frac_on_clean"], frac_on_all=o["frac_on_all"], p95_normal_clean=o["p95_normal_clean"],
                        p95_normal_straddle=o.get("p95_normal_straddle")))
    return res


def main():
    out = dict(spot=[], plants=[], preview=[])
    keys = ("verdict", "n_clean", "frac_clean", "frac_on_clean", "inplane_p95_clean", "inplane_bar_px", "p95_normal_clean", "offset_norm_px")
    for k in SPOT:
        p = PAIRS[k]; t0 = time.time()
        F, C, n_src, M = load(p)
        o = official(p, F, C, n_src, M)
        ref = S2[k]["lattice_clean"]
        same = {kk: (o[kk] == ref[kk]) if isinstance(o[kk], (str, bool, int)) else bool(np.isclose(o[kk], ref[kk], rtol=1e-9, atol=0))
                for kk in keys}
        ind = independent(p, F, C, n_src, M)
        ind1 = independent(p, F, C, n_src, M, seed=1)
        rec = dict(sample=k[0], segment=k[1], n_src=list(n_src), rerun={kk: o[kk] for kk in keys}, identical_to_s2=same,
                   n_valid_eroded=o["n_valid_eroded"], share_eroded_measured=o["n_clean"] / o["n_valid_eroded"],
                   bar_plain_pass=bool(o["inplane_p95_clean"] <= 0.003 and o["frac_on_clean"] >= 0.95),
                   independent_seed0=ind, independent_seed1=ind1, seconds=round(time.time() - t0, 1))
        print(json.dumps(rec), flush=True); out["spot"].append(rec)
        del F, C
    for k in PLANT:
        p = PAIRS[k]; F, C, n_src, M = load(p)
        for r in planted(p, F, C, n_src, M):
            r.update(sample=k[0], segment=k[1]); print(json.dumps(r), flush=True); out["plants"].append(r)
        del F, C
    for k in PREVIEW:
        p = PAIRS[k]
        if not cached(p):
            continue
        F, C, n_src, M = load(p)
        if n_src is None:
            continue
        o = official(p, F, C, n_src, M)
        rec = dict(sample=k[0], segment=k[1], label="APPLICATION PREVIEW (S1 pair; rule designed on S1)", **{kk: o[kk] for kk in keys})
        print(json.dumps(rec), flush=True); out["preview"].append(rec)
    # readable preview: first cached IN_FILES pair of 0009B and of 1667 (coarse is copy -> forward)
    fwd = {}
    for ln in open(W / "out" / "mesh_records.jsonl"):
        r = json.loads(ln); fwd[(r["sample"], r["segment"])] = (r.get("lines") or {}).get("M1_in_files")
    for s in ("PHerc0009B", "PHerc1667"):
        for k, p in PAIRS.items():
            if k[0] == s and fwd.get(k) == "IN_FILES" and cached(p):
                F, C, n_src, M = load(p)
                if n_src is None:
                    continue
                o = official(p, F, C, n_src, M)
                rec = dict(sample=k[0], segment=k[1], label="APPLICATION PREVIEW (readable pair)", coarse_is_copy=p["coarse"]["is_copy"],
                           **{kk: o[kk] for kk in keys})
                print(json.dumps(rec), flush=True); out["preview"].append(rec)
                break
    json.dump(out, open(HERE / "judge_check.json", "w"), indent=1, default=float)


def addendum():
    """Re-run only the independent check of the flatboi with seed 1, dropping non-finite tangent solves (degenerate cells)."""
    k = ("PHerc1667", "20260612121456"); p = PAIRS[k]; F, C, n_src, M = load(p)
    rec = dict(sample=k[0], segment=k[1], independent_seed1_finite=independent(p, F, C, n_src, M, seed=1))
    js = json.load(open(HERE / "judge_check.json")); js["addendum_flatboi_seed1"] = rec
    json.dump(js, open(HERE / "judge_check.json", "w"), indent=1, default=float); print(json.dumps(rec))


if __name__ == "__main__":
    addendum() if "--addendum" in sys.argv else main()
