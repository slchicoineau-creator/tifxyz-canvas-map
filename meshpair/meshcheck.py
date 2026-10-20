"""Mesh-route check of one fine/coarse pair -> one JSON record (stages M1-M3 of the prereg draft)."""
from __future__ import annotations

import numpy as np

from . import conventions as cv
from .exporter import compare_meshes, transform_tifxyz
from .fit import affine_point_deviation, free_affine_fit, inverse_map, linear_fit_axis, split_half_noise
from .tifxyz import Mesh

# Frozen decision lines (copied into PREREG_DRAFT.md; change only there, before any confirmatory run).
IN_FILES_MIN_PX = 3.0          # |route-I offset vs naive| at the canvas centre, coarse px: >= -> offset is in the files
NOT_IN_FILES_MAX_PX = 1.0      # <= -> the files agree with naive index pairing (offset was our assumption)
MODEL_MATCH_P95_PX = 0.25      # p95 |route-I p_C - exporter-model p_C| over sampled vertices, coarse px
REPRO_MEDIAN_VX = 0.01         # exporter replica vs published mesh, target-volume voxels
REPRO_P99_VX = 0.05
REPRO_MASK_AGREE = 0.999
MESH_NOISE_P95_PX = 0.05       # split-half /2 noise of route I, coarse px
ROUTE_I_MIN_OK = 0.95          # fraction of sampled fine vertices with a coincident coarse surface point


def to44(A):
    A = np.asarray(A, np.float64)
    M = np.eye(4)
    M[:3, :] = A[:3, :4]
    return M


def compose_coarse_to_fine(A_src_to_fine, A_src_to_coarse):
    """M = A_F o A_C^-1 (3x4). Either may be None (= identity: mesh is a plain copy of the source)."""
    F = to44(A_src_to_fine) if A_src_to_fine is not None else np.eye(4)
    C = to44(A_src_to_coarse) if A_src_to_coarse is not None else np.eye(4)
    return (F @ np.linalg.inv(C))[:3]


def check_pair(fine: Mesh, coarse: Mesh, px_fine: float, px_coarse: float, M_coarse_to_fine,
               source: Mesh | None = None, A_src_to_fine=None, A_src_to_coarse=None, fine_is_copy=False,
               coarse_is_copy=False, n_inv: int = 60000, n_free: int = 30000, do_free: bool = True) -> dict:
    gF = fine.geom(px_fine)
    gC = coarse.geom(px_coarse)
    k = px_fine / px_coarse
    vF, vC = fine.valid, coarse.valid
    rr, cc = np.nonzero(vF)
    box = (float(gF.grid_to_canvas(rr.min(), 0)), float(gF.grid_to_canvas(rr.max(), 0)),
           float(gF.grid_to_canvas(cc.min(), 1)), float(gF.grid_to_canvas(cc.max(), 1)))
    rec = dict(fine_shape=list(fine.points.shape[:2]), coarse_shape=list(coarse.points.shape[:2]),
               px_fine=px_fine, px_coarse=px_coarse, fine_valid=int(vF.sum()), coarse_valid=int(vC.sum()),
               fine_canvas_valid_box=box)
    rec["exporter_model"] = cv.offset_summary(gF, gC, k, valid_box=box)

    def naive(py, px):
        return cv.naive_canvas_map(py, k), cv.naive_canvas_map(px, k)

    M = np.asarray(M_coarse_to_fine, np.float64)
    inv = inverse_map(fine.points, vF, coarse.points, vC, gF, gC, M, naive, n=n_inv)
    ok = inv["ok"]
    py, px, cy, cx = inv["py"][ok], inv["px"][ok], inv["cy"][ok], inv["cx"][ok]
    tc = inv["tangential_converged"]
    ri = dict(n_sampled=int(ok.size), n_ok=int(ok.sum()), n_ok_first_start=inv["n_ok_first_start"],
              frac_ok=float(ok.mean()),
              normal_dist_all_converged_um=dict(
                  median=float(np.median(inv["normal_d"][tc]) * px_coarse) if tc.any() else float("nan"),
                  p95_abs=float(np.percentile(np.abs(inv["normal_d"][tc]), 95) * px_coarse) if tc.any() else float("nan")))
    if ok.sum() >= 100:
        ay, by, resy = linear_fit_axis(py, cy)
        ax, bx, resx = linear_fit_axis(px, cx)
        yc = 0.5 * (box[0] + box[1]); xc = 0.5 * (box[2] + box[3])
        ey = cy - cv.exact_canvas_map(py, gF, gC, 0)
        ex = cx - cv.exact_canvas_map(px, gF, gC, 1)
        e = np.hypot(ey, ex)
        off_y = (ay + by * yc) - k * yc
        off_x = (ax + bx * xc) - k * xc
        ri.update(ay=ay, by=by, ax=ax, bx=bx,
                  linear_resid_p95_px=float(np.percentile(np.hypot(resy, resx), 95)),
                  offset_vs_naive_centre_px=[float(off_y), float(off_x)],
                  offset_vs_naive_centre_um=[float(off_y * px_coarse), float(off_x * px_coarse)],
                  offset_vs_naive_centre_norm_px=float(np.hypot(off_y, off_x)),
                  offset_vs_naive_range_px=dict(
                      y=[float(np.min(cy - k * py)), float(np.max(cy - k * py))],
                      x=[float(np.min(cx - k * px)), float(np.max(cx - k * px))]),
                  vs_exporter_model_px=dict(median=float(np.median(e)), p95=float(np.percentile(e, 95)),
                                            max=float(e.max()), n_gt_0p25=int((e > 0.25).sum())),
                  normal_dist_um=dict(median=float(np.median(inv["normal_d"][ok]) * px_coarse),
                                      p95_abs=float(np.percentile(np.abs(inv["normal_d"][ok]), 95) * px_coarse)),
                  split_half_noise_px=split_half_noise(py, px, cy, cx))
    rec["route_I"] = ri

    if do_free and ok.sum() >= 100:
        fr = free_affine_fit(fine.points, vF, coarse.points, vC, gF, gC, (0.0, k, 0.0, k), n=n_free)
        sel = np.nonzero(vC)
        pts = coarse.points[sel].astype(np.float64)
        pts = pts[np.random.default_rng(2).choice(len(pts), min(20000, len(pts)), replace=False)]
        dev = affine_point_deviation(fr["affine"], M, pts)
        yc = 0.5 * (box[0] + box[1]); xc = 0.5 * (box[2] + box[3])
        rec["route_F"] = dict(ay=fr["ay"], by=fr["by"], ax=fr["ax"], bx=fr["bx"], n=fr["n"],
                              resid_median_um=fr["resid_median_vx"] * px_fine, resid_p95_um=fr["resid_p95_vx"] * px_fine,
                              offset_vs_naive_centre_px=[float(fr["ay"] + fr["by"] * yc - k * yc),
                                                         float(fr["ax"] + fr["bx"] * xc - k * xc)],
                              affine_vs_catalog_um=dict((kk, v * px_fine) for kk, v in dev.items()),
                              status=fr["status"], nfev=fr["nfev"])

    if source is not None:
        rep = {}
        for name, mesh, A, is_copy in (("fine", fine, A_src_to_fine, fine_is_copy),
                                       ("coarse", coarse, A_src_to_coarse, coarse_is_copy)):
            if is_copy:
                rep[name] = compare_meshes(source.points, mesh.points) | {"kind": "copy"}
                continue
            pred, info = transform_tifxyz(source.points, A)
            rep[name] = compare_meshes(pred, mesh.points) | {"kind": "transform", "replica_info": {
                kk: (list(v) if isinstance(v, tuple) else v) for kk, v in info.items()}}
        rec["reproduction"] = rep

    rec["lines"] = decision_lines(rec)
    return rec


def decision_lines(rec: dict) -> dict:
    ri = rec.get("route_I", {})
    out = {}
    if ri and ri.get("frac_ok", 0) < ROUTE_I_MIN_OK:
        out["M1_in_files"] = "VOID_ROUTE_I"
    elif "offset_vs_naive_centre_norm_px" in ri:
        o = ri["offset_vs_naive_centre_norm_px"]
        out["M1_in_files"] = "IN_FILES" if o >= IN_FILES_MIN_PX else ("NOT_IN_FILES" if o <= NOT_IN_FILES_MAX_PX else "MIXED")
        out["M2b_model_matches_route_I"] = bool(ri["vs_exporter_model_px"]["p95"] <= MODEL_MATCH_P95_PX)
        sh = ri["split_half_noise_px"]
        out["M3_mesh_noise_ok"] = bool(max(sh.values()) <= MESH_NOISE_P95_PX)
    rep = rec.get("reproduction")
    if rep:
        oks = []
        for v in rep.values():
            oks.append(bool(v.get("same_shape") and v.get("err_median_vx", 9) <= REPRO_MEDIAN_VX
                            and v.get("err_p99_vx", 9) <= REPRO_P99_VX and v.get("mask_agreement", 0) >= REPRO_MASK_AGREE))
        out["M2a_reproduced"] = all(oks)
    return out
