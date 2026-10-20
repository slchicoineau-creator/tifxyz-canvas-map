"""Phase-2 readout against the frozen bars of prereg/2026-10-11_partner_w00field.md (commit 8ca5789f).
Reads out/mesh_records.jsonl and out/ct_<seg>.jsonl; prints JSON summaries (RESULTS.md is written by hand from them)."""
import collections
import json
import sys
from pathlib import Path

import numpy as np

H = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(H))
from meshpair.ct import fit_field  # noqa: E402

# frozen bars (prereg)
FAM_IN, FAM_OUT, FAM_ROOT = 0.90, 0.90, 0.95
C2_PASS, C2_MARG = 0.67, 1.0
C3_RES, C3_RATIO = 2.0, 3.0
C1_TOL, C1_FRAC = 0.5, 0.95
C4_LAYERS = 1.0
C5_LOO = 2.0
MIN_REACH = 5


def mesh_summary():
    recs = {}
    for l in (H / "out/mesh_records.jsonl").read_text().splitlines():
        r = json.loads(l)
        recs[(r["sample"], r["segment"])] = r          # last record wins (none duplicated expected)
    fams = collections.defaultdict(list)
    for r in recs.values():
        fams[(r["sample"], r["px_fine"], r["px_coarse"])].append(r)
    out = {}
    for k, rs in sorted(fams.items()):
        lines = [r.get("lines", {}) for r in rs]
        m1 = collections.Counter(l.get("M1_in_files", "NONE") for l in lines)
        noisy = sum(1 for l in lines if l.get("M3_mesh_noise_ok") is False)
        nonvoid = [l for l in lines if l.get("M1_in_files") in ("IN_FILES", "NOT_IN_FILES", "MIXED")
                   and l.get("M3_mesh_noise_ok") is not False]
        n = len(nonvoid)
        fin = sum(l["M1_in_files"] == "IN_FILES" for l in nonvoid)
        fout = sum(l["M1_in_files"] == "NOT_IN_FILES" for l in nonvoid)
        if n == 0:
            v1 = "VOID (no non-void pair)"
        elif fin / n >= FAM_IN:
            v1 = "CONFIRMED IN FILES"
        elif fout / n >= FAM_OUT:
            v1 = "OUR ASSUMPTION"
        else:
            v1 = "MIXED"
        a = [l.get("M2a_reproduced") for l in lines]
        b = [l.get("M2b_model_matches_route_I") for l in lines]
        both = sum(1 for x, y in zip(a, b) if x is True and y is True)
        bonly = sum(1 for y in b if y is True)
        N = len(lines)
        b_eval = sum(1 for y in b if y is not None)
        if b_eval == 0:
            v2 = "NOT EXPLAINED by the frozen tree (M2b not evaluable: route I VOID on every pair)"
        elif both / N >= FAM_ROOT:
            v2 = "EXPORTER CONVENTION CONFIRMED"
        elif bonly / N >= FAM_ROOT:
            v2 = "PARTIAL"
        else:
            v2 = "NOT EXPLAINED"
        offs = [r["route_I"]["offset_vs_naive_centre_norm_px"] for r in rs if "offset_vs_naive_centre_norm_px" in r.get("route_I", {})]
        offs_um = [np.hypot(*r["route_I"]["offset_vs_naive_centre_um"]) for r in rs if "offset_vs_naive_centre_um" in r.get("route_I", {})]
        model = [np.hypot(r["exporter_model"]["y_offset_box_centre_um"], r["exporter_model"]["x_offset_box_centre_um"]) for r in rs]
        drift = [max(abs(r["exporter_model"]["y_b_minus_k"]), abs(r["exporter_model"]["x_b_minus_k"])) / r["exporter_model"]["k_naive"] for r in rs]
        rangeum = [max(abs(r["exporter_model"][f"{ax}_offset_box_{e}_um"]) for ax in "yx" for e in ("min", "max")) for r in rs]
        fracok = [r["route_I"].get("frac_ok", np.nan) for r in rs]
        vsm = [r["route_I"]["vs_exporter_model_px"]["p95"] for r in rs if "vs_exporter_model_px" in r.get("route_I", {})]
        rfdev = [r["route_F"]["affine_vs_catalog_um"]["max"] for r in rs if "route_F" in r]
        rf_vs_model = []
        for r in rs:
            if "route_F" in r:
                ey, ex = r["route_F"]["offset_vs_naive_centre_px"]
                rf_vs_model.append(float(np.hypot(ey - r["exporter_model"]["y_offset_box_centre_px"], ex - r["exporter_model"]["x_offset_box_centre_px"])))
        rep_err = [max(v.get("err_p99_vx", np.nan) for v in r["reproduction"].values()) for r in rs if "reproduction" in r]
        mask = [min(v.get("mask_agreement", np.nan) for v in r["reproduction"].values()) for r in rs if "reproduction" in r]
        q = lambda x: [round(float(np.nanmin(x)), 4), round(float(np.nanmedian(x)), 4), round(float(np.nanmax(x)), 4)] if len(x) else None
        out["/".join(map(str, k))] = dict(
            n=N, M1_counts=dict(m1), noisy=noisy, n_nonvoid=n, M1_family=v1, M2a_pass=sum(x is True for x in a),
            M2b_pass=bonly, M2b_evaluable=b_eval, M2_both=both, M2_family=v2,
            routeI_offset_centre_px_min_med_max=q(offs), routeI_offset_centre_um_min_med_max=q(offs_um),
            model_offset_centre_um_min_med_max=q(model), model_drift_rel_min_med_max=q(drift),
            model_offset_max_abs_over_valid_box_um=q(rangeum), routeI_frac_ok=q(fracok), routeI_vs_model_p95_px=q(vsm),
            routeF_vs_model_centre_px=q(rf_vs_model), routeF_affine_vs_catalog_max_um=q(rfdev),
            repro_err_p99_vx_max=q(rep_err), repro_mask_agree_min=q(mask))
    return out


def ct_summary(seg_file):
    tiles = [json.loads(l) for l in Path(seg_file).read_text().splitlines()]
    T = [t for t in tiles if t.get("kind") == "tile"]
    ex = {tuple(t["centre_c"]): t for t in T if t["map"] == "exact"}
    nv = {tuple(t["centre_c"]): t for t in T if t["map"] == "naive"}
    reach = [c for c, t in ex.items() if t["valid"]]
    res = dict(file=str(seg_file), n_tiles=len(ex), n_reachable=len(reach))
    if len(reach) < MIN_REACH:
        res["verdict"] = "VOID (< 5 reachable tiles)"
        res["exact_invalid_reasons"] = _reasons([ex[c] for c in ex])
        return res
    plants = [t for t in tiles if t.get("kind") == "plants"]
    ok_by = collections.Counter(); n_by = collections.Counter()
    for p in plants:
        for key, v in p["plants"].items():
            want = np.array(list(map(float, key.split(","))))
            n_by[key] += 1
            ok_by[key] += bool(np.abs(np.array(v["delta"]) - want).max() <= C1_TOL)
    c1 = {k: round(ok_by[k] / n_by[k], 3) for k in n_by}
    c1_pass = bool(n_by) and all(v >= C1_FRAC for v in c1.values()) and len(plants) == len(reach)
    res["C1_plant_recovery"] = c1
    res["C1"] = "PASS" if c1_pass else "FAIL"
    noise = np.array([ex[c]["noise_px"] for c in reach])
    p95 = float(np.percentile(noise, 95))
    res["C2_noise_p95_px"] = round(p95, 3)
    res["C2"] = "PASS" if p95 <= C2_PASS else ("MARGINAL" if p95 <= C2_MARG else "ABOVE 1.0 (bars x3 p95)")
    scale = 1.0 if p95 <= C2_MARG else 3 * p95 / 1.0
    resid = np.array([np.hypot(*ex[c]["map_correction_px"]) for c in reach])
    naive = np.array([np.hypot(*nv[c]["map_correction_px"]) for c in reach if c in nv])
    px_c = None
    res["exact_residual_px_median_p95"] = [round(float(np.median(resid)), 3), round(float(np.percentile(resid, 95)), 3)]
    res["exact_residual_vec_median_px"] = [round(float(np.median([ex[c]["map_correction_px"][i] for c in reach])), 3) for i in (0, 1)]
    res["naive_correction_px_median"] = round(float(np.median(naive)), 3) if naive.size else None
    res["naive_valid_tiles"] = int(sum(nv[c]["valid"] for c in reach if c in nv))
    c3_bar = C3_RES * (scale if p95 > C2_MARG else 1.0)
    res["C3"] = "EXPLAINED" if (np.median(resid) <= c3_bar and naive.size and np.median(naive) >= C3_RATIO * np.median(resid)) else "NOT EXPLAINED"
    dl = np.array([ex[c]["depth_offset_layers"] for c in reach])
    du = np.array([ex[c]["depth_offset_um"] for c in reach])
    res["C4_depth_median_layers"] = round(float(np.median(dl)), 3)
    res["C4_depth_median_um"] = round(float(np.median(du)), 2)
    res["C4_depth_iqr_um"] = [round(float(np.percentile(du, 25)), 2), round(float(np.percentile(du, 75)), 2)]
    res["C4_depth_on_edge_tiles"] = int(sum(ex[c]["depth_on_edge"] for c in reach))
    res["C4"] = "DEPTH OFFSET REAL" if abs(np.median(dl)) > C4_LAYERS else "centre-to-centre stands"
    ff = fit_field([ex[c] for c in reach])
    res["C5_field"] = {k: (round(v, 3) if isinstance(v, float) else v) for k, v in ff.items() if k in ("n", "loo_p95_px", "loo_median_px", "insample_p95_px")}
    c5_bar = C5_LOO * (scale if p95 > C2_MARG else 1.0)
    res["C5"] = ("ONE CORRECTION PER SEGMENT" if ff.get("loo_p95_px", 1e9) <= c5_bar else "LOCAL FIELD REQUIRED") if ff.get("n", 0) >= 4 else "n<4"
    res["r_median"] = round(float(np.median([ex[c]["r"] for c in reach])), 3)
    res["exact_invalid_reasons"] = _reasons([ex[c] for c in ex if not ex[c]["valid"]])
    return res


def _reasons(ts):
    c = collections.Counter()
    for t in ts:
        if t["valid"]:
            continue
        if t["r"] < 0.30: c["r<0.30"] += 1
        if t["runner_up"] > 0.70 * t["r"]: c["runner_up"] += 1
        if t["on_edge"]: c["on_edge"] += 1
        if t["coarse_valid_frac"] < 0.9 or t["fine_valid_frac"] < 0.9: c["valid_frac<0.9"] += 1
        if t.get("split_spread_px", 0) > 1.0: c["split_spread"] += 1
        if t.get("autocorr_at_shift") is not None and t["autocorr_at_shift"] > 0.30: c["autocorr"] += 1
    return dict(c)


if __name__ == "__main__":
    what = sys.argv[1]
    if what == "mesh":
        print(json.dumps(mesh_summary(), indent=1))
    else:
        print(json.dumps(ct_summary(what), indent=1))
