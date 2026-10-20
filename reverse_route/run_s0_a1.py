"""Amendment 1, synthetic stage S0-A1 (realistic model of the diagnosed cause; no real data).
Expected false-fail and planted-shift detection rates of the lattice-clean rule -> s0_a1.json."""
import json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..")); sys.path.insert(0, HERE)
from synth_a1 import CONFIGS, META, PXF, PXC, jittered_source, make_pair    # noqa: E402
from meshpair.reverse_lc import run_pair                                   # noqa: E402
from meshpair.reverse import reverse_route                                 # noqa: E402
from meshpair.tifxyz import Mesh                                           # noqa: E402

ZS = (-11000.0, 4000.0, 20000.0, 52000.0)      # coordinate magnitude cycle (fine coords up to ~2.7e5 vx)


def one(S, name, shift=(0.0, 0.0), orig=False):
    AF, AC = CONFIGS[name]
    F, C, M, fc, cc = make_pair(S, AF, AC, shift)
    Fm, Cm = Mesh(F, META), Mesh(C, META)
    r = run_pair(F, Fm.valid, C, Cm.valid, Fm.geom(PXF), Cm.geom(PXC), M, S.shape[:2], PXF, PXC, coarse_is_copy=cc, fine_is_copy=fc)
    if orig:
        o = reverse_route(F, Fm.valid, C, Cm.valid, Fm.geom(PXF), Cm.geom(PXC), M, PXF, PXC)
        r["orig_pair_pass"] = o["pair_pass"]; r["orig_frac_on"] = o["frac_on"]
    return r


def summ(rs, planted=False):
    v = [r["verdict"] for r in rs]
    d = dict(n=len(rs), n_pass=v.count("PASS"), n_fail=v.count("FAIL"), n_not_evaluable=v.count("NOT_EVALUABLE"),
             max_inplane_p95_over_bar=float(max(r["inplane_p95_clean"] / r["inplane_bar_px"] for r in rs)),
             max_p95_normal_clean=float(max(r["p95_normal_clean"] for r in rs)), min_frac_on_clean=float(min(r["frac_on_clean"] for r in rs)),
             median_frac_clean=float(np.median([r["frac_clean"] for r in rs])), max_bar_px=float(max(r["inplane_bar_px"] for r in rs)))
    if planted:
        d["detection_rate"] = v.count("FAIL") / len(rs)
        d["median_inplane_median_clean"] = float(np.median([r["inplane_median_clean"] for r in rs]))
    else:
        d["false_fail_rate"] = (len(rs) - v.count("PASS")) / len(rs)
    if "orig_pair_pass" in rs[0]:
        d["orig_rule_pass_rate"] = float(np.mean([r["orig_pair_pass"] for r in rs]))
        d["orig_median_frac_on"] = float(np.median([r["orig_frac_on"] for r in rs]))
    return d


if __name__ == "__main__":
    t0 = time.time(); out = dict(stage="S0-A1", seeds="0-99", jitter_source_vx=[0, 1, 3], z_shift_cycle=ZS, clean={}, planted={}, extreme={})
    for name in CONFIGS:
        for j in (0.0, 1.0, 3.0):
            rs = [one(jittered_source(s, jitter=j, z_shift=ZS[s % 4]), name, orig=True) for s in range(100)]
            out["clean"][f"{name}|jitter{j:g}"] = summ(rs); print(name, j, out["clean"][f"{name}|jitter{j:g}"], flush=True)
        for sh in (0.05, 0.25, 1.0):
            rs = []
            for s in range(100):
                d = [(sh, 0.0), (0.0, sh), (sh / np.sqrt(2), sh / np.sqrt(2))][s % 3]
                rs.append(one(jittered_source(s, jitter=1.0, z_shift=ZS[s % 4]), name, shift=d))
            out["planted"][f"{name}|shift{sh:g}"] = summ(rs, True); print(name, sh, out["planted"][f"{name}|shift{sh:g}"], flush=True)
    # extreme size: 1700 x 1700 source (fine ~6.8k x 6.8k = 46 M vertices, like the 1667 flatboi), fine coords ~6.5e4
    for name in ("both_ratio4", "coarse_copy"):
        rs = [one(jittered_source(s, rows=1700, cols=1700, spacing=9.0, jitter=1.0, z_shift=-11000.0), name) for s in range(3)]
        ps = [one(jittered_source(s, rows=1700, cols=1700, spacing=9.0, jitter=1.0, z_shift=-11000.0), name, shift=(0.05, 0.0)) for s in range(3)]
        out["extreme"][name] = dict(clean=summ(rs), planted_0p05=summ(ps, True)); print(name, out["extreme"][name], flush=True)
    allc = list(out["clean"].values()) + [v["clean"] for v in out["extreme"].values()]
    allp = list(out["planted"].values()) + [v["planted_0p05"] for v in out["extreme"].values()]
    out["overall"] = dict(n_clean_runs=sum(c["n"] for c in allc), false_fail=sum(c["n"] - c["n_pass"] for c in allc),
                          n_planted_runs=sum(p["n"] for p in allp), planted_detected=sum(p["n_fail"] for p in allp),
                          min_detection_rate=min(p["detection_rate"] for p in allp))
    out["seconds"] = round(time.time() - t0)
    json.dump(out, open(os.path.join(HERE, "s0_a1.json"), "w"), indent=1)
    print(json.dumps(out["overall"]))
