"""Amendment 1 scorer (prereg/2026-10-11_partner_w00field_amend1.md): amended C1 from out/ct_plants_a1_<seg>.jsonl
(plants (4,0),(0,8),(12,0),(0,-12); within 0.5 px on >= 95 % of reachable tiles, every plant, every reachable tile
measured). Only where it PASSES, C2-C5 are read from the pass-1 records with the frozen bars (analyze.ct_summary)."""
import glob, json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
from analyze import ct_summary, C1_TOL, C1_FRAC

H = Path(__file__).resolve().parents[1]


def c1_a1(path, n_reach):
    recs = [json.loads(l) for l in Path(path).read_text().splitlines() if l.strip()]
    ok, n = {}, {}
    for r in recs:
        for key, v in r["plants"].items():
            want = np.array(list(map(float, key.split(","))))
            n[key] = n.get(key, 0) + 1
            ok[key] = ok.get(key, 0) + bool(np.abs(np.array(v["delta"]) - want).max() <= C1_TOL)
    frac = {k: round(ok[k] / n[k], 3) for k in n}
    complete = len(recs) == n_reach and len(n) == 4
    return frac, bool(complete and all(v >= C1_FRAC for v in frac.values())), complete


def main():
    rows = []
    for f in sorted(glob.glob(str(H / "out/ct_plants_a1_*.jsonl"))):
        seg = Path(f).stem[len("ct_plants_a1_"):]
        s = ct_summary(str(H / f"out/ct_{seg}.jsonl"))
        frac, passed, complete = c1_a1(f, s["n_reachable"])
        row = dict(segment=seg, n_reachable=s["n_reachable"], C1_a1_recovery=frac, C1_a1_complete=complete,
                   C1_a1="PASS" if passed else "FAIL")
        if passed:
            for k in ("C2_noise_p95_px", "C2", "exact_residual_px_median_p95", "naive_correction_px_median", "C3",
                      "C4_depth_median_layers", "C4_depth_median_um", "C4", "C5_field", "C5", "r_median"):
                row[k] = s.get(k)
        else:
            row["note"] = "C1 failed (or incomplete): C2-C5 not read"
        rows.append(row)
        print(json.dumps(row))
    json.dump(rows, open(H / "out/ct_summary_a1.json", "w"), indent=1)


if __name__ == "__main__":
    main()
