"""Amendment 2 scorer: out/amend2/ct_plants_a1_<seg>.jsonl + out/amend2/ct_<seg>.jsonl, frozen bars (analyze.ct_summary)."""
import glob, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from analyze import ct_summary
from ct_table_a1 import c1_a1
H = Path(__file__).resolve().parents[1]; O = H / "out/amend2"
rows = []
for f in sorted(glob.glob(str(O / "ct_plants_a1_*.jsonl"))):
    seg = Path(f).stem[len("ct_plants_a1_"):]
    s = ct_summary(str(O / f"ct_{seg}.jsonl"))
    if "n_reachable" not in s or "C2_noise_p95_px" not in s:
        row = dict(segment=seg, n_tiles=s.get("n_tiles"), n_reachable=s.get("n_reachable"), verdict=s.get("verdict")); rows.append(row); print(json.dumps(row)); continue
    frac, passed, complete = c1_a1(f, s["n_reachable"])
    row = dict(segment=seg, n_tiles=s["n_tiles"], n_reachable=s["n_reachable"], C1_a1_recovery=frac, C1_a1_complete=complete, C1_a1="PASS" if passed else "FAIL")
    if passed:
        for k in ("C2_noise_p95_px", "C2", "exact_residual_px_median_p95", "naive_correction_px_median", "naive_valid_tiles", "C3", "C4_depth_median_layers", "C4_depth_median_um", "C4", "C5_field", "C5", "r_median"):
            row[k] = s.get(k)
    else:
        row["note"] = "C1 failed or incomplete: C2-C5 not read"
    rows.append(row); print(json.dumps(row))
json.dump(rows, open(O / "ct_summary_a2.json", "w"), indent=1)
