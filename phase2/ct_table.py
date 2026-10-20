import json, sys, glob
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from analyze import ct_summary
rows = []
for f in sorted(glob.glob(str(Path(__file__).parents[1] / "out/ct_2*.jsonl"))):
    s = ct_summary(f); rows.append(s)
    seg = Path(f).stem[3:]
    if "C1" not in s:
        print(seg, s.get("verdict"), s.get("n_tiles"), s.get("n_reachable")); continue
    print(f"| {seg} | {s['n_tiles']} / {s['n_reachable']} | {s['C1']} ({s['C1_plant_recovery']}) | {s['C2_noise_p95_px']} {s['C2']} | "
          f"{s['exact_residual_px_median_p95']} vec {s['exact_residual_vec_median_px']}; naive {s['naive_correction_px_median']} ({s['naive_valid_tiles']} valid) -> {s['C3']} | "
          f"{s['C4_depth_median_um']} um = {s['C4_depth_median_layers']} L, IQR {s['C4_depth_iqr_um']} -> {s['C4']} | {s['C5_field']['loo_p95_px']} -> {s['C5']} | r {s['r_median']} | inval {s['exact_invalid_reasons']} |")
json.dump(rows, open(Path(__file__).parents[1] / "out/ct_summary.json", "w"), indent=1)
