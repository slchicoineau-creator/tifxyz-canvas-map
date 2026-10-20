"""S0 of prereg 2026-10-13_partner_tifxyz_reverse_route: 100 seeds, synthetic only. Writes s0.json next to this file."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", "tests"))
from test_reverse import fine_is_source, run, A_coarse, META, PXF, PXC     # noqa: E402
from meshpair.exporter import transform_tifxyz                              # noqa: E402
from meshpair.reverse import reverse_route                                  # noqa: E402
from meshpair.synth import curved_sheet, demo_affine                        # noqa: E402
from meshpair.tifxyz import Mesh                                            # noqa: E402

out = dict(prereg="prereg/2026-10-13_partner_tifxyz_reverse_route.md", seeds=100)
p95s, fracs, inpl = [], [], []
for s in range(100):
    F, C, M = fine_is_source(s)
    r = run(F, C, M)
    p95s.append(r["p95_normal_coarse_vx"]); fracs.append(r["frac_on"]); inpl.append(r["inplane_p95_px"])
out["clean_fine_is_source"] = dict(max_p95_normal=max(p95s), median_p95_normal=float(np.median(p95s)), min_frac_on=min(fracs),
                                   max_inplane_p95=max(inpl), all_pass_pair_rule=bool(all(f >= .95 for f in fracs) and max(inpl) <= .003),
                                   gate_p95_le_0p01=bool(max(p95s) <= 0.01))
out["planted"] = {}
for shift in (0.05, 0.25, 1.0):
    fails, ip = 0, []
    for s in range(100):
        d = [(shift, 0.0), (0.0, shift), (shift / np.sqrt(2), shift / np.sqrt(2))][s % 3]
        F, C, M = fine_is_source(s, shift_px=d)
        r = run(F, C, M)
        fails += int(not r["pair_pass"]); ip.append(r["inplane_median_px"])
    out["planted"][str(shift)] = dict(fail_rate=fails / 100, median_inplane_px=float(np.median(ip)), gate_ge_0p95=bool(fails / 100 >= 0.95))
# both-transformed (descriptive)
bt = []
for s in range(20):
    S = curved_sheet(rows=110, cols=130, seed=s)
    AF, AC = demo_affine(), demo_affine(scale=1.027, rot_deg=-1.5, t=(6677.0, 6689.0, 624.0))
    Fp, _ = transform_tifxyz(S, AF); Cp, _ = transform_tifxyz(S, AC)
    M = (np.vstack([AF, [0, 0, 0, 1]]) @ np.linalg.inv(np.vstack([AC, [0, 0, 0, 1]])))[:3]
    r = run(Mesh(Fp, META), Mesh(Cp, META), M)
    bt.append((r["frac_on"], r["median_normal_coarse_vx"], r["p95_normal_coarse_vx"], r["inplane_p95_px"]))
bt = np.array(bt)
out["both_transformed_descriptive"] = dict(n=20, median_frac_on=float(np.median(bt[:, 0])), median_median_normal=float(np.median(bt[:, 1])),
                                           median_p95_normal=float(np.median(bt[:, 2])), median_inplane_p95=float(np.median(bt[:, 3])))
out["S0_PASS"] = bool(out["clean_fine_is_source"]["gate_p95_le_0p01"] and all(v["gate_ge_0p95"] for v in out["planted"].values()))
json.dump(out, open(os.path.join(HERE, "s0.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
