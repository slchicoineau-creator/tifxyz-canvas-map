"""Frozen scorer of amendment 1, stage S2 (prereg/2026-10-13_partner_tifxyz_reverse_amend1.md). Reads s2.jsonl -> summary.json."""
import collections, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
want = {(p["sample"], p["segment"]) for p in json.load(open(os.path.join(HERE, "s2_pairs.json")))["pairs"]}
R = {}
IN = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "s2.jsonl")
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "summary.json")
for ln in open(IN):
    r = json.loads(ln); R[(r["sample"], r["segment"])] = r
missing = sorted(want - set(R))
v = {k: (R[k].get("lattice_clean") or {}).get("verdict", "ERROR") if k in R else "MISSING" for k in want}
c = collections.Counter(v.values())
n = len(want)
n_pass, n_fail = c.get("PASS", 0), c.get("FAIL", 0)
verdict = "CONFIRMED" if n_pass >= 47 else ("REFUTED" if n_fail >= 6 else "INCONCLUSIVE")
fam = {}
for k in sorted(want):
    s = k[0]; lc = (R.get(k) or {}).get("lattice_clean") or {}
    f = fam.setdefault(s, dict(n=0, PASS=0, FAIL=0, NOT_EVALUABLE=0, other=0, in_files_among_pass=0, frac_clean=[], inplane_p95_clean=[],
                               frac_on_clean=[], offset_norm_px=[]))
    f["n"] += 1
    if v[k] in ("PASS", "FAIL", "NOT_EVALUABLE"):
        f[v[k]] += 1
    else:
        f["other"] += 1
    if v[k] == "PASS" and lc.get("in_files"):
        f["in_files_among_pass"] += 1
    for kk in ("frac_clean", "inplane_p95_clean", "frac_on_clean", "offset_norm_px"):
        if lc.get(kk) is not None and np.isfinite(lc.get(kk)):
            f[kk].append(lc[kk])
for s, f in fam.items():
    for kk in ("frac_clean", "inplane_p95_clean", "frac_on_clean", "offset_norm_px"):
        a = np.array(f[kk]); f[kk] = dict(min=float(a.min()), median=float(np.median(a)), max=float(a.max())) if a.size else None
out = dict(n=n, counts=dict(c), verdict=verdict, bars=dict(confirmed="PASS >= 47/52", refuted="FAIL >= 6/52"), missing=missing,
           families=fam, fails=[dict(sample=k[0], segment=k[1], **{kk: (R[k].get("lattice_clean") or {}).get(kk) for kk in
                    ("frac_clean", "frac_on_clean", "inplane_p95_clean", "inplane_bar_px", "median_normal_clean", "n_clean")})
                    for k in sorted(want) if v[k] in ("FAIL", "NOT_EVALUABLE", "ERROR")])
json.dump(out, open(OUT, "w"), indent=1)
print(json.dumps(out, indent=1))
