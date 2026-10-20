"""Frozen scorer of the amendment-1 APPLICATION to the 154 non-S2 catalog pairs. apply154.jsonl (+ s2.jsonl for the full-catalog table) -> summary.json."""
import collections, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
pairs = json.load(open(os.path.join(HERE, "apply154_pairs.json")))["pairs"]
setof = {(p["sample"], p["segment"]): p["set"] for p in pairs}
R = {}
for ln in open(os.path.join(HERE, "apply154.jsonl")):
    r = json.loads(ln); R[(r["sample"], r["segment"])] = r
S2 = {}
for ln in open(os.path.join(HERE, "..", "s2", "s2.jsonl")):
    r = json.loads(ln); S2[(r["sample"], r["segment"])] = r
def rec(r):
    lc = (r or {}).get("lattice_clean") or {}
    v = lc.get("verdict", "ERROR") if r else "MISSING"
    ipp = lc.get("inplane_p95_clean"); foc = lc.get("frac_on_clean"); nc = lc.get("n_clean")
    plain = bool(v in ("PASS", "FAIL") and ipp is not None and foc is not None and foc >= 0.95 and ipp <= 0.003)
    return dict(verdict=v, plain003_pass=plain, n_clean=nc, n_clean_off=(None if nc is None else nc - lc.get("n_clean_on", 0)),
                n_clean_nonconverged=(None if nc is None else nc - lc.get("n_clean_converged", 0)),
                frac_clean=lc.get("frac_clean"), frac_on_clean=foc, inplane_p95_clean=ipp, bar=lc.get("inplane_bar_px"),
                p95_normal_clean=lc.get("p95_normal_clean"), offset_norm_px=lc.get("offset_norm_px"), in_files=lc.get("in_files"),
                direction=lc.get("direction"), error=(r or {}).get("error"))
rows = []
for k in sorted(setof): rows.append(dict(sample=k[0], segment=k[1], set=setof[k], **rec(R.get(k))))
for k in sorted(S2): rows.append(dict(sample=k[0], segment=k[1], set="s2", **rec(S2[k])))
def agg(rs):
    c = collections.Counter(r["verdict"] for r in rs)
    ev = [r for r in rs if r["verdict"] in ("PASS", "FAIL")]
    return dict(n=len(rs), counts=dict(c), plain003_pass=sum(r["plain003_pass"] for r in rs),
                in_files_among_evaluable=sum(1 for r in ev if r["in_files"]),
                sum_n_clean=sum(r["n_clean"] or 0 for r in rs), sum_n_clean_off=sum(r["n_clean_off"] or 0 for r in ev),
                pairs_with_clean_off=sum(1 for r in ev if r["n_clean_off"]))
by = collections.defaultdict(list)
for r in rows: by[(r["set"], r["sample"])].append(r)
out = dict(apply154=agg([r for r in rows if r["set"] != "s2"]), s2=agg([r for r in rows if r["set"] == "s2"]), all206=agg(rows),
           by_set_sample={f"{a}/{b}": agg(v) for (a, b), v in sorted(by.items())},
           not_evaluable=[r for r in rows if r["verdict"] == "NOT_EVALUABLE"],
           fails=[r for r in rows if r["verdict"] == "FAIL"], errors=[r for r in rows if r["verdict"] not in ("PASS", "FAIL", "NOT_EVALUABLE")])
json.dump(out, open(os.path.join(HERE, "summary.json"), "w"), indent=1)
with open(os.path.join(HERE, "pairs_table.jsonl"), "w") as f:
    for r in rows: f.write(json.dumps(r) + "\n")
print(json.dumps({k: out[k] for k in ("apply154", "s2", "all206")}, indent=1))
