"""Amendment 2 queue (prereg fcb2012f): fresh C-set. Per segment: fetch -> Stage 0b -> ct pass 1 (no plants) -> plants (A1 plants) -> delete meshes.
Output out/amend2/. Reserves only for a Stage-0b failure (named in the prereg)."""
import json, shutil, subprocess, time
from pathlib import Path
H = Path(__file__).resolve().parents[1]
PY = str(H / ".venv/bin/python")
PREREG = "fcb2012f"
O = H / "out/amend2"; O.mkdir(parents=True, exist_ok=True)
P = {p["segment"]: p for p in json.loads((H / "out/pairs.json").read_text())["pairs"]}
QUEUE = [("PHerc0139", "20260220000000", ["20260127000000", "20260210000000"]), ("PHerc0139", "20260306000001", ["20260127000000", "20260210000000"]),
         ("PHerc0500P2", "20250825181922", ["20250825181859"]), ("PHerc0500P2", "20250825181907", ["20250825181859"]),
         ("PHerc0009B", "20250919132115", ["20250919125754"]), ("PHerc0009B", "20250919135433", ["20250919125754"]),
         ("PHerc0814", "20260226121302", []), ("PHerc0814", "20250926051122", ["20250925204843"])]
used = set()

def rm(seg):
    p = P[seg]; d = H / "cache/meshes" / p["sample"] / "segments" / p["fine"]["key"].split("/")[2]
    if d.exists(): shutil.rmtree(d)

def run(s, seg):
    jl = O / f"ct_{seg}.jsonl"
    if (O / f"ct_plants_a1_{seg}.jsonl").exists() and (O / f"done_{seg}").exists():
        print("skip done", seg, flush=True); return True
    subprocess.run([PY, str(H / "check.py"), "fetch", "--samples", s, "--segments", seg], check=True, capture_output=True)
    pf = O / f"stage0b_pyramid_{seg}.json"
    r = subprocess.run([PY, str(H / "check.py"), "pyramid", "--samples", s, "--segments", seg, "--out", str(pf)], capture_output=True, text=True)
    ok = r.returncode == 0 and pf.exists() and all(x["ok"] and x["nonzero_frac"] > 0.5 for x in json.loads(pf.read_text()))
    print("0b", s, seg, ok, flush=True)
    if not ok:
        rm(seg); return False
    if not jl.exists():
        t0 = time.time()
        with open(O / f"ct_{seg}.log", "w") as lf:
            rc = subprocess.run(["nice", PY, "-u", str(H / "check.py"), "ct", "--samples", s, "--segments", seg, "--prereg", PREREG, "--out", str(jl)], stdout=lf, stderr=subprocess.STDOUT).returncode
        print("ct", seg, "rc", rc, round(time.time() - t0), "s", flush=True)
        if rc: return True
    t0 = time.time()
    rc = subprocess.run(["nice", PY, "-u", str(H / "check.py"), "plants", "--samples", s, "--segments", seg, "--prereg", PREREG, "--out-dir", str(O)],
                        stdout=open(O / f"plants_{seg}.log", "a"), stderr=subprocess.STDOUT).returncode
    print("plants", seg, "rc", rc, round(time.time() - t0), "s", flush=True)
    if rc == 0: (O / f"done_{seg}").write_text("ok")
    rm(seg); return True

for s, seg, reserves in QUEUE:
    cands = [seg] + reserves
    for c in cands:
        if c in used: continue
        used.add(c)
        if run(s, c): break
        print("Stage 0b failed for", c, "-> reserve", flush=True)
print("A2_DONE", flush=True)
