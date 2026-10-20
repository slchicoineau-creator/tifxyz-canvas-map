"""Phase-2 CT queue for the frozen C-set (prereg 8ca5789f). Sequential, one process. Per segment: fetch meshes ->
Stage 0b pyramid check -> check.py ct --plants -> delete meshes. Selection of (iii)-(v) from M-route records only."""
import json, shutil, subprocess, sys, time
from pathlib import Path
H = Path(__file__).resolve().parents[1]
PY = str(H / ".venv/bin/python")
PREREG = "8ca5789f"
P = {(p["sample"], p["segment"]): p for p in json.loads((H / "out/pairs.json").read_text())["pairs"]}

def mrecs():
    out = {}
    for l in (H / "out/mesh_records.jsonl").read_text().splitlines():
        r = json.loads(l); out[(r["sample"], r["segment"])] = r
    return out

def wait_for(sample):
    n = sum(1 for k in P if k[0] == sample)
    while sum(1 for k in mrecs() if k[0] == sample) < n:
        time.sleep(30)

def top(sample, n, pred=lambda p: True):
    wait_for(sample)
    R = mrecs()
    c = [(R[k]["fine_valid"], k) for k in R if k[0] == sample and pred(P[k])]
    return [k for _, k in sorted(c, reverse=True)[:n]]

def run(k):
    p = P[k]; seg = k[1]
    out = H / f"out/ct_{seg}.jsonl"
    if out.exists():
        print("skip (exists)", k, flush=True); return
    subprocess.run([PY, str(H / "check.py"), "fetch", "--samples", k[0], "--segments", seg], check=True, capture_output=True)
    pyr = subprocess.run([PY, str(H / "check.py"), "pyramid", "--samples", k[0], "--segments", seg, "--out", str(H / f"out/stage0b_pyramid_{seg}.json")], capture_output=True, text=True)
    ok = pyr.returncode == 0 and all(r["ok"] for r in json.loads((H / f"out/stage0b_pyramid_{seg}.json").read_text()))
    print("0b", k, ok, pyr.stdout.strip()[-200:], flush=True)
    if ok:
        args = ["nice", PY, "-u", str(H / "check.py"), "ct", "--samples", k[0], "--segments", seg, "--prereg", PREREG, "--plants", "--out", str(out)]
        if k == ("PHerc0841", "20260220213127"):
            args += ["--exclude-fine-box", "10208,12992,9280,15776"]
        t0 = time.time()
        with open(H / f"out/ct_{seg}.log", "w") as lf:
            r = subprocess.run(args, stdout=lf, stderr=subprocess.STDOUT)
        print("ct", k, "rc", r.returncode, round(time.time() - t0), "s", flush=True)
    if k != ("PHerc0841", "20260220213127"):
        segdir = H / "cache/meshes" / k[0] / "segments" / p["fine"]["key"].split("/")[2]
        if segdir.exists():
            shutil.rmtree(segdir)

queue = [("PHerc0841", "20260220214732"), ("PHerc0841", "20260221022814")]
for k in queue:
    run(k)
for k in top("PHerc0139", 3):
    run(k)
a = top("PHerc0814", 1, lambda p: p["coarse"]["is_copy"]); b = top("PHerc0814", 1, lambda p: p["fine"]["is_copy"])
sel0814 = a + b if (a and b) else top("PHerc0814", 2)
print("0814 selection", sel0814, flush=True)
for k in sel0814:
    run(k)
for k in top("PHerc0500P2", 2):
    run(k)
print("CSET_DONE", flush=True)
