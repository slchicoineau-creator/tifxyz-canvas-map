"""Phase-2 driver for the frozen M-set (prereg/2026-10-11_partner_w00field.md, commit 8ca5789f).
Per segment: check.py fetch -> check.py mesh --prereg 8ca5789f -> delete that segment's meshes (manifest kept).
Does not change the frozen tool; resumable (skips segments already in the records file)."""
import json, shutil, subprocess, sys, time
from pathlib import Path
H = Path(__file__).resolve().parents[1]
PY = str(H / ".venv/bin/python")
PREREG = "8ca5789f"
REC = H / "out/mesh_records.jsonl"
P = json.loads((H / "out/pairs.json").read_text())["pairs"]
DEV = {("PHerc0841", "20260220213127")}
order = sorted([p for p in P if (p["sample"], p["segment"]) not in DEV], key=lambda p: (p["sample"] == "PHercParis4", p["sample"], p["segment"]))
done = set()
if REC.exists():
    for l in REC.read_text().splitlines():
        r = json.loads(l); done.add((r["sample"], r["segment"]))
for i, p in enumerate(order):
    if (p["sample"], p["segment"]) in done:
        continue
    t0 = time.time()
    f = subprocess.run([PY, str(H / "check.py"), "fetch", "--segments", p["segment"], "--samples", p["sample"]], capture_output=True, text=True)
    if f.returncode:
        print("FETCH_FAIL", p["sample"], p["segment"], f.stderr[-500:], flush=True); continue
    m = subprocess.run(["nice", PY, str(H / "check.py"), "mesh", "--segments", p["segment"], "--samples", p["sample"], "--prereg", PREREG, "--out", str(REC)], capture_output=True, text=True)
    print(i, p["sample"], p["long_id"], m.stdout.strip()[-300:], ("ERR " + m.stderr[-800:]) if m.returncode else "", round(time.time() - t0, 1), "s", flush=True)
    segdir = H / "cache/meshes" / p["sample"] / "segments" / p["fine"]["key"].split("/")[2]
    if segdir.exists() and m.returncode == 0:
        shutil.rmtree(segdir)
print("ALL_DONE", flush=True)
