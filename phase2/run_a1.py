"""Amendment 1 queue: plants-only on the 10 frozen C-set segments, sequential. fetch -> plants -> delete meshes (keep w00)."""
import json, shutil, subprocess, sys, time
from pathlib import Path
H = Path(__file__).resolve().parents[1]
PY = str(H / ".venv/bin/python")
P = {p["segment"]: p for p in json.loads((H / "out/pairs.json").read_text())["pairs"]}
SEGS = "20260220213127 20260220214732 20260221022814 20250611170257 20250628074500 20260226000000 20250928235953 20260206000002 20251226000000 20260227000000".split()
for seg in SEGS:
    p = P[seg]; s = p["sample"]
    subprocess.run([PY, str(H / "check.py"), "fetch", "--samples", s, "--segments", seg], check=True, capture_output=True)
    t0 = time.time()
    r = subprocess.run(["nice", PY, "-u", str(H / "check.py"), "plants", "--samples", s, "--segments", seg, "--prereg", "13df5696"])
    print("plants", s, seg, "rc", r.returncode, round(time.time() - t0), "s", flush=True)
    if seg != "20260220213127":
        d = H / "cache/meshes" / s / "segments" / p["fine"]["key"].split("/")[2]
        if d.exists(): shutil.rmtree(d)
print("A1_DONE", flush=True)
