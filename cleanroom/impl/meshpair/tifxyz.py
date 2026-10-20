# SPDX-License-Identifier: MIT
# Authors: Chris Thomas; clean-room implementation by Claude Code (Opus 5.5), implementer session 2026-10-18.
"""tifxyz mesh read / write / fetch (folder with meta.json, x.tif, y.tif, z.tif).

clean-room implementation from cleanroom/SPEC.md v1 (2026-10-18), Part C.
"""
from __future__ import annotations

import hashlib
import json
import os
import urllib.request
from pathlib import Path

import numpy as np
import tifffile

from .conventions import MeshGeom

S3_HTTP = "https://vesuvius-challenge-open-data.s3.amazonaws.com"
FILES = ("meta.json", "x.tif", "y.tif", "z.tif")


class Mesh:
    """A tifxyz mesh: points float32 (rows, cols, 3) with -1 = invalid, plus its meta.json dict."""

    def __init__(self, points, meta, path=""):
        self.points = points
        self.meta = meta
        self.path = path

    @property
    def valid(self):
        p = self.points
        return (p[..., 0] != p.dtype.type(-1)) & np.isfinite(p).all(axis=-1)

    def geom(self, px_um=float("nan")):
        sc = self.meta.get("scale", [0.05, 0.05])
        rows, cols = self.points.shape[:2]
        return MeshGeom(rows, cols, scale_x=float(sc[0]), scale_y=float(sc[1]), px_um=px_um)


def read_tifxyz(path, vc_load=False) -> Mesh:
    path = Path(path)
    meta = json.loads((path / "meta.json").read_text())
    pts = np.stack([tifffile.imread(path / f"{c}.tif").astype(np.float32) for c in "xyz"], axis=-1)
    if vc_load:
        pts[pts[..., 2] <= 0] = np.float32(-1.0)
    return Mesh(pts, meta, str(path))


def write_tifxyz(path, points, meta):
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    points = np.asarray(points)
    for i, c in enumerate("xyz"):
        tifffile.imwrite(path / f"{c}.tif", np.ascontiguousarray(points[..., i], dtype=np.float32))
    (path / "meta.json").write_text(json.dumps(meta, indent=1))


def fetch_tifxyz(s3_key_prefix, cache_root, base=S3_HTTP):
    prefix = str(s3_key_prefix).strip("/")
    local_dir = Path(cache_root) / prefix
    local_dir.mkdir(parents=True, exist_ok=True)
    shas = {}
    for name in FILES:
        dst = local_dir / name
        if not dst.exists() or dst.stat().st_size == 0:
            part = local_dir / (name + ".part")
            with urllib.request.urlopen(f"{base}/{prefix}/{name}", timeout=120) as r, open(part, "wb") as f:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
            os.replace(part, dst)
        h = hashlib.sha256()
        with open(dst, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        shas[name] = h.hexdigest()
    return local_dir, shas
