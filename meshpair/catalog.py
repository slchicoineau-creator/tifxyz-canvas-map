"""Enumerate fine/coarse tifxyz mesh pairs from the open-data catalog (s3://vesuvius-challenge-open-data/metadata.json).

Everything here comes from the catalog's own provenance fields: data[type == tifxyz-transformed] entries carry
parameters.target_volume and creation_info.{process, derived_from.source_path, provenance.parameters.{commit,
affine-matrix}}; data[type == layers-zarr] entries are the rendered surface volumes (num-slices, render commit).
"""
from __future__ import annotations

import gzip
import json
import math

FINE_RANGE = (2.0, 2.5)      # um
COARSE_RANGE = (7.5, 10.0)   # um
EXCLUDED_SAMPLES: set = set()   # samples to skip (none by default)


def load_catalog(path: str) -> dict:
    raw = open(path, "rb").read()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return json.loads(raw)


def _key(s3url: str | None) -> str | None:
    if not s3url:
        return None
    for pre in ("s3://vesuvius-challenge-open-data/", "s3://vesuvius-challenge/"):
        if s3url.startswith(pre):
            return s3url[len(pre):].rstrip("/")
    return None


def _affine(prov_params: dict):
    a = prov_params.get("affine-matrix")
    if not a:
        return None
    m = json.loads(a) if isinstance(a, str) else a
    return [list(map(float, r)) for r in m]


def segment_records(cat: dict):
    out = []
    for sname, s in cat["samples"].items():
        if sname in EXCLUDED_SAMPLES:
            continue
        vols = s.get("volumes", {}) or {}

        def vinfo(vid):
            v = vols.get(vid) or {}
            p = v.get("properties") or {}
            shape = p.get("shape")
            return dict(px=p.get("pixel_size_um"), shape=shape,
                        nvox=math.prod(shape) if shape else None, long_id=v.get("long_id"))

        vt = {}
        for e in (s.get("sample", {}).get("properties", {}) or {}).get("volume_transforms", []) or []:
            for t in e.get("transforms", []):
                vt[(e["from_volume_id"], t["to_volume_id"])] = t["matrix"]
        for segid, seg in (s.get("segments") or {}).items():
            meshes, svs = [], []
            for e in seg.get("data", []):
                ci = e.get("creation_info") or {}
                pp = (ci.get("provenance") or {}).get("parameters") or {}
                tv = (e.get("parameters") or {}).get("target_volume")
                path = e["origins"][0]["path"].rstrip("/") if e.get("origins") else None
                if e.get("type") == "tifxyz-transformed":
                    vi = vinfo(tv)
                    meshes.append(dict(target_volume=tv, px=vi["px"], vol_nvox=vi["nvox"], vol_long_id=vi["long_id"],
                                       key=path, process=ci.get("process"), commit=pp.get("commit"),
                                       affine_src_to_this=_affine(pp),
                                       source_key=_key((ci.get("derived_from") or {}).get("source_path")),
                                       is_copy=ci.get("process") == "copy-normalized-for-original"))
                elif e.get("type") == "layers-zarr":
                    svs.append(dict(target_volume=tv, key=path, num_slices=int(pp["num-slices"]) if pp.get("num-slices") else None,
                                    render_commit=pp.get("commit"), render_scale=pp.get("scale"),
                                    derived_from_key=_key((ci.get("derived_from") or {}).get("source_path"))))
            out.append(dict(sample=sname, segment=segid, long_id=seg.get("long_id"),
                            original_volume=seg.get("original_volume_id"), meshes=meshes, surface_volumes=svs,
                            volume_transforms={f"{a}->{b}": m for (a, b), m in vt.items()}))
    return out


def _in(px, rng):
    return px is not None and rng[0] <= px <= rng[1]


def pairs(cat: dict, fine_range=FINE_RANGE, coarse_range=COARSE_RANGE):
    """One record per segment that has >= 1 fine and >= 1 coarse mesh. 'primary' = the mesh on the largest volume
    (by voxel count) in each class; all combinations are listed under 'all_combos'."""
    res = []
    for seg in segment_records(cat):
        fine = [m for m in seg["meshes"] if _in(m["px"], fine_range)]
        coarse = [m for m in seg["meshes"] if _in(m["px"], coarse_range)]
        if not fine or not coarse:
            continue
        pick = lambda ms: sorted(ms, key=lambda m: (-(m["vol_nvox"] or 0), m["target_volume"]))[0]
        f, c = pick(fine), pick(coarse)
        sv = {v["target_volume"]: v for v in seg["surface_volumes"]}
        src = f["source_key"] if not f["is_copy"] else None
        if src is None and not c["is_copy"]:
            src = c["source_key"]
        res.append(dict(sample=seg["sample"], segment=seg["segment"], long_id=seg["long_id"],
                        original_volume=seg["original_volume"], fine=f, coarse=c,
                        source_key=f["source_key"] or c["source_key"],
                        fine_sv=sv.get(f["target_volume"]), coarse_sv=sv.get(c["target_volume"]),
                        all_combos=[[a["key"], b["key"]] for a in fine for b in coarse],
                        volume_transforms=seg["volume_transforms"]))
    return res


def coarse_to_fine_affine(pair: dict):
    """M (3x4) mapping coarse-volume coordinates to fine-volume coordinates, from the meshes' own provenance:
    M = A_src->fine o (A_src->coarse)^-1 (identity for a plain copy of the source)."""
    import numpy as np
    from .meshcheck import compose_coarse_to_fine
    f, c = pair["fine"], pair["coarse"]
    Af = None if f["is_copy"] else f["affine_src_to_this"]
    Ac = None if c["is_copy"] else c["affine_src_to_this"]
    if (not f["is_copy"] and Af is None) or (not c["is_copy"] and Ac is None):
        return None
    return np.asarray(compose_coarse_to_fine(Af, Ac))
