#!/usr/bin/env python
"""W00-FIELD / mesh-pair offset checker CLI (CPU only; meshes are small; surface volumes are read tile-wise).

  inventory  catalog -> pairs.json (every segment with a fine ~2.4 um and a coarse ~7.9-9.4 um tifxyz mesh)
  fetch      download the fine, coarse and source meshes of selected pairs (sha256 manifest)
  mesh       mesh route (M1-M3) for selected pairs -> JSONL
  ct         CT residual stage (C1-C4) on a tile lattice of one pair -> JSONL
  plants     Amendment 1: re-measure only the new C1 plants on the frozen reachable tiles (needs --prereg 13df5696)
  reverse    reverse route (coarse vertex -> fine surface; prereg 2026-10-13)
  reverse-lc lattice-clean route (amendment 1 of the 2026-10-13 prereg; source grid shape from the provenance source mesh)
  selftest   synthetic self-tests + plumbing on the development pair (PHerc0841 w00)

The authors' confirmatory runs used --require-pushed-prereg with --prereg <commit of the pushed pre-registration>;
both are optional for anyone else (the --prereg value is only recorded). Examples (from the repository root):
  python check.py inventory                                  # downloads metadata.json -> out/pairs.json (+ .tsv)
  python check.py mesh --segments 20260220214732 --out out/mesh_local.jsonl   # forward mesh route, one pair
  python check.py reverse-lc --samples PHerc0139 --out out/lc_0139.jsonl
  python check.py selftest --out out/selftest_local.json     # synthetic tests (pytest)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from meshpair import __version__                                    # noqa: E402
from meshpair import conventions as cv                             # noqa: E402
from meshpair.catalog import coarse_to_fine_affine, load_catalog, pairs  # noqa: E402
from meshpair.meshcheck import check_pair                          # noqa: E402
from meshpair.tifxyz import fetch_tifxyz, read_tifxyz              # noqa: E402

DEV_SEGMENTS = {("PHerc0841", "20260220213127")}     # w00: development / plumbing pair (seen before the prereg)
CATALOG_URL = "https://vesuvius-challenge-open-data.s3.amazonaws.com/metadata.json"


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


REQUIRE_PUSHED_PREREG = False   # set by --require-pushed-prereg (the authors' confirmatory-run guard)


def _guard(sel, prereg):
    """The authors ran every confirmatory measurement with this guard on: a --prereg commit had to be on a pushed
    branch of the repository holding the pre-registration. For everyone else the guard is off by default, and the
    --prereg value (if any) is only written into each record."""
    conf = [p for p in sel if (p["sample"], p["segment"]) not in DEV_SEGMENTS]
    if not conf or not REQUIRE_PUSHED_PREREG:
        return
    if not prereg:
        sys.exit(f"refusing: {len(conf)} confirmatory pair(s) selected and no --prereg <commit> given "
                 "(rule: no decisive measurement before the prereg is pushed)")
    ok = subprocess.run(["git", "-C", str(HERE), "branch", "-r", "--contains", prereg], capture_output=True, text=True)
    if ok.returncode != 0 or not ok.stdout.strip():
        sys.exit(f"refusing: prereg commit {prereg} is not on any pushed (remote) branch")


def _select(P, args):
    sel = P
    if args.samples:
        s = set(args.samples.split(","))
        sel = [p for p in sel if p["sample"] in s]
    if args.segments:
        spec = Path(args.segments).read_text() if Path(args.segments).is_file() else args.segments
        s = {x.strip() for x in spec.replace("\n", ",").split(",") if x.strip()}
        sel = [p for p in sel if p["segment"] in s or (p["long_id"] or "") in s]
    if getattr(args, "exclude_dev", False):
        sel = [p for p in sel if (p["sample"], p["segment"]) not in DEV_SEGMENTS]
    if getattr(args, "sample_n", None):
        rng = np.random.default_rng(args.seed)
        if len(sel) > args.sample_n:
            idx = sorted(rng.choice(len(sel), args.sample_n, replace=False))
            sel = [sel[i] for i in idx]
    return sel


def cmd_inventory(args):
    cat_path = args.catalog
    if not Path(cat_path).exists():
        import urllib.request
        Path(cat_path).parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(CATALOG_URL, cat_path)
    cat = load_catalog(cat_path)
    P = pairs(cat)
    out = dict(catalog=str(cat_path), catalog_sha256=_sha(cat_path), made=time.strftime("%Y-%m-%dT%H:%M:%S"),
               tool_version=__version__, n_pairs=len(P), pairs=P)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=1))
    by = {}
    for p in P:
        k = (p["sample"], p["fine"]["px"], p["coarse"]["px"])
        by[k] = by.get(k, 0) + 1
    for k, v in sorted(by.items()):
        print(v, *k)
    tsv = Path(args.out).with_suffix(".tsv")
    cols = ["sample", "segment", "long_id", "fine_px", "coarse_px", "fine_is_copy", "coarse_is_copy", "fine_key",
            "coarse_key", "source_key", "fine_sv_key", "coarse_sv_key", "exporter_commit_fine", "exporter_commit_coarse"]
    with open(tsv, "w", newline="\n") as fo:
        fo.write("\t".join(cols) + "\n")
        for p in P:
            fo.write("\t".join(str(x) for x in [
                p["sample"], p["segment"], p["long_id"], p["fine"]["px"], p["coarse"]["px"], p["fine"]["is_copy"],
                p["coarse"]["is_copy"], p["fine"]["key"], p["coarse"]["key"], p["source_key"],
                (p["fine_sv"] or {}).get("key"), (p["coarse_sv"] or {}).get("key"), p["fine"]["commit"],
                p["coarse"]["commit"]]) + "\n")
    print(len(P), "pairs ->", args.out, "and", tsv)


def _load_pairs(path):
    return json.loads(Path(path).read_text())["pairs"]


def cmd_fetch(args):
    P = _select(_load_pairs(args.pairs), args)
    man = {}
    for p in P:
        for key in {p["fine"]["key"], p["coarse"]["key"], p["source_key"]} - {None}:
            local, shas = fetch_tifxyz(key, args.cache)
            man[key] = shas
            print("ok", key)
    mp = Path(args.cache) / "manifest.json"
    old = json.loads(mp.read_text()) if mp.exists() else {}
    old.update(man)
    mp.write_text(json.dumps(old, indent=1))


def cmd_mesh(args):
    P = _select(_load_pairs(args.pairs), args)
    _guard(P, args.prereg)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    done = set()
    if Path(args.out).exists():
        for ln in open(Path(args.out)):
            r = json.loads(ln); done.add((r["sample"], r["segment"]))
    print(len(P), "pairs selected;", len(done), "already in", args.out, flush=True)
    with open(args.out, "a") as fo:
        for p in P:
            if (p["sample"], p["segment"]) in done:
                continue
            t0 = time.time()
            f, c = p["fine"], p["coarse"]
            for key in {f["key"], c["key"], p["source_key"]} - {None}:
                fetch_tifxyz(key, args.cache)
            F = read_tifxyz(Path(args.cache) / f["key"])
            C = read_tifxyz(Path(args.cache) / c["key"])
            M = coarse_to_fine_affine(p)
            rec = dict(sample=p["sample"], segment=p["segment"], long_id=p["long_id"], fine_key=f["key"],
                       coarse_key=c["key"], dev_pair=(p["sample"], p["segment"]) in DEV_SEGMENTS,
                       exporter_commits=[f["commit"], c["commit"]], prereg=args.prereg, tool_version=__version__)
            if M is None or not F.valid.any() or not C.valid.any():
                rec["error"] = "no affine in provenance" if M is None else "a mesh has no valid vertex"
                fo.write(json.dumps(rec) + "\n"); continue
            S = None
            if p["source_key"] and (Path(args.cache) / p["source_key"]).exists():
                S = read_tifxyz(Path(args.cache) / p["source_key"], vc_load=True)
            else:
                rec["warning"] = "no source mesh: M2a (replica) not run for this pair"
            rec.update(check_pair(F, C, f["px"], c["px"], M, source=S,
                                  A_src_to_fine=None if f["is_copy"] else f["affine_src_to_this"],
                                  A_src_to_coarse=None if c["is_copy"] else c["affine_src_to_this"],
                                  fine_is_copy=f["is_copy"], coarse_is_copy=c["is_copy"], do_free=not args.no_free))
            rec["seconds"] = round(time.time() - t0, 1)
            fo.write(json.dumps(rec, default=float) + "\n"); fo.flush()
            print(p["sample"], p["long_id"], rec.get("lines"), rec["seconds"], "s")


def cmd_reverse(args):
    """TIFXYZ-REVERSE (prereg 2026-10-13): reverse route per pair -> JSONL. Fetches missing fine/coarse meshes (public S3)."""
    from meshpair.reverse import reverse_route
    P = _select(_load_pairs(args.pairs), args)
    if args.forward_verdict:
        fwd = {}
        for ln in open(args.forward_records):
            r = json.loads(ln)
            fwd[(r["sample"], r["segment"])] = (r.get("lines") or {}).get("M1_in_files")
        P = [p for p in P if (fwd.get((p["sample"], p["segment"])) == "IN_FILES") == (args.forward_verdict == "readable")]
    _guard(P, args.prereg)
    done = set()
    if Path(args.out).exists():
        for ln in open(Path(args.out)):
            r = json.loads(ln); done.add((r["sample"], r["segment"]))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    print(len(P), "pairs selected;", len(done), "already in", args.out, flush=True)
    with open(args.out, "a") as fo:
        for p in P:
            if (p["sample"], p["segment"]) in done:
                continue
            t0 = time.time()
            f, c = p["fine"], p["coarse"]
            for key in (f["key"], c["key"]):
                fetch_tifxyz(key, args.cache)
            M = coarse_to_fine_affine(p)
            rec = dict(sample=p["sample"], segment=p["segment"], long_id=p["long_id"], fine_key=f["key"], coarse_key=c["key"],
                       fine_is_copy=f["is_copy"], coarse_is_copy=c["is_copy"], prereg=args.prereg, tool_version=__version__)
            if M is None:
                rec["error"] = "no affine in provenance"
            else:
                F = read_tifxyz(Path(args.cache) / f["key"])
                C = read_tifxyz(Path(args.cache) / c["key"])
                rec.update(fine_shape=list(F.points.shape[:2]), coarse_shape=list(C.points.shape[:2]))
                if not (F.valid.any() and C.valid.any()):
                    rec["error"] = "a mesh has no valid vertex"
                    F = C = None
            if M is not None and "error" not in rec:
                rec["reverse"] = reverse_route(F.points, F.valid, C.points, C.valid, F.geom(f["px"]), C.geom(c["px"]), M,
                                               f["px"], c["px"], n_max=100000, seed=0)
                del F, C
            rec["seconds"] = round(time.time() - t0, 1)
            fo.write(json.dumps(rec, default=float) + "\n"); fo.flush()
            rv = rec.get("reverse", {})
            print(p["sample"], p["segment"], "pass" if rv.get("pair_pass") else "FAIL", "frac_on=%.4f" % rv.get("frac_on", float("nan")),
                  "med_n=%.4g" % rv.get("median_normal_coarse_vx", float("nan")), "inpl_p95=%.4g" % rv.get("inplane_p95_px", float("nan")),
                  rec["seconds"], "s", flush=True)


def cmd_reverse_lc(args):
    """TIFXYZ-REVERSE amendment 1 (prereg 2026-10-13_partner_tifxyz_reverse_amend1): lattice-clean route per pair -> JSONL.
    Fetches the fine, coarse and source meshes (the source only for its grid shape)."""
    import tifffile
    from meshpair.reverse_lc import run_pair
    P = _select(_load_pairs(args.pairs), args)
    if args.forward_verdict:
        fwd = {}
        for ln in open(args.forward_records):
            r = json.loads(ln)
            fwd[(r["sample"], r["segment"])] = (r.get("lines") or {}).get("M1_in_files")
        P = [p for p in P if (fwd.get((p["sample"], p["segment"])) == "IN_FILES") == (args.forward_verdict == "readable")]
    _guard(P, args.prereg)
    done = set()
    if Path(args.out).exists():
        for ln in open(Path(args.out)):
            r = json.loads(ln); done.add((r["sample"], r["segment"]))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    print(len(P), "pairs selected;", len(done), "already in", args.out, flush=True)
    with open(args.out, "a") as fo:
        for p in P:
            if (p["sample"], p["segment"]) in done:
                continue
            t0 = time.time()
            f, c = p["fine"], p["coarse"]
            shas = {}
            for key in (f["key"], c["key"], p["source_key"]):
                shas[key] = fetch_tifxyz(key, args.cache)[1]
            M = coarse_to_fine_affine(p)
            rec = dict(sample=p["sample"], segment=p["segment"], long_id=p["long_id"], fine_key=f["key"], coarse_key=c["key"],
                       source_key=p["source_key"], fine_is_copy=f["is_copy"], coarse_is_copy=c["is_copy"], prereg=args.prereg,
                       tool_version=__version__, sha256=shas)
            if M is None:
                rec["error"] = "no affine in provenance"
            else:
                with tifffile.TiffFile(Path(args.cache) / p["source_key"] / "x.tif") as tf:
                    n_src = tuple(int(x) for x in tf.pages[0].shape[:2])
                F = read_tifxyz(Path(args.cache) / f["key"])
                C = read_tifxyz(Path(args.cache) / c["key"])
                rec.update(fine_shape=list(F.points.shape[:2]), coarse_shape=list(C.points.shape[:2]), source_shape=list(n_src))
                if f["is_copy"] and tuple(F.points.shape[:2]) != n_src:
                    rec["warning"] = "fine is a copy but its grid differs from the source grid"
                if c["is_copy"] and tuple(C.points.shape[:2]) != n_src:
                    rec["warning"] = "coarse is a copy but its grid differs from the source grid"
                if not (F.valid.any() and C.valid.any()):
                    rec["error"] = "a mesh has no valid vertex"
            if M is not None and "error" not in rec:
                rec["lattice_clean"] = run_pair(F.points, F.valid, C.points, C.valid, F.geom(f["px"]), C.geom(c["px"]), M, n_src,
                                                f["px"], c["px"], coarse_is_copy=c["is_copy"], fine_is_copy=f["is_copy"],
                                                n_max=100000, seed=0)
                del F, C
            rec["seconds"] = round(time.time() - t0, 1)
            fo.write(json.dumps(rec, default=float) + "\n"); fo.flush()
            rv = rec.get("lattice_clean", {})
            print(p["sample"], p["segment"], rv.get("verdict", rec.get("error")), rv.get("direction"),
                  "frac_clean=%.3f" % rv.get("frac_clean", float("nan")), "frac_on_clean=%.4f" % rv.get("frac_on_clean", float("nan")),
                  "inpl_p95=%.4g bar=%.4g" % (rv.get("inplane_p95_clean", float("nan")), rv.get("inplane_bar_px", float("nan"))),
                  rec["seconds"], "s", flush=True)


def cmd_ct(args):
    from meshpair.ct import CanvasMap, TileSpec, fit_field, measure_tile, remote_sv
    P = _select(_load_pairs(args.pairs), args)
    if len(P) != 1:
        sys.exit("ct: select exactly one pair")
    _guard(P, args.prereg)
    p = P[0]
    if not (p["fine_sv"] and p["coarse_sv"]):
        sys.exit("ct: pair has no surface volume for both meshes")
    F = read_tifxyz(Path(args.cache) / p["fine"]["key"])
    C = read_tifxyz(Path(args.cache) / p["coarse"]["key"])
    gF, gC = F.geom(p["fine"]["px"]), C.geom(p["coarse"]["px"])
    k = p["fine"]["px"] / p["coarse"]["px"]
    ay, by = cv.exact_linear_coeffs(gF, gC, 0)
    ax, bx = cv.exact_linear_coeffs(gF, gC, 1)
    maps = {"exact": CanvasMap(ay, by, ax, bx), "naive": CanvasMap(0.0, k, 0.0, k)}
    base = "https://vesuvius-challenge-open-data.s3.amazonaws.com/"
    fsv = remote_sv(base + p["fine_sv"]["key"], p["fine"]["px"], level=args.fine_level,
                    n_layers0=p["fine_sv"]["num_slices"])
    csv = remote_sv(base + p["coarse_sv"]["key"], p["coarse"]["px"], level=0, n_layers0=p["coarse_sv"]["num_slices"])
    spec = TileSpec(T=args.tile, S=args.search, depth_range_um=args.depth_range)
    # lattice of tile centres (coarse canvas) where the coarse mesh is valid around the whole tile
    vC = C.valid
    step = args.lattice
    tiles = []
    for yc in np.arange(spec.T, gC.canvas_h - spec.T, step):
        for xc in np.arange(spec.T, gC.canvas_w - spec.T, step):
            g0 = gC.canvas_to_grid(np.array([yc - spec.T, yc + spec.T]), 0)
            g1 = gC.canvas_to_grid(np.array([xc - spec.T, xc + spec.T]), 1)
            blk = vC[max(0, int(g0[0])):int(np.ceil(g0[1])) + 1, max(0, int(g1[0])):int(np.ceil(g1[1])) + 1]
            if blk.size and blk.mean() >= 0.98:
                tiles.append((float(yc), float(xc)))
    if args.max_tiles and len(tiles) > args.max_tiles:
        rng = np.random.default_rng(args.seed)
        tiles = [tiles[i] for i in sorted(rng.choice(len(tiles), args.max_tiles, replace=False))]
    if args.exclude_fine_box:
        y0f, y1f, x0f, x1f = map(float, args.exclude_fine_box.split(","))
        ey = cv.exact_canvas_map(np.array([y0f, y1f]), gF, gC, 0)
        ex = cv.exact_canvas_map(np.array([x0f, x1f]), gF, gC, 1)
        h = spec.T / 2
        tiles = [(yc, xc) for (yc, xc) in tiles
                 if not (yc + h > ey[0] and yc - h < ey[1] and xc + h > ex[0] and xc - h < ex[1])]
    if args.dry_run:
        print(json.dumps(dict(segment=p["segment"], n_tiles=len(tiles), fine_sv_shape=fsv.shape,
                              coarse_sv_shape=csv.shape, fine_layers0=fsv.n_layers0, coarse_layers=csv.n_layers0,
                              exact=maps["exact"].__dict__, first_tiles=tiles[:3]), default=float))
        return 0
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    seg = p["segment"]
    n_c = p["coarse_sv"]["num_slices"]
    recs = []
    with open(args.out, "a") as fo:
        # pass 1: unshifted maps on every lattice tile
        for (yc, xc) in tiles:
            for name in args.maps.split(","):
                r = measure_tile(fsv, csv, maps[name], (yc, xc), spec, n_c=n_c)
                r.update(sample=p["sample"], segment=seg, map=name, prereg=args.prereg, kind="tile")
                fo.write(json.dumps(r, default=float) + "\n"); fo.flush()
                recs.append(r)
                print(name, yc, xc, np.round(r["map_correction_px"], 2), round(r["depth_offset_um"], 1), r["valid"])
        # reachable tiles fixed on unshifted data (exact map) BEFORE any plant
        reach = [r["centre_c"] for r in recs if r["map"] == "exact" and r["valid"]]
        rp_path = Path(args.out).parent / f"reachable_{seg}.json"
        rp_path.write_text(json.dumps(dict(segment=seg, tiles=reach, spec=spec.__dict__), default=float))
        print("reachable", len(reach), "of", len(tiles), "sha256", _sha(rp_path))
        # pass 2: differential plants on reachable tiles only
        if args.plants:
            base = {tuple(r["centre_c"]): r for r in recs if r["map"] == "exact"}
            for c in reach:
                b0 = base[tuple(c)]
                plants = {}
                for pl in ((4, 0), (0, 8), (16, 0), (0, -16)):
                    rp = measure_tile(fsv, csv, maps["exact"].planted(*pl), tuple(c), spec, splits=False, n_c=n_c)
                    plants[f"{pl[0]},{pl[1]}"] = dict(
                        delta=[rp["map_correction_px"][0] - b0["map_correction_px"][0],
                               rp["map_correction_px"][1] - b0["map_correction_px"][1]], valid=rp["valid"])
                fo.write(json.dumps(dict(kind="plants", segment=seg, centre_c=c, plants=plants), default=float) + "\n")
                fo.flush()
    ex = [r for r in recs if r["map"] == "exact" and r["valid"]]
    if ex:
        print("field:", fit_field(ex))


A1_PLANTS = ((4, 0), (0, 8), (12, 0), (0, -12))
A1_PREREG = "13df5696"
A2_PREREG = "fcb2012f"  # amendment 2 (prereg/2026-10-13_partner_w00field_amend2.md): same plants command, fresh C-set


def _a1_inputs(out_dir, seg, prereg):
    """Input checks for `plants` (no network): prereg, frozen reachable list vs sha in the pass-1 log, base records.
    Returns (reachable_path, reach_centres, spec_dict, base_exact_records_by_centre)."""
    if prereg not in (A1_PREREG, A2_PREREG):
        sys.exit(f"refusing: plants (Amendment 1/2) requires --prereg {A1_PREREG} or {A2_PREREG}")
    out_dir = Path(out_dir)
    rp = out_dir / f"reachable_{seg}.json"
    logs = [out_dir / f"ct_{seg}.log"] + ([out_dir / "ct_w00.log"] if seg == "20260220213127" else [])
    log = next((l for l in logs if l.exists()), None)
    if log is None or not rp.exists():
        sys.exit(f"refusing: missing {rp.name} or pass-1 log for {seg}")
    want = [ln.split("sha256")[1].strip() for ln in log.read_text(errors="replace").splitlines()
            if ln.startswith("reachable") and "sha256" in ln]
    if not want or _sha(rp) != want[-1]:
        sys.exit(f"refusing: sha256 of {rp.name} ({_sha(rp)}) does not match the pass-1 log ({want[-1:] })")
    rd = json.loads(rp.read_text())
    jl = out_dir / f"ct_{seg}.jsonl"
    if not jl.exists():
        sys.exit(f"refusing: missing {jl.name}")
    base = {}
    for ln in jl.read_text().splitlines():
        r = json.loads(ln)
        if r.get("kind") == "tile" and r.get("map") == "exact":
            base[tuple(r["centre_c"])] = r
    reach = [tuple(c) for c in rd["tiles"]]
    miss = [c for c in reach if c not in base or not base[c]["valid"]]
    if miss:
        sys.exit(f"refusing: {len(miss)} reachable tile(s) without a valid base exact record")
    return rp, reach, rd["spec"], base


def cmd_plants(args):
    """Amendment 1 (prereg/2026-10-11_partner_w00field_amend1.md): re-measure ONLY the new C1 plants on the frozen
    reachable tiles. Base (unshifted, exact-map) records are read from ct_<seg>.jsonl, which is never modified."""
    from meshpair.ct import CanvasMap, TileSpec, measure_tile, remote_sv
    P = _select(_load_pairs(args.pairs), args)
    if len(P) != 1:
        sys.exit("plants: select exactly one pair")
    p = P[0]
    seg = p["segment"]
    out_dir = Path(args.out_dir)
    _, reach, spec_d, base = _a1_inputs(out_dir, seg, args.prereg)
    _guard(P, args.prereg)
    spec = TileSpec(T=args.tile, S=args.search, depth_range_um=args.depth_range)
    if json.loads(json.dumps(spec.__dict__, default=float)) != spec_d:
        sys.exit("refusing: TileSpec differs from the one stored with the frozen reachable list")
    if not (p["fine_sv"] and p["coarse_sv"]):
        sys.exit("plants: pair has no surface volume for both meshes")
    F = read_tifxyz(Path(args.cache) / p["fine"]["key"])
    C = read_tifxyz(Path(args.cache) / p["coarse"]["key"])
    gF, gC = F.geom(p["fine"]["px"]), C.geom(p["coarse"]["px"])
    ay, by = cv.exact_linear_coeffs(gF, gC, 0)
    ax, bx = cv.exact_linear_coeffs(gF, gC, 1)
    exact = CanvasMap(ay, by, ax, bx)
    url = "https://vesuvius-challenge-open-data.s3.amazonaws.com/"
    fsv = remote_sv(url + p["fine_sv"]["key"], p["fine"]["px"], level=args.fine_level, n_layers0=p["fine_sv"]["num_slices"])
    csv = remote_sv(url + p["coarse_sv"]["key"], p["coarse"]["px"], level=0, n_layers0=p["coarse_sv"]["num_slices"])
    n_c = p["coarse_sv"]["num_slices"]
    outp = out_dir / f"ct_plants_a1_{seg}.jsonl"
    done = set()
    if outp.exists():
        done = {tuple(json.loads(l)["centre_c"]) for l in outp.read_text().splitlines() if l.strip()}
    with open(outp, "a") as fo:
        for c in reach:
            if c in done:
                continue
            b0 = base[c]
            plants = {}
            for pl in A1_PLANTS:
                rp_ = measure_tile(fsv, csv, exact.planted(*pl), c, spec, splits=False, n_c=n_c)
                plants[f"{pl[0]},{pl[1]}"] = dict(
                    delta=[rp_["map_correction_px"][0] - b0["map_correction_px"][0],
                           rp_["map_correction_px"][1] - b0["map_correction_px"][1]], valid=rp_["valid"])
            fo.write(json.dumps(dict(kind="plants_a1", segment=seg, centre_c=list(c), plants=plants,
                                     prereg=args.prereg, tool_version=__version__), default=float) + "\n")
            fo.flush()
            print(seg, c, {k: np.round(v["delta"], 2).tolist() for k, v in plants.items()}, flush=True)


def cmd_pyramid(args):
    """Stage 0b: level-l of the fine surface volume == 2^l block mean of level 0 (within 1 grey level), all layers."""
    from meshpair.ct import remote_sv
    from meshpair.synth import block_mean
    P = _select(_load_pairs(args.pairs), args)
    base = "https://vesuvius-challenge-open-data.s3.amazonaws.com/"
    out = []
    for p in P:
        if not p["fine_sv"]:
            continue
        url = base + p["fine_sv"]["key"]
        l = args.fine_level
        f = 2 ** l
        s0, sl = remote_sv(url, p["fine"]["px"], 0), remote_sv(url, p["fine"]["px"], l)
        H, W = sl.shape[1], sl.shape[2]
        z = sl.shape[0] // 2
        y, x = H // 2 // 8 * 8, W // 2 // 8 * 8          # centre block (I/O check only, no measurement)
        a = sl.read(z - 3, z + 3, y, y + 40, x, x + 40)
        b = block_mean(s0.read(z - 3, z + 3, y * f, (y + 40) * f, x * f, (x + 40) * f), f)
        nz = (b != 0).mean()
        rec = dict(segment=p["segment"], key=p["fine_sv"]["key"], level=l, layers_level=sl.shape[0],
                   layers0=s0.shape[0], max_abs_diff=float(np.abs(a - b).max()), nonzero_frac=float(nz),
                   ok=bool(sl.shape[0] == s0.shape[0] and np.abs(a - b).max() <= 1.0 and nz > 0.5))
        out.append(rec); print(json.dumps(rec))
    if args.out:
        Path(args.out).write_text(json.dumps(out, indent=1))


def cmd_selftest(args):
    import pytest
    out = dict(tool_version=__version__, when=time.strftime("%Y-%m-%dT%H:%M:%S"))
    rc = pytest.main(["-q", str(HERE / "tests")])
    out["synthetic_pytest_exit"] = int(rc)
    if args.plumbing_mesh:
        Pdir, Ndir = args.plumbing_mesh
        A = np.array(json.loads(args.affine)) if args.affine else np.array(W00_AFFINE)
        Fm, Cm = read_tifxyz(Pdir), read_tifxyz(Ndir)
        rec = check_pair(Fm, Cm, 2.403, 9.366, A, source=read_tifxyz(Ndir, vc_load=True), A_src_to_fine=A,
                         coarse_is_copy=True)
        out["plumbing_mesh_w00"] = rec
        print(json.dumps(rec["lines"]), "route_I offset vs naive at centre (px):",
              np.round(rec["route_I"]["offset_vs_naive_centre_px"], 3))
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out, indent=1, default=float))
    print("selftest ->", args.out, "pytest exit", rc)
    return rc


# catalog affine PHerc0841 20250821151531 (9.366 um) -> 20260319124803 (2.403 um), as published in metadata.json
W00_AFFINE = [[3.8989743209135876, -0.01255505722529341, -0.0018339025635139032, 343.9127006249706],
              [0.010672926356536081, 3.8982252200029714, 0.0022465882056711857, 214.7521699852929],
              [7.354390744793846e-05, -0.001896739717178478, 3.897597266482002, -49844.690396924205]]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    a = sp.add_parser("inventory"); a.add_argument("--catalog", default=str(HERE / "cache" / "metadata.json"))
    a.add_argument("--out", default=str(HERE / "out" / "pairs.json")); a.set_defaults(fn=cmd_inventory)

    def common(x):
        x.add_argument("--pairs", default=str(HERE / "out" / "pairs.json"))
        x.add_argument("--cache", default=str(HERE / "cache" / "meshes"))
        x.add_argument("--samples"); x.add_argument("--segments")
        x.add_argument("--exclude-dev", action="store_true")
        x.add_argument("--sample-n", type=int); x.add_argument("--seed", type=int, default=20261011)
        x.add_argument("--prereg", help="pre-registration commit; written into each record")
        x.add_argument("--require-pushed-prereg", action="store_true",
                       help="refuse non-development pairs unless --prereg is on a pushed branch of this git repository")
    b = sp.add_parser("fetch"); common(b); b.set_defaults(fn=cmd_fetch)
    c = sp.add_parser("mesh"); common(c); c.add_argument("--out", default=str(HERE / "out" / "local" / "mesh.jsonl"))
    c.add_argument("--no-free", action="store_true"); c.set_defaults(fn=cmd_mesh)
    r = sp.add_parser("reverse"); common(r); r.add_argument("--out", default=str(HERE / "out" / "local" / "reverse.jsonl"))
    r.add_argument("--forward-verdict", choices=["readable", "void"], help="filter by the stored forward M1 verdict")
    r.add_argument("--forward-records", default=str(HERE / "out" / "mesh_records.jsonl")); r.set_defaults(fn=cmd_reverse)
    q = sp.add_parser("reverse-lc"); common(q); q.add_argument("--out", default=str(HERE / "out" / "local" / "reverse_lc.jsonl"))
    q.add_argument("--forward-verdict", choices=["readable", "void"], help="filter by the stored forward M1 verdict")
    q.add_argument("--forward-records", default=str(HERE / "out" / "mesh_records.jsonl")); q.set_defaults(fn=cmd_reverse_lc)
    d = sp.add_parser("ct"); common(d); d.add_argument("--out", default=str(HERE / "out" / "local" / "ct_tiles.jsonl"))
    d.add_argument("--fine-level", type=int, default=2); d.add_argument("--tile", type=int, default=213)
    d.add_argument("--search", type=int, default=16); d.add_argument("--depth-range", type=float, default=28.0)
    d.add_argument("--lattice", type=float, default=534.0, help="tile spacing, coarse px (534 ~ 5 mm at 9.4 um)")
    d.add_argument("--max-tiles", type=int); d.add_argument("--maps", default="exact,naive")
    d.add_argument("--plants", action="store_true")
    d.add_argument("--dry-run", action="store_true", help="list tiles + SV shapes only (no measurement)")
    d.add_argument("--exclude-fine-box", help="y0,y1,x0,x1 in FINE canvas px; tiles overlapping it are skipped")
    d.set_defaults(fn=cmd_ct)
    h = sp.add_parser("plants"); common(h); h.add_argument("--out-dir", default=str(HERE / "out"))
    h.add_argument("--fine-level", type=int, default=2); h.add_argument("--tile", type=int, default=213)
    h.add_argument("--search", type=int, default=16); h.add_argument("--depth-range", type=float, default=28.0)
    h.set_defaults(fn=cmd_plants)
    g = sp.add_parser("pyramid"); common(g); g.add_argument("--fine-level", type=int, default=2)
    g.add_argument("--out"); g.set_defaults(fn=cmd_pyramid)
    e = sp.add_parser("selftest"); e.add_argument("--plumbing-mesh", nargs=2, metavar=("FINE_DIR", "COARSE_DIR"))
    e.add_argument("--affine"); e.add_argument("--out", default=str(HERE / "out" / "local" / "selftest.json"))
    e.set_defaults(fn=cmd_selftest)
    args = ap.parse_args()
    global REQUIRE_PUSHED_PREREG
    REQUIRE_PUSHED_PREREG = bool(getattr(args, "require_pushed_prereg", False))
    rc = args.fn(args)
    sys.exit(int(rc or 0))


if __name__ == "__main__":
    main()
