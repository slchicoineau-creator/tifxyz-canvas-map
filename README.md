# tifxyz-canvas-map

Map pixels between the **fine (~2.2-2.4 µm)** and **coarse (~7.9-9.4 µm)** surface volumes of the same Vesuvius
Challenge segment, and check every published mesh pair against the rule that relates them.

> **Author's note.** We started looking into this when we tried to move ink labels and predictions between the 2.4 µm and 9.4 µm surface volumes of the same segment and noticed that multiplying pixel coordinates by the voxel-size ratio landed about 85 µm away from the same papyrus on PHerc0841 w00, even though the 3-D meshes agreed. Working with my teammate Chris Thomas, Claude Code traced the offset to the half-pixel resampling and grid-sizing convention in `volume-cartographer`'s mesh exporter and renderer. We trust this rule because our Python replica rebuilds all 206 published transformed meshes to within 0.008 voxel, passes 52 of 52 held-out mesh pairs (plus 151 of 154 applied pairs, 0 failures), and cuts the CT surface-volume mismatch on 7 fresh pre-registered segments from 6.3–14.4 coarse pixels down to 0.48–1.46 coarse pixels. Note that this tool only maps 2-D canvas pixels between the two surface volumes—it does not fix the underlying 3-D scan-to-scan affine transform, which still leaves that 0.5–1.5 coarse-pixel residual.

## Why this matters

Cross-frame training and evaluation depend on moving ink labels, model predictions, and tile coordinates accurately between the fine (2.2–2.4 µm) and coarse (7.9–9.4 µm) surface volumes of the same segment. Because the fine and coarse canvases do not share a 2-D pixel origin, and their scales do not follow the nominal voxel-size ratio, naive index scaling puts coordinates in the wrong place by roughly 0.1 mm on the 9.4 µm and 8.64 µm scrolls, and by 0.3–2.3 mm on the 7.91 µm scrolls (where index pairing is unusable). This repository gives a one-call Python function and a lookup table for all 207 published pairs so teams can map 2-D pixels to the exact same mesh point, separating the 2-D canvas convention from the remaining 3-D scan-to-scan transform residual.

## What it is, in one paragraph

The organizers publish, per segment, a mesh on the fine scan and a mesh on the coarse scan, and a surface volume
rendered along each. The two meshes come from one source mesh and agree in 3-D. Their **2-D pixel grids do not line up
by the voxel-size ratio**: pairing fine pixel `p` with coarse pixel `p × (fine µm / coarse µm)` lands tens of µm to
millimetres away from the same papyrus (median 56-132 µm at the canvas centre on the 9.4 µm scrolls, 89-112 µm on the
8.64 µm scrolls, 0.3-2.3 mm on the 7.91 µm scrolls). The cause is a convention in the exporter and the renderer, not
an error in the transform: per axis the same mesh point sits at

```
p_coarse = (p_fine + 10.5) × n_coarse / n_fine − 10.5        (grid step s = 0.05, render --scale 1)
```

where `n` are the two meshes' grid sizes on that axis. This repository gives that map as a function, a table of it for
all 207 published pairs, and the checker that tested it against the published files.

## Quick start

```bash
pip install -e .            # or: pip install numpy scipy tifffile imagecodecs numcodecs
```

```python
from meshpair.canvasmap import fine_to_coarse_px, coarse_to_fine_px, warp_canvas_image

fine   = "PHerc0139/segments/<seg>/mesh/<id>-on-<fine volume>.tifxyz"     # local copies of the two mesh folders
coarse = "PHerc0139/segments/<seg>/mesh/<id>-on-<coarse volume>.tifxyz"   # the surface volumes were rendered from

yc, xc = fine_to_coarse_px(1200.0, 3400.0, fine, coarse)                 # level-0 pixels
yc2, xc2 = fine_to_coarse_px(300.0, 850.0, fine, coarse, fine_level=2)   # fine pixel given at pyramid level 2
labels_on_coarse = warp_canvas_image(fine_label_mask, fine, coarse, "fine_to_coarse", order=0)  # 0 = nearest for masks
```

- From the command line: `python -m meshpair.canvasmap FINE_DIR COARSE_DIR --fine-px 1200 3400`.
- Without the meshes: `results/pair_maps.tsv` lists `p_coarse = a + b × p_fine` (`y_a, y_b, x_a, x_b`, level 0) for
  every pair in the catalog of 2026-10-11 (`metadata.json` sha256 `7ec8dcaa…abdd3`). A later re-export of a mesh changes
  its row; the helper run on the mesh folders is always current.
- Pyramid levels: level `l` pixel `i` is centred on level-0 pixel `2^l · i + (2^l − 1)/2`.
- What the map does **not** fix: the residual of the scan-to-scan affine itself. After the map, the two surface
  volumes still differ by about 0.5-1.5 coarse px (≈ 4-14 µm) on the segments we measured (CT route below).

## Where the convention comes from (villa `volume-cartographer`)

1. `vc_transform_geom` → `transformSurfacePoints` (`core/src/AffineTransform.cpp`) maps the source mesh with the
   affine and regrids it with `QuadSurface::resample(INTER_LINEAR)`, with the new grid size
   `round(n · measured spacing ratio)` (it follows the affine's measured scale, not the nominal voxel-size ratio).
2. `resamplePointsLinearPreservingInvalids` (`core/src/QuadSurface.cpp`) samples half-pixel centres:
   `src = (dst + 0.5) · n_src/n_dst − 0.5`.
3. `vc_render_tifxyz` (`computeCanvasOrigin` + `QuadSurface::gen`) maps canvas pixel `p` to grid `s · (p + 0.5)`.

Two meshes made from one source this way satisfy the formula above. The constant part is −7.3 to −7.8 coarse px at the
canvas origin on every scroll; the slope part (grid size ratio vs voxel ratio, 0.03-1.23 % by scroll) dominates on
large canvases and on the 7.91 µm scrolls. The exporter commits on the published meshes (d9a472ea, 650bd566, 77627ab1)
and the render commits (530584554c, 650bd566, 862bd172, ae2ecce9, bc28564d) have the same code for these steps, and
villa main has not changed it as of 2026-10-18. `meshpair/exporter.py` and `meshpair/conventions.py` are clean-room Python
implementations of that behaviour, written from a functional specification (see Licence and `cleanroom/`).

## What was checked (catalog of 2026-10-11, 207 pairs on 8 scrolls)

All measurements are mesh-only on CPU except the CT route. Each was pre-registered (rules and bars committed and pushed
before the run, in the authors' private repository); the full run records are in `out/` and `reverse_route/`.
Scrubbed copies of the pre-registration files (named `prereg/…` in the records) are in `preregs/`, with an index of their
commits and times (`preregs/README.md`).

**1. The exporter replica rebuilds the published meshes.** A Python copy of steps 1-2 rebuilds every published
transformed mesh of the 206 catalog pairs from `mesh/intermediate/tifxyz_normalized` and the catalog affine: same grid
shape, identical valid mask, vertex error p99 ≤ 0.008 voxel (`out/RESULTS.md`, M2a).

**2. Forward mesh route (no convention assumed).** Fine vertices are projected through the catalog affine onto the
coarse mesh surface. Where the coarse mesh is a plain copy of the source (70 pairs: all of PHerc0841, 0009B and 0343P,
19 of 20 PHerc1667, and part of 0814 and Paris4), the measured map matches the formula to p95 ≤ 0.003 coarse px (median 0.0007), and the
offset against index pairing is in the files on 70 / 70. The route is VOID (by its pre-registered on-surface bar) on the
136 pairs where the coarse mesh is itself a resample (`out/RESULTS.md`).

**3. Lattice-clean route.** On those 136 pairs both meshes are resamples of one source mesh whose bilinear surface has
a crease along every source grid line; vertices whose cell straddles a source line sit off the other mesh by
construction (`reverse_route/DIAGNOSIS.md`, a post-hoc diagnosis on 6 pairs). The route therefore measures only
lattice-clean vertices.
- **Held-out test:** 52 / 52 pairs pass (bar ≥ 47), 50 / 52 at the plain 0.003 px in-plane bar; the other 2 are
  explained by the exporter's float32 source index (`reverse_route/s2/RESULT.md`).
- **Application (not a test):** the frozen rule on the other 154 pairs gives 151 PASS, 3 NOT_EVALUABLE (PHerc0500P2,
  fewer than 100 clean vertices), 0 FAIL; 150 / 151 at the plain bar (`reverse_route/apply154/RESULT.md`). The rule
  was designed after seeing those pairs, so these 154 carry no independent evidence weight.
- **Offsets:** every evaluable pair (203) has the offset against index pairing in its files; 192 / 203 exceed 3 coarse
  px at the valid-box centre, 11 are 0.92-2.9 px there (larger towards the edges).
- **Limits (state them with the claim):** the lattice selection was tested on held-out pairs at fine/source ratio 4
  only (42 PHercParis4 pairs); ratio 1.95 (most of PHerc0500P2) has no held-out pair, and the forward direction has one
  (PHerc1667, passing through the float32 floor). The rule is blind to errors confined to straddling cells and to
  defects on fewer than ~5 % of the clean vertices; measured off-surface clean vertices: 9 in 5 of 206 pairs.

**4. CT route (surface volumes, not just meshes).** Tiles of the two surface volumes are correlated after mapping by
the formula or by index pairing (`out/RESULTS.md`, `out/amend2/RESULT.md`).
- The first 10-segment run is **exploratory** (its positive control was redesigned after the values were seen).
- On a fresh, pre-registered set of 7 measurable segments the readings **repeat**: the residual after the formula is
  0.48-1.46 coarse px (median per segment) against 6.3-14.4 px under index pairing; depth stays within ±0.53 coarse
  layers (one PHerc0009B segment +0.92); one affine correction per segment suffices on PHerc0500P2 and PHerc0814
  (coarse-copy tracings), while **a single affine is borderline on PHerc0139** (leave-one-out p95 2.05-2.24 px against a
  2.0 px bar).
- Not covered by the fresh set: PHerc0841 (exploratory only), PHerc0009B (a first reading, not a repeat), PHerc0814
  fine-copy tracings, PHerc0343P, 1667 and Paris4 (no CT run).

**Offset against index pairing at the valid-box centre** (µm, min / median / max over segments; `results/pair_maps.tsv`):

| scroll | fine / coarse µm | pairs | offset (µm) |
|---|---|---|---|
| PHerc0841 | 2.403 / 9.366 | 2 (+ w00) | 102 / 106 / 109 |
| PHerc0009B | 2.401 / 8.64 | 18 | 97 / 112 / 124 |
| PHerc0343P | 2.215 / 8.64 | 8 | 85 / 89 / 100 |
| PHerc0139 | 2.399 / 9.362 | 38 | 17 / 69 / 115 |
| PHerc0500P2 | 2.215 / 9.362 | 46 | 9 / 56 / 133 |
| PHerc0814 | 2.399 / 9.362 | 19 | 39 / 132 / 711 |
| PHerc1667 | 2.399 / 7.91 | 20 | 191 / 304 / 5,954 |
| PHercParis4 | 2.4 / 7.91 | 55 | 667 / 2,342 / 3,060 |

On the 7.91 µm scrolls index pairing is not usable; the slope term dominates on their large canvases.

## Reproduce

```bash
python check.py selftest                                     # synthetic / unit tests, ~1 min, no download
python check.py inventory                                    # metadata.json -> out/pairs.json (207 pairs)
python check.py mesh --segments 20260220214732               # forward route, one pair
python check.py reverse-lc --samples PHerc0139               # lattice-clean route, 38 pairs
python check.py fetch --segments 20250825181922              # meshes for the CT example (ct does not fetch them)
python check.py ct --segments 20250825181922                  # CT route, one segment (streams CT)
python scripts/make_pair_maps.py                             # results/pair_maps.tsv from the stored records
python reverse_route/apply154/score_apply154.py              # the 206-pair table from the stored records
```

- Meshes are fetched from the public bucket into `cache/meshes/` (4 small files each; all 206 pairs ≈ 16 GB with
  PHercParis4, 4.3 GB without). The CT route streams surface-volume chunks into `cache/`.
- New runs write to `out/local/` by default, so the published records are never touched; a command skips pairs
  already in its output file. `--prereg` only labels records. `--require-pushed-prereg` restores the authors' rule (refuse a run unless the
  pre-registration commit is pushed).
- What changed against the code that made the records: `CHANGES.md`. Packaging regression check:
  `results/REGRESSION_2026-10-18.md`.

## Files

| path | what |
|---|---|
| `meshpair/canvasmap.py` | the public helper |
| `meshpair/conventions.py`, `meshpair/exporter.py` | renderer / exporter conventions; exporter replica |
| `meshpair/meshcheck.py`, `fit.py`, `reverse.py`, `reverse_lc.py`, `ct.py` | the routes |
| `check.py` | CLI |
| `out/RESULTS.md`, `out/amend2/RESULT.md`, `reverse_route/**/RESULT.md`, `DESIGN.md` | run reports and design |
| `out/*.jsonl`, `reverse_route/**/*.jsonl` | raw per-pair / per-tile records |
| `results/pair_maps.tsv` | the map for all 207 pairs |

## Related work

- villa `tifxyz_label_transfer` (#1457 by erdpx, merged 2026-08-14; fixed in #1740): moves a label image from one tifxyz
  canvas to another through the two meshes, including fine (2.4 µm) to coarse (9 µm). This repo gives the per-pair
  pixel map as a table and checks it; it does not replace that tool.
- villa #1727 (ge-al) and `tifxyz_canvas.py` (#1699, merged 2026-09-18): the canvas-size formula and the first report
  that fine and coarse canvases differ per axis.
- scrollashlab `sv-canvas-1727` / `xframe` (2026-10-07): explains the per-axis regrid behind #1727's scales and rebuilds
  all 902 published meshes; this repo covers 207 pairs.
- campbelltron1 `vesuvius-frame-mesh` (2026-10-12; see #1845, #1944): picks the published mesh for a volume and checks
  canvas size.
- TAUIL-Abd-Elilah, villa #1912: the PHerc0139 offset and scale measurements credited below.

## Credits

- Code and measurements: Chris Thomas (@Baz248); the formula was first stated publicly in his villa #1912 comment
  (6116916788, 2026-10-12). Release packaging: Sebastien Chicoineau (@slchicoineau-creator). `meshpair/zio.py` is the
  team's reader from its spiral-fit audit package (MIT, Sebastien Chicoineau and Chris Thomas).
- The Vesuvius Challenge organizers: the scans, meshes, surface volumes and catalog (open data, CC BY-NC 4.0).
- `volume-cartographer` (EduceLab and ScrollPrize/villa contributors, GPLv3), whose exporter and renderer behaviour this
  repository follows; it is linked, not copied (see Licence).
- jrudolph (#1912) noted that the 2.4 → 9.4 µm transform itself can be improved; that is a separate 3-D question, and
  this repository does not change the transform.
- TAUIL-Abd-Elilah, villa issue #1912: the PHerc0139 measurements (−8, −7 px offset; 0.3-0.4 % scale) that this
  convention fits, and the separate transform-accuracy numbers it does not replace.
- Related open villa PR #2112 (scrollashlab): another frame field proposed for `meta.json` (`original_volume_downscale`).

## AI use

Prepared with help from Claude Code: it found the exporter/renderer pixel convention in the code, built the Python exporter replica, ran the pre-registered mesh and CT checks, built the release copy (the pixel-map helper, tests and fixes) and drafted the documentation. Chris Thomas (@Baz248) directed and reviewed the checks, and Sebastien Chicoineau (@slchicoineau-creator) reviewed the results and the release.

## Licence

Code: MIT (`LICENSE`). The modules that follow `volume-cartographer`'s pixel conventions (`meshpair/conventions.py`,
`exporter.py`, `tifxyz.py`, `slicecentre.py`) are a clean-room implementation: one person wrote a functional
specification (`cleanroom/SPEC.md`), and a separate implementer who never saw the earlier GPL-derived code or the
C++ wrote them from it; all 61 tests pass, all 207 stored pair maps are identical, and on all 207 pairs the new code gives
bit-identical `mesh` and `reverse-lc` records to the code before the rewrite (`results/PARITY_CLEANROOM.md`). Numbers and images
derived from the Vesuvius Challenge scans: CC BY-NC 4.0 (https://creativecommons.org/licenses/by-nc/4.0/), with
attribution to the Vesuvius Challenge.
