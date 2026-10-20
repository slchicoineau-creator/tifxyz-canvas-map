# W00-FIELD + mesh-pair offset checker: design (phase 1 of 2, 2026-10-11)

**Status:** built and self-tested on CPU. No decisive measurement has been run on the pre-registered set. The
development pair is PHerc0841 w00. Everything already seen on it is listed in section 1 and counts as known before the
prereg.
**Code:** this repository (`check.py` CLI, `meshpair/` package, `tests/`). The original venv was
(uv, py3.12, numpy 2.5.4, scipy 1.18.1, tifffile 2026.9.20, numcodecs 0.17.0, pytest 9.1.1). CPU only.

## 1. What phase 1 already found (development pair w00; disclosed, not confirmatory)

The open caveat was: "First confirm it is in the organizers' files, not our mapping assumption". On w00 the
answer is that the offset **is in the organizers' files**, and it comes from **the exporter's grid convention**. It is
not a scan-registration error, and it is not something our scripts introduced.

- **Provenance (catalog `metadata.json`, sha256 `7ec8dcaa…abdd3`):**
  - The 9.366 µm mesh is a plain copy of `intermediate/tifxyz_normalized` (`copy-normalized-for-original`).
  - The 2.403 µm mesh is `transform-tifxyz-normalized` (villa `d9a472ea`) applied to that same mesh, with the catalog's
    affine `20250821151531 → 20260319124803`.
  - Both surface volumes are `vc_render_tifxyz` (villa `530584554c`, `--scale 1`).
- **Mechanism (villa source; the same code at exporter commits 77627ab1, 650bd566 and d9a472ea, and on main d63d3e968):**
  - `vc_transform_geom` calls `transformSurfacePoints` (`core/src/AffineTransform.cpp`). That function maps every point
    with the affine, measures the median neighbour spacing before and after (`measureGridAxisSpacing`), then calls
    `QuadSurface::resample(ratio_x, ratio_y, INTER_LINEAR)` (`core/src/QuadSurface.cpp`).
  - The resample sets `n_new = round(n_old * ratio)` and samples **half-pixel centres**:
    `src = (dst + 0.5) * n_old / n_new - 0.5` (`resamplePointsLinearPreservingInvalids`).
  - The renderer maps canvas px to grid as `grid = s (p + 0.5)` with s = 0.05. This comes from `computeCanvasOrigin`
    plus `QuadSurface::gen`, and it is the same at all 5 render commits in the catalog.
  - The same mesh vertex therefore sits at **p_C = (p_F + 10.5) · ρ − 10.5** per axis, with ρ = n_C / n_F (211/823 in
    y, 238/928 in x). A script that pairs by index assumes p_C = p_F · k, with k = 2.403/9.366.
  - The difference is a constant ≈ −7.81 coarse px (−73 µm) from the half-pixel rule, plus a drift p_F (ρ − k) from
    the rounded grid size (−1.9e-4 per px in y, −1.0e-4 in x).
- **Exact reproduction:**
  - The Python replica of that code (`meshpair/exporter.py`) rebuilds the published 2.403 µm w00 mesh from the
    published 9.366 µm mesh and the catalog affine.
  - Grid size 823 × 928 matches; 211 × 3.89819 = 822.5 rounds up, so the rounding is part of the result.
  - The valid mask is identical (0 mismatches out of 552,561 vertices). Vertex error: median 0.0010 voxel, max 0.0044
    voxel, which is float32 rounding.
- **Model-free mesh routes agree:**
  - Route I gives the same answer with no convention assumed. It inverts every fine vertex through the catalog affine
    onto the coarse bilinear surface.
    - The offset against naive pairing at the valid-box centre is (−9.35, −8.74) coarse px = (−87.5, −81.9) µm.
    - It differs from the formula by p95 0.0004 px.
    - The on-surface fraction is 99.997 %, and the split-half noise is 0.0008 px.
  - Route F fits a free affine and the index map together. It recovers the same coefficients, and its affine matches
    the catalog affine to 1e-4 µm.
- **Fit with the earlier CT numbers:**
  - REGISTRATION.md measured CT-to-CT (8.66, 8.49) px "array convention". In canvas convention that is about
    (8.72, 9.46).
  - At strip R the exporter model predicts (9.6, 8.7). A first plumbing tile of `ct.py` (below) gives a naive-map
    offset of (8.62, 9.56). The exact-map residual on that tile is (+1.41, −0.48) px ≈ 13 µm.
  - That residual is what is left for the scan-to-scan affine. It is the same size as TAUIL's 12.8 µm median for the
    transform error on PHerc0139 (villa #1912).
  - TAUIL's PHerc0139 canvas offset ("median −8, −7 native px") matches the −7.8 px constant of the half-pixel rule.
  - His "0.3-0.4 % scale" is larger than w00's rounding drift (0.04-0.07 %). Phase 2 checks whether it comes from size
    rounding or from the affine.
- **Disclosed plumbing CT tile on strip R (seen before the prereg):**
  - One 2 × 2 mm tile at the strip-R centre, read from the cached `R_sv24.npy` / `R_sv94.npy`.
  - Naive map: correction (−8.62, −9.56) px.
  - Exact map: residual (+1.41, −0.48) px, depth offset +10.5 µm (1.12 coarse layers), r 0.84, split-half noise
    0.08-0.09 px.
  - Both tiles were marked invalid by the original runner-up rule: the plain maximum outside 3 px was 0.60-0.62 against
    a bar of 0.588. That maximum was the main peak's own skirt.
  - The rule was changed, before any confirmatory data, to "the highest distinct local maximum outside 3 px" (the
    `ct._peak` docstring has the reason). The tile was not re-run under the new rule.
  - Strip R is excluded from every confirmatory CT line.
- **I/O plumbing:**
  - The remote level-0 read of the w00 2.4 µm surface volume equals the cached crop.
  - Level 2 has 109 layers (no depth downsampling) and equals the 4 × 4 block mean of level 0 to within 1 grey level.
    This confirms the pyramid-centre convention that `ct.py` assumes (level-l pixel i is centred on 2^l·i + (2^l − 1)/2).

So the open question is no longer "is it real on w00". It is: (a) is it the same exporter convention on every
public fine/coarse pair, including pairs where both meshes were transformed from a common source (PHerc0139, 0500P2,
Paris4); (b) how large is the remaining scan-to-scan residual (shift field and depth) once the exact map is used; and
(c) what is the noise floor of that CT measurement.

## 2. Public segment pairs (read-only listing of the catalog, 2026-10-11)

There are 207 segments with a fine (2.0-2.5 µm) and a coarse (7.5-10 µm) tifxyz mesh. The full list with S3 keys for
the fine, coarse and source meshes and both surface volumes is in `out/pairs.tsv` (and `out/pairs.json`), made by
`check.py inventory` from `cache/metadata.json`. Bucket: `s3://vesuvius-challenge-open-data/` (CC BY-NC 4.0,
no auth).

| Sample | n | fine mesh `-on-` | coarse mesh `-on-` | how they were made | SVs (fine/coarse) |
|---|---|---|---|---|---|
| PHerc0841 | 3 | `20260319124803-2.403um` | `20250821151531-9.366um` | coarse = copy of source; fine = transform | 3/3 |
| PHerc0139 | 38 | `20260102150214-2.399um` | `20250728140407-9.362um` | both transforms of the 2.399 source | 38/38 |
| PHerc0814 | 19 | `20260309142202-2.399um` | `20250804134230-9.362um` | 12 traced at 9.4 (coarse = copy), 7 traced at 2.4 (fine = copy) | 19/19 |
| PHerc0500P2 | 46 | `20250526151718-2.215um` | `20250820143440-9.362um` | both transforms of the 4.317 source | 46/46 |
| PHerc0009B | 18 | `20250820154339-2.401um` | `20250521125136-8.64um` | coarse = copy; fine = transform | 18/18 |
| PHerc0343P | 8 | `20260304131111-2.215um` | `20250521134555-8.64um` | coarse = copy; fine = transform | 8/8 |
| PHerc1667 | 20 | `20251217075048-2.399um` | `20231117161658-7.91um` | coarse = copy; fine = transform | 19/0 |
| PHercParis4 | 55 | `20260411134726-2.4um` | `20230205180739-7.91um` | mostly both transforms of a 2.4 source | 54/0 |

- **Path pattern:** `<sample>/segments/<seg-long-id>/mesh/<segid>-on-<volume>.tifxyz/{x,y,z}.tif,meta.json`. The source
  is `<…>/mesh/intermediate/tifxyz_normalized/`. Surface volumes are under `<…>/surface-volumes/<…>.zarr`.
- **Example (w00):**
  - `PHerc0841/segments/20260220213127-w00/mesh/20260220213127-on-20260319124803-2.403um.tifxyz`
  - `…/20260220213127-on-20250821151531-9.366um.tifxyz`
  - `…/surface-volumes/2.403um-0.22m-77keV-volume-20260319124803.zarr`
  - `…/surface-volumes/9.366um-1.2m-113keV-volume-20250821151531.zarr`
- **Strict 2.4 ↔ 9.4 pairs:** 106 (0139, 0814, 0841 and 0500P2 at 2.215 µm).
- **Other pairs:** 26 at 8.64 µm and 75 at 7.91 µm (no coarse surface volume for those).
- **Download sizes (x+y+z, measured by HEAD):** about 16.0 GB for all 621 meshes including sources. PHercParis4 is
  11.7 GB of that. All other samples together are 4.3 GB (PHerc0841 0.03 GB). Meshes only; no volume downloads beyond
  CT tiles.
- **Exporter commits on these meshes:** d9a472ea, 650bd566 and 77627ab1. All three have byte-identical
  transform/resample code.
- **Render commits:** 530584554c, 650bd566, 862bd172, ae2ecce9 and bc28564d. All five have the same canvas convention.

## 3. What is measured

**Mesh route (every pair; CPU seconds per pair; `check.py mesh`, `meshpair/meshcheck.py`):**
- **M0 exporter model:** the exact canvas map from grid sizes and the two conventions. Reports the offset against naive
  index pairing at the canvas origin, the valid-box centre and its extremes, in px and µm, plus the drift b − k.
- **M1 route I (model-free; uses only the published meshes and the published affine):**
  - The correspondence field p_C(p_F) for 60k fine vertices. Gauss-Newton projection onto the coarse bilinear surface,
    with whole-cell restarts (no convention assumed).
  - Per-axis linear fit and its residual.
  - Offset against naive at the centre: **this is the index offset**.
  - Normal distance: **this is the mesh-only depth offset**. It is about 0 by construction when one mesh is an affine
    image of the other's surface. A non-zero value would itself be a finding.
  - On-surface fraction.
- **M2 root cause:**
  - (a) Replica reproduction of each transformed mesh from the published source mesh and its provenance affine: shape,
    vertex error, valid-mask agreement.
  - (b) |route I − exporter model| p95.
- **M3 noise floor of route I:** fit (a, b) on left/right and top/bottom halves; whole-mesh noise = |half_a − half_b| / 2.
- **Route F (descriptive):**
  - A free 3-D affine and the per-axis index map fitted together. It is identified only through sheet curvature.
  - Reported as a cross-check: offset against route I, and its affine against the catalog affine (µm).

**CT route (pairs with both surface volumes; `check.py ct`, `meshpair/ct.py`):**
- **Tiles:** 2 × 2 mm tiles (213 coarse px at 9.4 µm) on a 5 mm lattice, where the coarse mesh is ≥ 98 % valid
  around the tile.
- **Fine data:** the fine surface volume is read at pyramid level 2 (≈ 9.6 µm in-plane, all layers). It is blurred to
  the coarse voxel (σ = 0.42 coarse px), averaged over each coarse layer's depth footprint (4 sub-samples), and
  resampled onto the coarse pixel grid through the map under test.
- **Per tile:**
  - Joint search over depth δ (±28 µm in fine-voxel steps) and in-plane shift (±16 coarse px), using exact FFT NCC of
    DoG(1, 6) stacks of the 8 central coarse layers.
  - Sub-pixel parabolic peak, in-plane and in depth.
- **Outputs per tile:**
  - `map_correction_px` (true p_C − map p_C): the **residual shift field** under the exact map, or the full offset
    under the naive map.
  - `depth_offset_um`: the **depth offset** (the coarse layer at depth d shows fine content at d + δ).
  - r, distinct runner-up, edge flags, valid fractions, split-half estimates.
- **Field:** a 6-parameter affine fit of the residual against tile position, with leave-one-out residuals
  (one correction per segment, or a local field).
- **Positive control (differential plants):** the same tile re-measured with the map moved by (4, 0), (0, 8), (16, 0)
  and (0, −16) px. The recovered Δcorrection must match the plant. Plants are differential, so a tile's own residual
  cannot cancel them.

## 4. Noise floor

- **Mesh route:** the floor is float32 rounding of the coordinates (about 1e-3 voxel). The split-half /2 rule measures
  it per pair. The bar is 0.05 coarse px, about 50× above w00's 0.0008.
- **CT route:**
  - Each reachable tile is measured again on two disjoint depth sub-slabs (lower 4 vs upper 4 of the 8 coarse layers,
    each with its own depth-matched fine slab) and on left/right halves.
  - Each half carries half the signal, so the whole-tile noise is n = |ŝ_a − ŝ_b| / 2 (Euclidean). Per tile,
    noise = min(n_z, n_h).
  - The floor is p95 over reachable tiles.
- **Reachable tiles** are fixed on unshifted data before any plant: valid peak under the exact map, ≥ 90 % non-zero
  voxels in both blocks, inside the lattice margin. They are written to `out/reachable_<seg>.json` with a sha256.
  Every rate is over reachable tiles only (the lesson from REG-PREC-0814's control ceiling).
- **Self-test result (synthetic):** the floor is 0.002-0.004 px at r 0.93-0.99. Real tiles will be noisier; the w00
  plumbing tile gave 0.08-0.09 px.

## 5. Decision lines (frozen in PREREG_DRAFT.md)

- **"In the organizers' files" (M1):** route I's offset against naive at the valid-box centre is ≥ 3 coarse px
  (≈ 28 µm) on ≥ 90 % of a family's non-void pairs. A family is sample × fine px × coarse px.
  - **"Our assumption only":** ≤ 1 px on ≥ 90 % (the files agree with index pairing).
  - **Otherwise:** MIXED, reported per pair.
- **"Root cause in the exporter" (M2):** on ≥ 95 % of a family's pairs, both of these hold:
  - the replica reproduces each transformed mesh (same shape; median ≤ 0.01 voxel; p99 ≤ 0.05 voxel; mask
    agreement ≥ 99.9 %);
  - route I matches the exporter formula (p95 ≤ 0.25 px).

  If (b) holds and (a) fails, the convention is confirmed but the replica differs in detail (for example a newer
  exporter build). That is reported as PARTIAL.
- **What it is, and what it is not, called:**
  - The 3-D geometry is right (route I normal distance about 0; route F affine = catalog affine).
  - The two canvases simply do not share an origin or a scale. This is an undocumented convention with a
    size-rounding drift, not a broken transform.
  - The upstream ask, if confirmed, is documentation and/or writing the canvas map (a, b per axis) into `meta.json`,
    with the 4-line formula. A resampler change is the organizers' call.
  - We claim no bug. The scan-to-scan affine error is a separate quantity, measured by the CT route.
- **CT lines:** they cover the residual against the exact map (EXPLAINED / NOT EXPLAINED), the depth offset
  (|median| > 1 coarse layer → consumers must apply it), the noise floor (≤ 0.67 px PASS, 0.67-1.0 MARGINAL,
  > 1.0 → bars become 3 × p95) and the one-affine field (LOO p95 ≤ 2 px). Exact numbers are in the prereg.

## 6. Compute and data

- **Compute:** CPU only, no GPU claim. The mesh route takes about 1-3 s per pair (w00: 1.3 s on one CPU core) → about
  10 min for 206 pairs.
- **CT:** about 15 s per tile per map at level 0. Level 2 is faster. Roughly 150-250 tiles → about 1-2 h, niced, one
  process.
- **Downloads:**
  - Meshes: 16 GB if all of Paris4 is included, 4.3 GB without. Rule 4 (ask before > 30 GB) is not triggered.
    Meshes go to `cache/meshes/` (git-ignored) and can be deleted after the run (sha256 manifest kept).
  - CT tiles: < 2 GB of chunks through `meshpair/zio.py`'s cache.

(Section 7, open questions for the authors before the confirmatory run, is omitted from this release; they were
answered in the pre-registration.)
