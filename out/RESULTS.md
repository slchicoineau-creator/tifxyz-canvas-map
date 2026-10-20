# W00-FIELD + mesh-pair offset checker: confirmatory results (phase 2, 2026-10-11)

- **Prereg:** `prereg/2026-10-11_partner_w00field.md`, commit 8ca5789f in the authors' private repository, pushed before any
  confirmatory run.
- **Tool:** the working tree equals 8ca5789f for `work/w00field` (checked with `git diff --quiet`). The only change is
  the environment fix listed under Deviations.
- **Catalog:** sha256 7ec8dcaa…abdd3. `pairs.tsv` regenerated identical to the committed one.
- **Machine:** CPU only, niced, one measuring process at a time. No GPU claim. Meshes were deleted after their records
  were written (only the dev pair w00 stays cached).
- **Raw records:**
  - `out/mesh_records.jsonl` (206 records)
  - `out/ct_<segment>.jsonl` (10 segments)
  - `out/reachable_<segment>.json` (sha256 printed in the CT logs)
  - `out/stage0*.json`
- **Summaries:**
  - `out/mesh_summary.json`, `out/ct_summary.json`, `out/ct_table.txt`
  - Made by `phase2/analyze.py` and `phase2/ct_table.py`, which apply the frozen bars mechanically.

## Headline
1. **The offset is in the organizers' files.**
   - It holds on all 6 families where the mesh route was readable (PHerc0841, 0009B, 0343P, 0814, 1667, Paris4).
   - On 4 of those families (0841, 0009B, 0343P, 1667) the exporter-convention root cause is CONFIRMED under the
     frozen bars.
   - Across all 206 segments, a replica of villa's `transformSurfacePoints` → `QuadSurface::resample` rebuilds every
     published transformed mesh exactly from `tifxyz_normalized` and the catalog affine: same grid shape, identical
     valid mask, vertex error p99 ≤ 0.008 voxel.
2. **Mesh route VOID on 2 families, by a design flaw.** On 0139 and 0500P2 (and on 7/19 of 0814 and 44/55 of Paris4)
   route I is VOID: < 95 % of fine vertices lie on the coarse surface.
   - Cause: in those pairs the coarse mesh is a down-sampled resample of a common source.
   - So the fine vertices sit off the coarse bilinear chords. Normal distance median ≈ 0.2 µm, p95 up to ~6 µm, against
     the 0.05 coarse-voxel (~0.5 µm) on-surface bar.
   - The bar was written with copy-type pairs in mind. This is reported, not fixed (see Deviations).
3. **The CT route FAILED its positive control (C1) on all 10 C-set segments, by a design error in the prereg.**
   - The 16 px plants equal the ±16 px search radius. 472 of the 474 failed 16 px plant recoveries needed a shift
     beyond the search window (base residual + plant > 15.5 px), so the peak was clamped at the edge.
   - The 4 px and 8 px plants were recovered on 100 % of tiles (0 failures out of 1,506 each).
   - Per the frozen rule, C2-C5 are **not read**. Their values are shown below as descriptive only.

## Stage 0
| stage | result |
|---|---|
| 0a synthetic self-test at the tool commit | PASS: 15/15 pytest (`out/stage0a_selftest.json`) |
| 0b pyramid convention | PASS on 9/10 C-set surface volumes: level 2 = 4×4 block mean within 1 grey level, all layers. For PHerc0814 `20250928235953` the tool's centre block was empty (0 % non-zero); see Deviation 3 |
| 0c provenance | PASS: all 207 pairs have their affine and source key in provenance; no pair VOID for lack of an affine |

## Mesh route (M-set = 206 segments, w00 excluded as the development pair)
- **Bars (frozen):**
  - **M1:** route-I offset against naive at the valid-box centre: ≥ 3 px → IN_FILES, ≤ 1 px → NOT_IN_FILES.
    VOID_ROUTE_I if fewer than 95 % of vertices are on-surface.
  - **Family M1:** ≥ 90 % of non-void pairs IN_FILES → CONFIRMED IN FILES.
  - **M2a:** replica matches: same shape, median ≤ 0.01 vx, p99 ≤ 0.05 vx, mask ≥ 0.999.
  - **M2b:** route I vs formula p95 ≤ 0.25 px.
  - **Family M2:** ≥ 95 % of pairs pass M2a + M2b → CONFIRMED; ≥ 95 % pass M2b only → PARTIAL; otherwise NOT EXPLAINED.
  - **M3:** split-half noise ≤ 0.05 px.

| family (fine/coarse µm) | n | M1 counts | M1 family verdict | M2a pass | M2b pass / evaluable | M2 family verdict |
|---|---|---|---|---|---|---|
| PHerc0841 2.403/9.366 | 2 | IN 2 | **CONFIRMED IN FILES** | 2/2 | 2/2 | **EXPORTER CONVENTION CONFIRMED** |
| PHerc0009B 2.401/8.64 | 18 | IN 18 | **CONFIRMED IN FILES** | 18/18 | 18/18 | **EXPORTER CONVENTION CONFIRMED** |
| PHerc0343P 2.215/8.64 | 8 | IN 8 | **CONFIRMED IN FILES** | 8/8 | 8/8 | **EXPORTER CONVENTION CONFIRMED** |
| PHerc1667 2.399/7.91 | 20 | IN 19, VOID 1 | **CONFIRMED IN FILES** (19/19 non-void) | 20/20 | 19/19 | **EXPORTER CONVENTION CONFIRMED** (19/20 = 0.95, exactly at the bar) |
| PHerc0814 2.399/9.362 | 19 | IN 12, VOID 7 | **CONFIRMED IN FILES** (12/12 non-void) | 19/19 | 12/12 | NOT EXPLAINED by the frozen tree (12/19 < 0.95; the 7 failures are VOID, not mismatches) |
| PHercParis4 2.4/7.91 | 55 | IN 11, VOID 44 | **CONFIRMED IN FILES** (11/11 non-void) | 55/55 | 11/11 | NOT EXPLAINED by the frozen tree (11/55; the 44 failures are VOID) |
| PHerc0139 2.399/9.362 | 38 | VOID 38 | **VOID** (no non-void pair) | 38/38 | 0/0 | NOT EXPLAINED by the frozen tree (M2b not evaluable) |
| PHerc0500P2 2.215/9.362 | 46 | VOID 46 | **VOID** (no non-void pair) | 46/46 | 0/0 | NOT EXPLAINED by the frozen tree (M2b not evaluable) |

M3 (noise): no pair was NOISY. The split-half noise on readable pairs was ≤ 0.002 px.

On all 70 readable (IN_FILES) pairs, route I matches the exporter formula to p95 ≤ 0.003 coarse px (median 0.0007).

**VOID pattern (diagnosed):** route I is VOID exactly where the coarse mesh is not a plain copy of the source.

| family | pairs where coarse is a copy | pairs where fine is a copy or both are transformed |
|---|---|---|
| 0139 | 0 | 38/38 VOID |
| 0500P2 | 0 | 46/46 VOID |
| 0814 | 12, all readable | 7/7 VOID |
| Paris4 | 11, all readable | 44/44 VOID |
| 1667 | 19 of 20 readable (1 VOID) | none |

- The one exception is PHerc1667 `w011…flatboi_straightened_v4`. It has 46.8 M fine vertices and an offset of
  ~6 mm against naive, and route I fails to converge there. It does not reach on-surface points from the naive start
  guess with ±2-cell steps and ±1-cell restarts. Its normal distances are ~0.
- On the VOID pairs, descriptive numbers still agree with the exporter formula:
  - route I on its on-surface subset matches the formula to p95 0.02-0.05 px (0139), median 0.10 px (0500P2);
  - route F (free affine) matches it to a median 0.013 px (0139);
  - M2a reproduces every mesh exactly.

  So nothing contradicts the convention there, but the frozen tree cannot confirm it.

**Offset magnitudes (exporter model at the valid-box centre, Euclidean, µm; min / median / max over segments):**

| family | min / median / max | max per-axis over the valid box (median) | relative scale drift (b − k)/k |
|---|---|---|---|
| 0841 | 102 / 106 / 109 | 83 | 0.03 % |
| 0009B | 97 / 112 / 124 | 110 | 0.08 % |
| 0343P | 85 / 89 / 100 | 71 | 0.10 % |
| 0139 | 17 / 69 / 115 | 197 (max 301) | 0.39 % |
| 0500P2 | 9 / 56 / 133 | 111 (max 230) | 0.57 % |
| 0814 | 39 / 132 / 711 | 270 (max 1,494) | 0.33 % |
| 1667 | 191 / 304 / 5,954 | 586 (max 11,961) | 0.83 % |
| Paris4 | 667 / 2,342 / 3,060 | 4,483 | 1.23 % |

- The offset has two parts.
  - The constant half-pixel part is −7.3 to −7.8 coarse px at the canvas origin on every family.
  - The slope part comes from grid sizes following the affine's measured scale, not the nominal voxel-size ratio.
    It dominates on large canvases and on the 7.91 µm pairs.
- On PHerc0139 the drift (0.2-0.5 %) matches TAUIL's "0.3-0.4 %" in villa #1912.
- **Route F affine against catalog:** median ≤ 0.6 µm on 0009B/0343P/0841/1667/0139. On some 0500P2/0814/Paris4 pairs it
  reaches tens of µm. This is descriptive: route F is identified only through curvature.

## CT route (C-set = 10 segments; strip R excluded on w00)
- **Selection (from M-route vertex counts only):**
  - 0139: w054 `20260206000002`, w055 `20251226000000`, w052 `20260227000000`.
  - 0814: `20250928235953` (traced at 9.4 µm), `20260226000000` (traced at 2.4 µm).
  - 0500P2: `20250628074500`, `20250611170257`.
  - Plus PHerc0841 w00, `20260220214732`, `20260221022814`.
- **Bars:**
  - C1: each plant recovered within 0.5 px on ≥ 95 % of reachable tiles. FAIL means C2-C5 are not read.
  - C2: noise p95 ≤ 0.67 px.
  - C3: median exact residual ≤ 2 px and naive ≥ 3 × residual.
  - C4: |median depth| > 1 coarse layer.
  - C5: LOO p95 ≤ 2 px.

| segment | tiles / reachable | **C1 (frozen verdict)** 16,0 / 0,−16 recovery | C2 noise p95 (desc.) | exact residual median / p95 px (desc.) | naive median px | C4 depth median (desc.) | C5 LOO p95 px (desc.) |
|---|---|---|---|---|---|---|---|
| 0841 w00 (minus strip R) | 31 / 30 | **FAIL** 0.23 / 0.73 | 0.145 | 1.35 / 2.21 | 11.8 | +2.2 µm (0.24 L) | 1.51 |
| 0841 `…214732` | 40 / 38 | **FAIL** 0.26 / 0.63 | 0.086 | 1.38 / 2.11 | 10.3 | +1.0 µm (0.11 L) | 1.11 |
| 0841 `…022814` | 35 / 32 | **FAIL** 0.13 / 0.75 | 0.117 | 1.45 / 2.10 | 10.9 | +1.7 µm (0.18 L) | 1.52 |
| 0139 w054 | 222 / 205 | **FAIL** 0.76 / 0.70 | 0.106 | 1.07 / 2.72 | 13.7* | +5.1 µm (0.55 L) | 2.23 |
| 0139 w055 | 208 / 189 | **FAIL** 0.70 / 0.71 | 0.115 | 1.10 / 2.68 | 13.0* | +3.9 µm (0.41 L) | 2.17 |
| 0139 w052 | 209 / 197 | **FAIL** 0.72 / 0.71 | 0.108 | 1.04 / 2.72 | 12.3* | +4.6 µm (0.49 L) | 2.22 |
| 0814 `…235953` | 103 / 102 | **FAIL** 0.97 / 0.90 | 0.090 | 0.77 / 2.40 | 13.8* | −3.4 µm (−0.37 L) | 1.59 |
| 0814 `…000000` (2um_try2) | 8 / 8 | **FAIL** 1.00 / 0.875 | 0.070 | 0.45 / 0.98 | 10.1 | −3.5 µm (−0.38 L) | 0.65 |
| 0500P2 `…074500` | 37 / 33 | **FAIL** 0.88 / 0.55 | 0.102 | 0.69 / 1.05 | 6.0 | −1.6 µm (−0.17 L) | 0.76 |
| 0500P2 `…170257` | 24 / 21 | **FAIL** 0.81 / 0.62 | 0.092 | 0.57 / 1.07 | 5.6 | −1.7 µm (−0.19 L) | 0.98 |

- **Plants:** 4 px and 8 px plants were recovered on 100 % of reachable tiles in every segment.
- **Starred naive values:** the naive offset exceeds the ±16 px search on many tiles (naive valid 60-65 % on 0139, 8/102 on
  0814 `…235953`), so the naive median is capped there.
- **Descriptive only (C1 failed, so no verdict is drawn):**
  - In every segment the exact-map residual is ~1 px (0.45-1.45 px ≈ 4-14 µm), against 6-14 px under index pairing.
  - 0139's ~10 µm sits in the range of TAUIL's 12.8 µm transform-error median.
  - The noise floors (0.07-0.15 px) are well under 0.67 px.
  - The depth medians are all within ±0.6 coarse layer, so no segment comes near the 1-layer line.
  - The strip-R glimpse on w00 (1.12 layers, one tile) is not repeated anywhere on the rest of w00 (median +0.24 layer,
    IQR −0.9..+4.6 µm).
  - The 0841 residual is systematic in +y (~+1.1-1.2 px ≈ 11 µm on all three segments).
  - The 0139 fields would fail C5's 2 px line (2.17-2.23) and need a local field. 0841, 0814 and 0500P2 would pass it.

## Deviations (all reported, none silently fixed)
1. **Environment fix during the M-run.**
   - The PHerc0139 tifxyz files are LZW-compressed. tifffile needs `imagecodecs`, so 4 segments crashed with
     `ValueError: LZW requires imagecodecs` before any record was written.
   - I installed `imagecodecs==2026.8.16` into `work/w00field/.venv` (pinned below the newest release on purpose) and
     resumed. The 4 segments were re-run cleanly. No code or bar changed.
   - Phase 1's venv list should add imagecodecs.
2. **Route-I on-surface bar (design flaw, not changed).** See "VOID pattern". A post-hoc fix needs an amendment, for
   example: accept pairs whose |normal| is within the chord sagitta of the coarse grid, or invert through the source
   mesh. That is for the main session; nothing was re-scored.
3. **Stage 0b on PHerc0814 `20250928235953` (judgement call, disclosed).**
   - The tool's centre block was empty. The tool required > 50 % non-zero, which is not in the prereg text; the prereg's
     own conditions (difference ≤ 1 grey level, all layers) held vacuously.
   - I then ran a non-vacuous check on a data block (level-2 position (1224, 12824)): difference 1.0, 109/109 layers →
     PASS.
   - I ran that segment's CT afterwards with the same frozen settings. Its C1 FAILs like every other segment.
4. **C1 design error (not changed).** Plants of 16 px with a ±16 px search cannot be recovered once the base residual
   points the same way. 472/474 failures are of this kind. A valid control needs an amendment, for example S = 24 or
   plants of 4/8/12 px, with a plants-only re-run on the frozen reachable tiles. C2-C5 above would then become readable.
   I did not run it.
5. **Ops (no effect on numbers).**
   - A stray wait-loop from my own `pgrep -f` self-match was found and killed. It never launched a duplicate run.
   - Total CT wall time was ~4 h, against the prereg's 1-2 h estimate. The 0139 and 0814 canvases are large: about
     200 tiles per 0139 segment.

## What this supports (within the frozen tree)
- **Supported:** the open caveat (is the offset in the files, or in our mapping assumption?) is resolved. The index offset is in the organizers' published meshes, CONFIRMED
  on 6 families.
- **Supported:** its root cause is the exporter's half-pixel / rounded-size grid convention, CONFIRMED on 4 families
  (0841, 0009B, 0343P, 1667).
  - On the other 4 families the tree reads NOT EXPLAINED only because route I is VOID where the coarse mesh is a
    downsample. No mismatch was found anywhere, and the replica reproduces all 206.
- **Not supported:** no CT verdict (C2-C5) is supported, because C1 failed by design. Nothing about the scan-to-scan
  residual or depth can be claimed from this run.


## Amendment 1 (post-sight; Judge to rule)
- **Status:** prereg `2026-10-11_partner_w00field_amend1.md` (13df5696), tool commit d113f600. Not blind: C2-C5 values were printed
  descriptively in the section above before the amendment. Whether the verdicts below count as confirmatory or exploratory is the
  Judge's ruling. No bar, map, search radius, spec or reachable list was changed; the pass-1 records were not rewritten.
- **What was re-measured:** only the plants (4,0), (0,8), (12,0), (0,-12) px, on the frozen reachable tiles (855 in all, sha256 of each
  list checked against the pass-1 logs). Records: `out/ct_plants_a1_<seg>.jsonl`; scorer `phase2/ct_table_a1.py`
  -> `out/ct_summary_a1.json`; run log `out/a1.log`. All 10 runs rc 0.
- **Amended C1: PASS on all 10 segments.** Every plant was recovered within 0.5 px on 100 % of reachable tiles (all 4 plants,
  every segment, record count = reachable count). C2-C5 are therefore read from the pass-1 records with the frozen bars:

| segment | reach. | C1 (A1) | C2 noise p95 px (<= 0.67) | C3 exact median / p95 px; naive median | C3 | C4 depth median | C4 | C5 LOO p95 px (<= 2) | C5 |
|---|---|---|---|---|---|---|---|---|---|
| 0841 w00 (minus strip R) | 30 | PASS | 0.145 PASS | 1.35 / 2.21; 11.8 | EXPLAINED | +2.2 um (0.24 L) | stands | 1.51 | ONE CORRECTION |
| 0841 `...214732` | 38 | PASS | 0.086 PASS | 1.38 / 2.11; 10.3 | EXPLAINED | +1.0 um (0.11 L) | stands | 1.11 | ONE CORRECTION |
| 0841 `...022814` | 32 | PASS | 0.117 PASS | 1.45 / 2.10; 10.9 | EXPLAINED | +1.7 um (0.18 L) | stands | 1.52 | ONE CORRECTION |
| 0139 w054 `...206000002` | 205 | PASS | 0.106 PASS | 1.07 / 2.72; 13.7* | EXPLAINED | +5.1 um (0.55 L) | stands | 2.23 | LOCAL FIELD REQUIRED |
| 0139 w055 `...251226000000` | 189 | PASS | 0.115 PASS | 1.10 / 2.68; 13.0* | EXPLAINED | +3.9 um (0.41 L) | stands | 2.17 | LOCAL FIELD REQUIRED |
| 0139 w052 `...260227000000` | 197 | PASS | 0.108 PASS | 1.04 / 2.72; 12.3* | EXPLAINED | +4.6 um (0.49 L) | stands | 2.22 | LOCAL FIELD REQUIRED |
| 0814 `...235953` | 102 | PASS | 0.090 PASS | 0.77 / 2.40; 13.8* | EXPLAINED | -3.4 um (-0.37 L) | stands | 1.59 | ONE CORRECTION |
| 0814 `...260226000000` | 8 | PASS | 0.070 PASS | 0.45 / 0.98; 10.1 | EXPLAINED | -3.5 um (-0.38 L) | stands | 0.65 | ONE CORRECTION |
| 0500P2 `...074500` | 33 | PASS | 0.102 PASS | 0.69 / 1.05; 6.0 | EXPLAINED | -1.6 um (-0.17 L) | stands | 0.76 | ONE CORRECTION |
| 0500P2 `...170257` | 21 | PASS | 0.092 PASS | 0.57 / 1.07; 5.6 | EXPLAINED | -1.7 um (-0.19 L) | stands | 0.98 | ONE CORRECTION |

- **Reading:** C3 EXPLAINED on 10/10 (exact-map residual about 0.5-1.5 px median, naive pairing 5.6-13.8 px); C4 "centre-to-centre
  stands" on 10/10 (all depth medians within +-0.55 coarse layer); C5 ONE CORRECTION on 7/10, LOCAL FIELD REQUIRED on the three
  PHerc0139 segments (LOO p95 2.17-2.23 px against the 2 px bar). C2 passes everywhere (noise p95 0.07-0.15 px).
- *Naive median is capped by the +-16 px search on 0139 and 0814 `...235953` (many naive tiles invalid or on the edge), as before;
  the C3 ratio bar is met with margin anyway.
- **Deviation 3 stands:** the 0814 `...235953` Stage 0b pyramid check had an empty centre block in phase 2. Its plants were run
  anyway because it is on the frozen list; treat its row with that caveat.
- **Deviation (runner):** the queue script `phase2/run_a1.py` (fetch, plants, delete meshes; w00 kept) is new; it only calls the
  committed `check.py plants`. Meshes were deleted after each segment (cache 12 MB).
