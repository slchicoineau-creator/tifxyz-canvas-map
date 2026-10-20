# Pre-registration (pushed 2026-10-11 before any confirmatory run): W00-FIELD + mesh-pair offset checker (2.4 ↔ 7.9/8.6/9.4 µm)

- **Date:** 2026-10-11 draft.
  - The main session pushes this file to `prereg/` before any confirmatory run.
  - The tool refuses confirmatory pairs without `--prereg <pushed commit>` (`check.py _guard`).
- **Owner:** Chris's side ([machine name removed]), team contribution (Chris Thomas + Sebastien Chicoineau).
- [Team-process lines removed: how this test relates to the team's internal task list. It was run by hand and keeps
  an earlier internal design (/2 split noise, reachable tiles first, differential plants, one-affine-or-local).]
- **Design and code:** `work/w00field/DESIGN.md`; tool `work/w00field/check.py` + `meshpair/` (version 0.1.0); upstream
  code read at ScrollPrize/villa main d63d3e968 plus the commits named in the catalog provenance.
- **Data:**
  - Public open data only: `s3://vesuvius-challenge-open-data/`, catalog `metadata.json` with sha256
    `7ec8dcaa3000b9716724c26a1991ee416491d8a6b3615238ca90ae11540abdd3` (fetched 2026-10-11; a copy is kept git-ignored
    in `work/w00field/cache/`).
  - Inputs are meshes and surface-volume tiles only. No labels, ink maps, keys, site lists, ratings or [an internal file class removed].
  - [Internal data names removed.] CPU only.

## Development pair and what was seen before this prereg (disclosed)
- **PHerc0841 w00 (`20260220213127`) is the development pair.** Seen on it:
  - the exporter replica reproduces its 2.403 µm mesh exactly;
  - route I/F offsets (−9.35, −8.74) px at the valid-box centre;
  - one CT plumbing tile at the strip-R centre: naive (−8.62, −9.56) px; exact-map residual (+1.41, −0.48) px; depth
    +10.5 µm = 1.12 coarse layers; marked invalid by the then runner-up rule.
- **Changes made after that tile:**
  - The runner-up rule was changed to "highest distinct local maximum outside 3 px". The old rule counted the main
    peak's own skirt.
  - The changed tile was not re-measured.
- **Exclusions that follow:**
  - w00 is excluded from the mesh-route family counts.
  - Strip R (fine canvas y 10208-12992, x 9280-15776) is excluded from every confirmatory CT line.
  - Earlier team measurements on strips R/R2 (REGISTRATION.md, FACTS) are prior knowledge, not evidence here.

## Questions
1. **In the files?** Is the (u, v) canvas offset between the fine and coarse meshes of the same segment present in the
   organizers' published meshes for every public pair, or only an artefact of index pairing (p_C = p_F · k)?
2. **Root cause?** Is it fully explained by the exporter's grid convention?
   - The convention: half-pixel resampling with a rounded grid size in `transformSurfacePoints` → `QuadSurface::resample`,
     plus the renderer's `grid = s (p + 0.5)`.
   - Explained means the exact map p_C = (p_F + 10.5) · (n_C/n_F) − 10.5 per axis (for s = 0.05, R = 1).
3. **What is left?** Once the exact map is used:
   - the residual in-plane shift field and the depth offset between the two surface volumes (scan-to-scan affine
     error);
   - their noise floor;
   - whether one affine correction per segment is enough.

## Sets (frozen now)
- **M-set (mesh route):** all 207 segments that `check.py inventory` lists from the frozen catalog
  (`out/pairs.tsv`: fine 2.0-2.5 µm, coarse 7.5-10 µm, primary pair = the mesh on the largest volume of each class),
  minus the development pair w00. That leaves 206.
  - Families are sample × (fine px, coarse px), 8 families.
  - All Paris4 pairs are included (11.7 GB of meshes). Meshes are deleted after their record is written; the sha256
    manifest is kept.
- **C-set (CT route; surface volumes for both meshes required):**
  - (i) w00, whole-canvas 5 mm lattice minus strip R (31 tiles by the dry run);
  - (ii) PHerc0841 `20260220214732` and `20260221022814`, whole lattice;
  - (iii) the 3 PHerc0139 segments with the most valid fine vertices in their M-route record;
  - (iv) the 2 PHerc0814 segments with the most valid fine vertices, one traced at 9.4 µm and one traced at 2.4 µm
    if both kinds exist;
  - (v) the 2 PHerc0500P2 segments with the most valid fine vertices.
  - The choice is made only from M-route records (vertex counts). No CT output is used to choose.
- **Tile settings:** 213 coarse px (≈ 2 mm), lattice 534 px (≈ 5 mm), coarse mesh ≥ 98 % valid around the tile,
  search ±16 px, depth ±28 µm, DoG(1, 6), fine surface volume at pyramid level 2, 8 central coarse layers.

## Stage 0 (VOID if it fails; no science read before it passes)
- **0a. Synthetic self-test:** `check.py selftest` → all pytest cases pass at the tool commit used for the run.
  - Covers the replica against an independent float64 half-pixel resampler.
  - Covers route I recovering the exporter formula without being told it.
  - Covers routes I and F recovering planted index shifts of 0.2-0.6 cells within 0.02 / 0.05 px, and the replica
    flagging a planted file.
  - Covers CT recovering planted shifts within 0.2 px and depth within 0.2 coarse layer, at level 0 and at level 2.
  - Covers differential plants within 0.5 px, the naive map returning the exporter offset, and the runner-up rules.
- **0b. Pyramid convention check:** for every surface volume used in the C-set, one 6 × 40 × 40 block of level 2 must
  equal the 4 × 4 block mean of level 0 within 1 grey level, and level 2 must have all layers. If not, CT for that
  scroll is VOID until the reader handles its pyramid.
- **0c. Provenance:**
  - Each pair's affine comes from its own provenance (`affine-matrix`).
  - A pair with no affine is VOID for M1/M2 and listed.

## Mesh route bars (per pair, then per family)
- **M1 (in the files):** route I offset against naive at the fine valid-box centre, as a Euclidean norm in coarse px:
  - ≥ 3.0 → IN_FILES;
  - ≤ 1.0 → NOT_IN_FILES;
  - otherwise → MIXED;
  - VOID_ROUTE_I if fewer than 95 % of the 60k sampled fine vertices find a coincident coarse surface point
    (|normal| ≤ 0.05 coarse vx, tangential < 1e-3).
  - **Family verdict:**
    - CONFIRMED IN FILES if ≥ 90 % of non-void pairs are IN_FILES;
    - OUR ASSUMPTION if ≥ 90 % are NOT_IN_FILES;
    - otherwise MIXED (listed per pair).
- **M2a (replica):** every transformed mesh of the pair is rebuilt from the published source mesh and its provenance
  affine. It passes when all of these hold:
  - same grid shape;
  - vertex error median ≤ 0.01 and p99 ≤ 0.05 target-volume voxels;
  - valid-mask agreement ≥ 0.999.

  Copies are compared directly with the source.
- **M2b (formula):** p95 |route I p_C − exact-map p_C| ≤ 0.25 coarse px.
- **Family root-cause verdict:**
  - EXPORTER CONVENTION CONFIRMED if ≥ 95 % of the family's pairs pass M2a and M2b;
  - PARTIAL if ≥ 95 % pass M2b but not M2a (the convention holds, the replica differs in detail);
  - NOT EXPLAINED otherwise. The failing pairs are listed with their numbers, and the cause is diagnosed (no bar
    change).
- **M3 (noise floor):** route I split-half /2 noise (left/right, top/bottom) ≤ 0.05 px. A pair above that is flagged
  NOISY and is left out of M1 counts.
- **Descriptive only:** route F coefficients, route F affine against the catalog affine (µm), mesh-only normal
  distance, offset range across the canvas (drift).

## CT route bars (per segment of the C-set)
- **C0 (reachable tiles):** fixed on unshifted data under the exact map, before any plant. A tile is reachable when all
  of these hold:
  - valid peak: r ≥ 0.30;
  - distinct runner-up ≤ 0.70 r;
  - not on the search edge;
  - split spread ≤ 1.0 px;
  - autocorrelation at the shift ≤ 0.30 when |s| ≥ 3 px;
  - ≥ 90 % non-zero voxels in both blocks.

  The list is written to `reachable_<seg>.json` with its sha256. If fewer than 5 reachable tiles remain, that segment
  is VOID for the CT lines.
- **C1 (positive control):** differential plants of (4, 0), (0, 8), (16, 0) and (0, −16) px must each be recovered
  within 0.5 px on ≥ 95 % of reachable tiles. If not, the CT route FAILS for that segment and its C2-C5 are not read.
- **C2 (noise floor):** p95 over reachable tiles of min(n_z, n_h), with n = |ŝ_a − ŝ_b| / 2.
  - ≤ 0.67 px → PASS.
  - 0.67-1.0 → MARGINAL (reported).
  - > 1.0 → every px bar below becomes 3 × that p95. Stated now; no new prereg.
- **C3 (does the exporter map explain the CT offset?):**
  - EXPLAINED if both hold:
    - median |exact-map residual| ≤ 2.0 px;
    - median |naive-map correction| ≥ 3 × median |exact-map residual|.
  - NOT EXPLAINED otherwise.
  - The residual is reported in µm next to TAUIL's PHerc0139 transform error (12.8 µm median, villa #1912).
- **C4 (depth):** median depth offset over reachable tiles.
  - |median| > 1.0 coarse layer → DEPTH OFFSET REAL. Consumers that pair the two surface volumes per pixel
    ([an internal test name] and any 2.4 ↔ 9.4 comparison on that segment) must apply it.
  - Otherwise the frozen depth correspondence (centre layer to centre layer) stands.
  - For w00, the strip-R glimpse (1.12 layers, one tile, invalid at the time) is prior information only.
- **C5 (one affine?):** 6-parameter affine of the exact-map residual against tile centre, with leave-one-out
  residuals.
  - LOO p95 ≤ 2.0 px → ONE CORRECTION PER SEGMENT suffices.
  - Otherwise LOCAL FIELD REQUIRED. Consumers outside the measured tiles must register locally.
- **Descriptive only:** per-tile field (pushed as CT-only numbers), naive-map offsets against the exporter-model
  prediction at each tile, and hp-filter variants.

## Overall readout (what each outcome means)
- **IN FILES + EXPORTER CONVENTION CONFIRMED:**
  - The [internal register reference removed] claim stands: the offset is in the published meshes on every family where that holds, and its
    root cause is the exporter's half-pixel / rounded-size grid convention.
  - The deliverable is the checker plus the exact canvas map, documented (proposed upstream ask: document it or write
    a, b per axis into `meta.json`).
  - We do not call it a bug. The 3-D geometry is right.
- **OUR ASSUMPTION:** the team's index pairing caused the offset. [internal register reference removed] becomes internal only (fix our scripts).
- **MIXED or NOT EXPLAINED:** report per family. No upstream claim for those families.
- **CT EXPLAINED with a small residual:** the residual is the scan-to-scan affine error, reported as a field. It is a
  separate quantity from the index offset.

## Not claimed
Nothing about ink, sheet identity, legibility or any reader. Nothing about scrolls outside the listed families. No
claim that the organizers' transform is "wrong". The index offset is a grid convention; the affine residual is
measured, not judged.

## Compute rules
- CPU only, niced, one process. No GPU claim.
- Meshes about 16 GB in total, deleted after use. CT chunks < 2 GB through `lib/zio.py`'s cache under
  `work/w00field/cache/` (git-ignored).
- Mesh route about 10 min; CT about 1-2 h (`nohup`, log polled).

## Release
Results go to `work/w00field/out/` (records JSONL; small, pushed) and a short `RESULTS.md`. Any outward note (#1912
comment or its own issue, small tool repo) is drafted only. Both partners OK it, and a human writes the commentary and
posts it (villa CONTRIBUTING). [Internal data name removed.]
