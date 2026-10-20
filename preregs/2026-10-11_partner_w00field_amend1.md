# Amendment 1 to `prereg/2026-10-11_partner_w00field.md` (8ca5789f): the CT positive control (C1)

- Date: 2026-10-11, written after the frozen phase-2 run (6bb9723c) and pushed before the re-run below.
- **Disclosed:** C2-C5 were not scored as verdicts, but their values were printed as descriptive numbers in
  `work/w00field/out/RESULTS.md` (6bb9723c), and both partners' sessions have seen them. So this amendment is **not blind**
  to the outcome. Sebastien's Judge decides whether a CT verdict read after this amendment counts as confirmatory or only
  as exploratory. Until he rules, it is reported as "post-sight amendment".

## Why
C1 planted shifts of (16, 0) and (0, −16) px while the CT search is ±16 px. Once the tile's base residual points the
same way, the planted peak lies past the search window and is clamped at the edge: 472 of the 474 failed 16 px
recoveries are of this kind. The 4 px and 8 px plants were recovered on 100 % of reachable tiles (0 failures out of 1,506
each). C1 failed because of how the control was designed, not because the measurement failed.

## Change (the only one)
- C1 plants become **(4, 0), (0, 8), (12, 0) and (0, −12) px**. Everything else in C1 is unchanged: each plant
  recovered within 0.5 px on ≥ 95 % of reachable tiles, and FAIL means that segment's C2-C5 are not read.
- One option chosen, none tried: the ±16 px search radius, the reachable-tile lists (`reachable_<seg>.json`, sha256 in
  the CT logs), the maps, the sampling and the C0, C2, C3, C4 and C5 bars stay exactly as frozen. C2-C5 are read from
  the pass-1 records already on file (`out/ct_<seg>.jsonl`). Nothing is re-measured except the plants.
- Re-run: plants only, on the frozen reachable tiles of the same 10 C-set segments, with the tool at the commit that
  adds the plants-only command. That command reads the base (unshifted, exact-map) records from the existing JSONL and
  does not rewrite them. It writes `out/ct_plants_a1_<seg>.jsonl`.
- The mesh route and its verdicts are untouched. Route I's on-surface rule (VOID on pairs whose coarse mesh is a
  downsample) is **not** amended: the frozen tree's "NOT EXPLAINED" on 0139, 0500P2, and parts of 0814 and Paris4
  stands, reported next to the M2a replica result.
