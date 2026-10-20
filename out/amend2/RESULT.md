# W00-FIELD amendment 2: CT route on a fresh C-set (confirmatory repeat)

- **Prereg:** `prereg/2026-10-13_partner_w00field_amend2.md`, commit fcb2012f, pushed before any run. **Tool:** commit 46906121
  (`check.py plants` accepts fcb2012f; `phase2/run_a2.py`, `phase2/ct_table_a2.py` new). Bars, tiles, maps, search and C0 as in 8ca5789f;
  C1 plants (4,0),(0,8),(12,0),(0,-12) as in 13df5696. CPU only, one run, no crash, no resume, no replacement segment used.
- **Records:** `ct_<seg>.jsonl/.log`, `reachable_<seg>.json`, `ct_plants_a1_<seg>.jsonl`, `stage0b_pyramid_<seg>.json`, `ct_summary_a2.json`, `a2.log`.
- Cache after the run: `work/cache` 130 GB (< 150 GB).

## Fresh set (confirmatory rows)
| family | segment | tiles / reachable | C1 | C2 noise p95 (<= 0.67) | C3 exact median / p95 px; naive median | C4 depth median | C5 LOO p95 px (bar 2) |
|---|---|---|---|---|---|---|---|
| 0139 | `20260220000000` | 204 / 186 | PASS | 0.094 | 1.05 / 2.74; 8.3 -> EXPLAINED | +4.2 um (0.45 L), stands | **2.242 -> LOCAL FIELD REQUIRED** |
| 0139 | `20260306000001` | 195 / 176 | PASS | 0.116 | 0.96 / 2.42; 14.4 -> EXPLAINED | +5.0 um (0.53 L), stands | **2.046 -> LOCAL FIELD REQUIRED** |
| 0814 (coarse is copy) | `20250926051122` | 77 / 77 | PASS | 0.083 | 0.84 / 2.69; 12.4 (11 naive-valid) -> EXPLAINED | -1.0 um (-0.11 L), stands | 1.766 -> ONE CORRECTION |
| 0814 (fine is copy) | `20260226121302` | 0 / 0 | n/a | n/a | VOID (< 5 reachable; no lattice tile has a >= 98 % valid coarse block), no replacement | | |
| 0500P2 | `20250825181922` | 34 / 33 | PASS | 0.064 | 0.55 / 1.05; 7.6 -> EXPLAINED | -1.4 um (-0.14 L), stands | 0.70 -> ONE CORRECTION |
| 0500P2 | `20250825181907` | 35 / 33 | PASS | 0.050 | 0.48 / 1.04; 6.3 -> EXPLAINED | -1.2 um (-0.13 L), stands | 0.65 -> ONE CORRECTION |
| 0009B (first CT reading) | `20250919132115` | 115 / 97 | PASS | 0.081 | 1.46 / 2.81; 11.8 -> EXPLAINED | +7.9 um (0.92 L), stands (bar 1.0 L) | 1.42 -> ONE CORRECTION |
| 0009B (first CT reading) | `20250919135433` | 108 / 71 | PASS | 0.086 | 1.04 / 1.54; 13.6 -> EXPLAINED | -1.5 um (-0.17 L), stands | 1.13 -> ONE CORRECTION |

Amended C1: every plant recovered within 0.5 px on 100 % of reachable tiles on all 7 measured segments. Stage 0b passed on all 8.

## Repeat verdicts (rule frozen in the prereg: every fresh C1-PASS segment of the family agrees with the pass-1 family reading)
| verdict | family | pass-1 (exploratory) | fresh | repeats? |
|---|---|---|---|---|
| C3 EXPLAINED | 0139 | 3/3 | 2/2 | **repeats** |
| C3 EXPLAINED | 0814 | 2/2 | 1/1 (fine-copy segment VOID) | **repeats** (coarse-copy kind only) |
| C3 EXPLAINED | 0500P2 | 2/2 | 2/2 | **repeats** |
| C3 EXPLAINED | 0841 | 3/3 | none exists | not tested (stays exploratory) |
| C4 centre-to-centre stands | 0139, 0814 (coarse-copy), 0500P2 | 7/7 | 5/5 | **repeats** |
| C4 centre-to-centre stands | 0841 | 3/3 | none | not tested |
| C5 ONE CORRECTION | 0814 (coarse-copy) | 2/2 | 1/1 (1.77 vs bar 2) | **repeats** (one segment, 0.23 px under the bar) |
| C5 ONE CORRECTION | 0500P2 | 2/2 | 2/2 (0.65, 0.70) | **repeats** |
| C5 LOCAL FIELD REQUIRED | 0139 | 3/3 (2.17-2.23) | 2/2 (2.242, 2.046) | **repeats**, see flag |
| C2 noise <= 0.67 | all | 10/10 | 7/7 (0.05-0.12) | **repeats** |
| 0009B, all lines | 0009B | no CT | 2/2 give C3 EXPLAINED, C4 stands, C5 ONE CORRECTION | first-and-only reading under frozen bars; not a repeat |

## Flag: C5 on PHerc0139 (LOCAL FIELD REQUIRED)
- Fresh LOO p95: 2.242 px and 2.046 px against the 2.0 bar. Pass-1 values: 2.17, 2.22, 2.23. All five sit within 0.25 px above the bar.
- It repeats on the rule as frozen (both fresh segments exceed 2.0), but **`20260306000001` clears the bar by only 0.046 px (its in-sample p95 is
  exactly 2.00)**. The 0139 verdict is "LOCAL FIELD REQUIRED" as a bar outcome; it is not evidence of a large field. Median LOO is 1.05-1.14 px.
  A reader should treat it as "a single affine is borderline on 0139", not as a clear field. The bar was frozen, not tuned.
- Other points near a line: 0009B `20250919132115` depth +0.92 layers (bar 1.0, so "stands" by 0.08 L); 0814 `20250926051122` LOO 1.77.

## Wording
- Repeated on the fresh set: C2, C3, C4 (all measured families), C5 per family as in the table. These are the confirmatory CT readings. They cover
  0139, 0814 (coarse-copy tracing only), 0500P2, and, as a first reading, 0009B.
- Not covered by any fresh segment: PHerc0841 (all its CT rows stay exploratory), 0814 fine-copy tracings (the one fresh candidate is VOID), and
  everything on Paris4/1667/0343P (no CT run). The pass-1 rows (`../RESULTS.md` Amendment 1 section) are the exploratory discovery set.
