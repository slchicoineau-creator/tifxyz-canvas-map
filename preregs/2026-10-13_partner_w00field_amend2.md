# Amendment 2 to `prereg/2026-10-11_partner_w00field.md` (8ca5789f): confirmatory CT route on a fresh C-set

- Date: 2026-10-13. Written and pushed BEFORE any measurement on the segments below. Author: Chris's Claude (partner).
- Route: Sebastien's Judge ruling `reviews/2026-10-11_seb_judge_w00field_amend1.md` §5. The Judge ruled that the CT verdicts
  C2-C5 of pass 1 (6bb9723c, with amendment 1, 13df5696) are **exploratory** (the C1 gate was changed after C2-C5 had been
  printed). This amendment freezes a repeat on segments that no CT measurement has touched.
- No CT output on any segment below exists or has been looked at. The segments were chosen from M-route vertex counts only
  (`out/mesh_records.jsonl`, `fine_valid`) and `out/pairs.json` (surface-volume presence, copy type).

## Frozen (nothing new is invented)
- **Tool, tiles, maps, search (±16 px), tile 213, lattice 534, depth ±28 µm, level 2, DoG(1,6), reachable rule C0:** exactly as
  in 8ca5789f. Reachable tiles are fixed on unshifted exact-map data before any plant.
- **C1 (amended, as in 13df5696):** plants (4, 0), (0, 8), (12, 0), (0, −12) px; each recovered within 0.5 px on ≥ 95 % of
  reachable tiles, every reachable tile measured. FAIL means that segment's C2-C5 are not read.
- **C2-C5 bars: unchanged from 8ca5789f** (C2 noise p95 ≤ 0.67 px; C3 median exact residual ≤ 2 px and naive ≥ 3× residual;
  C4 |median depth| > 1 coarse layer = real; C5 LOO p95 ≤ 2 px = ONE CORRECTION, else LOCAL FIELD REQUIRED). The scoring code
  is the committed `phase2/analyze.py::ct_summary`.
- **Order of work:** pass 1 (unshifted exact + naive maps; no plants) -> reachable list written with sha256 -> plants
  (`check.py plants`, amended plants only). Output folder `work/w00field/out/amend2/`.

## Fresh C-set (chosen by the pass-1 rule: top valid-fine-vertex count among segments not used in pass 1, surface volumes for both meshes present)
| family | segments (rank by `fine_valid`) | slot |
|---|---|---|
| PHerc0139 | `20260220000000` (2,515,696), `20260306000001` (2,486,640) | next 2 after pass 1's top 3 |
| PHerc0814 | `20250926051122` (3,374,847, coarse is a copy), `20260226121302` (13,742, fine is a copy) | next-ranked of each pass-1 kind |
| PHerc0500P2 | `20250825181922` (587,006), `20250825181907` (576,630) | next 2 after pass 1's top 2 |
| PHerc0009B | `20250919132115` (1,240,139), `20250919135433` (1,195,257) | top 2; both have surface volumes for both meshes; no pass-1 CT |

- **PHerc0841:** all 3 pass-1 segments (w00 + 2) are the whole eligible set (2 M-route segments + w00). No fresh 0841 segment
  exists, so the 0841 CT rows stay exploratory.
- **Stage 0b (pyramid check) fails or has an empty centre block:** that segment is VOID for CT and is reported. Pre-named
  reserves, used only for a Stage-0 failure (never for a CT result): 0139 `20260127000000`, `20260210000000`; 0814 coarse-copy
  `20250925204843`; 0500P2 `20250825181859`; 0009B `20250919125754`.
- **< 5 reachable tiles:** VOID for CT, reported, no replacement.

## Repeat rule (what "repeats" means; frozen now)
A pass-1 family reading is: C3 EXPLAINED (all families); C4 centre-to-centre stands (all families); C5 LOCAL FIELD REQUIRED
(0139) and ONE CORRECTION (0814, 0500P2). A verdict **repeats** for a family only if **every fresh segment of that family with
C1 PASS** gives the same verdict (at least one such segment). A split (one yes, one no) is "not repeated". A family with no
usable fresh segment is "not tested". PHerc0009B has no pass-1 CT: its fresh verdicts are first-and-only readings under bars
frozen before the run; they are reported as such, not as repeats.
Judge's wording rules: "confirmed" and "PASS" are used for no CT verdict that does not repeat. The pass-1 rows are reported next
to the fresh rows as the exploratory discovery set. **C5 on PHerc0139 (pass-1 LOO p95 2.17-2.23 px against the 2 px bar) is flagged
as the most fragile line** and is reported on its own.

## Cost and limits
CPU only, niced, one process at a time. Estimated 6-9 h wall (0139 ~2,500 s per segment in pass 1; 0814 `20250926051122` is large).
Downloads go to the existing cache; meshes are deleted after each segment; `work/cache` stays < 150 GB. Run once; a crash is
resumed, never re-scored with changed settings.

## Tool change (disclosed, before the run)
`check.py plants` accepts this amendment's prereg commit in addition to 13df5696 (one line); new files: `phase2/run_a2.py`
(queue) and `phase2/ct_table_a2.py` (scorer). The mesh route code under `meshpair/` is not touched.
