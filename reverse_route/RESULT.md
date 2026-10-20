# TIFXYZ-REVERSE result (prereg `prereg/2026-10-13_partner_tifxyz_reverse_route.md`, commit ed1e660c)

CPU only, public meshes. Code: `meshpair/reverse.py`, `check.py reverse`, `tests/test_reverse.py` (39 tests pass; the 19 old ones untouched).
Forward-route results and old VOID verdicts are unchanged.

| gate | result | verdict |
|---|---|---|
| S0 synthetic (100 seeds) | clean fine-is-source: max p95 normal 0.00038 coarse vx (bar 0.01). Planted 0.05 / 0.25 / 1 px in-plane shifts FAIL the pair rule in 100 % / 100 % / 100 % of 100 trials (bar 95 %) | PASS |
| S1 0139 (38) + 0500P2 (46) | pair rule (>= 95 % of vertices within 0.05 coarse vx AND in-plane p95 <= 0.003 px) holds on **1 / 84** pairs (0139 0/38, 0500P2 1/46); bar >= 76 | PASS bar missed |
| S1 KILL line (0139 median reverse distance >= 0.1 vx) | median over pairs of the per-pair median = **0.00027** coarse vx | not triggered |
| => S1 | between the two lines | **INCONCLUSIVE** (as the prereg defines it) |
| S1-readable (70 IN_FILES pairs) | **70 / 70** keep IN_FILES under the reverse route (offset >= 3 px) | gate met |
| S2 (0814, Paris4 VOID, 1667) | not run (S1 not PASS) | not run |

Numbers (S1, per-pair medians over pairs): 0139 frac on-surface 0.855 (min 0.823), in-plane p95 0.0114 px, pair p95 normal 0.118 vx.
0500P2 frac 0.909 (min 0.59), in-plane p95 0.021 px, p95 normal 0.071 vx. Only 1 pair of 84 has in-plane p95 <= 0.003, and 1 has frac >= 0.95.

## What it says
- The hypothesis holds at the centre: the median coarse vertex lies 0.0003 coarse vx (0139) / 0.003 (0500P2) from the fine surface, i.e. float noise,
  where the forward route had VOIDed (it started from the wrong side and compared a 3.9x finer lattice to the coarse chords). So the "common source" premise is supported for
  the typical vertex on 0139/0500P2.
- The frozen pair rule is not met, because 13-15 % of vertices (0139) sit 0.05-0.3 coarse vx off the fine bilinear surface and the in-plane p95 is 0.01-0.02 px.
  Post-hoc description (`diag.py`, not a gate, 6 pairs): these off-surface vertices are interior (2.5 % within 3 cells of the mask edge vs 3.9 % of all),
  with a smooth tail (q80 0.03, q90 0.07, q99 0.26 vx). That looks like a lattice effect (a fine bilinear cell that straddles two source cells is not the
  source surface between its corners), not a convention error. This is a hypothesis, untested.
- The pair rule is also not calibrated on the readable pairs: 0 / 70 pass it (frac on-surface 0.24-0.43), although all 70 keep IN_FILES. There the coarse mesh is the source and its
  vertices are the corners of the coarse cells that the 3.9x fine lattice chords across. So the 0.05 vx line is tighter than the lattice allows on real meshes, while S0 (smooth synthetic sheets) could not show that.
- Headline target (>= 195 / 206 verdicts): not reached by this rule. By the frozen bars a "whole catalog" claim does not follow from the reverse route; the 134 VOIDs stay VOID.

## Decision for a human (not taken here)
A new prereg would need a lattice-aware bar, for example a normal-distance bound from the fine cell sagitta computed per vertex, or a bar on the median and the
in-plane median instead of p95 tails. It must be fixed before looking at more of these numbers; the S1 numbers above are now seen, so any such bar is post-hoc
and the verdict stays an amendment, not a pass. Not done.

Files: `s0.json`, `s1.jsonl`, `s1_summary.json`, `s1_readable.jsonl`, logs/times, `run_s0.py`, `summarize.py`, `diag.py`.

## Amendment 1 (added 2026-10-13; post-hoc design, held-out validation)
- **Diagnosis:** `DIAGNOSIS.md`. The off-surface vertices come from source-lattice creases (fine cells that straddle a source grid line). Smoothing, affine,
  mask and float explanations were ruled out.
- **Rule:** `prereg/2026-10-13_partner_tifxyz_reverse_amend1.md` (8a4c6e6c). Same bars, lattice-clean vertices only, plus a float32 floor.
  - Synthetic: 0 / 1206 false fail, 1206 / 1206 planted shifts caught.
- **S2 (held out):** 52 / 52 PASS, which is **CONFIRMED** (`s2/RESULT.md`). The float floor decided 2 pairs; their excess is the organizers' float32 resampler index.
- **The S1 verdict above is unchanged: INCONCLUSIVE.**
