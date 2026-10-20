# Pre-registration: APPLICATION of the frozen TIFXYZ-REVERSE amendment-1 rule to the other 154 catalog pairs (CPU only; public meshes; no CT, no key, no rating)

- **Date:** 2026-10-13. Owner: [internal lane name removed] (Chris, [machine name removed]). [Internal register reference removed.]
- **Basis:** judge review `reviews/2026-10-13_partner_judge_tifxyz_reverse_amend1.md` (PASS WITH CONDITIONS, section 6 conditions 1-5). This file meets condition 1.
- **This is an APPLICATION, not a test.**
  - The rule was designed after S1 was seen (84 pair numbers seen; 6 pairs used for the diagnosis). The coarse-copy forward direction was added after the readable pairs failed the original rule.
  - The results have no independent evidence weight. S1's original verdict stays **INCONCLUSIVE as filed**. The run is never called a re-score or a test of S1.
  - The forward-route VOID verdicts stay as filed. The #1912 drafts stay human-posted.

## 1. Rule: frozen, no changes
- The rule, bars, projection, sampling and seeds are exactly those of `prereg/2026-10-13_partner_tifxyz_reverse_amend1.md`, section 2 (commit **8a4c6e6c**, pushed 07:30:08 -0500).
- Code sha1 at this registration (unchanged since 8a4c6e6c; `git diff 8a4c6e6c HEAD` touches no code):
  - `work/w00field/meshpair/reverse_lc.py` 1461ffe0b6b644db6d8d65456ceb3c1495bb0623
  - `work/w00field/check.py` 070d3ed2e179a29fab4b026055a6269e42ace487
- Not allowed: changing any number, bar, direction rule or seed; re-running pairs; rescue analyses. A crash is reported as ERROR. One retry only for a network/download error, logged.

## 2. The 154 pairs, defined by rule
- Every catalog pair that is not in S2: the S1 set (84: PHerc0139 38, PHerc0500P2 46) and the S1-readable set (70: 0814 12, 0343P 8, 0841 2, 1667 19, Paris4 11, 0009B 18).
- With the 52 S2 pairs these are all 206 catalog pairs. The list is in `work/w00field/reverse_route/apply154/apply154_pairs.json` (no overlap with `s2/s2_pairs.json`, checked).
- Includes the dev pair PHerc0841 20260220213127 (seen before; the existing guard treats it as non-confirmatory). It is listed, not hidden.

## 3. Command and outputs
- `cd work/w00field && nice .venv/bin/python check.py reverse-lc --segments <segments.txt> --prereg <this commit> --out reverse_route/apply154/apply154.jsonl`
  (the segment list is `apply154/segments.txt`; the run must report exactly 154 pairs selected, or it is stopped).
- Scorer (frozen with this commit): `reverse_route/apply154/score_apply154.py` -> `summary.json`, `pairs_table.jsonl`; then `RESULT.md` (full-catalog table 84 + 70 + 52 = 206).
- Per pair: PASS / FAIL / NOT_EVALUABLE (+ ERROR), offset vs naive and IN_FILES class (>= 3 px), as in amendment section 2.6-2.7.
- Also reported (judge condition 1): the count of pairs that pass at the plain 0.003 px bar next to the floor-bar count; the off-surface clean-vertex counts per pair (n_clean minus n_clean_on) and their total.

## 4. Expectations (stated before running)
- Most S1 0500P2 pairs have ratio 1.95 and ~20 % clean vertices. Three of them have only 55, 189 and 362 eroded coarse vertices, so about 3 NOT_EVALUABLE (< 100 clean) are expected. Coverage is therefore not 206/206 in advance.
- Preview on 4 cached pairs (judge): PASS on 0139, 0500P2, 0009B and 1667.
- No pass threshold and no verdict gate: nothing here is a hypothesis test. Every FAIL and NOT_EVALUABLE pair is listed with its numbers. A FAIL is reported as a FAIL and not explained away in the headline.

## 5. How the result may be worded ([internal register reference removed] claim; judge section 6 item 4)
> Mesh-only check of the published 2.4/9.4 um (and 7.91 um) mesh pairs: on lattice-clean vertices the fine and coarse meshes agree with the exporter formula (common source + half-pixel / rounded-size convention) to float32 resolution. Rule validated on 52 held-out pairs (52/52; 50/52 at the plain 0.003 px bar, the other 2 explained by the exporter's float32 index), then applied to the other 154 pairs: N PASS, N NOT_EVALUABLE (too few clean vertices), N FAIL. Offsets vs naive index pairing are in the files on every evaluable pair.

- N are filled from the run. **Not acceptable:** "confirmed on all 206 pairs"; "S1 now passes"; any statement about geometry inside straddling cells; "no local defects" without the 5 % tolerance and the measured off counts.

## 6. Thin spots (to be stated wherever the claim is written)
- The lattice selection is validated on 42 both-transformed pairs, all Paris4, all at fine/source ratio 4. Held-out validation of ratio-1.95 geometry (most of 0500P2) is **zero**.
- The forward (coarse-copy) direction is validated on **one** held-out pair (the 1667 flatboi), which passes only through the float floor. The 9 fine-copy pairs of S2 are 100 % clean and also pass the original rule, so they do not test the selection.
- Blind spot C3: an error phase-locked to the straddling cells is invisible to the rule (only the lattice effect itself has that structure). Tolerance C2: local defects on less than ~5 % of the clean vertices are invisible to the 0.95 bar. The off counts are therefore quoted, not only the bar.
- The power is against a wrong convention, a wrong affine, or a different/edited source. It says nothing about geometry inside straddling cells.
- The diagnosis used 6 hand-chosen S1 pairs and is post hoc.

## 7. Ops
- CPU only, niced, one process; no GPU, [internal GPU-sharing note removed]. Meshes are cached in `work/w00field/cache/meshes/` (15 GB now); the total cache stays < 150 GB.
