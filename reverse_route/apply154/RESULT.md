# TIFXYZ-REVERSE amendment 1: APPLICATION to the other 154 catalog pairs

Prereg `prereg/2026-10-13_partner_tifxyz_reverse_apply154.md` (commit b8930ac6, pushed before the run). Rule, code and bars frozen at 8a4c6e6c
(`reverse_lc.py` sha1 1461ffe0..., `check.py` sha1 070d3ed2...). CPU only, public meshes, 154 pairs selected as registered, 0 errors, 0 retries, about 6 min.
**This is an application, not a test.** The rule was designed on S1 data and the readable pairs. These results carry no independent evidence weight.
**S1's original verdict stays INCONCLUSIVE as filed.** Forward-route VOID verdicts stay as filed. Files: `apply154.jsonl` (raw), `pairs_table.jsonl` (206 rows), `summary.json`, `score_apply154.py`.

## Full-catalog table (206 pairs)
| set | sample | n | PASS | FAIL | NOT_EVAL | PASS at plain 0.003 px | IN_FILES (offset >= 3 px) | clean vertices | off-surface clean vertices |
|---|---|---|---|---|---|---|---|---|---|
| S1 (reverse) | PHerc0139 | 38 | 38 | 0 | 0 | 38 | 35 | 2,059,235 | 0 |
| S1 (reverse) | PHerc0500P2 | 46 | 43 | 0 | 3 | 43 | 35 of 43 | 224,044 | 7 (3 pairs) |
| readable (forward) | PHerc0009B | 18 | 18 | 0 | 0 | 18 | 18 | 1,800,000 | 0 |
| readable (forward) | PHerc0343P | 8 | 8 | 0 | 0 | 8 | 8 | 800,000 | 0 |
| readable (forward) | PHerc0814 | 12 | 12 | 0 | 0 | 11 | 12 | 1,200,000 | 1 (1 pair) |
| readable (forward) | PHerc0841 | 2 | 2 | 0 | 0 | 2 | 2 | 200,000 | 0 |
| readable (forward) | PHerc1667 | 19 | 19 | 0 | 0 | 19 | 19 | 1,900,000 | 0 |
| readable (forward) | PHercParis4 | 11 | 11 | 0 | 0 | 11 | 11 | 1,100,000 | 0 |
| **application subtotal** | | **154** | **151** | **0** | **3** | **150** | **140 of 151** | 9,283,279 | 8 (4 pairs) |
| S2 (held out, 8a4c6e6c) | 0814 7, Paris4 44, 1667 1 | 52 | 52 | 0 | 0 | 50 | 52 | 2,675,306 | 1 (1 pair) |
| **all 206** | | **206** | **203** | **0** | **3** | **200** | **192 of 203** | 11,958,585 | 9 (5 pairs) |

- Direction: S1 84 reverse; readable 70 forward (coarse is a plain copy). The 0841 development pair 20260220213127 (w00) is not among the 206; the 2 PHerc0841 pairs are 20260220214732 and 20260221022814 (release correction, see CHANGES.md).
- **FAIL pairs: none.**
- **NOT_EVALUABLE (3, as expected in the prereg):** PHerc0500P2 20250919184428 (55 clean), 20250919184429 (24 clean), 20250919184430 (23 clean). Their clean vertices that were measured were all on the surface, but the rule needs >= 100.
- Smallest evaluable 0500P2 PASS has 128 clean vertices. Thin power there: 0500P2 as a whole contributes only 224,044 clean vertices over 43 pairs.
- Plain-0.003 count: 150 of 151 application PASS pairs would also pass at the plain bar. The one exception is PHerc0814 20250928235953 (in-plane p95 0.00304 px, floor bar 0.0059). Largest p95/bar ratio over all application PASS pairs 0.62.
- Off-surface clean vertices (on = |d| <= 0.05 coarse vx; counts of clean vertices not on): 0500P2 20250609204953 1 of 2,328; 20250611171318 2 of 2,619; 20250716055230 4 of 5,809; 0814 20250929040153 1 of 100,000; S2 Paris4 20260603145540 1 of 56,150. Lowest frac_on_clean in the application is 0.99924.
- Offsets: all 151 evaluable pairs have an offset against naive pairing in the files. 140 are >= 3 px (IN_FILES). 11 PASS pairs are below the line (0.92-2.9 px; 8 PHerc0500P2, 3 PHerc0139). Reported, not gated.

## Thin spots and limits (state them wherever the claim is written)
- Held-out validation of ratio-1.95 geometry (most of 0500P2) is zero; the forward (coarse-copy) direction has one held-out pair (the 1667 flatboi, passing through the float floor). The selection itself is validated on 42 Paris4 both-transformed pairs at ratio 4 only. The 9 fine-copy pairs of S2 also pass the original rule, so they do not test the lattice selection.
- Blind spot C3: an error phase-locked to the straddling cells is invisible. Tolerance C2: defects on under ~5 % of the clean vertices are invisible to the 0.95 bar; the off counts above are the measured values.
- The power is against a wrong convention, wrong affine, or a different/edited source. It says nothing about geometry inside straddling cells. The diagnosis used 6 hand-chosen S1 pairs and is post hoc.
- Not acceptable wording: "confirmed on all 206 pairs"; "S1 now passes"; "no local defects" without the 5 % tolerance and the counts.

## Summary wording (the form accepted by an independent review, a separate AI session that did not do the work, filled with these numbers)
> Mesh-only check of the published 2.4/9.4 um (and 7.91 um) mesh pairs: on lattice-clean vertices the fine and coarse meshes agree with the exporter formula (common source + half-pixel / rounded-size convention) to float32 resolution. Rule validated on 52 held-out pairs (52/52; 50/52 at the plain 0.003 px bar, the other 2 explained by the exporter's float32 index), then applied to the other 154 pairs: 151 PASS, 3 NOT_EVALUABLE (too few clean vertices; PHerc0500P2), 0 FAIL. Offsets vs naive index pairing are in the files on every evaluable pair.

Caveat line to keep beside it: lattice selection validated on ratio-4 pairs only; ratio-1.95 (0500P2) and the forward direction are applications with no or one held-out pair; 5 of 206 pairs have a total of 9 off-surface clean vertices; local defects under ~5 % of clean vertices are invisible to the verdict.
