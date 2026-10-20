# TIFXYZ-REVERSE amendment 1: S2 result (held-out pairs; prereg `prereg/2026-10-13_partner_tifxyz_reverse_amend1.md`, commit 8a4c6e6c)

CPU only, public meshes. Command as frozen: `check.py reverse-lc --samples PHerc0814,PHercParis4,PHerc1667 --forward-verdict void --prereg 8a4c6e6c`.
Total run time 22 min. No errors, no retries. Scorer `score_s2.py` writes `summary.json`. Records are in `s2.jsonl`, the log in `s2.log`.
S1 stays **INCONCLUSIVE as filed**. This rule was designed after S1 was seen; this run is its validation on held-out data.

| gate | result | verdict |
|---|---|---|
| S2: 52 pairs (0814 7, Paris4 44, 1667 flatboi 1), lattice-clean rule | **52 / 52 PASS**, 0 FAIL, 0 NOT_EVALUABLE (bar: >= 47 PASS) | **CONFIRMED** |

| group | n | clean share | frac_on_clean (min) | p95 normal, clean (max, coarse vx) | in-plane p95, clean (px) | bar (px) | offset vs naive (px) | IN_FILES |
|---|---|---|---|---|---|---|---|---|
| fine is copy (0814 7 + Paris4 2), reverse | 9 | 1.00 | 1.000 | 0.0003 | 0.00007-0.0012 | 0.003-0.0048 | 5.1-134.8 | 9/9 |
| both transformed (Paris4), reverse | 42 | 0.55-0.59 | 1.000 | 0.0005 | 0.0010-0.0033 | 0.0063-0.0153 | 119-387 | 42/42 |
| coarse is copy (1667 flatboi, 46.8 M fine vertices), forward | 1 | 1.00 | 1.000 | 0.0003 | 0.0112 | 0.0217 | 753 (~6.0 mm) | 1/1 |

## Read it with these caveats
- **The float floor decided 2 of 52 pairs.** With the plain 0.003 px bar, 50 / 52 would pass. The two exceptions are:
  - the 1667 flatboi, in-plane p95 0.0112 px. It has a 9204-column source grid, so float32 index resolution is 0.0195 px;
  - Paris4 `20260623153216`, in-plane p95 0.00326 px.

  A post-verdict descriptive check (`post_float_check.py` -> `post_float_check.json`, not a gate) predicts each vertex's residual from the exporter's
  float32 source index minus the float64 formula, over the clean vertices:
  - corr along x 0.976 / 0.977;
  - p95 after removing the predicted float32 index error 0.00057 / 0.00092 px.

  So on these two the residual above 0.003 px is the organizers' float32 resampler. It is not a convention mismatch.
- **Fine-copy pairs.** On the 9 fine-copy pairs every vertex is lattice-clean, so the frozen original rule also holds there (frac on all vertices 1.000,
  in-plane p95 <= 0.0012). On the 42 both-transformed pairs the original rule would fail as on 0139 (frac on all vertices 0.82-0.83). That is the lattice effect the
  diagnosis predicted, now seen on held-out data.
- **Flatboi offset.** The flatboi's ~6 mm offset matches the exporter-model value in `out/RESULTS.md` ("an offset of ~6 mm against naive"). The forward route could not
  start there from the naive guess; the formula start reaches it.
- **Scope.** CONFIRMED covers these 52 pairs only. It says the published fine and coarse meshes relate by the exporter formula, i.e. a common source plus the
  half-pixel / rounded-size convention, and that the offset is IN_FILES on 52 / 52. It does not re-score S1 (84 pairs, still INCONCLUSIVE as filed).
  Running the now-validated rule on S1 and on the 70 readable pairs would be an application, not a test. Whether to do that, and whether a summary may then
  claim the whole catalog, was left to an independent review (a separate AI session that did not do the work; it ruled: application, labelled; see `../apply154/RESULT.md`).
