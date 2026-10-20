# Pre-registration: TIFXYZ-REVERSE. Does the reverse route (coarse vertex -> fine surface) turn the 84 VOID pairs of PHerc0139 and PHerc0500P2 into verdicts, and does it keep the 70 readable verdicts? (CPU only; public catalog meshes only; no CT volume, no key, no rating.)

- **Date:** 2026-10-13. Committed and pushed BEFORE any real-mesh run of the new code. The only things run before this
  commit: the existing self-tests (19 passed, `work/w00field/tests`, 10 s) and HEAD requests for file sizes.
- **Owner:** [internal lane name removed] (Chris, [machine name removed]). [Internal register and brainstorm-note references removed.] The vet reports the rule "fine is a copy, or both
  transformed" as the VOID cause (`work/w00field/out/RESULTS.md` l.71-79) and the amendment "invert through the source mesh"
  (l.158-159) as never done.
- **Data:** public tifxyz meshes (x/y/z.tif + meta.json) from `s3://vesuvius-challenge-open-data/` through the code's own
  fetch path (`meshpair.tifxyz.fetch_tifxyz`, listed in `work/w00field/out/pairs.json`, catalog sha `7ec8dcaa...`). Cache under
  `work/w00field/cache/meshes/` (git-ignored; measured size of S1 meshes ~1.3 GB, the 70 readable ~4 GB; total cache < 150 GB).
  No surface volume, no CT. [Internal data name removed.]
- **Compute:** CPU only. No GPU, [internal GPU-sharing note removed].

## 1. Question and hypothesis

The forward route I (fine vertex -> coarse surface) was VOID on 0139 (38/38) and 0500P2 (46/46), and on 0814 7/7, Paris4 44/44.
Hypothesis H: on these pairs both meshes are (re)samples of the same source surface. Then the catalog's coarse vertices lie on
the fine bilinear surface to float noise (the fine surface is a 3.9x finer resample, so the chord error of the fine cells is the
only expected deviation), and the position along the sheet follows the exporter formula
p_C = (p_F + 10.5) * rho - 10.5 (rho = n_C / n_F per axis). The reverse route measures this directly.

## 2. The reverse route (frozen definition)

Per pair (fine mesh F in the fine volume, coarse mesh C in the coarse volume, catalog affine M = coarse-volume -> fine-volume
exactly as `coarse_to_fine_affine` gives it, no change):
1. Coarse vertices: all valid coarse vertices after a 2-cell binary erosion of the valid mask (same erosion as route I). If more than
   100,000, a seeded subsample (numpy `default_rng(0)`).
2. Y = M applied to the coarse vertex (fine-volume voxels).
3. Start: the fine canvas position from the formula, `exact_canvas_map_inv` (p_F = (p_C + 10.5)/rho - 10.5), converted to a fine
   grid position. Then 12 Gauss-Newton steps on the fine bilinear patch (same step rule and clip as route I `_gn_project`; own
   bilinear that gathers the 4 corners, so no full-array float64 copy). No restart cells (a start that does not converge counts as a failure).
4. Measures per vertex:
   - (a) normal distance to the fine surface, in COARSE voxels (fine-voxel normal distance x px_fine / px_coarse);
   - converged = tangential residual < 1e-3 fine vx and all 4 corners valid; a vertex is "on-surface" iff converged and |(a)| <= 0.05 coarse vx;
   - (b) in-plane residual vs the formula: the found fine grid position -> fine canvas -> exact map -> coarse canvas, minus the vertex's
     own coarse canvas position; norm in coarse px (over converged vertices).
5. Per-pair: `frac_on` = on-surface / sampled; `inplane_p95` = p95 of (b); `median_normal` = median |(a)| (converged).
   **Pair PASS iff frac_on >= 0.95 AND inplane_p95 <= 0.003 coarse px.**
6. Offset vs naive (for the readable pairs): linear fit p_C = a + b p_F over the on-surface vertices (all converged if < 100), offset at the centre of the fine
   valid box vs the naive map k * p (k = px_fine/px_coarse), norm in coarse px. IN_FILES iff >= 3.0 px (the existing line).

Code: new module `work/w00field/meshpair/reverse.py` + `check.py reverse` subcommand. The forward route (`fit.py`, `meshcheck.py`)
and the stored forward results are not changed. The old VOID verdicts stay as filed in `out/RESULTS.md` and `out/mesh_records.jsonl`.

## 3. Gates (as the vet states them)

- **S0 (synthetic, `meshpair/synth.py` `curved_sheet` + `planted_resample`, new self-tests):** a fine-is-source pair (source = fine
  sheet, spacing 20 vx; coarse = `transform_tifxyz(source, A_C)` with a ~1/3.9 affine) in 100 seeds (0-99).
  - Reverse normal distance p95 (pooled over the pair's vertices) must be <= 0.01 coarse vx on every seed.
  - Planted in-plane shifts of 0.05, 0.25 and 1 coarse px (coarse vertices sampled at a shifted source index along y, x or both, 100 trials each, seeds
    0-99) must FAIL the pair rule in >= 95 % of the trials per shift size.
  - Also reported, not gated: a both-transformed synthetic (source -> A_F and A_C), the same measures.
  - S0 must pass before any real mesh is read. If it fails, stop and report; the real data is not touched.
- **S1 (real, CPU):** PHerc0139 (38 pairs) + PHerc0500P2 (46 pairs), the primary pair of each segment as in `out/pairs.json`.
  - PASS iff on >= 90 % of the 84 pairs (>= 76) the pair rule holds (frac_on >= 0.95 and inplane_p95 <= 0.003).
  - KILL iff the median over 0139 pairs of `median_normal` is >= 0.1 coarse vx.
  - Between the two (neither PASS nor KILL): INCONCLUSIVE, reported with the distribution.
  - **S1-readable:** the 70 pairs whose forward M1 verdict is IN_FILES (list from `out/mesh_records.jsonl`: 0841 2 + the 0841 dev pair is excluded, 0009B 18, 0343P 8, 0814 12, Paris4 11, 1667 19 = 70 after
    the count is read from the file; the count read is the one that stands). Each must give IN_FILES (offset >= 3 px) under the reverse route: 70/70
    required ("keep their verdicts"). Their pair-rule pass rate is reported, not gated.
- **S2 (real; only if S1 PASS and the S1-readable gate passes):** the 7 PHerc0814 fine-is-copy pairs, the 44 Paris4 VOID pairs (~11 GB) and
  the PHerc1667 pairs seeded from the formula (the flatboi segment that route I could not start). Target: the headline >= 195/206 pair-rule PASS
  verdicts over all 206 catalog pairs (the 84 + 70 + these). Reported honestly if lower.

Not allowed: changing 0.05, 0.003, 0.95, 0.9, 0.1, 3.0 or the sampling after seeing S1; picking pairs; re-running with other seeds
to rescue a miss. A post-hoc idea goes in a separate amendment commit, labelled as such.

## 4. What stays as filed
Forward route I verdicts (VOID on 134 pairs counting 0814/Paris4), M2a/M2b, the CT stage, `out/RESULTS.md`: unchanged. The reverse route
adds a column; it does not rewrite a verdict. The #1912 drafts stay human-posted.

## 5. Seed and outputs
Seed 0 for vertex subsampling; synthetic seeds 0-99. Results: `work/w00field/reverse_route/` (`s0.json`, `s1.jsonl`, `s1_summary.json`,
`RESULT.md`).
