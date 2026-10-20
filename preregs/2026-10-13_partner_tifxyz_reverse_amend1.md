# Pre-registration amendment 1 to TIFXYZ-REVERSE: the lattice-clean rule, validated only on the unseen S2 pairs (CPU only; public catalog meshes; no CT, no key, no rating)

- **Date:** 2026-10-13. This is an amendment to `prereg/2026-10-13_partner_tifxyz_reverse_route.md` (commit ed1e660c). It is committed and
  pushed BEFORE any S2 mesh is run through either reverse rule.
- **Owner:** [internal lane name removed] (Chris, [machine name removed]). [Internal register reference removed.]
- **Status of S1: unchanged. S1 stays INCONCLUSIVE as filed** (`work/w00field/reverse_route/RESULT.md`).
  - **This rule was designed AFTER seeing the S1 numbers.** It is therefore not a re-score of S1.
  - It is validated only on held-out data: the S2 pairs (section 4). The reverse route has never been run on any S2 mesh, and no reverse-route number from S2 has been looked at.
  - Disclosure: the S2 pairs were part of the earlier forward run of 2026-10-11 (`work/w00field/out/RESULTS.md`).
    - Their forward VOID verdicts are known, and so is the fact that the replica reproduces their meshes (M2a).
    - Nothing else about them is known.

## 1. Why the frozen rule failed (diagnosis, `work/w00field/reverse_route/DIAGNOSIS.md`)
- **Data used:** 6 S1 pairs, their public source meshes, and synthetic data only.
- **Mechanism:** both published meshes are exporter resamples of one source mesh S.
  - S is bilinear per source cell, so it has a crease along each source grid line.
  - A fine cell that lies inside one source cell reproduces the source surface exactly.
  - A fine cell that straddles a source line chords across the crease.
- **Distinguishing checks:**
  - The coarse vertices lie on the source surface itself: |d| p99 <= 0.0005 coarse vx.
  - 99.9-100 % of the off-surface vertices sit in straddling cells. The off rate in the other cells is 0-0.08 %.
  - The exact source chain predicts the measured distance per vertex: residual p99 <= 0.0011 vx, correlation 0.94-1.000.
- **Rivals ruled out, each by its own check:**
  - re-meshing or smoothing: the float32 replica equals the published meshes;
  - an imprecise affine: |Y - A_F(S(s_C))| p99 <= 0.0017 vx;
  - masked cells: the off vertices are as far from invalid cells as all the others;
  - float rounding: ulp = 0.0005 vx, against off values of 0.05-0.49, and the off vertices are phase-locked to the lattice.
- **Synthetic check:** a jittered synthetic source reproduces the failure with no tuning (`synth_a1.py`).
  - Ratio 4 with ~1 vx jitter gives frac on ~0.85.
  - Coarse-copy gives ~0.24.

## 2. The amended pair rule (frozen; code `work/w00field/meshpair/reverse_lc.py`, CLI `check.py reverse-lc`)
The bars are those of the original prereg. Only WHICH vertices are measured changes, plus a float-resolution floor on the in-plane bar.
1. **Inputs:**
   - the fine and coarse meshes and the catalog M, exactly as in the original prereg;
   - the grid shape (rows, cols) of the provenance source mesh `source_key` (`.../intermediate/tifxyz_normalized`). Only its shape is read, from x.tif.
2. **Direction:**
   - **Reverse** (coarse vertex onto the fine surface) unless the coarse mesh is a plain copy of the source.
   - **Forward** when the coarse mesh is a plain copy of the source: fine vertex mapped by M^-1 onto the coarse surface, starting from the
     exporter formula position. In a coarse-copy pair every reverse query sits on a source node, so no reverse query is lattice-clean,
     while every coarse cell is a source cell.
3. **Sampling and projection:**
   - sampling is unchanged: the query mesh's valid mask is eroded by 2 cells, then a 100,000 subsample (seed 0);
   - projection is unchanged: start at the exporter formula position, 12 Gauss-Newton steps, converged = tangential residual < 1e-3 reference
     voxels and all 4 corners valid.
4. **Lattice-clean vertex:** the reference cell [floor(g), floor(g)+1] holding the vertex's FORMULA position contains no source grid line strictly inside, on both axes.
   - Per axis, the cell's source interval is [s(j), s(j+1)] with s(g) = (g + 0.5) n_src / n_ref - 0.5, the exporter relation.
   - A copy mesh has n_ref = n_src, so all its cells are clean.
5. **Per-vertex measures:** as in the original prereg.
   - Normal distance in coarse voxels; "on" = converged and |d| <= 0.05.
   - In-plane residual of the found position against the formula, in coarse canvas px.
6. **Pair verdict** (over the clean vertices):
   - **NOT_EVALUABLE** if fewer than 100 clean vertices.
   - **PASS** iff `frac_on_clean` (on and clean, divided by clean) >= 0.95 AND `inplane_p95_clean` <= bar.
   - **FAIL** otherwise.
   - bar = max(0.003 px, u_coord + u_index), with:
     - u_coord = ulp32(largest |fine coordinate|) x px_f/px_c + ulp32(largest |coarse coordinate|);
     - u_index = (number of the two meshes made by the float32 resampler) x max over axes of ulp32(n_src) x (n_C / n_src) / s_C.
     - Reason, from first principles: the files store float32 coordinates, and the resampler computes the source index in float32.
       - On clean vertices of meshes with fine coordinates above ~8e4 vx, or source grids of ~1700 cells, these alone exceed 0.003 px in synthetic.
       - Where the files allow it, the bar stays 0.003.
7. **Offset against naive:** reported per pair (IN_FILES iff >= 3 px, the existing line). It is computed from a linear fit over the clean on-surface vertices at the
   fine valid-box centre. It is reported, not gated.

## 3. S0-A1: synthetic expected rates (`work/w00field/reverse_route/run_s0_a1.py` -> `s0_a1.json`; run before this commit, synthetic only)
- **Model:** the exporter replica (float32) applied to a curved synthetic source sheet with isotropic per-vertex jitter of 0, 1 or 3 source vx.
  - Four catalog configurations:
    - both transformed, ratio 4 (as 0139 and Paris4);
    - both transformed, ratio 1.95 (as most of 0500P2);
    - fine is copy (as the 0814 and Paris4 copies);
    - coarse is copy (as 1667 and the readable pairs).
  - Coordinate magnitudes cycle up to fine coordinates of ~2.7e5 vx.
  - 100 seeds per cell, plus an extreme-size case: 1700 x 1700 source, fine 6800 x 6800 = 46 M vertices, like the 1667 flatboi, 3 seeds each.
- **Results** (`s0_a1.json`, 603 s):

| configuration | clean runs PASS (jitter 0 / 1 / 3) | max inplane_p95 / bar | max p95 normal (clean, vx) | median share of clean vertices | frozen original rule pass rate (jitter 0 / 1 / 3) |
|---|---|---|---|---|---|
| both, ratio 4 | 100 / 100 / 100 | 0.55 | 0.0024 | 0.59 | 0.25 / 0 / 0 |
| both, ratio 1.95 | 100 / 100 / 100 | 0.53 | 0.0009 | 0.20 | 0.75 / 0 / 0 |
| fine is copy | 100 / 100 / 100 | 0.40 | 0.0005 | 1.00 | 1 / 1 / 1 |
| coarse is copy (forward) | 100 / 100 / 100 | 0.39 | 0.0013 | 1.00 | 0 / 0 / 0 |
| extreme size, both ratio 4 / coarse copy | 3 / 3 each | 0.50 / 0.63 | 0.0018 / 0.0010 | 0.57 / 1.00 | not run |

  - Planted in-plane shifts of 0.05 / 0.25 / 1 coarse px (jitter 1, 100 trials per configuration per size, plus 3 + 3 extreme-size trials at 0.05):
    1206 / 1206 FAIL. The measured median in-plane residual equals the planted shift to 1e-4 px.
  - The bar ranged from 0.003 to 0.016 px. It is above 0.003 only for fine coordinates > ~6.5e4 vx or large grids.
- **Expected rates for S2:**
  - false FAIL / NOT_EVALUABLE on a genuine common-source pair with the exporter convention: 0 / 1206 in synthetic (95 % upper bound ~0.25 %); the largest clean p95 used 63 % of its bar;
  - a planted in-plane convention error of 0.05 / 0.25 / 1 coarse px FAILs in 100 % (1206 / 1206) of the trials;
  - a shift smaller than the float floor (<= ~0.003-0.016 px) is not detectable by construction.
- **Limit:** the synthetic sheet is smoother than real traced meshes in its large-scale shape. The rule's validity does not depend on that
  (clean cells are exact for any source), but real meshes can have folds or self-overlaps that synthetic does not.

## 4. S2 (held out; real; CPU)
- **Set:** the 52 pairs listed in `work/w00field/reverse_route/s2/s2_pairs.json`: every pair whose stored forward verdict is VOID_ROUTE_I in
  PHerc0814 (7, fine is copy), PHercParis4 (44: 42 both transformed + 2 fine is copy) and PHerc1667 (1, the flatboi `20260612121456`,
  coarse is copy, forward direction).
- **Command:**
  `.venv/bin/python check.py reverse-lc --samples PHerc0814,PHercParis4,PHerc1667 --forward-verdict void --prereg <this commit> --out reverse_route/s2/s2.jsonl`
  Then `reverse_route/s2/score_s2.py` (frozen in this commit) writes `s2/summary.json`.
- **Gates:**
  - **CONFIRMED** iff >= 47 / 52 pairs PASS (90 %, the original S1 fraction).
  - **REFUTED** iff >= 6 / 52 pairs FAIL.
  - **INCONCLUSIVE** otherwise (e.g. many NOT_EVALUABLE).
  - Every FAIL and NOT_EVALUABLE pair is listed with its numbers.
- **What CONFIRMED supports:** the lattice-clean rule is validated on held-out pairs, and on those 52 pairs the published fine and coarse meshes relate
  by the exporter formula (a common source plus the half-pixel/rounded-size convention).
  - It does NOT change S1's filed verdict.
  - Applying the validated rule to S1 and to the readable pairs would be an application, labelled as such. That is a decision for the Judge / a human, and it is not part of this amendment.
- **Not allowed after this commit:** changing any number in section 2, the S2 list or the scorer; re-running pairs; choosing seeds. A crash on a pair is reported as ERROR
  (counts against PASS); one retry is allowed only for a network/download error, and it is logged.

## 5. Ops
- CPU only, niced, one process. No GPU, [internal GPU-sharing note removed].
- Downloads go to `work/w00field/cache/meshes/` (git-ignored, the tool's cache), about 12-13 GB (Paris4 ~11.7 GB, the flatboi ~0.6 GB). Caches stay < 150 GB in total
  (`work/cache` 77 GB + `work/w00field/cache` 3.9 GB before the run).
- Outputs go to `work/w00field/reverse_route/s2/`: `s2.jsonl`, `s2.log`, `summary.json`, `RESULT.md`.
