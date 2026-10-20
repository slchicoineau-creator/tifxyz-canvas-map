# TIFXYZ-REVERSE: why 13-15 % of coarse vertices sit off the fine surface (step-1 diagnosis, post-hoc, 2026-10-13)

Post-hoc and descriptive: the S1 numbers were seen first (`RESULT.md`), and S1 stays **INCONCLUSIVE as filed**.
Data used: 6 S1 pairs (0139 x2, 0500P2 x4, chosen to cover the worst, the best and the smallest S1 cases) and their public
source meshes (`.../mesh/intermediate/tifxyz_normalized`, the provenance source of both published meshes), plus synthetic data.
No S2 mesh was opened. Script `diagnose.py`, records `diag/diagnose_s1.jsonl` (16 s CPU).

## Verdict: cause found. It is the source lattice (H-LAT). The four rivals are each ruled out by a separate check.
Both published meshes are exporter resamples of one source mesh S. S is bilinear per source cell, so its surface has a
crease (kink) along every source grid line. The fine mesh is a 1.95x or 4x finer float32 resample, and at ratio 4 its vertices sit at
source index -0.375 + k/4, so no fine vertex lies on a source line. A fine cell that lies wholly inside one source cell
reproduces the source surface exactly. A fine cell that straddles a source grid line chords across the crease. A coarse vertex (a
bilinear sample of S) that falls in such a cell is off the fine chord by the crease sagitta. That sagitta grows with
the source mesh's own jitter, not with any convention error.

| check (6 pairs) | what each explanation predicts | observed (range over pairs) |
|---|---|---|
| (i) project the coarse vertices onto the SOURCE surface A_F(S) (no resample, so no lattice) | H-LAT: on it to float noise | \|d\| p99 0.00024-0.00050 coarse vx; 100 % within 0.05 |
| (ii) off-surface (\|d\| > 0.05) vertices vs straddling cells | H-LAT: only in straddling cells | 99.9-100 % of off vertices are in straddling cells; off rate in non-straddling cells 0-0.08 %; non-straddling \|d\| p99 0.00018-0.0008 |
| (iii) exact chain prediction e = A_F(bilerp_S(s_C)) - bilerp_F(formula position), normal part, vs measured d | H-LAT: predicts d per vertex | residual p99 0.00027-0.0011 coarse vx; corr 0.94-1.000 |
| H-SMO (fine/coarse re-meshed or smoothed after the transform) | replica(S) differs from the published mesh | float32 replica = published: fine p99 0.0022-0.0028 vx, coarse p99 <= 0.001 vx, masks identical, so ruled out |
| H-AFF (catalog affine imprecise) | smooth field \|Y - A_F(bilerp_S(s_C))\| far above float | p99 0.0005-0.0017 coarse vx, so ruled out |
| H-MSK (masked-cell interpolation) | off vertices next to invalid fine corners | median distance to invalid: off 13-230 cells vs all 16-229 (same), so ruled out |
| H-FLT (float rounding of stored coords) | \|d\| of a few ulp, no phase dependence | ulp = 0.0005 coarse vx vs off values of 0.05-0.49, all phase-locked to straddling cells, so ruled out |

Per pair (coarse vx; "straddle" = the found fine cell contains a source grid line):

| pair | fine/source ratio | frac off | share of off in straddle | off rate non-straddle | \|d\| p99 non-straddle | \|d\| p99 straddle | in-plane p95 non-straddle / straddle (px) |
|---|---|---|---|---|---|---|---|
| 0500P2 20250611171318 | 1.95 | 0.19 | 0.999 | 0.0008 | 0.00034 | 0.31 | 0.0009 / 0.27 |
| 0500P2 20250716055231 | 1.95 | 0.075 | 0.999 | 0.0002 | 0.00030 | 0.15 | 0.0005 / 0.023 |
| 0500P2 20250825181917 | 3.90 | 0.071 | 1 | 0 | 0.00034 | 0.18 | 0.0005 / 0.009 |
| 0500P2 20250919184429 | 1.95 | 0.37 | 1 | 0 | 0.00024 | 0.49 | 0.0003 / 0.050 |
| 0139 20260311000000 | 4.00 | 0.16 | 1 | 0 | 0.00045 | 0.36 | 0.0006 / 0.026 |
| 0139 20260206000002 | 4.00 | 0.18 | 1 | 0 | 0.00052 | 0.45 | 0.0006 / 0.042 |

Supporting observations:
- The 7 0500P2 pairs with ratio 3.9 are exactly the S1 pairs with the smallest in-plane p95 (0.005). At ratio 1.95, the
  straddling share is 78 % (vs 44 % at ratio 4), so the in-plane tail is larger.
- A model built only from the fine mesh (second differences D_j reveal each crease: chord error = -min(t|D_j|, (1-t)|D_{j+1}|))
  predicts d at ratio 4 (corr 0.9998). At ratio 1.95 it predicts less well (corr 0.76-0.92), because two creases fall
  within 2 fine cells. So a per-vertex sagitta tolerance taken from the fine mesh alone is not exact at ratio ~2. The amendment
  does not use one; it uses the source grid shape instead (below).
- The readable pairs (coarse = source copy) fail the frozen rule for the same reason, and more often. Each coarse vertex IS a source node, and
  a node sits at the crossing of two creases. So every reverse query falls in a cell that straddles on both axes (frac on 0.24-0.43).
- Synthetic reproduction (`synth_a1.py`, jittered source, exporter replica). The frozen rule passes when the source has no jitter, which
  is what S0 tested. It fails once the source has ~1 vx jitter: ratio 4 gives frac on 0.85, coarse-copy gives 0.24. The model was not tuned to S1.

## Consequence for a rule (used in amendment 1)
Inside a fine cell that holds no source grid line, theory says the fine surface IS the source surface. There the frozen bars (0.05 vx, 0.95, 0.003 px)
apply as designed. So the rule keeps the bars and changes WHICH vertices are measured: only coarse vertices whose formula position
falls in a lattice-clean cell. The source grid shape comes from the provenance source mesh. When the coarse mesh is the source
copy, the rule runs in the forward direction instead (fine vertex onto the coarse surface, from the formula start), because there
every coarse cell is a source cell. One more term is added from first principles: the in-plane bar is max(0.003 px, the float32 resolution
of the files). That resolution is one ulp of the largest stored coordinate of each mesh, plus one ulp of the largest float32 source
index per resample. It is needed because synthetic meshes with fine coordinates above ~8e4 vx, or with source grids of ~1700 cells, exceed
0.003 px on clean vertices alone. Details and expected rates: `prereg/2026-10-13_partner_tifxyz_reverse_amend1.md`, `s0_a1.json`.
