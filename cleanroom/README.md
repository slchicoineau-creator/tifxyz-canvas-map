# Clean-room rewrite of the VC3D-derived modules (card 02)

Why: the prize needs a permissive licence (MIT). `meshpair/conventions.py`, `meshpair/exporter.py`, the vc_load rule in
`meshpair/tifxyz.py` and the slice centre in `meshpair/ct.py` re-implement conventions of GPLv3 volume-cartographer.
Chris's condition (review correction `bb7d203fe`, 2026-10-18): rewrite them, get **bit-identical** results to his
`11b478606` on all 207 pairs before release, and Chris stays named as author.

## Two steps, two people/sessions

| Step | Who | When | Output |
|---|---|---|---|
| 1. Specification ("dirty room") | Claude Code (Opus 5.5), Sebastien's session, approved by Sebastien | 2026-10-18 | `SPEC.md` (this folder). The writer read the GPL modules and turned them into facts, formulas and test vectors, not code. |
| 2. Implementation ("clean room") | a separate implementer (a fresh Claude session or a person) that has **never** seen the GPL sources | after step 1 | new `conventions.py`, `exporter.py`, `tifxyz.py`, and a slice-centre helper for `ct.py`, MIT-licensed |
| 3. Integration + parity | someone other than the implementer (may see the old code) | after step 2 | swap the modules in, run the tests and the parity checks below |

## Rules for the implementer (step 2)

You may read ONLY:
1. `cleanroom/SPEC.md`;
2. the public test files `tests/test_meshpair.py`, `tests/test_release.py`, `tests/test_reverse.py`,
   `tests/test_reverse_lc.py`, `tests/test_plants_a1.py`;
3. MIT-licensed villa code (github.com/ScrollPrize/villa), **excluding** `volume-cartographer/` (GPLv3) in full;
4. data files the tests read: `out/mesh_records.jsonl`, `out/dev_w00_mesh.jsonl`, `results/pair_maps.tsv`;
5. public documentation of Python, NumPy, SciPy and tifffile.

You must NOT open (not even to "check one thing"):
- volume-cartographer sources (any copy, any commit);
- this package's `meshpair/conventions.py`, `meshpair/exporter.py`, `meshpair/tifxyz.py`, `meshpair/ct.py`;
- files that restate the same conventions: `meshpair/synth.py` (`planted_resample`, `synth_sv_pair`),
  `meshpair/reverse_lc.py` (its exporter-relation helper), `meshpair/canvasmap.py` (docstring), `DESIGN.md` (section 2),
  `reverse_route/judge/judge_check.py`, `reverse_route/synth_a1.py`, `CHANGES.md`;
- the git history of any of these files, other branches of the authors' private repository, that repository outside
  this folder, and earlier sessions' notes about them.

You may **run** the rest of the package (the tests import `synth.py`, `meshcheck.py`, `fit.py`, `canvasmap.py`,
`check.py`) without reading it. Work in a fresh copy of the package in which the four GPL modules are deleted; write your
modules there. If SPEC.md is unclear, write the question down and ask the spec writer through the person running the
process; the answer comes back as a dated addition to SPEC.md, never as code.

Put at the top of each new file: `SPDX-License-Identifier: MIT`, the authors (Chris Thomas, plus the implementer), and
"clean-room implementation from cleanroom/SPEC.md v1 (2026-10-18)".

## Acceptance (step 3)

1. All existing tests pass unchanged (`python -m pytest -q tests`).
2. Every acceptance vector in SPEC.md §7 matches exactly (SHA-256 and `repr` equality).
3. The full-record parity of SPEC.md §7.1: 207 / 207 records reproduce `exporter_model` exactly.
4. Chris's condition: the integrator re-runs the mesh stage on all 207 pairs with the new modules and diffs against
   the records made with `11b478606`; any numeric difference blocks the MIT release.
5. Only then: change `LICENSE` to MIT and the README licence note; record who did each step here.

## Record

- 2026-10-18: SPEC.md v1 written (step 1). Reference code `11b478606` (work/w00field/meshpair) checked numerically
  identical to `submission/tifxyz-canvas-map/meshpair` (only an import path in ct.py differs). Not yet implemented.
- 2026-10-18: the 4 files that restated the conventions now call the clean modules where that is bit-identical:
  `meshpair/synth.py` (`planted_resample` index rule -> `exporter_grid_map`; `synth_sv_pair` -> `slice_centre`,
  `layer_depth_um`), `meshpair/reverse_lc.py` (`source_index` -> `exporter_grid_map`). `reverse_route/synth_a1.py`
  already used the clean exporter; its `* 0.05` canvas-px -> cell factor is kept (no clean function).
  `reverse_route/judge/judge_check.py` kept its own `((g + 0.5) * n_src) / n_ref` order (the clean order differs by
  1 ulp on 27 % of indices and flips clean flags on 4 real shape pairs; comment added). Also kept (not covered by the
  clean modules): `planted_resample`'s float64 bilinear body + x-only validity, `block_mean`, `float_floor_px` and
  judge scale factors, `own_bilerp`. Checked: outputs bit-identical before/after (1131 cases), 61 tests pass, spec
  vectors FAILS: []. `judge_check.py` could not be re-run here (no `out/pairs.json`, no mesh cache).
- 2026-10-19, step 3 item 4 (integrator: Claude Code, Sebastien's session): full-catalog run on Windows, 621 mesh files
  fetched fresh. New modules vs the release before the rewrite (`7077a5f21`, which Chris showed bit-identical to
  `11b478606` on Linux): `mesh` 207 / 207 and `reverse-lc` 207 / 207 records bit-identical in every field.
  Against the published Linux records, the differences come from the platform and are the same for the old code
  (`results/PARITY_CLEANROOM.md`). Open: Chris either runs the new code on Linux or accepts the two links.
  Note on order: `LICENSE` was changed to MIT in `709aecec4` (2026-10-18), before item 4 finished, so not in the
  order of item 5. The repository stays private until Chris closes item 4.
- 2026-10-19, item 4 closed (Chris's side, Linux): the new tree `aa9212d3e` gives `mesh` 207 / 207 and `reverse-lc`
  207 / 207 records identical in every field (`seconds`, `tool_version`, `prereg` and `warning` ignored) to the original `11b478606` run, and 207 / 207 + 206 / 206 to the
  published records (`results/PARITY_CLEANROOM.md`).
