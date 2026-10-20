# Clean-room implementer log (card 02, step 2)

Implementer: Claude Code (Opus 5.5), a fresh subagent session that had never seen the GPL sources. 2026-10-18.

## What I read
- `cleanroom/README.md`, `cleanroom/SPEC.md` (v1).
- Tests: `tests/test_meshpair.py`, `tests/test_release.py`, `tests/test_reverse.py`, `tests/test_reverse_lc.py`
  (imports only for `tests/test_plants_a1.py`).
- Data: the first 3 lines of `out/mesh_records.jsonl` (field names); `pyproject.toml`.
- Public NumPy / tifffile behaviour from memory. No villa code was needed.

## What I did NOT open
`meshpair/conventions.py`, `exporter.py`, `tifxyz.py`, `ct.py`, `synth.py`, `reverse_lc.py`, `canvasmap.py`, `DESIGN.md`,
`CHANGES.md`, `reverse_route/judge/judge_check.py`, `reverse_route/synth_a1.py`, any volume-cartographer code, any git
history. In the scratch copy the three GPL modules were deleted by path and `ct.py` was moved (unread) to a quarantine
folder. Tests ran with `--tb=line`, so no source lines of unread files were printed.

## Output (`cleanroom/impl/`)
| File | SPEC part |
|---|---|
| `meshpair/conventions.py` | A |
| `meshpair/exporter.py` | B |
| `meshpair/tifxyz.py` | C |
| `meshpair/slicecentre.py` | D: `slice_centre(N)`, `layer_depth_um(k, N, px, delta=0)`, `depth_to_layer(d, px, N0)`, `coarse_layers(n_c, count=None)`. New module name; the integrator wires it into `ct.py` (SPEC D says the integrator does this). |
| `acceptance/check_spec_vectors.py` | (extra, for step 3) every SPEC §7 vector + the §7.1 full-record parity; run `python check_spec_vectors.py <package root>`; exit 1 on any mismatch; expects 15 keys for the no-px_um case (gap 1) |

Scratch copies holding the old `ct.py` were deleted after the runs.

Each module has the SPDX MIT header, the authors (Chris Thomas + this implementer) and the clean-room line.

## Results (scratch copy, NumPy 2.4.3, SciPy 1.17.1, Windows x64)
- Without `ct.py` (the 4 modules replaced, `ct.py` removed): `test_release`, `test_reverse`, `test_reverse_lc`,
  `test_plants_a1`: **46 / 46 pass**. `test_meshpair.py` cannot be collected without `ct.py` (it imports `CanvasMap`,
  `TileSpec`, ... from it).
- With the old `ct.py` put back unread next to the new modules: **61 / 61 pass** (`python -m pytest -q tests`).
  Note: this run used the OLD `ct.py`, so its slice-centre code, not `slicecentre.py`, ran in the `test_ct_*` tests;
  `slicecentre.py` is checked only by the SPEC D vectors. The integrator must rerun the tests after wiring it in.
- SPEC §7 vectors: all match exactly (coefficients, canvas sizes, canvas_to_grid order check, maps, offset_summary,
  every exporter SHA-256, counts, ratios, shapes, AF-edge, AF-m1, CM, CM2, vc_load, geom, coarse_layers), except the
  key-count line noted below (a SPEC error, not an implementation difference).
- **207-pair parity: yes.** `test_pair_maps_table_reproduces_records` passes (difference 0.0), and the stronger §7.1
  full-record check gives 207 / 207 `exporter_model` dicts identical (same keys, same order, every value `==`);
  split (f32, 0.05) 72, (f32, f32) 126, (0.05, f32) 9, as SPEC says.

## Spec gaps / questions for the spec writer
1. **Key count in A.3 and §7.1 is off by one.** SPEC says 22 keys (14 without px_um). The key list it gives
   (`k_naive` + 7 per axis + 8 `_um`) is 23 (15 without px_um); the §7.1 listing has 23 entries and the 207 stored
   records match my 23-key dict exactly. My `k=0.25` / no-px_um call returns 15 keys; the "14" in §7.1 should read 15.
2. SPEC D gives no module or function names for the slice-centre helper; I chose `meshpair/slicecentre.py` with the
   names above. `ct.py` itself was not rewritten (SPEC scope says the integrator wires the helper).
3. `Mesh.geom()` with no `"scale"` key: SPEC says use `[0.05, 0.05]`; implemented with `meta.get("scale", ...)`.
4. B.6: SPEC does not say what to do if `A` is not None but the 3x3 part is singular (μ = 0 → division by zero in the
   uniformity test). Not covered; NumPy gives a warning and the uniform branch is not taken.
5. Non-float32 input to the resampler follows NumPy promotion (SPEC says only float32 occurs); not tested.
