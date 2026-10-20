# Packaging regression check (2026-10-18): release 0.2.0 vs the published records (0.1.0)

**Not a new result.** It checks that the release edits in `CHANGES.md` do not change what the code measures.

- **Pair list:** `check.py inventory` on the live catalog of 2026-10-18 (`metadata.json` sha256 `55aa3d21…44e3`, newer
  than the 2026-10-11 catalog of the records) gives 207 pairs, line for line identical to `out/pairs.tsv`
  (keys, voxel sizes, copy flags, surface volumes, exporter commits).
- **Reruns:** `check.py mesh` and `check.py reverse-lc` with release 0.2.0 on 5 public pairs, one per kind:
  PHerc0841 20260220214732 (coarse is a copy, forward), PHerc0139 20260206000002 (both resampled, reverse),
  PHercParis4 20260603185441 (held-out S2, reverse), PHerc1667 20251212185248 and PHerc0009B 20250919131352 (forward).
  Windows, Python 3.12, numpy 2.5.4, scipy 1.18.1 (the records were made on Linux with the same numpy / scipy versions).
- **Compared** field by field with `out/mesh_records.jsonl` and `reverse_route/{s2/s2.jsonl, apply154/apply154.jsonl}`,
  ignoring only `seconds`, `tool_version`, `prereg` and the new `warning` field.

| route | verdicts | identical fields | largest difference |
|---|---|---|---|
| mesh (M1-M3, M2a) | 5 / 5 identical | exporter model (the map), valid box, replica result | route F (free affine fit) 9.0e-8 absolute, 3.1e-5 relative; route I 2.9e-11; route F used 6 solver evaluations instead of 5 on one pair |
| lattice-clean | 5 / 5 identical (all PASS, same direction) | mesh sha256 of all 15 files, counts, bars | 8.4e-12 absolute |

The differences are floating-point (operating system / BLAS); no verdict, count or map coefficient moved.

- Also: `reverse_route/apply154/score_apply154.py` regenerates `summary.json` and `pairs_table.jsonl` identical to the
  stored files (apart from line ends on Windows); `scripts/make_pair_maps.py` reproduces the per-scroll offset table of
  `out/RESULTS.md` exactly; `tests/test_release.py` checks that all 207 stored maps equal the conventions code to 0.0.
- Files: `regression_2026-10-18_summary.txt` (the comparison above), `regression_2026-10-18_detail.txt` (first
  differences per pair).
