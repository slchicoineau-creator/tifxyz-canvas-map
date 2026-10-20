# Changes in this release (0.2.0) against the code that produced the published records (0.1.0)

Every record under `out/` and `reverse_route/` was made with version 0.1.0, exported from the authors' private
repository at commit 11b478606 (2026-10-13; `check.py` at 46906121, `meshpair/reverse_lc.py` frozen at 8a4c6e6c,
`meshpair/ct.py` and `meshpair/catalog.py` at 8ca5789f). Those records are unchanged. The edits below make the code run
outside that repository and add a public helper; none changes a measurement on the inputs the records used.
Each line names the outside-review finding it answers (R1-R7) where there is one.

## Code
1. **R1, bundled reader.** `meshpair/ct.py` imported `zio` from a `lib/` folder three levels up (the private
   repository). The same reader now ships as `meshpair/zio.py` and is imported relatively. The CT route, `plants` and
   `pyramid` crashed without it.
2. **R2, pre-registration guard opt-in.** `check.py _guard` refused every non-development pair unless `--prereg` named a
   commit on a pushed branch of the repository two levels up. The authors' confirmatory runs used that guard. It is now
   off unless `--require-pushed-prereg` is given (and it then checks this repository). `--prereg` is still written into
   each record.
3. **R3, `--segments` takes a file.** It also accepts a path to a file of comma- or newline-separated ids. (The
   published run passed the content of `reverse_route/apply154/segments.txt`; the log shows 154 pairs selected.)
4. **R4, `mesh` command.** It now fetches the fine, coarse and source meshes itself (as `reverse` and `reverse-lc`
   always did), skips pairs already in the output file, and writes `warning` if the catalog gives no source mesh (M2a not run).
5. **R6, empty meshes.** `mesh`, `reverse` and `reverse-lc` write an `error` record for a mesh with no valid vertex
   instead of crashing. No published mesh has none.
6. **R7, float floor and render scale.** `reverse_lc.float_floor_px` converted the float32 index resolution to canvas
   pixels with `1 / s`; the canvas map uses `R / s`. It now multiplies by the render scale `R`. Every published run used
   `R = 1`, so no bar changed.
7. **Catalog.** The sample exclusion list is empty (`EXCLUDED_SAMPLES`); it named two samples that have no segments in
   the catalog, so the pair list is unchanged (207 pairs, identical to `out/pairs.tsv`, checked 2026-10-18).
8. **`inventory` writes `pairs.tsv` with `\n` line ends on every OS** (it wrote `\r\n` on Windows).
9. **New: `meshpair/canvasmap.py`** (public helper: `fine_to_coarse_px`, `coarse_to_fine_px`, `linear_coeffs`,
   `warp_canvas_image`, `mesh_geom`; `python -m meshpair.canvasmap FINE COARSE`). It reads each mesh's stored `scale`
   (published meshes store 0.05 either as a double or as float32 0.0500000007).
10. **New: `scripts/make_pair_maps.py` -> `results/pair_maps.tsv`** (the map of all 207 pairs from the stored records).
11. **New tests:** `tests/test_release.py` (8 tests). Total 61 pass (53 original + 8).
12. **Default outputs** of `mesh`, `reverse`, `reverse-lc`, `ct` and `selftest` moved to `out/local/` (they pointed at
    the published records, and a run appended to them). `--segments` also strips spaces. The CT chunk cache of
    `meshpair/zio.py` is `<repository>/cache` (it was the private repository's `cache/`); no measurement effect.
13. **Clean-room rewrite for the MIT licence** of the modules that follow `volume-cartographer` behaviour
    (`conventions.py`, `exporter.py`, `tifxyz.py`, new `slicecentre.py` used by `ct.py`): see `cleanroom/README.md`
    and the full-catalog parity check in `results/PARITY_CLEANROOM.md`.

## Documents (release redactions and corrections)
- `DESIGN.md`: the owner/team-process paragraph, one line about excluded samples, machine names and section 7 (open
  questions for the authors, answered in the pre-registration) are removed; paths made relative. Nothing technical changed.
- `out/RESULTS.md`, `reverse_route/RESULT.md`, `reverse_route/s2/RESULT.md`, `reverse_route/apply154/RESULT.md`:
  references to the private repository's branch and internal tracking documents are reworded.
- `out/ct_summary.json`, `out/mset.log`: absolute paths on the authors' machine made relative.
- **Correction** in `reverse_route/apply154/RESULT.md`: the caveat line said "9 of 206 pairs have a total of 9
  off-surface clean vertices"; the table above it gives 9 vertices in **5** pairs. Fixed to 5.
- **Correction** in `reverse_route/apply154/RESULT.md`: it said the PHerc0841 development pair (w00, 20260220213127) is
  among the 2 PHerc0841 pairs; it is not (they are 20260220214732 and 20260221022814; w00 is excluded from all 206).
- Removed: the pre-registration draft and the internal upstream-comment drafts (and the line in `out/RESULTS.md` that
  pointed at them); an internal item id in `meshpair/ct.py` and `DESIGN.md`; private paths in
  `reverse_route/judge/judge_check.py`. "Independent judge" reworded to say it was a separate AI review session.

- Added: `preregs/`, scrubbed copies of the six pre-registrations the records cite (only team-internal references
  removed, each replaced by a bracketed note), with an index of commits and times (`preregs/README.md`).

## Regression check (packaging only, not a new result)
See `results/REGRESSION_2026-10-18.md`.
