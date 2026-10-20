# Full-catalog parity of the clean-room rewrite (2026-10-18/19)

Condition for the MIT licence (Chris Thomas, review correction `bb7d203fe`): the rewritten modules must give
bit-identical results to his original code `11b478606` on all 207 pairs. This file records the check, done in two links.

## Link 1 (Chris, Linux, 2026-10-18): release 0.2.0 before the rewrite = original

Chris's review of the release copy at `7077a5f21` (`reviews/2026-10-18_partner_judge_tifxyz_release.md` in the
authors' private repository, parity table) ran the original `11b478606` and that release head to head on Linux
(Python 3.12, numpy 2.5.4, scipy 1.18.1, the environment that made the records):

- `mesh`, whole catalog + w00: 207 / 207 records identical, and 207 / 207 identical to the published
  `out/mesh_records.jsonl` + `out/dev_w00_mesh.jsonl`;
- `reverse-lc`, whole catalog: 207 / 207 identical, and 206 / 206 identical to the published `reverse_route/s2/s2.jsonl`
  + `reverse_route/apply154/apply154.jsonl`, including the sha256 of every fine, coarse and source mesh file.

## Link 2 (this check, Windows, same PC, same inputs): clean-room code = release 0.2.0 before the rewrite

- **Old:** `check.py` + `meshpair/` exactly as in `7077a5f21` (14 files, byte-compared with the commit).
- **New:** the current tree, with the clean-room `conventions.py`, `exporter.py`, `tifxyz.py`, `slicecentre.py` and the
  modules changed to call them (`synth.py`, `reverse_lc.py`, `ct.py`; see `cleanroom/README.md`).
- **Inputs:** 207 pairs from the live catalog (sha256 55aa3d21…44e3; pair list identical to `out/pairs.tsv`); all 621 mesh
  files downloaded fresh from the open-data bucket (15 GB) into one cache that both versions read.
- **Runs:** `check.py mesh` then `check.py reverse-lc` with `--pairs`, `--cache` and `--out` set, old and new in parallel;
  all four runs exit 0 (finished 2026-10-19 00:15, CPU only). Windows 11, Python 3.12.11, numpy 2.5.4,
  scipy 1.18.1.
- **Comparison:** record by record, every field, floats by exact equality, ignoring only `seconds`, `tool_version`,
  `prereg` and `warning` (the same rule as Link 1).

| route | pairs | new vs old: records with every field bit-identical | other differences |
|---|---|---|---|
| `mesh` (M1-M3, M2a, routes F and I, exporter model, replica) | 207 | **207 / 207** | none (0 missing, 0 extra, 0 structural) |
| `reverse-lc` (lattice-clean, incl. mesh sha256) | 207 | **207 / 207** | none (0 missing, 0 extra, 0 structural) |

Verdicts of the new code: `mesh` M1 136 VOID_ROUTE_I, 71 IN_FILES (70 + w00), M2a reproduced 207 / 207; `reverse-lc`
204 PASS, 3 NOT_EVALUABLE, 0 FAIL (= 52 held-out + 151 applied + w00, which is not one of the 206).

## Windows against the published (Linux) records, for context

Old and new are bit-identical to each other, so the differences below come from the platform. They are the same for
both versions.

- Identical on 207 / 207: the verdict lines (M1-M3, M2a), `exporter_model` (the pixel map), grid shapes and the M2a
  replica result. `reverse-lc`: 206 / 206 verdicts and every sha256 identical. Its floats differ by at most 3.6e-10.
- The least-squares fits of `mesh` routes F and I differ:
  - optimizer stop code (`route_F.status`) on 30 pairs, all VOID_ROUTE_I;
  - function-evaluation counts (`nfev`) on 61 pairs;
  - converged-vertex counts (`route_I.n_ok`) by 1-3 on 18 pairs, 9 of them IN_FILES;
  - floats: at most 1.3e-7 on IN_FILES pairs, and at most 1.9e-3 µm on VOID pairs (`normal_dist_all_converged_um.p95_abs`).
- The 9e-8 in `REGRESSION_2026-10-18.md` was measured on 5 pairs only. Across the full catalog, the platform
  differences in routes F and I are the ones listed above.

## What this shows, and what is left

- The clean-room code gives the same bits as the release code on every pair and every field, on the same machine.
  Link 1 shows the release code gives the same bits as the original on Linux. Together, the rewrite changes no
  result.
- **Done on Linux too (2026-10-19, Chris's side):** the new tree (`aa9212d3e`), run on the machine and in the environment
  that made the records, with all 621 meshes fetched fresh. Results: `mesh` 207 / 207 and `reverse-lc` 207 / 207 records
  identical in every field (`seconds`, `tool_version`, `prereg` and `warning` ignored, as in Link 1) (maximum difference 0) to the original `11b478606` run. Against the published records:
  207 / 207 and 206 / 206. Chris's condition is met as written, and step 4 of `cleanroom/README.md` is closed.
- The `ct` route was not re-run (it streams CT). Its only code change is that the old expressions `(n - 1) / 2.0` and
  `coarse_layers` now come from `slicecentre.py`. Checked directly: `slice_centre(n)` equals `(n - 1) / 2.0` bit for bit,
  and `coarse_layers` equals the old function, for every n from 1 to 399 and every layer count (80,199 cases, 0 mismatches).
- Not re-run: `reverse_route/judge/judge_check.py` (it kept its own index order; see `cleanroom/README.md`).
