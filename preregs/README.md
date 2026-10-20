# Pre-registrations (scrubbed public copies)

Each file here is a public copy of a pre-registration or amendment (the records name it `prereg/<file name>`). Chris Thomas's side committed it to the authors'
private repository; the records under `out/` and `reverse_route/` name it by its commit (`"prereg"` field). Each copy is
the file as it stood at that commit. Each of these files has only that one commit.

**What was changed in the copies:** only team-internal references were removed: internal task ids and lists, the
team's lane name, a private note path, internal test and data names, the machine name and a note on its GPU-sharing
tool. Each removal is replaced by a visible
`[bracketed note]`, and nothing else was edited.

**Paths and names inside the copies:**
- `work/w00field/` is the root of this repository (`check.py`, `meshpair/`, `out/`, `reverse_route/`, `phase2/`).
- `reviews/...` and `lib/zio.py` point into the private repository. `lib/zio.py` ships here as `meshpair/zio.py`.
- "Judge" and "the vet" are separate AI review sessions that did not do the work. Their rulings are summarised in the
  result files.
- The commit times are in that private repository. The authors can show the commits on request.

| file | private commit | committed (UTC) | records that cite it | note |
|---|---|---|---|---|
| [2026-10-11_partner_w00field.md](2026-10-11_partner_w00field.md) | `8ca5789f227e1d1c230932cfb9649b6253fef46a` | 2026-10-11T13:11:08Z | mesh route (`out/mesh_records.jsonl`), CT pass 1 | The file states it was pushed before any confirmatory run. What was seen on the development pair w00 before it is disclosed in the file. |
| [2026-10-11_partner_w00field_amend1.md](2026-10-11_partner_w00field_amend1.md) | `13df569677a86fb20e261cc932bf06d5fc36f3e7` | 2026-10-11T16:19:02Z | CT pass-1 plants (`out/ct_plants_a1_*.jsonl`) | **Post-sight, disclosed in the file.** The CT C1 control changed after the C2-C5 values had been printed, so the review ruled the pass-1 CT verdicts exploratory (see amendment 2). |
| [2026-10-13_partner_w00field_amend2.md](2026-10-13_partner_w00field_amend2.md) | `fcb2012f49624c1872d1d2cdb82a9d22619e5785` | 2026-10-13T20:33:47Z | CT repeat on 7 fresh segments (`out/amend2/`) | The file states it was written and pushed before any measurement on the fresh segments. |
| [2026-10-13_partner_tifxyz_reverse_route.md](2026-10-13_partner_tifxyz_reverse_route.md) | `ed1e660c4480aac22c0acdeb2d0c9d133dadd567` | 2026-10-13T11:43:58Z | reverse route S0/S1 (`reverse_route/s0.json`, `s1*.json*`) | The file states it was pushed before any real-mesh run of the new code. S1 stayed INCONCLUSIVE as filed. |
| [2026-10-13_partner_tifxyz_reverse_amend1.md](2026-10-13_partner_tifxyz_reverse_amend1.md) | `8a4c6e6c07cab72772f2d4833afdf8824dd64127` | 2026-10-13T12:30:07Z | lattice-clean rule, 52 held-out S2 pairs (`reverse_route/s2/`) | The rule was designed after S1 was seen (disclosed). It was tested only on the S2 pairs, pushed before any S2 mesh was run. |
| [2026-10-13_partner_tifxyz_reverse_apply154.md](2026-10-13_partner_tifxyz_reverse_apply154.md) | `b8930ac6cc14b3ac894a3a52c8485fe875403ed5` | 2026-10-13T13:02:37Z | application to the other 154 pairs (`reverse_route/apply154/`) | **An application, not a test** (the file says so). The rule, code and scorer were frozen before the run. |
