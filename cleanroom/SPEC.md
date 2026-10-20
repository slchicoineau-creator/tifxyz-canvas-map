# Clean-room functional specification: tifxyz pixel conventions (card 02, `meshpair`)

Version 1, 2026-10-18. Reference behaviour = Chris Thomas's code at commit `11b478606` (paths `work/w00field/meshpair/`),
which is numerically identical to the release copy `submission/tifxyz-canvas-map/meshpair/` (diffed 2026-10-18: the only
difference is an import path in `ct.py`, not in any number).

This document describes **behaviour**: names, signatures, inputs and outputs, data types, formulas, rounding, and test
vectors. It is written so that someone who has never read the GPL sources (volume-cartographer C++, or this package's
`conventions.py`, `exporter.py`, `tifxyz.py`, `ct.py`) can write modules that give **bit-identical** results. Read
`README.md` in this folder first for the rules.

Scope (four parts):

| Part | Module to write | What it is |
|---|---|---|
| A | `meshpair/conventions.py` | renderer canvas <-> grid map, exporter half-pixel index map, the exact fine->coarse canvas map, offset summary |
| B | `meshpair/exporter.py` | Python model of the mesh transfer step (affine in float32, grid-spacing measure, size rule, half-pixel bilinear resample that keeps invalids), plus a mesh comparison |
| C | `meshpair/tifxyz.py` | tifxyz read/write/fetch, including the "vc_load" loader rule |
| D | slice centre (used by `meshpair/ct.py`) | which surface-volume layer sits at depth 0, and the default coarse-layer choice |

---

## 0. Notation and number rules (apply everywhere)

- `fl64(e)` = evaluate `e` in IEEE-754 binary64, rounding to nearest after **every** operation, in the order written.
  `fl32(e)` = the same in binary32. All parentheses in this document are significant: they give the order of operations.
  Where no parentheses are written, the order is left to right (`a + b - c` means `(a + b) - c`; `a * b / c` means
  `(a * b) / c`).
- "float32" / "float64" = NumPy `float32` / `float64`. Integer → float conversions are exact for all sizes that occur.
- A mesh point array has shape `(rows, cols, 3)`. Component 0 = x, 1 = y, 2 = z. Axis 0 of the array = rows = **y of the
  canvas**, axis 1 = cols = **x of the canvas**. The value -1 marks an invalid vertex.
- meta.json `"scale"` = `[s_x, s_y]`. **`s_x` pairs with cols (axis 1), `s_y` with rows (axis 0).** Published files hold
  0.05 either as the double 0.05 or as the float32 value 0.05000000074505806 (the double nearest to float32(0.05)). The
  code always uses the stored value as a Python float; it never rounds it to 0.05.
- In this spec, "axis" arguments are integers: 0 = y/rows, 1 = x/cols.
- Byte layout used for all SHA-256 vectors in §7: the array as float32, C order (row-major), little-endian, shape
  `(rows, cols, 3)`, raw bytes with no header.

Two different integer rounding rules occur. Do not mix them up:

| Where | Rule |
|---|---|
| canvas size (`canvas_w`, `canvas_h`) | `lround(x)` = `floor(x + 0.5)` for x >= 0, `-floor(-x + 0.5)` for x < 0; computed in float64 (round half away from zero, evaluated as written; note this is not identical to C `lround` for x = 0.49999999999999994, and the replica's rule is the one to copy) |
| resample grid size (`transform_tifxyz`) | NumPy `round` of a **float32** value = round half to **even**, then `int` (e.g. 4.5 → 4, 7.5 → 8) |

---

## A. `meshpair/conventions.py`

### A.1 What the conventions are (facts)

Renderer (volume-cartographer `vc_render_tifxyz`, render scale R = `--scale`, default 1): a mesh with `n` grid vertices
along an axis and grid step `s` (meta scale for that axis) is drawn on a full-resolution canvas of
`W = lround(n * R / s)` pixels along that axis. Canvas pixel `p` (0-based, pixel centres at integers) sits at the
fractional grid index

    g(p) = s * (p - 0.5 W + 0.5) / R + n / 2           (centre-anchored map)

With R = 1 and W = n / s exactly, this reduces to `g = s (p + 0.5)`, i.e. `p = g / s - 0.5`.

Exporter (`vc_transform_geom` → resample): when a mesh with `n_src` vertices along an axis is resampled to `n_dst`
vertices, destination vertex `d` takes its value from source index `(d + 0.5) * (n_src / n_dst) - 0.5` (half-pixel
centres). Two meshes made from one source by such resamples therefore satisfy, per axis,

    g_C = (g_F + 0.5) * (n_C / n_F) - 0.5               (exporter_grid_map, fine F -> coarse C)

The exact canvas map of a fine pixel to the coarse canvas = fine canvas → fine grid → coarse grid → coarse canvas. It is
affine per axis, `p_C = a + b p_F`. With s = 0.05 and R = 1 on both meshes it is `p_C = (p_F + 10.5) * n_C / n_F - 10.5`.

### A.2 `MeshGeom`

An immutable, hashable dataclass (frozen; equality and hash by field values; normal dataclass repr). Fields in this
**positional order**, with defaults:

| # | field | type | default |
|---|---|---|---|
| 1 | `rows` | int | (required) |
| 2 | `cols` | int | (required) |
| 3 | `scale_x` | float | 0.05 |
| 4 | `scale_y` | float | 0.05 |
| 5 | `render_scale` | float | 1.0 |
| 6 | `px_um` | float | NaN (voxel size in µm of the volume the mesh lives in) |

Callers construct it positionally (`MeshGeom(926, 738, s, s)`) and by keyword (`MeshGeom(rows=..., cols=..., px_um=...)`),
and test `isinstance(x, MeshGeom)`.

Read-only properties (Python int):
- `canvas_w = max(1, lround(fl64((cols * render_scale) / scale_x)))`
- `canvas_h = max(1, lround(fl64((rows * render_scale) / scale_y)))`

Per-axis parameters: axis 0 → `(n, s, W) = (rows, scale_y, canvas_h)`; axis 1 → `(cols, scale_x, canvas_w)`.
`R = render_scale`.

Methods (input `p` or `g`: a Python scalar or array-like; it is converted to a float64 NumPy array first; the result is
a float64 NumPy array of the same shape, a 0-d array for a scalar input):

- `canvas_to_grid(p, axis)` = `fl64( ((s * ((p - 0.5*W) + 0.5)) / R) + (n / 2.0) )`
  (inner order: `0.5*W` first, subtract it from p, add 0.5, multiply by s, divide by R; separately `n / 2.0`; add.)
- `grid_to_canvas(g, axis)` = `fl64( ((((g - n/2.0) * R) / s) + 0.5*W) - 0.5 )`

These are exact inverses only in real arithmetic; in float64 they round-trip within `allclose` (test below).

### A.3 Free functions (all float64; inputs converted to float64 arrays; outputs float64 arrays unless stated)

- `exporter_grid_map(g_f, n_f, n_c)` = `fl64( (g_f + 0.5) * q - 0.5 )` with `q = fl64(n_c / n_f)` computed first as a
  true division of the two integers. No clamping.
- `exporter_grid_map_inv(g_c, n_f, n_c)` = `fl64( (g_c + 0.5) * q' - 0.5 )` with `q' = fl64(n_f / n_c)`.
- `exact_canvas_map(p_f, fine, coarse, axis)`: take `n_f`, `n_c` = rows (axis 0) or cols (axis 1) of `fine` / `coarse`;
  compute `g_f = fine.canvas_to_grid(p_f, axis)`, then `g_c = exporter_grid_map(g_f, n_f, n_c)`, then return
  `coarse.grid_to_canvas(g_c, axis)`. Each step uses its own mesh's s, W, R.
- `exact_canvas_map_inv(p_c, fine, coarse, axis)`: `g_c = coarse.canvas_to_grid(p_c, axis)`,
  `g_f = exporter_grid_map_inv(g_c, n_f, n_c)`, return `fine.grid_to_canvas(g_f, axis)`.
- `naive_canvas_map(p_f, k)` = `fl64(p_f * k)` (what index-pairing scripts did).
- `exact_linear_coeffs(fine, coarse, axis)` → `(a, b)` as two Python floats, obtained **by evaluating the map, not from a
  closed form**: `p0 = exact_canvas_map(0.0, ...)`, `p1 = exact_canvas_map(1.0, ...)`, `a = p0`, `b = fl64(p1 - p0)`.
  The 207-pair test compares these with **zero tolerance**; an analytic `b = (n_C/n_F)(s_F/s_C)` does not match bits.
- `offset_summary(fine, coarse, k=None, valid_box=None)` → dict (key order matters: records are JSON).
  - `k` defaults to `fl64(fine.px_um / coarse.px_um)`; if given, it is used and stored unchanged.
  - `valid_box = (y0, y1, x0, x1)` in fine canvas px; default `(0, fine.canvas_h - 1, 0, fine.canvas_w - 1)`.
  - Output keys, in this order: `k_naive` (= k). Then for axis 0 with prefix `y` (bounds `lo, hi = y0, y1`) and then axis
    1 with prefix `x` (bounds `x0, x1`):
    - `(a, b) = exact_linear_coeffs(fine, coarse, axis)`; `ends = [lo, hi]` as a float64 array;
      `off = fl64( (a + b*ends) - k*ends )` element-wise;
    - `<p>_a` = a; `<p>_b` = b; `<p>_b_minus_k` = `fl64(b - k)`; `<p>_offset_at_origin_px` = a;
      `<p>_offset_box_min_px` = min(off); `<p>_offset_box_max_px` = max(off);
      `<p>_offset_box_centre_px` = `fl64( (a + b*mid) - k*mid )` with `mid = fl64(0.5 * (lo + hi))`.
      All values as Python floats.
  - Then, only if `coarse.px_um` is finite: for every key present so far whose name ends in `_px`, in the order they
    appear, append key `<name without "_px">_um` = `fl64(value * coarse.px_um)`. That gives, in order:
    `y_offset_at_origin_um, y_offset_box_min_um, y_offset_box_max_um, y_offset_box_centre_um, x_offset_at_origin_um,
    x_offset_box_min_um, x_offset_box_max_um, x_offset_box_centre_um` (22 keys in total; 14 if px_um is NaN).

Private helpers in the reference (`_lround`, `MeshGeom._axis`) are not imported anywhere else; names are free.

---

## B. `meshpair/exporter.py`

Module constant `F32` = the NumPy `float32` type (imported by `synth.py`, `tests/test_reverse.py`,
`reverse_route/synth_a1.py`; they use it as `F32(-1)`, `F32(x)` and as a dtype).

### B.1 Two validity predicates (they differ on purpose)

- **Point validity** (used by the affine step, the spacing measure, and `compare_meshes`): a vertex is valid when its
  x component is not equal to -1 **and** all three components are finite. A vertex with only y = -1 or only z = -1 is
  therefore **valid** here.
- **Sample validity** (used only inside the resampler): valid when **none** of the three components equals -1 **and**
  all three are finite.

Comparisons with -1 are done against float32 -1 on the array's own dtype.

### B.2 `_median_cpp` rule (median used by the spacing measure)

Sort ascending (NumPy sort, float64). With `n` values and 0-based index: if n = 0 → NaN. Take `m = v[n // 2]` (the upper
middle). If n is even, `m = fl64(0.5 * (m + v[n//2 - 1]))` (upper + lower, then times 0.5). Return a Python float.
Examples: [3, 1, 2, 4] → 2.5; [3, 1, 2] → 2.0.

### B.3 `measure_grid_axis_spacing(points)` → `(med_h, med_v, n_h, n_v)`

`med_h` = median distance between horizontal neighbours (same row, columns c and c-1, i.e. along x / cols);
`med_v` = median distance between vertical neighbours (rows r and r-1, along y / rows); `n_h`, `n_v` = how many
distances went into each median (Python ints).

- If rows < 2 or cols < 2: return `(NaN, NaN, 0, 0)`.
- Sample window:
  - if rows < 20 or cols < 20: rows `r = 1, 2, ..., rows-1` and cols `c = 1, ..., cols-1` (step 1);
  - else: `r0 = max(1, int(fl64(rows * 0.1)))`, `r1 = int(fl64(rows * 0.9))`, the same for cols with `c0, c1`
    (multiply by the double literals 0.1 / 0.9, then truncate toward zero); then `r1 = max(r0 + 1, r1)`,
    `c1 = max(c0 + 1, c1)`; step 4: `r = r0, r0+4, ... < r1`, `c = c0, c0+4, ... < c1`.
- Every (r, c) pair of the two index lists is a sample (outer product).
- Horizontal distance for a sample: Euclidean distance between `P[r, c]` and `P[r, c-1]`, used only if both are valid
  (point validity). Vertical: between `P[r, c]` and `P[r-1, c]`, if both valid.
- Distance arithmetic: convert both points to float64, subtract component-wise, then the NumPy 2-norm along the last
  axis, which is `sqrt( fl64((dx*dx + dy*dy) + dz*dz) )`.
- Keep only distances that are finite and > 0. Medians by B.2. Order in which samples are collected does not matter
  (the median sorts).

### B.4 `apply_affine_f32(points, A, scale_before=1.0, scale_after=1.0)` → new float32 array

`A` is a 3×4 affine `[L | t]` (row i: `L[i,0..2]`, `t[i]`). The input is not modified.

1. `out` = copy of `points` as float32.
2. Valid set `V` = point validity (B.1) on `out`.
3. For each valid vertex: `p = fl32(point * float32(scale_before))` (component-wise float32 product).
4. `q = p (as float64) times L-transpose, plus t`: computed with NumPy's float64 matrix product (`@`) of the N×3 array of
   valid points with the transpose of the 3×3 linear part, then the length-3 translation added in float64. (Bit note: a
   hand-written dot product, or one using fused multiply-add, can differ in the last float64 bit, which very rarely
   changes the float32 result. Parity is defined against NumPy's float64 matmul on the same machine/library.)
5. `q = fl32( float32(q) * float32(scale_after) )` (cast to float32 first, then a float32 product).
6. Any transformed vertex with a non-finite component, or whose new x equals -1, becomes `(-1, -1, -1)`.
7. Write the transformed values back into the valid positions. **Invalid input vertices are copied unchanged** (they are
   not normalised to -1; e.g. `(NaN, 2, 3)` stays `(NaN, 2, 3)`, `(-1, 0, 0)` stays `(-1, 0, 0)`). A vertex `(5, -1, 7)` is
   valid by B.1 and is transformed (with A = 2·I: `(10, -2, 14)`).

### B.5 `resample_linear_preserving_invalids(points, new_h, new_w)` → float32 array `(new_h, new_w, 3)`

Source size `rows, cols` = `points.shape[:2]`. The output starts filled with float32 -1.

Per axis (x shown; y is the same with rows/new_h):
1. `x_scale = fl32(float32(cols) / float32(new_w))` (float32 division).
2. For destination column `d = 0 .. new_w-1` (d as float32): `sx = fl32( ((d + 0.5) * x_scale) - 0.5 )`, every operation in
   float32, in that order.
3. Clamp `sx` to `[0, cols-1]` (float32 bounds).
4. `x0 = floor(sx)` as an integer; `x1 = min(x0 + 1, cols - 1)`; weight `fx = fl32(sx - float32(x0))`.
   (At the clamped far edge `x0 = cols-1`, `x1 = x0`, `fx = 0`.)

Value at destination (row i, col j), with corners `p00 = P[y0_i, x0_j]`, `p01 = P[y0_i, x1_j]`, `p10 = P[y1_i, x0_j]`,
`p11 = P[y1_i, x1_j]`, and `one = float32 1`:

    top = (p00 * (one - fx)) + (p01 * fx)
    bot = (p10 * (one - fx)) + (p11 * fx)
    res = (top * (one - fy)) + (bot * fy)                then cast to float32

For float32 input every operation is a float32 operation (`one - fx` is rounded first). Do **not** use the algebraically
equal form `p00 + fx (p01 - p00)`: it gives different bits. (For non-float32 input the reference follows NumPy type
promotion of the input dtype with float32; only float32 input occurs in practice.)

A destination vertex is written only if **all four** corners pass **sample validity** (B.1); otherwise it stays
`(-1, -1, -1)`. Corners are taken from the clamped indices, so no out-of-range invalidation happens.

### B.6 `transform_tifxyz(points, A, scale_before=1.0, scale_after=1.0)` → `(out, info)`

1. `base = measure_grid_axis_spacing(points)` (on the input as given).
2. `out = apply_affine_f32(points, A_or_identity, scale_before, scale_after)`; if `A is None` the identity `[I | 0]` is used.
3. `after = measure_grid_axis_spacing(out)`.
4. `ratio_x = after.med_h / base.med_h` (float64) if `base.n_h > 0` and `after.n_h > 0` and `base.med_h > 0`, else NaN.
   `ratio_y` the same with `med_v` / `n_v`. (x pairs with the horizontal = cols spacing.)
5. Fallback scale: `fb = scale_before * scale_after` (float64). If `A` is not None: take the singular values `σ` of the 3×3
   linear part (float64, NumPy `linalg.svd` with `compute_uv=False`), `μ = mean(σ)`; if
   `max_i |σ_i - μ| / μ <= 1e-4` (inclusive) then `fb = fb * μ`. (A non-uniform linear part leaves `fb` unchanged.)
6. If `ratio_x` is not (finite and > 0) → `ratio_x = fb`; same for `ratio_y`.
7. `info` = dict with keys in this order: `spacing_before` = `(base.med_h, base.med_v)` (a 2-tuple),
   `spacing_after` = `(after.med_h, after.med_v)`, `ratio_x`, `ratio_y`, `resampled` = False.
8. Resample decision:
   - if `|ratio_x - 1.0| > 1e-4` **or** `|ratio_y - 1.0| > 1e-4` (float64, strict):
     - `fx = float32(ratio_x)`, `fy = float32(ratio_y)`;
     - if `|fx - 1| < float32(0.001)` **and** `|fy - 1| < float32(0.001)` (float32 arithmetic, strict): return
       `(out, info)` with no resample (the affine output is returned, `resampled` False, no `new_shape` key);
     - else `new_w = max(1, int(round_half_even(fl32(float32(cols) * fx))))`,
       `new_h = max(1, int(round_half_even(fl32(float32(rows) * fy))))` (cols/rows of `out`), then
       `out = resample_linear_preserving_invalids(out, new_h, new_w)`, and set `info["resampled"] = True`,
       `info["new_shape"] = (new_h, new_w)` (a tuple of Python ints, added as the last key).
   - else: no resample.
   - A resample happens even if the new shape equals the old one (vector T-win2 in §7).
9. Return `(out, info)`; `out` is float32.

Downstream code turns tuples in `info` into lists for JSON, and reads `info["new_shape"]`, `info["ratio_x"]`,
`info["ratio_y"]`, `info["resampled"]`.

### B.7 `compare_meshes(pred, pub)` → dict

Keys in order: `pred_shape` = list of the first two dims of pred; `pub_shape` likewise; `same_shape` = whether the
**full** shapes (including the last dim) are equal. If not equal, return now (3 keys). Otherwise, with point validity
`vp` (pred), `vq` (pub), `both = vp & vq`, and `e` = float64 Euclidean distances (B.3 arithmetic) between pred and pub at
the `both` vertices, append:

`n_valid_pub` = count(vq), `n_valid_pred` = count(vp), `n_both` = count(both), `mask_mismatch` = count(vp XOR vq),
`mask_agreement` = `1 - mask_mismatch / max(1, count(vp OR vq))`, `err_median_vx` = NumPy median of e,
`err_p99_vx` = NumPy `percentile(e, 99)` (default linear interpolation), `err_max_vx` = max of e. The three error values
are NaN if `e` is empty. Counts are Python ints, the rest Python floats.

Callers merge extra keys into this dict (`| {"kind": ...}`), so it must be a plain `dict`.

---

## C. `meshpair/tifxyz.py`

A tifxyz mesh is a folder with `meta.json`, `x.tif`, `y.tif`, `z.tif` (each a single-page 2-D float image of shape
`(rows, cols)`).

- Constants: `S3_HTTP = "https://vesuvius-challenge-open-data.s3.amazonaws.com"`;
  `FILES = ("meta.json", "x.tif", "y.tif", "z.tif")` (this order).
- `class Mesh(points, meta, path="")`: attributes `points` (float32 `(rows, cols, 3)`, invalid = -1), `meta` (dict),
  `path` (str).
  - property `valid` → bool array `(rows, cols)`: point validity of B.1 (x != -1 and all components finite).
  - `geom(px_um=NaN)` → `MeshGeom(rows, cols, scale_x = float(meta["scale"][0]), scale_y = float(meta["scale"][1]),
    px_um = px_um)`, render scale left at its default 1.0. If meta has no `"scale"`, use `[0.05, 0.05]`.
- `read_tifxyz(path, vc_load=False)` → `Mesh`: parse `meta.json` as JSON text; read `x.tif`, `y.tif`, `z.tif` with
  `tifffile.imread`, convert each to float32, stack on a new last axis in the order x, y, z. `Mesh.path` = `str(path)`.
  - **vc_load rule** (VC3D's loader, `load_quad_from_tifxyz`): when `vc_load=True`, every vertex whose **z <= 0** gets all
    three components set to -1.0. A NaN z is not `<= 0` and is left unchanged; z = -1 and z = 0 both become (-1,-1,-1).
    Applied after the float32 stack.
- `write_tifxyz(path, points, meta)`: create the folder (parents, no error if it exists); write each component
  `points[..., i]` as a C-contiguous float32 2-D image with `tifffile.imwrite` to `x.tif`, `y.tif`, `z.tif`; write
  `meta.json` as `json.dumps(meta, indent=1)`. (Tests write meshes with this and read them back through `canvasmap`,
  which reads the scale from meta.json and the shape from the first page of `x.tif`.)
- `fetch_tifxyz(s3_key_prefix, cache_root, base=S3_HTTP)` → `(local_dir: Path, shas: dict)`: strip leading/trailing "/"
  from the prefix; `local_dir = cache_root / prefix` (created). For each name in `FILES`: if the local file is missing or
  empty, download `"{base}/{prefix}/{name}"` (timeout 120 s, read in 1 MiB chunks) into `<name>.part` next to it, then
  atomically rename to `<name>`. Then `shas[name]` = SHA-256 hex digest of the file's bytes (also for files already
  cached). Nothing downloaded is executed.

---

## D. Slice centre (used by `meshpair/ct.py`)

Fact (renderer): a surface volume with `N` layers (slices) along the mesh normal, one level-0 voxel per layer, has its
layer at the mesh surface at the fractional index `c(N) = (N - 1) / 2` (float64; exact for every integer N). Layer `k`
therefore sits at signed depth `fl64((k - c(N)) * px)` µm along the normal, with `px` the level-0 voxel size of that
volume. For an even N the surface falls halfway between layers N/2 - 1 and N/2.

Uses inside `ct.py` (the ct.py integrator wires these; the implementer provides the helper(s)):
- coarse layer `k` of an `n_c`-layer coarse volume is compared at depth `fl64(((k - c(n_c)) * px_c) + delta_um)`;
- the fine volume's fractional level-0 layer for a depth `d` is `fl64((d / px_f) + c(n_f0))`, with `n_f0` the fine
  volume's layer count at level 0;
- the depth range that the fine block must cover uses `(min(ks) - c(n_c)) * px_c` and `(max(ks) - c(n_c)) * px_c`.

`coarse_layers(n_c, count=None)` → sorted list of ints: the default `count` is 8 if `n_c` is even, else 7. Rank all
layers `k = 0 .. n_c-1` by distance `|k - c(n_c)|`, ties broken by the **smaller k first**; take the first `count`
(all of them if `count >= n_c`); return them sorted ascending.

Vectors: `coarse_layers(12)` = [2..9]; `(7)` = [0..6]; `(13)` = [3..9]; `(12, 3)` = [4, 5, 6]; `(5, 9)` = [0..4];
`(1)` = [0]; `(2)` = [0, 1].

(The pyramid-level rule "level-l pixel i is centred on level-0 pixel `2^l i + (2^l - 1)/2`" used by `ct.py` and
`canvasmap.py` is this project's own block-mean assumption, not part of this spec.)

---

## 5. Interface contract (who imports what)

| Importer | Names it needs |
|---|---|
| `meshpair/canvasmap.py` | `MeshGeom`, `exact_canvas_map`, `exact_canvas_map_inv`, `exact_linear_coeffs` |
| `meshpair/meshcheck.py` | `conventions` as a module (`offset_summary`, `naive_canvas_map`, `exact_canvas_map`); `compare_meshes`, `transform_tifxyz`; `Mesh` |
| `meshpair/fit.py`, `reverse.py`, `reverse_lc.py` | `MeshGeom`, `canvas_to_grid`, `grid_to_canvas`, `exact_canvas_map`, `exact_canvas_map_inv` |
| `meshpair/synth.py` | `F32`, `apply_affine_f32`, `resample_linear_preserving_invalids`, `transform_tifxyz` |
| `check.py` | `conventions` (`exact_linear_coeffs`, `exact_canvas_map`), `read_tifxyz(..., vc_load=True)`, `fetch_tifxyz` |
| `meshpair/ct.py` | `coarse_layers` and the slice-centre rule of §D |
| tests | see §6 |

---

## 6. Behaviours the existing tests depend on

Public tests: `tests/test_meshpair.py`, `tests/test_release.py`, `tests/test_reverse.py`, `tests/test_reverse_lc.py`,
`tests/test_plants_a1.py` (the last one does not touch these modules).

| Test | Behaviour it checks |
|---|---|
| `test_canvas_grid_roundtrip_and_renderer_rule` | `grid_to_canvas(canvas_to_grid(p))` ≈ p (allclose) on both axes for `MeshGeom(211, 238)`, p = [0, 17.3, 4759]; with s = 0.05, R = 1, `canvas_to_grid(p, 1)` ≈ `0.05 (p + 0.5)` |
| `test_exact_map_inverse_and_w00_numbers` | inverse map round-trips (allclose) for F = (823, 928), C = (211, 238); `exact_linear_coeffs(F, C, 0)`: b = 211/823 within 1e-12, a = 10.5·211/823 - 10.5 within 1e-9 |
| `test_resampler_matches_independent_float64_halfpixel_bilinear` | B.5 on a 60×70 sheet → 233×271 equals a float64 half-pixel, clamped bilinear within 2e-3 |
| `test_transform_size_rule_and_invalid_propagation` | `new_shape` = round(float32(n) · float32(ratio)); identity affine → `resampled` False and output **array_equal** to input; invalid count grows by the 4-corner rule |
| `test_route_I_*`, `test_reproduction_flags_a_planted_file` (via `meshcheck.check_pair`) | `offset_summary` keys `y_offset_box_centre_px` / `x_offset_box_centre_px`; `exact_canvas_map`; `transform_tifxyz` reproduces a mesh it made itself (M2a: median error ≤ 0.01, p99 ≤ 0.05, mask agreement ≥ 0.999) and `compare_meshes` keys `same_shape`, `err_median_vx`, `err_p99_vx`, `mask_agreement`; `info["new_shape"]`; `apply_affine_f32` |
| `test_ct_*` | `ct.py` with §D (coarse layers, depth of layer k); planted depth recovered within 0.2 coarse layer |
| `tests/test_reverse.py` (all) | `F32`, `transform_tifxyz` (`info["new_shape"]`), `apply_affine_f32`; `Mesh.valid`, `Mesh.geom` |
| `tests/test_reverse_lc.py` (all; uses `reverse_route/synth_a1.py`) | `transform_tifxyz`, `apply_affine_f32`, `F32`; `offset_summary(gF, gC, k=..., valid_box=...)` and `grid_to_canvas` on scalar ints |
| `test_helper_matches_published_formula` | `write_tifxyz`; map = (p + 10.5)·n_C/n_F - 10.5 within 1e-9 per axis |
| `test_helper_round_trip_and_levels` | inverse within 1e-8 with a float32-0.05 fine scale |
| `test_helper_reads_float32_scale` | `MeshGeom(926, 738, F32_005, F32_005)` positional; `scale_x == float(float32(0.05))` read from meta; `linear_coeffs` (via `exact_linear_coeffs`) **exactly equal** (`==`) to a second evaluation |
| `test_warp_moves_content_to_mapped_place` | `canvas_h`, `canvas_w`, maps both ways |
| `test_pair_maps_table_reproduces_records` | for all 207 records in `out/mesh_records.jsonl` + `out/dev_w00_mesh.jsonl`, `exact_linear_coeffs` on `MeshGeom(fine_shape, sf, sf)` / `MeshGeom(coarse_shape, sc, sc)` for some sf, sc ∈ {0.05, float32 0.05} gives `y_a, y_b, x_a, x_b` with **difference exactly 0.0** |
| `test_canvasmap_cli` | the CLI prints `linear_coeffs` as JSON |

---

## 7. Acceptance vectors (computed 2026-10-18 with the reference code; NumPy 2.4.3, SciPy 1.17.1, Windows x64)

Float values are Python `repr` (shortest round-trip). SHA-256 per §0 byte layout.

### 7.1 Conventions

`exact_linear_coeffs` for fine (823, 928) → coarse (211, 238), R = 1, returned as `[y_a, y_b, x_a, x_b]`:

| s_fine | s_coarse | y_a | y_b | x_a | x_b |
|---|---|---|---|---|---|
| 0.05 | 0.05 | -7.808019441069064 | 0.2563791008501539 | -7.80711206896558 | 0.25646551724184974 |
| 0.05 | f32 0.05 | -7.807987890721506 | 0.2563790970298214 | -7.807076495318142 | 0.256465513420153 |
| f32 0.05 | 0.05 | -7.808050880609244 | 0.25637910467094116 | -7.807147531818373 | 0.25646552106309173 |
| f32 0.05 | f32 0.05 | -7.8080193302612315 | 0.25637910085060867 | -7.807111958170481 | 0.256465517241395 |

(f32 0.05 = 0.05000000074505806.) Canvas size of (823, 928): `canvas_h = 16460`, `canvas_w = 18560` with either scale.

`MeshGeom(211, 238)` (R = 1, s = 0.05; canvas 4220 × 4760), p = [0.0, 17.3, 4759.0]:
- `canvas_to_grid(p, 1)` = [0.024999999999991473, 0.8900000000000006, 237.97500000000002]
- `canvas_to_grid(p, 0)` = [0.024999999999991473, 0.8900000000000006, 237.975]  (the last element differs only through W; this is an order-of-operations check)

`exact_canvas_map(5000.0, MeshGeom(823, 928), MeshGeom(211, 238), 0)` = 1274.0874848116646;
`exact_canvas_map_inv(1270.0, same, axis 1)` = 4982.373949579833.

Canvas sizes: `MeshGeom(10, 10, 0.05, 0.05, render_scale=0.5).canvas_w` = 100; `MeshGeom(1, 1, 0.05, 0.05,
render_scale=0.0).canvas_w` = 1 (the max(1, ·) floor); `MeshGeom(7, 5, 0.04, 0.04).canvas_w` = 125.
`lround`: 2.5 → 3, -2.5 → -3, 3.5 → 4.

`offset_summary(MeshGeom(823, 928, px_um=2.403), MeshGeom(211, 238, px_um=9.366))` (default k and box):
```
k_naive 0.25656630365150546
y_a -7.808019441069064            y_b 0.2563791008501539          y_b_minus_k -0.0001872028013515381
y_offset_at_origin_px -7.808019441069064   y_offset_box_min_px -10.889190348513694
y_offset_box_max_px -7.808019441069064     y_offset_box_centre_px -9.348604894791606
x_a -7.80711206896558             x_b 0.25646551724184974         x_b_minus_k -0.00010078640965571761
x_offset_at_origin_px -7.80711206896558    x_offset_box_min_px -9.677607045766308
x_offset_box_max_px -7.80711206896558      x_offset_box_centre_px -8.742359557365944
y_offset_at_origin_um -73.12991008505284   y_offset_box_min_um -101.98815680417925
y_offset_box_max_um -73.12991008505284     y_offset_box_centre_um -87.55903344461818
x_offset_at_origin_um -73.12141163793162   x_offset_box_min_um -90.64046759064723
x_offset_box_max_um -73.12141163793162     x_offset_box_centre_um -81.88093961428943
```
Same meshes without px_um and `k=0.25`: 14 keys, `y_b_minus_k` 0.006379100850153918, `y_offset_box_max_px`
97.18560145161427, `y_offset_box_centre_px` 44.68879100527238, `x_offset_box_max_px` 112.18642242252372,
`x_offset_box_centre_px` 52.18965517677907, no `_um` keys.

**Full-record parity (stronger than the shipped test; checked 2026-10-18):** for each of the 207 records, build
`MeshGeom(*fine_shape, sf, sf, px_um=px_fine)` and `MeshGeom(*coarse_shape, sc, sc, px_um=px_coarse)` and call
`offset_summary(gF, gC, px_fine / px_coarse, valid_box=tuple(fine_canvas_valid_box))`. For every record one (sf, sc)
combination reproduces the stored `exporter_model` dict **exactly** (same keys, same order, every value `==`): 207 / 207.
Split: (f32, 0.05) 72 records, (f32, f32) 126, (0.05, f32) 9.

### 7.2 Exporter

Inputs come from the repo's own (non-GPL) `meshpair/synth.py` (`curved_sheet`, `demo_affine`), which the implementer
**runs** but need not read; check the input SHA first.

| id | recipe | expected |
|---|---|---|
| S80 | `curved_sheet(rows=80, cols=90, seed=1)` | input SHA `d640295bbb43a0cf86688a5f2274031c76c5744b3b7fd96d20891941eacac4c9` |
| SP1 | `measure_grid_axis_spacing(S80)` | (20.00825044360449, 20.00935154363681, 271, 270) |
| T1 | `transform_tifxyz(S80, demo_affine())` | out (312, 351) SHA `1df4347843d55fe734a2e47f48c02f63dd6fd77bfca21b4033ad2ea2ef806b08`, 11165 vertices with x = -1; spacing_after (77.99673724418912, 77.97294373053788); ratio_x 3.898228756383859; ratio_y 3.89682512001915; resampled True; new_shape (312, 351) |
| T-id | `transform_tifxyz(S80, [I \| 0])` | output identical to S80 (same SHA); ratios 1.0, 1.0; resampled False; no new_shape |
| T-win | `transform_tifxyz(S80, [1.0005·I \| 0])` | ratio_x 1.000502824758029, ratio_y 1.000527017154908; resampled False (float32 0.001 window), shape (80, 90), SHA `aa8d590a495b330fd85464bbb7ef9be4395d9dcf2e31e1ce9aa806b8d7088ec1` |
| T-win2 | `transform_tifxyz(S80, [1.002·I \| 0])` | ratio_x 1.0019980692007522, ratio_y 1.0020020491101451; resampled True, new_shape (80, 90) (same size, still resampled); 641 invalid; SHA `eff60163cd8b709b19f2d14791e7269837bd47c13db7353b56e93d414e3ae1e8` |
| S60 | `curved_sheet(rows=60, cols=70, holes=0, seed=3)`, then rows 0-1 set to row 2, cols 0-3 set to col 4 (as in the resampler test) | input SHA `1c9f80ffcacaf50a7081557d737fb04635f332a0465b7366aa107c05c57b7335` |
| R1 | `resample_linear_preserving_invalids(S60, 233, 271)` | SHA `b0c4a8d3f3022d31c51b364e7339964ba5bf2f817091f048b6cf12c98d801109` |
| R2 | `resample_linear_preserving_invalids(S80, 233, 271)` | SHA `fe29e885c128671f071f47c00b161040b8ef1b8297f93c78b8efb3bd453616c3`, 6376 invalid |

Formula-defined vector V1 (no synth needed): grid 6 × 7, for row r and col c (float64, then cast to float32):
`x = c·20 + 0.1·r·r`, `y = r·20 + 0.05·c`, `z = 1000 + 0.3·r·c`; then vertex (3, 4) set to (-1, -1, -1) and the y component
of vertex (1, 2) set to -1. Input SHA `5b8b07e9a68c09a31681ebd3e964b1bee7548791cb627ad79ec17770c4d0b14e`.
A = `[[3.9, 0, 0, 10.25], [0, 3.9, 0, -3.5], [0, 0, 3.9, 7.0]]`.

| id | recipe | expected |
|---|---|---|
| SP-V1 | `measure_grid_axis_spacing(V1)` | (20.02816521491387, 20.037216479929768, 28, 28) |
| AF1 | `apply_affine_f32(V1, A, 2.0, 0.5)` | SHA `a60dc5dac8cea5e87d02503454750df836d50c5e75fa152573224be1d432cb1c`; vertex (0,0) = [5.125, -1.75, 3903.5]; (5,6) = [482.875, 389.4200134277344, 3938.60009765625]; (1,2) (y was -1, still valid) = [161.51499938964844, -5.650000095367432, 3905.83984375] |
| TV1 | `transform_tifxyz(V1, A)` | out (23, 27), SHA `7492b7162611290571e302dc63c92a104a129914f034b37401709ce411e8cdf2`, 56 invalid; ratio_x 3.8999999112942025, ratio_y 3.899999914065268; new_shape (23, 27); vertex (0,0) = [10.25, -3.5, 3907.0], (22,26) = [488.0, 387.6700134277344, 3942.10009765625] |
| TL | `transform_tifxyz(L, [1.5·I \| 0])`, L = one row of 3 vertices (0,0,0), (20,0,0), (40,0,0) float32 | spacing NaN (rows < 2) → fallback 1.5 on both axes (uniform SVD); new_shape **(2, 4)** (1.5 → 2; 3 × 1.5 = 4.5 → **4**, half to even); both rows = x [0.0, 18.75, 41.25, 60.0], y = z = 0 |
| TL-nu | same L with `diag(1.5, 1.5, 1.6)` | non-uniform SVD → fallback 1.0 → ratios 1.0, resampled False, shape (1, 3) |
| AF-edge | `apply_affine_f32` on one row [(1,2,3), (NaN,2,3), (5,-1,7), (-1,0,0)] with `[2·I \| 0]` | [(2,4,6), (NaN,2,3), (10,-2,14), (-1,0,0)] |
| AF-m1 | `apply_affine_f32([[(2,0,0)]], [I \| (-3,0,0)])` | [(-1,-1,-1)] (new x = -1 → invalid) |
| CM | `compare_meshes(S80, S80[:-1])` | `{pred_shape: [80, 90], pub_shape: [79, 90], same_shape: False}` |
| CM2 | `compare_meshes(X, X)` for T-win2's output | n_valid 6559 (pub, pred, both), mask_mismatch 0, mask_agreement 1.0, errors 0.0 |

### 7.3 tifxyz

`read_tifxyz(dir, vc_load=True)` on a mesh whose z values include 0, -1, 5 and NaN: the vertices with z = 0 and z = -1
become (-1, -1, -1); z = 5 and z = NaN are unchanged. `Mesh.geom()` of a mesh with `"scale": [0.04, 0.05]` and shape
(r, c) gives `scale_x = 0.04` (cols), `scale_y = 0.05` (rows).

---

## 8. Things that cannot be pinned down without naming a library primitive

These are described by naming the NumPy primitive rather than by an independent formula, because their last bit depends
on the library's internal reduction order:
1. The float64 affine product (B.4 step 4): NumPy `@` (BLAS) on an N×3 by 3×3 product.
2. The 3-term Euclidean norm (B.3, B.7): NumPy `linalg.norm` along the last axis (sum of squares in order x, y, z, then
   sqrt).
3. The singular values (B.6 step 5): NumPy `linalg.svd(compute_uv=False)`; only matters at the 1e-4 threshold.
4. `median` and `percentile` in `compare_meshes`: NumPy defaults.

A re-implementation that uses these same NumPy calls is bit-identical on the reference machine; a different BLAS can
move (1) by one float64 ulp in rare cases. The 207-pair parity (§7.1) uses only Part A, which has no such dependency.
