"""CT residual stage: per-tile in-plane shift and depth offset between a fine and a coarse surface volume of one segment,
after a given canvas map (exact exporter map, or the naive index map for comparison).

Conventions
- Canvas px of each surface volume = full-resolution canvas of its own mesh (see conventions.py).
- Coarse layer k sits at depth (k - (n_C - 1)/2) * px_C along the mesh normal; fine layer j at (j - (n_F - 1)/2) * px_F
  (vc_render_tifxyz: center = 0.5 * (numSlices - 1), one level-0 voxel per slice).
- Pyramid level l of a surface volume is assumed to be a 2^l block mean in-plane: level-l pixel i is centred on
  level-0 pixel 2^l * i + (2^l - 1) / 2 (checked by the plumbing test 'level' in the self-test; DESIGN.md 4.3).
- Result sign: NCC peak s means coarse pixel q shows the fine content the map puts at q + s, so the map's correction is
  -s (true p_C = map(p_F) - s). depth_offset_um = delta: the coarse layer at depth d shows the fine content at d + delta.
- Noise floor (split-half check): the tile is re-measured on two disjoint depth sub-slabs and on left/right
  halves; each half carries half the signal, so the whole-tile noise is n = |s_a - s_b| / 2.
- Validity of a peak: r >= 0.30, runner-up (highest distinct local maximum outside a 3 px radius) <= 0.70 r, peak not on the search edge, split spread
  (max of |s_a - s_b|) <= 1.0 px; the autocorrelation-at-shift rule (<= 0.30) is applied only when |s| >= 3 px
  (it guards against periodic aliases of a large shift; a ~1 px residual sits inside the autocorrelation core).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.signal import fftconvolve

R_MIN = 0.30
RUNNER_UP_MAX_FRAC = 0.70
SPLIT_SPREAD_MAX_PX = 1.0
AUTOCORR_MAX = 0.30
AUTOCORR_MIN_SHIFT_PX = 3.0


@dataclass
class SV:
    """A surface volume accessor at one pyramid level.
    read(z0, z1, y0, y1, x0, x1) returns array[z0:z1, y0:y1, x0:x1] at this level (out-of-range filled with 0)."""
    read: callable
    shape: tuple            # (n_layers_at_level, H_level, W_level)
    px_um: float            # level-0 voxel size
    level: int = 0
    n_layers0: int | None = None   # layers at level 0 (defaults to shape[0])

    @property
    def f(self) -> int:
        return 2 ** self.level

    @property
    def zf(self) -> int:
        n0 = self.n_layers0 or self.shape[0]
        return max(1, int(round(n0 / self.shape[0])))


def array_sv(arr: np.ndarray, px_um: float, level: int = 0, n_layers0: int | None = None, origin=(0, 0),
             full_shape=None) -> SV:
    """origin = (y, x) level-l canvas index of arr[:, 0, 0] (for cached crops of a larger surface volume)."""
    oy, ox = origin

    def read(z0, z1, y0, y1, x0, x1):
        y0, y1, x0, x1 = y0 - oy, y1 - oy, x0 - ox, x1 - ox
        out = np.zeros((z1 - z0, y1 - y0, x1 - x0), np.float32)
        Z, H, W = arr.shape
        zs, ys, xs = max(z0, 0), max(y0, 0), max(x0, 0)
        ze, ye, xe = min(z1, Z), min(y1, H), min(x1, W)
        if ze > zs and ye > ys and xe > xs:
            out[zs - z0:ze - z0, ys - y0:ye - y0, xs - x0:xe - x0] = arr[zs:ze, ys:ye, xs:xe]
        return out
    return SV(read=read, shape=full_shape or arr.shape, px_um=px_um, level=level, n_layers0=n_layers0)


def remote_sv(url: str, px_um: float, level: int = 0, n_layers0: int | None = None, workers: int = 16) -> SV:
    from .zio import RemoteArray         # bundled public open-data reader with a local chunk cache
    ra = RemoteArray(url, level)

    def read(z0, z1, y0, y1, x0, x1):
        Z, H, W = ra.shape
        out = np.zeros((z1 - z0, y1 - y0, x1 - x0), np.float32)
        zs, ys, xs = max(z0, 0), max(y0, 0), max(x0, 0)
        ze, ye, xe = min(z1, Z), min(y1, H), min(x1, W)
        if ze > zs and ye > ys and xe > xs:
            out[zs - z0:ze - z0, ys - y0:ye - y0, xs - x0:xe - x0] = ra.read(zs, ze, ys, ye, xs, xe, workers=workers)
        return out
    return SV(read=read, shape=ra.shape, px_um=px_um, level=level, n_layers0=n_layers0)


@dataclass
class CanvasMap:
    """Per-axis affine fine-canvas -> coarse-canvas map p_c = a + b p_f (y and x), plus an optional plant (coarse px)
    added on the coarse side (used for planted controls)."""
    ay: float
    by: float
    ax: float
    bx: float
    plant: tuple = (0.0, 0.0)

    def inv(self, cy, cx):
        return (cy + self.plant[0] - self.ay) / self.by, (cx + self.plant[1] - self.ax) / self.bx

    def planted(self, dy, dx):
        return CanvasMap(self.ay, self.by, self.ax, self.bx, (self.plant[0] + dy, self.plant[1] + dx))


def _dog(a, s1, s2):
    return gaussian_filter(a, s1) - gaussian_filter(a, s2)


from .slicecentre import coarse_layers, slice_centre  # clean-room helpers (cleanroom/SPEC.md Part D)


def _ncc_stack(tmpl: np.ndarray, img: np.ndarray):
    """Masked-free NCC of a (L, T, T) template against every in-plane placement inside img (L, T+2S, T+2S).
    Returns r of shape (2S+1, 2S+1); r[S, S] is zero shift."""
    L, T, Tw = tmpl.shape
    t = tmpl - tmpl.mean()
    tn = np.sqrt((t * t).sum())
    num = np.zeros((img.shape[1] - T + 1, img.shape[2] - Tw + 1))
    s1 = np.zeros_like(num)
    s2 = np.zeros_like(num)
    ones = np.ones((T, Tw))
    for l in range(L):
        num += fftconvolve(img[l], t[l, ::-1, ::-1], mode="valid")
        s1 += fftconvolve(img[l], ones, mode="valid")
        s2 += fftconvolve(img[l] * img[l], ones, mode="valid")
    N = L * T * Tw
    var = np.maximum(s2 - s1 * s1 / N, 1e-12)
    return num / (tn * np.sqrt(var) + 1e-12)


def _subpix(r, i, j):
    def par(a, b, c):
        d = a - 2 * b + c
        return 0.0 if d >= 0 else 0.5 * (a - c) / d
    dy = par(r[i - 1, j], r[i, j], r[i + 1, j]) if 0 < i < r.shape[0] - 1 else 0.0
    dx = par(r[i, j - 1], r[i, j], r[i, j + 1]) if 0 < j < r.shape[1] - 1 else 0.0
    return i + dy, j + dx


def _peak(r, S, excl=3):
    """Best integer peak + parabolic sub-pixel. Two runner-ups outside an excl-px radius:
    runner_up_any = plain maximum there (includes the main peak's own skirt);
    runner_up = highest LOCAL maximum (8-neighbourhood) there, i.e. a distinct second peak. Validity uses runner_up."""
    i, j = np.unravel_index(np.argmax(r), r.shape)
    yy, xx = np.mgrid[:r.shape[0], :r.shape[1]]
    far = (yy - i) ** 2 + (xx - j) ** 2 > excl ** 2
    ru_any = float(r[far].max()) if far.any() else -1.0
    rp = np.pad(r, 1, constant_values=-np.inf)
    nb = np.max(np.stack([rp[1 + a:rp.shape[0] - 1 + a, 1 + b:rp.shape[1] - 1 + b]
                          for a in (-1, 0, 1) for b in (-1, 0, 1) if (a, b) != (0, 0)]), axis=0)
    lmax = (r > nb) & far
    ru = float(r[lmax].max()) if lmax.any() else -1.0
    si, sj = _subpix(r, i, j)
    return dict(r=float(r[i, j]), shift=(si - S, sj - S), int_shift=(int(i - S), int(j - S)), runner_up=ru,
                runner_up_any=ru_any, on_edge=bool(i in (0, r.shape[0] - 1) or j in (0, r.shape[1] - 1)))


@dataclass
class TileSpec:
    T: int = 213                 # coarse px (2 mm at 9.4 um)
    S: int = 16                  # in-plane search radius, coarse px
    depth_range_um: float = 28.0  # +- depth search (3 coarse layers at 9.4 um)
    depth_step_um: float | None = None   # default: fine level-0 voxel
    dog: tuple = (1.0, 6.0)      # coarse px
    n_coarse_layers: int | None = None
    subsamples: int = 4          # depth sub-samples per coarse layer (box PSF along the normal)
    min_valid_frac: float = 0.90  # fraction of non-zero voxels required in both blocks


def _build_fine_on_coarse(fine: SV, cmap: CanvasMap, y0c, x0c, n, ks, n_c, px_c, delta_um, spec: TileSpec, cache):
    """Fine content resampled on the coarse pixel grid [y0c, y0c+n) x [x0c, x0c+n) for coarse layers ks, at depth
    offset delta_um. cache holds the blurred fine block (shared across deltas)."""
    f = fine.f
    if "blk" not in cache:
        yy = np.array([y0c, y0c + n - 1], float)
        xx = np.array([x0c, x0c + n - 1], float)
        fy, fx = cmap.inv(yy, xx)                       # fine level-0 canvas
        pad = 8 * f + 4
        ly0 = int(np.floor((min(fy) - pad - (f - 1) / 2) / f)); ly1 = int(np.ceil((max(fy) + pad) / f)) + 1
        lx0 = int(np.floor((min(fx) - pad - (f - 1) / 2) / f)); lx1 = int(np.ceil((max(fx) + pad) / f)) + 1
        n_f0 = fine.n_layers0 or fine.shape[0]
        c_f0 = slice_centre(n_f0)
        dmin = (min(ks) - slice_centre(n_c)) * px_c - spec.depth_range_um - px_c
        dmax = (max(ks) - slice_centre(n_c)) * px_c + spec.depth_range_um + px_c
        j0 = int(np.floor(dmin / fine.px_um + c_f0)) - 2
        j1 = int(np.ceil(dmax / fine.px_um + c_f0)) + 3
        zf = fine.zf
        lz0, lz1 = max(0, j0 // zf), min(fine.shape[0], j1 // zf + 1)
        blk = fine.read(lz0, lz1, ly0, ly1, lx0, lx1).astype(np.float32)
        sig = 0.42 * px_c / (fine.px_um * f)
        if sig > 0.3:
            blk = gaussian_filter(blk, (0, sig, sig))
        cache.update(blk=blk, lz0=lz0, ly0=ly0, lx0=lx0, c_f0=c_f0, zf=zf,
                     valid=float((blk != 0).mean()))
    blk, lz0, ly0, lx0, c_f0, zf = (cache[k] for k in ("blk", "lz0", "ly0", "lx0", "c_f0", "zf"))
    cy = np.arange(y0c, y0c + n, dtype=np.float64)
    cx = np.arange(x0c, x0c + n, dtype=np.float64)
    fy, _ = cmap.inv(cy, np.zeros_like(cy))
    _, fx = cmap.inv(np.zeros_like(cx), cx)
    ly = (fy - (f - 1) / 2.0) / f - ly0                 # fractional level-l index in blk
    lx = (fx - (f - 1) / 2.0) / f - lx0
    out = np.zeros((len(ks), n, n), np.float32)
    LY, LX = np.meshgrid(ly, lx, indexing="ij")
    m = spec.subsamples
    for li, k in enumerate(ks):
        d = (k - slice_centre(n_c)) * px_c + delta_um
        acc = np.zeros((blk.shape[1], blk.shape[2]), np.float32)
        for q in range(m):
            dq = d + px_c * ((q + 0.5) / m - 0.5)
            j = dq / fine.px_um + c_f0                     # fractional level-0 layer
            jl = (j - (zf - 1) / 2.0) / zf - lz0            # fractional layer in blk
            j0 = int(np.floor(jl)); w = jl - j0
            j0 = min(max(j0, 0), blk.shape[0] - 2)
            acc += (1 - w) * blk[j0] + w * blk[j0 + 1]
        acc /= m
        out[li] = map_coordinates(acc, [LY, LX], order=1, mode="constant", cval=0.0)
    return out


def measure_tile(fine: SV, coarse: SV, cmap: CanvasMap, centre_c, spec: TileSpec = TileSpec(), splits: bool = True,
                 n_c: int | None = None):
    """Measure one tile centred at coarse canvas (yc, xc). Returns a dict (see module doc for signs)."""
    T, S = spec.T, spec.S
    n_c = n_c or coarse.shape[0]
    px_c = coarse.px_um
    ks = coarse_layers(n_c, spec.n_coarse_layers)
    y0 = int(round(centre_c[0] - T / 2)); x0 = int(round(centre_c[1] - T / 2))
    C = coarse.read(min(ks), max(ks) + 1, y0, y0 + T, x0, x0 + T).astype(np.float32)
    C = C[[k - min(ks) for k in ks]]
    cvalid = float((C != 0).mean())
    step = spec.depth_step_um or fine.px_um
    deltas = np.arange(-spec.depth_range_um, spec.depth_range_um + 1e-9, step)
    cache = {}
    Cf = np.stack([_dog(c, *spec.dog) for c in C])
    best = None
    prof = []
    rs = {}
    for di, dl in enumerate(deltas):
        G = _build_fine_on_coarse(fine, cmap, y0 - S, x0 - S, T + 2 * S, ks, n_c, px_c, dl, spec, cache)
        Gf = np.stack([_dog(g, *spec.dog) for g in G])
        r = _ncc_stack(Cf, Gf)
        rs[di] = (r, Gf)
        pk = _peak(r, S)
        prof.append(pk["r"])
        if best is None or pk["r"] > best[1]["r"]:
            best = (di, pk)
    di, pk = best
    prof = np.array(prof)
    if 0 < di < len(deltas) - 1:
        a, b, c = prof[di - 1], prof[di], prof[di + 1]
        d2 = a - 2 * b + c
        ddi = 0.0 if d2 >= 0 else 0.5 * (a - c) / d2
    else:
        ddi = 0.0
    delta = float(deltas[di] + ddi * step)
    res = dict(centre_c=[float(centre_c[0]), float(centre_c[1])], box_c=[y0, y0 + T, x0, x0 + T],
               r=pk["r"], ncc_shift_px=[float(pk["shift"][0]), float(pk["shift"][1])],
               map_correction_px=[float(-pk["shift"][0]), float(-pk["shift"][1])],
               runner_up=pk["runner_up"], runner_up_any=pk["runner_up_any"], on_edge=pk["on_edge"], depth_offset_um=delta,
               depth_offset_layers=delta / px_c, depth_on_edge=bool(di in (0, len(deltas) - 1)),
               coarse_valid_frac=cvalid, fine_valid_frac=cache.get("valid", float("nan")),
               depth_profile_r=[float(x) for x in prof])
    r_best, Gf = rs[di]
    if splits:
        nl = len(ks)
        lo, hi = list(range(nl // 2)), list(range(nl - nl // 2, nl))
        sz = [_peak(_ncc_stack(Cf[sel], Gf[sel]), S)["shift"] for sel in (lo, hi)]
        h = T // 2
        sh = []
        for xs in ((0, h), (T - h, T)):
            tm = Cf[:, :, xs[0]:xs[1]]
            im = Gf[:, :, xs[0]:xs[1] + 2 * S]
            sh.append(_peak(_ncc_stack(tm, im), S)["shift"])
        dz = np.subtract(sz[0], sz[1]); dh = np.subtract(sh[0], sh[1])
        res.update(split_z=[list(map(float, s)) for s in sz], split_h=[list(map(float, s)) for s in sh],
                   n_z=float(np.linalg.norm(dz) / 2), n_h=float(np.linalg.norm(dh) / 2),
                   split_spread_px=float(max(np.abs(dz).max(), np.abs(dh).max())))
        res["noise_px"] = min(res["n_z"], res["n_h"])
    shift_norm = float(np.hypot(*pk["shift"]))
    ac = None
    if shift_norm >= AUTOCORR_MIN_SHIFT_PX:
        Cext = np.pad(Cf, ((0, 0), (S, S), (S, S)))
        ra = _ncc_stack(Cf, Cext)
        iy, ix = pk["int_shift"]
        ac = float(ra[S + iy, S + ix])
    res["autocorr_at_shift"] = ac
    valid = (pk["r"] >= R_MIN and pk["runner_up"] <= RUNNER_UP_MAX_FRAC * pk["r"] and not pk["on_edge"]
             and cvalid >= spec.min_valid_frac and res["fine_valid_frac"] >= spec.min_valid_frac
             and (not splits or res["split_spread_px"] <= SPLIT_SPREAD_MAX_PX)
             and (ac is None or ac <= AUTOCORR_MAX))
    res["valid"] = bool(valid)
    return res


def fit_field(tiles: list[dict], key: str = "map_correction_px"):
    """6-parameter affine fit of the per-tile correction vs tile centre, with leave-one-out residuals (coarse px)."""
    P = np.array([t["centre_c"] for t in tiles], float)
    V = np.array([t[key] for t in tiles], float)
    n = len(P)
    if n < 4:
        return dict(n=n)
    A = np.c_[np.ones(n), P]
    coef, *_ = np.linalg.lstsq(A, V, rcond=None)
    loo = []
    for i in range(n):
        m = np.arange(n) != i
        c, *_ = np.linalg.lstsq(A[m], V[m], rcond=None)
        loo.append(np.linalg.norm(V[i] - A[i] @ c))
    loo = np.array(loo)
    return dict(n=n, coef_y=coef[:, 0].tolist(), coef_x=coef[:, 1].tolist(), loo_p95_px=float(np.percentile(loo, 95)),
                loo_median_px=float(np.median(loo)), insample_p95_px=float(np.percentile(np.linalg.norm(V - A @ coef, axis=1), 95)))
