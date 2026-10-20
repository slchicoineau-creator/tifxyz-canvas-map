"""Map pixels between the fine and coarse surface volumes of one segment.

The organizers publish, per segment, a fine (~2.2-2.4 um) and a coarse (~7.9-9.4 um) tifxyz mesh and a surface volume
rendered along each. Both meshes come from one source mesh through `vc_transform_geom` (half-pixel resample, grid size
from the affine's measured scale) and are rendered by `vc_render_tifxyz`. The same papyrus point therefore sits at,
per axis (y = rows, x = cols),

    p_coarse = (p_fine + 10.5) * n_coarse / n_fine - 10.5      (grid step s = 0.05, render --scale 1;
                                                                  n = the two meshes' grid sizes on that axis)

and NOT at p_fine * (fine voxel / coarse voxel). This module gives that map from the two mesh folders.

What it does not do: it does not correct the scan-to-scan transform itself. After this map, the two surface volumes
still differ by the residual of the organizers' affine (about 0.5-1.5 coarse px on the segments we measured).

Pyramid levels: level l of a surface volume is a 2^l block mean in-plane, so level-l pixel i is centred on level-0
pixel 2^l * i + (2^l - 1) / 2 (the same convention meshpair.ct uses).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .conventions import MeshGeom, exact_canvas_map, exact_canvas_map_inv, exact_linear_coeffs

__all__ = ["mesh_geom", "fine_to_coarse_px", "coarse_to_fine_px", "linear_coeffs", "warp_canvas_image"]


def mesh_geom(mesh, render_scale: float = 1.0) -> MeshGeom:
    """MeshGeom of a tifxyz folder (reads meta.json's "scale" and x.tif's shape only), or pass a MeshGeom through.

    The scale is read, not assumed: published meshes store 0.05 either as a double or as float32 (0.0500000007),
    and the map uses the stored value."""
    if isinstance(mesh, MeshGeom):
        return mesh
    import tifffile
    path = Path(mesh)
    meta = json.loads((path / "meta.json").read_text())
    sx, sy = meta.get("scale", [0.05, 0.05])
    with tifffile.TiffFile(path / "x.tif") as tf:
        rows, cols = (int(v) for v in tf.pages[0].shape[:2])
    return MeshGeom(rows=rows, cols=cols, scale_x=float(sx), scale_y=float(sy), render_scale=render_scale)


def _to_level0(p, level: int):
    f = 2 ** level
    return np.asarray(p, dtype=np.float64) * f + (f - 1) / 2.0


def _from_level0(p0, level: int):
    f = 2 ** level
    return (np.asarray(p0, dtype=np.float64) - (f - 1) / 2.0) / f


def fine_to_coarse_px(y, x, fine, coarse, fine_level: int = 0, coarse_level: int = 0):
    """Fine surface-volume pixel (y, x) -> coarse surface-volume pixel of the same mesh point.

    fine, coarse: tifxyz folders (the `mesh/<id>-on-<volume>.tifxyz` the surface volumes were rendered from) or
    MeshGeom objects. *_level: pyramid level of the surface volume the pixel indices refer to (0 = full resolution)."""
    gF, gC = mesh_geom(fine), mesh_geom(coarse)
    yc = exact_canvas_map(_to_level0(y, fine_level), gF, gC, 0)
    xc = exact_canvas_map(_to_level0(x, fine_level), gF, gC, 1)
    return _from_level0(yc, coarse_level), _from_level0(xc, coarse_level)


def coarse_to_fine_px(y, x, fine, coarse, fine_level: int = 0, coarse_level: int = 0):
    """Inverse of fine_to_coarse_px."""
    gF, gC = mesh_geom(fine), mesh_geom(coarse)
    yf = exact_canvas_map_inv(_to_level0(y, coarse_level), gF, gC, 0)
    xf = exact_canvas_map_inv(_to_level0(x, coarse_level), gF, gC, 1)
    return _from_level0(yf, fine_level), _from_level0(xf, fine_level)


def linear_coeffs(fine, coarse):
    """Level-0 map as p_coarse = a + b * p_fine per axis: returns dict(y_a, y_b, x_a, x_b)."""
    gF, gC = mesh_geom(fine), mesh_geom(coarse)
    ya, yb = exact_linear_coeffs(gF, gC, 0)
    xa, xb = exact_linear_coeffs(gF, gC, 1)
    return dict(y_a=ya, y_b=yb, x_a=xa, x_b=xb)


def warp_canvas_image(img, fine, coarse, direction: str = "fine_to_coarse", out_shape=None,
                      src_level: int = 0, dst_level: int = 0, order: int = 1, cval: float = 0.0):
    """Resample a 2-D image (label mask, ink prediction, one surface-volume layer) from one canvas onto the other.

    direction: "fine_to_coarse" (img is on the fine canvas) or "coarse_to_fine".
    out_shape: (H, W) of the destination canvas at dst_level; default = the destination mesh's full canvas at that level.
    order: 0 for label masks (nearest), 1 for continuous maps (bilinear). Pixels that fall outside img get cval.
    Memory: the index grids cover the whole destination canvas in float64. For coarse -> fine on large segments use
    dst_level=2 (or pass out_shape); a 926 x 738 fine mesh at level 0 is ~273 M pixels, about 10 GB."""
    from scipy.ndimage import map_coordinates
    gF, gC = mesh_geom(fine), mesh_geom(coarse)
    dst = gC if direction == "fine_to_coarse" else gF
    if out_shape is None:
        f = 2 ** dst_level
        out_shape = (int(np.ceil(dst.canvas_h / f)), int(np.ceil(dst.canvas_w / f)))
    yy, xx = np.meshgrid(np.arange(out_shape[0], dtype=np.float64), np.arange(out_shape[1], dtype=np.float64),
                         indexing="ij")
    if direction == "fine_to_coarse":
        sy, sx = coarse_to_fine_px(yy, xx, gF, gC, fine_level=src_level, coarse_level=dst_level)
    elif direction == "coarse_to_fine":
        sy, sx = fine_to_coarse_px(yy, xx, gF, gC, fine_level=dst_level, coarse_level=src_level)
    else:
        raise ValueError("direction must be 'fine_to_coarse' or 'coarse_to_fine'")
    img = np.asarray(img)
    return map_coordinates(img.astype(np.float32 if order else img.dtype), [sy, sx], order=order, mode="constant",
                           cval=cval)


def _cli(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="Map fine <-> coarse surface-volume pixels of one segment.")
    ap.add_argument("fine", help="fine tifxyz folder"); ap.add_argument("coarse", help="coarse tifxyz folder")
    ap.add_argument("--fine-px", nargs=2, type=float, metavar=("Y", "X"), help="fine level-0 pixel to map")
    a = ap.parse_args(argv)
    c = linear_coeffs(a.fine, a.coarse)
    print(json.dumps(c, indent=1))
    if a.fine_px:
        print("coarse px:", [float(v) for v in fine_to_coarse_px(*a.fine_px, a.fine, a.coarse)])


if __name__ == "__main__":      # python -m meshpair.canvasmap FINE_DIR COARSE_DIR [--fine-px Y X]
    _cli()
