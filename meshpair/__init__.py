"""meshpair: mesh-only checker for the (u, v) index relation between the organizers' fine (~2.4 um) and coarse
(~7.9-9.4 um) tifxyz meshes and surface volumes of the same segment, plus an optional CT residual stage.

Modules
- canvasmap:   public helper: map pixels / warp images between a segment's fine and coarse surface volumes.
- conventions: canvas <-> grid conventions of vc_render_tifxyz and the exporter's half-pixel index map.
- exporter:    Python replica of vc_transform_geom (transformSurfacePoints + QuadSurface::resample).
- tifxyz:      read / fetch tifxyz meshes.
- catalog:     enumerate fine/coarse mesh pairs from the open-data metadata.json.
- fit:         model-free index correspondence (inverse map through the published affine; free-affine fit).
- ct:          CT residual shift / depth offset per tile, split-half noise floor, planted controls.
- synth:       synthetic meshes and surface volumes for self-tests.
- zio:         small remote Zarr v2 reader with a local chunk cache (CT route only).
"""
__version__ = "0.2.0"   # 0.1.0 = the version every published record was made with; see CHANGES.md
