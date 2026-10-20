"""Build results/pair_maps.tsv: the fine -> coarse canvas map of every published mesh pair, from the stored records.

Per pair and axis (y = rows, x = cols): p_coarse = a + b * p_fine, level-0 canvas pixels, render --scale 1.
Sources: out/mesh_records.jsonl (206 catalog pairs), out/dev_w00_mesh.jsonl (PHerc0841 w00, the development pair),
reverse_route/apply154/pairs_table.jsonl (lattice-clean verdict of the 206 catalog pairs).
Catalog: metadata.json sha256 7ec8dcaa...abdd3 (2026-10-11). A re-export of a mesh changes its row; the helper
(meshpair.canvasmap) run on the mesh folders is always current.

    python scripts/make_pair_maps.py            # writes results/pair_maps.tsv and prints the per-family offsets
"""
import json
import math
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def jl(p):
    return [json.loads(l) for l in (ROOT / p).read_text().splitlines() if l.strip()]


def main():
    recs = jl("out/mesh_records.jsonl") + jl("out/dev_w00_mesh.jsonl")
    lc = {(r["sample"], r["segment"]): r for r in jl("reverse_route/apply154/pairs_table.jsonl")}
    cols = ["sample", "segment", "long_id", "fine_um", "coarse_um", "fine_rows", "fine_cols", "coarse_rows", "coarse_cols",
            "y_a", "y_b", "x_a", "x_b", "naive_k", "offset_centre_um", "offset_centre_y_um", "offset_centre_x_um",
            "forward_M1", "lattice_clean_verdict", "fine_key", "coarse_key"]
    rows, fam = [], {}
    for r in recs:
        em = r["exporter_model"]
        oy, ox = em["y_offset_box_centre_um"], em["x_offset_box_centre_um"]
        off = math.hypot(oy, ox)
        dev = bool(r.get("dev_pair"))
        v = lc.get((r["sample"], r["segment"]), {}).get("verdict", "dev pair (not in the catalog runs)" if dev else "")
        rows.append([r["sample"], r["segment"], r["long_id"], r["px_fine"], r["px_coarse"], *r["fine_shape"][:2],
                     *r["coarse_shape"][:2], f'{em["y_a"]:.9f}', f'{em["y_b"]:.12f}', f'{em["x_a"]:.9f}',
                     f'{em["x_b"]:.12f}', f'{em["k_naive"]:.12f}', f"{off:.1f}", f"{oy:.1f}", f"{ox:.1f}",
                     (r.get("lines") or {}).get("M1_in_files", ""), v, r["fine_key"], r["coarse_key"]])
        if not dev:
            fam.setdefault(r["sample"], []).append(off)
    out = ROOT / "results" / "pair_maps.tsv"
    out.parent.mkdir(exist_ok=True)
    with open(out, "w", newline="\n") as fo:
        fo.write("\t".join(cols) + "\n")
        for row in rows:
            fo.write("\t".join(str(x) for x in row) + "\n")
    print(len(rows), "rows ->", out.relative_to(ROOT))
    print("offset vs naive index pairing at the valid-box centre (um), min / median / max over segments:")
    for s, v in sorted(fam.items()):
        print(f"  {s:12s} n={len(v):3d}  {min(v):7.0f} / {statistics.median(v):7.0f} / {max(v):7.0f}")


if __name__ == "__main__":
    main()
