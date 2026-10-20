"""Self-tests of the lattice-clean route (amendment 1; synthetic only)."""
import os, sys
import numpy as np
import pytest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "reverse_route"))
from meshpair.reverse_lc import clean_axis, run_pair        # noqa: E402
from meshpair.reverse import reverse_route                  # noqa: E402
from meshpair.tifxyz import Mesh                            # noqa: E402
from synth_a1 import CONFIGS, META, PXF, PXC, jittered_source, make_pair   # noqa: E402


def _run(S, name, shift=(0.0, 0.0)):
    AF, AC = CONFIGS[name]
    F, C, M, fc, cc = make_pair(S, AF, AC, shift)
    Fm, Cm = Mesh(F, META), Mesh(C, META)
    return run_pair(F, Fm.valid, C, Cm.valid, Fm.geom(PXF), Cm.geom(PXC), M, S.shape[:2], PXF, PXC,
                    coarse_is_copy=cc, fine_is_copy=fc), (F, C, M, Fm, Cm)


def test_clean_axis_copy_and_ratio4():
    g = np.arange(0, 50) + 0.3
    assert clean_axis(g, 100, 100).all()                      # copy: cells = source cells
    c = clean_axis(np.arange(0, 400) + 0.5, 400, 100)         # ratio 4: one cell in 4 straddles
    assert abs(c.mean() - 0.75) < 0.01


@pytest.mark.parametrize("name", list(CONFIGS))
def test_clean_pass(name):
    r, _ = _run(jittered_source(3, jitter=2.0), name)
    assert r["verdict"] == "PASS", r


def test_original_rule_fails_on_jittered_ratio4():
    """The diagnosed cause reproduced: jittered source, both transformed, ratio 4 -> the frozen reverse rule fails."""
    r, (F, C, M, Fm, Cm) = _run(jittered_source(4, jitter=2.0), "both_ratio4")
    o = reverse_route(F, Fm.valid, C, Cm.valid, Fm.geom(PXF), Cm.geom(PXC), M, PXF, PXC)
    assert not o["pair_pass"] and r["verdict"] == "PASS"


@pytest.mark.parametrize("name", list(CONFIGS))
def test_planted_005_fails(name):
    r, _ = _run(jittered_source(5, jitter=1.0), name, shift=(0.05, 0.0))
    assert r["verdict"] == "FAIL", r


@pytest.mark.parametrize("name", list(CONFIGS))
def test_offset_matches_exporter_model(name):
    """The measured offset vs naive at the fine valid-box centre equals the exporter model's (conventions.offset_summary)."""
    from meshpair import conventions as cv
    r, (F, C, M, Fm, Cm) = _run(jittered_source(6, jitter=1.0), name)
    gF, gC = Fm.geom(PXF), Cm.geom(PXC)
    rr, cc = np.nonzero(Fm.valid)
    box = (gF.grid_to_canvas(rr.min(), 0), gF.grid_to_canvas(rr.max(), 0), gF.grid_to_canvas(cc.min(), 1), gF.grid_to_canvas(cc.max(), 1))
    m = cv.offset_summary(gF, gC, k=PXF / PXC, valid_box=box)
    oy, ox = r["offset_vs_naive_centre_px"]
    assert abs(oy - m["y_offset_box_centre_px"]) < 0.02 and abs(ox - m["x_offset_box_centre_px"]) < 0.02, (r["offset_vs_naive_centre_px"], m)
