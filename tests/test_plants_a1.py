"""Input checks of `check.py plants` (Amendment 1); no network, no scroll data."""
import hashlib, json, os, sys
import pytest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import check  # noqa: E402


def _mk(tmp_path, sha=None, valid=True):
    seg = "S1"
    rp = tmp_path / f"reachable_{seg}.json"
    rp.write_text(json.dumps(dict(segment=seg, tiles=[[1.0, 2.0]], spec={})))
    good = hashlib.sha256(rp.read_bytes()).hexdigest()
    (tmp_path / f"ct_{seg}.log").write_text(f"reachable 1 of 1 sha256 {sha or good}\n")
    rec = dict(kind="tile", map="exact", centre_c=[1.0, 2.0], valid=valid)
    jl = tmp_path / f"ct_{seg}.jsonl"
    jl.write_text(json.dumps(rec) + "\n")
    return seg, jl


def test_requires_prereg(tmp_path):
    seg, _ = _mk(tmp_path)
    for bad in (None, "8ca5789f"):
        with pytest.raises(SystemExit, match="requires --prereg 13df5696"):
            check._a1_inputs(tmp_path, seg, bad)


def test_sha_mismatch_refused(tmp_path):
    seg, _ = _mk(tmp_path, sha="0" * 64)
    with pytest.raises(SystemExit, match="does not match"):
        check._a1_inputs(tmp_path, seg, "13df5696")


def test_ok_and_base_file_untouched(tmp_path):
    seg, jl = _mk(tmp_path)
    before = jl.read_bytes()
    rp, reach, spec, base = check._a1_inputs(tmp_path, seg, "13df5696")
    assert reach == [(1.0, 2.0)] and (1.0, 2.0) in base
    assert jl.read_bytes() == before


def test_invalid_base_refused(tmp_path):
    seg, _ = _mk(tmp_path, valid=False)
    with pytest.raises(SystemExit, match="without a valid base"):
        check._a1_inputs(tmp_path, seg, "13df5696")
