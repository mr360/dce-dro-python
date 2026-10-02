"""Unit tests for the FTV validation tools.  Run from the repository root:  python -m pytest validation/tests"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "slicer"))

import reference_ftv as ref  # noqa: E402
from ground_truth import expected_outcome, on_threshold_mask  # noqa: E402
from phantom import build_phantom, to_slicer_order  # noqa: E402
from validate_phantom import reference_masks, roi_presets  # noqa: E402


# --------------------------------------------------------------------------------------------- helpers
def voxels(s0, s1, s2, nt=3):
    """4D volume [t, k, j, i] with one 3x3x3 block per (s0, s1, s2) triple along i."""
    s0, s1, s2 = (np.atleast_1d(np.asarray(x, float)) for x in (s0, s1, s2))
    n = s0.size
    vol = np.zeros((nt, 3, 3, 3 * n))
    for idx in range(n):
        sl = (slice(None), slice(None), slice(3 * idx, 3 * idx + 3))
        vol[0][sl] = s0[idx]
        vol[1][sl] = s1[idx]
        vol[2:][(slice(None),) + sl] = s2[idx]
    return vol


def full_roi(vol):
    _, nk, nj, ni = vol.shape
    return {"IJKmin": (0, 0, 0), "IJKmax": (ni, nj, nk)}


def run(vol, **kw):
    kw.setdefault("background_threshold", 0.0)
    kw.setdefault("pe_threshold", 0.0)
    return ref.process(vol, full_roi(vol), pre=0, early=1, late=vol.shape[0] - 1, **kw)


@pytest.fixture(scope="module")
def phantom():
    return build_phantom(n_slices=6)


# ------------------------------------------------------------------------------------ module equations
def test_pe_and_ser_formulas():
    res = run(voxels(100.0, 180.0, 200.0))
    assert np.allclose(res.pe_raw, 100 * 80 / (100 + ref.EPSILON))
    assert np.allclose(res.ser_raw, 80 / (100 + ref.EPSILON))


def test_ser_interval_rule_is_lower_exclusive_upper_inclusive():
    # SER = 0.9 exactly (S1-S0 = 90, Sn-S0 = 100), slightly above, and 3.0 (upper limit, included)
    s0 = 100.0
    vol = voxels([s0, s0, s0, s0], [190, 190.001, 400, 401], [200 - ref.EPSILON] * 2 + [200 - ref.EPSILON] * 2)
    res = run(vol)
    labels = [int(res.ser_map[1, 1, 3 * i + 1]) for i in range(4)]
    # 0.9 -> bin 1 (0 < SER <= 0.9); 0.90001 -> bin 2; 3.0 -> bin 3 (1.3 < SER <= 3.0); 3.01 -> non-SER
    assert labels == [1, 2, 3, 0]
    assert res.counts["FTV0"] == 3 * 27


def test_negative_ser_is_non_ser():
    res = run(voxels(100.0, 150.0, 90.0))   # signal falls below baseline late -> SER < 0
    assert res.counts["FTV0"] == 0
    assert res.counts["ETV"] == 27           # still enhancing early


def test_pe_threshold_is_inclusive_but_epsilon_lowers_pe():
    res = run(voxels(100.0, 170.0, 200.0), pe_threshold=70.0)
    # PE = 100*70/(100+1e-6) = 69.99999930 < 70: excluded because of EPSILON in the denominator
    assert res.counts["FTV0"] == 0
    res = run(voxels(128.0, 128.0 * 1.75, 300.0), pe_threshold=75.0)
    assert res.counts["FTV0"] == 0           # same effect for any S0
    res = run(voxels(100.0, 170.001, 200.0), pe_threshold=70.0)
    assert res.counts["FTV0"] == 27


def test_background_threshold_uses_95th_percentile_of_box():
    vol = voxels([10.0, 100.0], [20.0, 200.0], [30.0, 300.0])
    res = run(vol, background_threshold=60.0)
    assert res.background_threshold_value == pytest.approx(0.6 * np.percentile(vol[0], 95))
    assert res.counts["FTV0"] == 27          # dim block (S0 = 10 < 60) removed


def test_etv_ignores_segment_background_and_pe_threshold():
    vol = voxels([10.0, 100.0], [11.0, 200.0], [12.0, 300.0])
    seg = np.zeros(vol.shape[1:], dtype=np.uint8)
    seg[..., 3:] = 1                         # segment = second block only
    res = ref.process(vol, full_roi(vol), segment_mask=seg, pre=0, early=1, late=2,
                      pe_threshold=70.0, background_threshold=60.0)
    assert res.counts["FTV0"] == 27
    assert res.counts["ETV"] == 54


def test_washout_delta_equals_inverse_ser():
    res = run(voxels([100.0] * 3, [200.0] * 3, [300.0, 190.0, 150.0]), use_ser_map=False)
    centre = [res.wo_delta[1, 1, 3 * i + 1] for i in range(3)]
    ser = [res.ser_raw[1, 1, 3 * i + 1] for i in range(3)]
    assert np.allclose(centre, [100 * (1 / s - 1) for s in ser], rtol=1e-6)
    assert [int(res.washout_map[1, 1, 3 * i + 1]) for i in range(3)] == [1, 2, 3]  # persistent/plateau/washout


def test_single_threshold_bins_and_ftvth():
    b = ref.ser_bins(display_ser_range=False, ser_threshold=1.4)
    assert b.intervals == [0.0, 1.4, np.round(1.4 * 13 / 9, 2), 3.0]
    assert b.ftv_threshold == 1.4
    b = ref.ser_bins(display_ser_range=True)
    assert b.intervals == [0.0, 0.9, 1.3, 3.0] and b.ftv_threshold == 0.9
    b = ref.ser_bins(display_ser_range=False, ser_threshold=0.0)
    assert b.intervals == [0.0, 3.0] and b.ftv_threshold == 0


def test_thin_refbox_fails_like_the_module():
    vol = np.ones((3, 2, 6, 6))
    vol[1] = 2.0
    with pytest.raises(ValueError, match="3 voxels"):
        run(vol)


# ------------------------------------------------------------------------------------------ phantom
def test_notebook_cells_found_and_geometry(phantom):
    assert phantom.signal.shape == (528, 528, 6, 8)
    assert (phantom.pre_index, phantom.early_index, phantom.late_index) == (0, 1, 7)


def test_float_phantom_recovers_encoded_values(phantom):
    """On the noiseless floating-point signal the module equations return the encoded PE and SER."""
    vol = to_slicer_order(phantom.signal)
    roi = roi_presets(phantom)["whole_phantom"]
    res = ref.process(vol, roi, pre=0, early=1, late=7, pe_threshold=0.0, background_threshold=0.0)
    sl = ref.roi_slices(roi)
    pe_t = np.broadcast_to(phantom.target_pe, (6,) + phantom.target_pe.shape)[sl]
    ser_t = np.broadcast_to(phantom.target_ser, (6,) + phantom.target_ser.shape)[sl]
    enh = pe_t > 0
    assert np.allclose(res.pe_raw[enh], pe_t[enh], rtol=1e-6)
    assert np.allclose(res.ser_raw[enh], ser_t[enh], rtol=1e-5)   # notebook adds 1e-6 in s2_ser


@pytest.mark.parametrize("kw", [{}, {"pe_threshold": 50.0}, {"display_ser_range": False, "ser_threshold": 1.0},
                                {"use_ser_map": False, "pe_threshold": 50.0}])
def test_dicom_phantom_matches_ground_truth_off_threshold(phantom, kw):
    vol = to_slicer_order(phantom.saved_uint16)
    roi = roi_presets(phantom)["breasts"]
    res = ref.process(vol, roi, pre=0, early=1, late=7, **kw)
    exp = expected_outcome(phantom, roi, **kw)
    off = ~on_threshold_mask(phantom, roi, **kw)
    label = res.ser_map if res.params["use_ser_map"] else res.washout_map
    assert np.array_equal(exp["label_map"][off], label[off].astype(int))
    for key, mask in reference_masks(res).items():
        assert np.array_equal(exp["masks"][key][off], mask[off]), key


# --------------------------------------------------------------------------------- compare script
def fake_export(tmp_path, phantom, corrupt=False):
    """Write an export folder as export_plugin_results.py would, using the reference as the 'module'."""
    roi = roi_presets(phantom)["breasts"]
    vol = to_slicer_order(phantom.saved_uint16)
    sl = ref.roi_slices(roi)
    box = vol[(slice(None),) + sl]
    label = np.ones(box.shape[1:], dtype=np.uint8)
    times = np.arange(8) * 64576.512
    res = ref.process(box, {"IJKmin": (0, 0, 0), "IJKmax": box.shape[:0:-1]}, segment_mask=label,
                      pre=0, early=1, late=7, spacing_mm=(0.68, 0.68, 1.5), timepoints_ms=times)
    (tmp_path / "tables").mkdir()
    np.save(tmp_path / "input_box.npy", box)
    np.save(tmp_path / "label_box.npy", label)
    ser_map = res.ser_map.copy()
    if corrupt:
        ser_map[0, 0, 0] += 1
    np.savez_compressed(tmp_path / "maps_box.npz", PE=res.pe, SER=ser_map, Washout=res.washout_map)
    with open(tmp_path / "tables" / "ser_table.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Ranges", "Volume (cm3)", "Distribution (%)"])
        w.writerows([[n, repr(float(v)), repr(float(d))] for n, v, d in res.bin_table])
    with open(tmp_path / "tables" / "summary_table.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Parameter", "Value", "Units"])
        w.writerows([[n, repr(float(v)), u] for n, v, u in res.summary_table])
    with open(tmp_path / "tables" / "tic_table.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Timepoint [min]", "PE (%)", "Linear Fit"])
        w.writerows([[repr(float(x)) for x in r] for r in res.tic_table])
    params = dict(pre=0, early=1, late=7, pe_threshold=70.0, background_threshold=60.0, display_ser_range=True,
                  ser_threshold=1.4, use_ser_map=True, ser_bins=[0.0, 0.9, 1.3, 3.0],
                  roi_ijk_min=list(roi["IJKmin"]), roi_ijk_max=list(roi["IJKmax"]), spacing_mm=[0.68, 0.68, 1.5],
                  label_shape_matches_volume=True, timepoints_ms=list(times), injection_time_ms=0.0)
    (tmp_path / "parameters.json").write_text(json.dumps(params))
    return tmp_path


def test_compare_script_passes_on_identical_output(tmp_path, phantom):
    import compare_with_plugin
    assert compare_with_plugin.run(fake_export(tmp_path, phantom))


def test_compare_script_detects_a_single_voxel_difference(tmp_path, phantom):
    import compare_with_plugin
    assert not compare_with_plugin.run(fake_export(tmp_path, phantom, corrupt=True))
