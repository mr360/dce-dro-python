"""
Expected (ground-truth) FTV outcome of the phantom, derived only from the PE / SER values the notebook
encoded in each region and the module's documented decision rules. No signal intensities are used except
the noiseless pre-contrast signal S0, which the background threshold rule needs.

This is deliberately independent of reference_ftv.py: agreement between the two shows that the module's
equations recover the encoded values; any disagreement is reported region by region (validate_phantom.py).
"""
from __future__ import annotations

import numpy as np

from phantom import Phantom
from reference_ftv import CAD_LEVEL, SER_UPPER_THRESHOLD, roi_slices, ser_bins


def broadcast_kji(map_rc: np.ndarray, n_slices: int) -> np.ndarray:
    """[row, col] -> [k, j, i] (identical slices)."""
    return np.broadcast_to(map_rc, (n_slices,) + map_rc.shape)


def expected_outcome(ph: Phantom, roi_ijk: dict, pe_threshold=70.0, background_threshold=60.0,
                     display_ser_range=True, ser_threshold=1.4, use_ser_map=True) -> dict:
    nz = ph.signal.shape[2]
    sl = roi_slices(roi_ijk)
    pe = broadcast_kji(ph.target_pe, nz)[sl]
    ser = broadcast_kji(ph.target_ser, nz)[sl]
    s0 = broadcast_kji(ph.s0, nz)[sl]
    bins = ser_bins(display_ser_range, ser_threshold)

    bkg = background_threshold / 100.0 * np.percentile(s0, 95)
    keep = (s0 >= bkg) & (pe >= pe_threshold)
    valid_ser = (ser > 0) & (ser <= SER_UPPER_THRESHOLD)

    label = np.zeros(ser.shape, dtype=int)
    for idx, (lb, ub) in enumerate(zip(bins.lower, bins.upper)):
        label[keep & (ser > lb) & (ser <= ub)] = idx + 1

    with np.errstate(divide="ignore"):
        delta = np.where(ser > 0, 100.0 * (1.0 / ser - 1.0), np.nan)  # (Sn-S1)/(S1-S0) = 1/SER - 1
    cad = np.zeros(ser.shape, dtype=int)
    cad[keep & (delta > CAD_LEVEL)] = 1                     # Persistent
    cad[keep & (np.abs(delta) <= CAD_LEVEL)] = 2            # Plateau
    cad[keep & (delta < -CAD_LEVEL)] = 3                    # Washout

    if use_ser_map:
        ftv0 = keep & valid_ser
        ftvth = keep & valid_ser & (ser > bins.ftv_threshold)
        label_map, legend = label, bins.legend
    else:
        ftv0 = cad > 0
        ftvth = (cad == 2) | (cad == 3)
        label_map, legend = cad, ["non CAD", "Persistent", "Plateau", "Washout"]
    etv = pe > 0
    return {
        "background_threshold_value": float(bkg),
        "counts": {"FTV0": int(ftv0.sum()), "FTVthresh": int(ftvth.sum()), "ETV": int(etv.sum())},
        "bin_counts": {name: int((label_map == v).sum()) for v, name in enumerate(legend) if v > 0},
        "label_map": label_map,
        "masks": {"FTV0": ftv0, "FTVthresh": ftvth, "ETV": etv},
    }


def on_threshold_mask(ph, roi_ijk: dict, pe_threshold=70.0, display_ser_range=True, ser_threshold=1.4,
                      use_ser_map=True, **_) -> np.ndarray:
    """Voxels whose encoded PE or SER lies exactly on one of the module's decision thresholds.

    For these voxels the outcome is decided by numerical round-off (the notebook's +1e-6 in s2_ser, the
    module's EPSILON, and uint16 truncation in the DICOM export), so ground truth is undefined.
    """
    nz = ph.signal.shape[2]
    sl = roi_slices(roi_ijk)
    pe = broadcast_kji(ph.target_pe, nz)[sl]
    ser = broadcast_kji(ph.target_ser, nz)[sl]
    bins = ser_bins(display_ser_range, ser_threshold)
    if use_ser_map:
        edges = set(bins.intervals) | {bins.ftv_threshold}
    else:
        edges = {1.0 / (1.0 + CAD_LEVEL / 100.0), 1.0 / (1.0 - CAD_LEVEL / 100.0)}
    edges = {e for e in edges if e > 0}
    on_ser = np.isin(ser, sorted(edges)) & (pe >= pe_threshold)
    on_pe = (pe == pe_threshold) & (pe > 0)
    return on_ser | on_pe
