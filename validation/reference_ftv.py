"""
Slicer-free reference implementation of the computations performed by the SeQ-DCEMRI 3D Slicer module
(parametricDCEMRI, quantification/quantification.py, quantificationLogic.process).

Purpose: make every number reported by the module traceable. Each step below is a transcription of the
module's numpy code; the comment ``[qN]`` refers to line N of quantification.py at the audited commit
(AUDITED_COMMIT). Behaviour is reproduced as-is, including quirks, so that the output can be compared
voxel-by-voxel with the module (see slicer/compare_with_plugin.py). Nothing here "fixes" the module.

Array convention (same as the module): 4D input is [time, k(slice), j(row), i(col)], i.e. the stack of
slicer.util.arrayFromVolume() for each frame [q1424-1440]. ROI boxes are given in IJK voxel indices with
an exclusive upper bound, as returned by quantificationLogic.getBoxROIIJKCoordinates [q1517-1540].
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import signal

AUDITED_COMMIT = "7490646dfb019dbdfb9d8789054e5333d17837f7"

EPSILON = 1.0e-6               # [q1234] quantificationLogic.EPSILON
SER_UPPER_THRESHOLD = 3.0      # [q267]  SER above this is "non-SER"
SER_DELTA_FACTOR = 4.0 / 9.0   # [q272]  width of the 2nd interval in single-threshold mode
CAD_LEVEL = 10.0               # [q1146-1149] Persistent / Plateau / Washout cut-off (%)
CAD_LEGEND = ["non CAD", "Persistent", "Plateau", "Washout"]  # [q1138-1144]
ROI_ELLIPSOID_SCALE = np.pi / 6.0  # [q1887]


@dataclass
class SerBins:
    intervals: list           # interval edges, bins are lb < SER <= ub
    ftv_threshold: float      # SERthreshold used for "FTVth (SER>SER_th)"
    legend: list              # names of the label values 0..n (0 = 'non SER')

    @property
    def lower(self):
        return self.intervals[:-1]

    @property
    def upper(self):
        return self.intervals[1:]


def ser_bins(display_ser_range: bool = True, ser_threshold: float = 1.4,
             ser_upper: float | None = SER_UPPER_THRESHOLD) -> SerBins:
    """Mirror of quantificationWidget.setSERColourMapDict(update=True, serUpperThreshold=3.0) [q1152-1212],
    as called from onApplyButton [q656]."""
    if display_ser_range:                                   # [q1166-1173] pre-defined range
        thr = 0.9
        intervals = [0.00, 0.90, 1.30]
    elif ser_threshold == 0.0:                              # [q1174-1179]
        thr = 0
        intervals = [0.0]
    else:                                                   # [q1180-1189] single threshold
        upper = (1.0 + SER_DELTA_FACTOR) * ser_threshold
        thr = ser_threshold
        intervals = [0.0, np.round(ser_threshold, 2), np.round(upper, 2)]
    if ser_upper is not None and ser_upper > max(intervals):  # [q1191-1193]
        intervals.append(ser_upper)
    legend = ["non SER"] + [f"{lb:.2f} < SER ≤ {ub:.2f}" for lb, ub in zip(intervals[:-1], intervals[1:])]
    return SerBins(intervals=intervals, ftv_threshold=thr, legend=legend)


def roi_slices(roi_ijk: dict) -> tuple:
    """IJK box -> numpy slices in [k, j, i] order, as used throughout process() (e.g. [q1687])."""
    lo, hi = roi_ijk["IJKmin"], roi_ijk["IJKmax"]
    return (slice(lo[2], hi[2]), slice(lo[1], hi[1]), slice(lo[0], hi[0]))


def default_timepoints_ms(nt: int) -> np.ndarray:
    """Fallback time axis when the sequence has no MultiVolume.FrameLabels attribute [q842-847].
    Note: nt points from 0 to nt minutes inclusive, i.e. spacing nt/(nt-1) min, not exactly 1 min."""
    return np.linspace(0, nt * 60.0 * 1.0e3, num=nt, endpoint=True)


@dataclass
class Result:
    params: dict
    bins: SerBins
    # Maps on the cropped grid [k, j, i]
    st0: np.ndarray
    pe_raw: np.ndarray            # 100*(S1-S0)/S0 before any masking (also used for ETV)
    pe: np.ndarray                # PE after masking (0 outside base_mask) -> "PE" output map
    ser_raw: np.ndarray           # (S1-S0)/(Sn-S0) before clipping/masking
    ser: np.ndarray               # SER after clipping to [0, 3] and masking -> used for FTV and peak SER
    base_mask: np.ndarray
    ser_map: np.ndarray           # label 0..n -> "SER" output map
    wo_delta: np.ndarray          # 100*(Sn-S1)/(S1-S0)
    washout_map: np.ndarray       # label 0..3 -> "Washout" output map
    background_threshold_value: float
    counts: dict                  # voxel counts of FTV0, FTVthresh, ETV
    volumes_cm3: dict
    bin_table: list               # rows of the "SER Table" (or CAD labels) as produced by the module
    summary_table: list           # rows of the "Summary Table"
    tic_table: np.ndarray         # "TIC Table": time [min], mean PE(t) [%], linear fit
    roi_stats: dict = field(default_factory=dict)


def process(volume4d: np.ndarray,
            roi_ijk: dict,
            segment_mask: np.ndarray | None = None,
            omit_boxes_ijk: tuple = (),
            pre: int = 0,
            early: int = 1,
            late: int = -1,
            pe_threshold: float = 70.0,
            background_threshold: float = 60.0,
            display_ser_range: bool = True,
            ser_threshold: float = 1.4,
            use_ser_map: bool = True,
            spacing_mm: tuple = (1.0, 1.0, 1.0),
            timepoints_ms: np.ndarray | None = None,
            injection_time_ms: float = 0.0) -> Result:
    """Reference transcription of quantificationLogic.process [q1573-2033].

    volume4d      [t, k, j, i] signal intensities (any numeric dtype; converted to float64 as in [q1433])
    roi_ijk       {'IJKmin': (i, j, k), 'IJKmax': (i, j, k)} exclusive upper bound (RefBox markup)
    segment_mask  [k, j, i] binary mask of the selected segment, full volume; None/empty -> RefBox used
    spacing_mm    voxel size (i, j, k) in mm; only used for volumes in cm3
    """
    if pre == early:                                                    # [q1608-1609]
        raise ValueError("Pre Contrast Index Cannot be the same as the Early Post Contrast")

    vol = np.asarray(volume4d, dtype=float)                             # [q1433-1438] float64 copy
    nt, nz, ny, nx = vol.shape
    bins = ser_bins(display_ser_range, ser_threshold)
    sl = roi_slices(roi_ijk)
    voxel_cm3 = float(np.prod(spacing_mm)) / 1000.0

    # --- Analysis mask: selected segment, or the RefBox minus omit regions if the segment is empty [q1663-1677]
    label = None if segment_mask is None else np.asarray(segment_mask).astype(np.uint8)
    if label is None or not label.any():
        label = np.zeros((nz, ny, nx), dtype=np.uint8)
        label[sl] = 1
        for omit in omit_boxes_ijk:
            label[roi_slices(omit)] = 0
    full_label = label

    # --- Crop everything to the RefBox [q1682-1688]
    vol = vol[(slice(None),) + sl]
    label = label[sl]

    # --- Pre-contrast image and background threshold [q1691-1698]
    st0 = vol[pre]
    bkg_value = (background_threshold / 100.0) * np.percentile(st0, 95)   # 95th pct over the whole box
    base_mask = (st0 >= bkg_value) & label                              # bool & uint8 -> uint8 (as in module)

    # --- Percentage enhancement [q1700-1709]
    st_minus_st0 = vol - st0
    st1_minus_st0 = st_minus_st0[early]
    stn_minus_st0 = st_minus_st0[late]
    pe_raw = 100 * st1_minus_st0 / (st0 + EPSILON)
    etv_pe = 100 * st1_minus_st0 / (st0 + EPSILON)                      # [q1706] ETVthreshold, never masked
    base_mask &= (pe_raw >= pe_threshold)
    pe = np.where(base_mask, pe_raw, 0)

    # --- Signal enhancement ratio [q1718-1723]
    ser_raw = st1_minus_st0 / (stn_minus_st0 + EPSILON)
    ser = ser_raw.copy()
    ser[ser < 0.0] = 0.0
    ser[ser > SER_UPPER_THRESHOLD] = 0.0
    base_mask &= (ser >= 0.0)                                           # always true after clipping: no-op
    ser = np.where(base_mask, ser, 0)
    seg_points = np.where(base_mask)                                    # [q1737]

    # --- SER label map, lb < SER <= ub [q1739-1751]
    ser_map = np.zeros_like(ser)
    for idx, (lb, ub) in enumerate(zip(bins.lower, bins.upper)):
        ser_map[(ser > lb) & (ser <= ub)] = idx + 1
    if len(bins.upper) == 0:
        ser_map[ser > 0.0] = 1
    ser_map *= base_mask

    # --- Washout (CAD-like) map [q1769-1784]. 'non CAD' leaves zeros (the module re-assigns 0 to the voxels of
    #     the last SER interval through a reused variable, which has no effect on a zero-initialised map).
    wo_delta = 100.0 * (stn_minus_st0 - st1_minus_st0) / (st1_minus_st0 + EPSILON)
    washout_map = np.zeros_like(wo_delta)
    for idx, name in enumerate(CAD_LEGEND):
        if name == "Plateau":
            washout_map[np.abs(wo_delta) <= CAD_LEVEL] = idx
        elif name == "Persistent":
            washout_map[wo_delta > CAD_LEVEL] = idx
        elif name == "Washout":
            washout_map[wo_delta < -CAD_LEVEL] = idx
    washout_map *= base_mask

    # --- Peak SER / PE: mean over 3x3x3 neighbourhoods on a non-overlapping grid, then max [q1812-1824]
    nrows, ncols, ndepth = st0.shape
    kernel = np.ones((3, 3, 3))
    kernel /= kernel.sum()
    mean_ser = signal.convolve(ser, kernel, mode="same")
    mean_pe = signal.convolve(pe, kernel, mode="same")
    grid = (slice(1, nrows - np.mod(nrows, 3), 3), slice(1, ncols - np.mod(ncols, 3), 3),
            slice(1, ndepth - np.mod(ndepth, 3), 3))
    if mean_ser[grid].size == 0:
        # The module raises the same numpy ValueError here: the RefBox must be >= 3 voxels along every axis
        raise ValueError("RefBox is thinner than 3 voxels along at least one axis; the module fails at [q1823]")
    peak_ser = mean_ser[grid].max()
    peak_pe = mean_pe[grid].max()

    # --- FTV / ETV voxel sets [q1826-1846]
    if use_ser_map:
        map_volumes = {"FTV0": ser > 0.0,
                       "FTVthresh": ser > bins.ftv_threshold,
                       "ETV": etv_pe > 0.0}
    else:
        map_volumes = {"FTV0": washout_map > 0.0,
                       "FTVthresh": (washout_map == CAD_LEGEND.index("Plateau")) |
                                    (washout_map == CAD_LEGEND.index("Washout")),
                       "ETV": etv_pe > 0}
    counts = {k: int(np.count_nonzero(v)) for k, v in map_volumes.items()}
    # Segment Statistics 'volume_cm3' = voxel count x voxel volume [q1849-1861, q1959-1961]
    volumes = {k: c * voxel_cm3 for k, c in counts.items()}

    # --- Time-intensity curve and enhancement summaries over the voxels in base_mask [q1622-1625, q1866-1874]
    if timepoints_ms is None:
        timepoints_ms = default_timepoints_ms(nt)
    tic = np.full((nt, 3), np.nan)
    tic[:, 0] = np.asarray(timepoints_ms, dtype=float) / (1000 * 60)
    uptake = 100 * st_minus_st0 / (st0 + EPSILON)
    for t in range(nt):
        tic[t, 1] = uptake[t][seg_points].mean() if seg_points[0].size else np.nan
    max_enh = np.max(uptake, axis=0)
    delta_enh = (uptake[late] - uptake[early])[seg_points]
    first_pass = uptake[early][seg_points]
    with np.errstate(all="ignore"):
        lin = np.polyfit(tic[1:, 0], tic[1:, 1], 1)
    tic[1:, 2] = np.polyval(lin, tic[1:, 0])
    slope = lin[0]

    # --- ROI statistics of the selected segment (whole segment, not cropped) [q1877-1890]
    roi_vox = int(full_label.sum())
    roi_volume_cm3 = roi_vox * voxel_cm3
    kk, jj, ii = np.nonzero(full_label)
    extent_mm = ([(ii.max() - ii.min() + 1) * spacing_mm[0], (jj.max() - jj.min() + 1) * spacing_mm[1],
                  (kk.max() - kk.min() + 1) * spacing_mm[2]] if roi_vox else [0.0, 0.0, 0.0])
    roi_stats = {"voxel_count": roi_vox,
                 "segment_volume_cm3": roi_volume_cm3,
                 "reported_roi_volume_cm3": roi_volume_cm3 * ROI_ELLIPSOID_SCALE,
                 # Oriented bounding box diameter; exact only for box-shaped, axis-aligned segments
                 "longest_axis_mm_approx": float(max(extent_mm))}

    summary = [  # [q1901-1941]
        ("Peak SER", np.round(peak_ser, 3), "[]"),
        ("Peak PE", np.round(peak_pe, 3), "%"),
        ("PE Threshold", pe_threshold, "%"),
        ("SER Upper Threshold", bins.ftv_threshold, "[]"),   # label in module; value is the FTV SER threshold
        ("ROI longest axis", roi_stats["longest_axis_mm_approx"], "mm"),
        ("ROI Volume", roi_stats["reported_roi_volume_cm3"], "cm3"),
        ("Bolus injection time", injection_time_ms / (1000 * 60), "min"),
        ("Early Phase Time", tic[early, 0], "min"),
        ("Late Phase Time", tic[late, 0], "min"),
        ("Maximum Enhancement", max_enh[seg_points].mean() if seg_points[0].size else np.nan, "%"),
        ("Delta Enhancement", delta_enh.mean() if delta_enh.size else np.nan, "%"),
        ("First Pass Enhancement", first_pass.mean() if first_pass.size else np.nan, "%"),
        ("Enhancement Slope", slope, "[]"),
    ]

    # --- Distribution table [q1957-2006]. Bin volumes come from the label map imported as segments.
    label_map, legend = (ser_map, bins.legend) if use_ser_map else (washout_map, CAD_LEGEND)
    ftv0 = counts["FTV0"]
    bin_table = []
    for value, name in enumerate(legend):
        if value == 0:
            continue
        n = int(np.count_nonzero(label_map == value))
        if n == 0:
            continue  # empty labels produce no segment, hence no row
        bin_table.append((name, np.round(n * voxel_cm3, 3), np.round(100 * n / ftv0, 2) if ftv0 else np.nan))
    bin_table.append(("FTVth (SER>SER_th)", np.round(volumes["FTVthresh"], 3),
                      np.round(100 * counts["FTVthresh"] / ftv0, 2) if ftv0 else np.nan))
    bin_table.append(("FTV (SER>0)", np.round(volumes["FTV0"], 3), 100.0 if ftv0 else np.nan))
    bin_table.append(("ETV (Enhanced Tumour Volume)", np.round(volumes["ETV"], 3), 100.0))

    params = dict(pre=pre, early=early, late=late, pe_threshold=pe_threshold,
                  background_threshold=background_threshold, display_ser_range=display_ser_range,
                  ser_threshold=ser_threshold, use_ser_map=use_ser_map, roi_ijk=roi_ijk,
                  spacing_mm=tuple(spacing_mm))
    return Result(params=params, bins=bins, st0=st0, pe_raw=pe_raw, pe=pe, ser_raw=ser_raw, ser=ser,
                  base_mask=base_mask.astype(bool), ser_map=ser_map, wo_delta=wo_delta, washout_map=washout_map,
                  background_threshold_value=float(bkg_value), counts=counts, volumes_cm3=volumes,
                  bin_table=bin_table, summary_table=summary, tic_table=tic, roi_stats=roi_stats)
