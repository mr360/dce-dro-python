"""
Validate the SeQ-DCEMRI FTV computations against the FTV digital phantom.

For every scenario three results are compared:
  expected   ground truth from the encoded PE / SER values (ground_truth.py)
  float      reference implementation (reference_ftv.py) on the notebook's floating-point signal
  dicom      reference implementation on the uint16 pixel values written to DICOM (what 3D Slicer reads)

Differences are explained region by region. Run:
    python validate_phantom.py                   # full phantom (120 slices), writes results/
    python validate_phantom.py --slices 6        # quick run
    python validate_phantom.py --dicom-dir DIR   # also check that DIR holds exactly this phantom
"""
from __future__ import annotations

import argparse
import csv
import datetime
import platform
import subprocess
from pathlib import Path

import numpy as np
import scipy

import reference_ftv as ref
from ground_truth import broadcast_kji, expected_outcome, on_threshold_mask
from phantom import build_phantom, load_dicom_series, to_slicer_order

HERE = Path(__file__).resolve().parent


# ----------------------------------------------------------------------------------------------- ROI boxes
def bbox_ijk(mask_rc: np.ndarray, n_slices: int, margin: int = 0) -> dict:
    rows, cols = np.nonzero(mask_rc)
    return {"IJKmin": (int(cols.min()) - margin, int(rows.min()) - margin, 0),
            "IJKmax": (int(cols.max()) + 1 + margin, int(rows.max()) + 1 + margin, n_slices)}


def roi_presets(ph) -> dict:
    nz = ph.signal.shape[2]
    names = ph.tissue_name_map()
    breast = np.isin(names, ["breast_fat", "breast_glandular"])
    ncols = ph.t1.shape[1]
    left = breast & (np.arange(ncols)[None, :] < ncols // 2)
    right = breast & (np.arange(ncols)[None, :] >= ncols // 2)
    shape = (ph.t1.shape[1], ph.t1.shape[0], nz)  # (i, j, k) sizes
    # Default RefBox of the module: fitted to the volume, then size / 4, centred [q927-931] (approximate)
    centre = [(n - 1) / 2.0 for n in shape]
    default = {"IJKmin": tuple(int(round(c - n / 8.0)) for c, n in zip(centre, shape)),
               "IJKmax": tuple(int(round(c + n / 8.0)) + 1 for c, n in zip(centre, shape))}
    return {
        "breasts": bbox_ijk(breast, nz),
        "left_breast": bbox_ijk(left, nz),
        "right_breast": bbox_ijk(right, nz),
        "whole_phantom": bbox_ijk(ph.target_pe > 0, nz),
        "module_default_box": default,
    }


# ------------------------------------------------------------------------------------------------ helpers
def run_reference(vol, ph, roi, **kw):
    sx, sy, sz = ph.resolution_mm[1], ph.resolution_mm[0], ph.resolution_mm[2]
    nt = vol.shape[0]
    times = np.arange(nt) * ph.frame_interval_s * 1000.0  # TriggerTime written by the notebook (ms)
    return ref.process(vol, roi, pre=kw.pop("pre", ph.pre_index), early=kw.pop("early", ph.early_index),
                       late=kw.pop("late", ph.late_index), spacing_mm=(sx, sy, sz), timepoints_ms=times, **kw)


def region_breakdown(ph, roi, res_dicom, res_float, exp, edge) -> list[dict]:
    nz = ph.signal.shape[2]
    sl = ref.roi_slices(roi)
    names = broadcast_kji(ph.tissue_name_map(), nz)[sl]
    pe_t = broadcast_kji(ph.target_pe, nz)[sl]
    ser_t = broadcast_kji(ph.target_ser, nz)[sl]
    lbl_d = res_dicom.ser_map if res_dicom.params["use_ser_map"] else res_dicom.washout_map
    lbl_f = res_float.ser_map if res_float.params["use_ser_map"] else res_float.washout_map
    legend = res_dicom.bins.legend if res_dicom.params["use_ser_map"] else ref.CAD_LEGEND
    rows = []
    keys = sorted({(n, p, s) for n, p, s in zip(names.ravel(), pe_t.ravel(), ser_t.ravel())},
                  key=lambda x: (x[0], x[1], x[2]))
    for name, p, s in keys:
        m = (names == name) & (pe_t == p) & (ser_t == s)
        def dist(lbl):
            vals, cnt = np.unique(lbl[m].astype(int), return_counts=True)
            return "; ".join(f"{legend[v] if v else 'excluded'}: {c}" for v, c in zip(vals, cnt))
        rows.append({
            "tissue": name, "target_PE": p, "target_SER": s, "voxels": int(m.sum()),
            "S0_float": float(broadcast_kji(ph.s0, nz)[sl][m][0]),
            "S0_dicom": float(res_dicom.st0[m][0]),
            "PE_dicom_min": float(res_dicom.pe_raw[m].min()), "PE_dicom_max": float(res_dicom.pe_raw[m].max()),
            "SER_dicom_min": float(res_dicom.ser_raw[m].min()), "SER_dicom_max": float(res_dicom.ser_raw[m].max()),
            "SER_float": float(res_float.ser_raw[m][0]),
            "expected": dist(exp["label_map"]), "float": dist(lbl_f), "dicom": dist(lbl_d),
            "on_threshold": bool(edge[m].any()),
            "agrees": bool(np.array_equal(exp["label_map"][m], lbl_d[m].astype(int))),
        })
    return rows


def reference_masks(res) -> dict:
    """Voxel sets behind FTV0 / FTVthresh / ETV in the reference result (same rules as reference_ftv.process)."""
    if res.params["use_ser_map"]:
        return {"FTV0": res.ser > 0.0, "FTVthresh": res.ser > res.bins.ftv_threshold, "ETV": res.pe_raw > 0.0}
    w = res.washout_map
    return {"FTV0": w > 0.0, "FTVthresh": (w == 2) | (w == 3), "ETV": res.pe_raw > 0.0}


def git_commit(path: Path) -> str:
    try:
        return subprocess.check_output(["git", "-C", str(path), "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def fmt(x, nd=3):
    if isinstance(x, (float, np.floating)):
        return f"{x:.{nd}f}"
    return str(x)


def md_table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(fmt(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def write_csv(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


# -------------------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--slices", type=int, default=None, help="number of slices (default: notebook value, 120)")
    ap.add_argument("--out", type=Path, default=HERE / "results")
    ap.add_argument("--dicom-dir", type=Path, default=None, help="folder with the phantom DICOM files")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    ph = build_phantom(n_slices=args.slices)
    nz = ph.signal.shape[2]
    vol_dicom = to_slicer_order(ph.saved_uint16)
    vol_float = to_slicer_order(ph.signal)
    voxel_cm3 = ph.voxel_volume_mm3 / 1000.0
    rois = roi_presets(ph)
    report = []
    R = report.append

    R("# Validation of SeQ-DCEMRI FTV computations with the FTV digital phantom\n")
    R(f"Generated {datetime.datetime.now().isoformat(timespec='seconds')} by `validation/validate_phantom.py`.\n")
    R("## Provenance\n")
    R(md_table(["Item", "Value"], [
        ("Phantom notebook SHA-256", ph.notebook_sha256),
        ("dce-dro-python commit", git_commit(HERE.parent)),
        ("parametricDCEMRI commit audited", ref.AUDITED_COMMIT),
        ("Python / numpy / scipy", f"{platform.python_version()} / {np.__version__} / {scipy.__version__}"),
        ("Phantom size [row, col, slice, time]", "×".join(map(str, ph.signal.shape))),
        ("Voxel size (mm)", " × ".join(map(str, ph.resolution_mm)) + f" = {ph.voxel_volume_mm3:.4f} mm³"),
        ("Pre / early / late frame", f"{ph.pre_index} / {ph.early_index} / {ph.late_index}"),
        ("Frame interval written to DICOM", f"{ph.frame_interval_s:.3f} s"),
    ]))

    if args.dicom_dir:
        dcm = load_dicom_series(args.dicom_dir)
        same = dcm.shape == vol_dicom.shape and np.array_equal(dcm, vol_dicom)
        R(f"\n**DICOM check** ({args.dicom_dir}): shape {dcm.shape}, "
          f"{'identical to the regenerated phantom' if same else 'DIFFERENT from the regenerated phantom'}.\n")

    R("\n## ROI boxes used (IJK, exclusive upper bound)\n")
    R(md_table(["ROI", "IJK min", "IJK max"], [(k, v["IJKmin"], v["IJKmax"]) for k, v in rois.items()]))

    # ---------------------------------------------------------------------------- scenario comparison
    scenarios = [
        ("A  defaults, preset SER range", "breasts", {}),
        ("B  PE threshold 0", "breasts", {"pe_threshold": 0.0}),
        ("C  PE threshold 50", "whole_phantom", {"pe_threshold": 50.0}),
        ("D  background 0, PE 0", "whole_phantom", {"pe_threshold": 0.0, "background_threshold": 0.0}),
        ("E  single SER threshold 1.0", "breasts", {"display_ser_range": False, "ser_threshold": 1.0}),
        ("F  single SER threshold 1.25", "breasts", {"pe_threshold": 50.0, "display_ser_range": False,
                                                     "ser_threshold": 1.25}),
        ("G  Washout (CAD) map", "breasts", {"use_ser_map": False, "pe_threshold": 50.0}),
        ("H  left breast only", "left_breast", {"pe_threshold": 50.0}),
        ("I  right breast only", "right_breast", {"pe_threshold": 50.0}),
    ]
    if nz >= 24:  # the default box is a quarter of the volume and must be >= 3 slices thick
        scenarios.append(("J  module default RefBox", "module_default_box", {}))
    rows_csv, rows_md, breakdowns = [], [], {}
    for name, roi_name, kw in scenarios:
        roi = rois[roi_name]
        exp = expected_outcome(ph, roi, **kw)
        rd = run_reference(vol_dicom, ph, roi, **kw)
        rf = run_reference(vol_float, ph, roi, **kw)
        for key in ("FTV0", "FTVthresh", "ETV"):
            rows_csv.append({"scenario": name, "roi": roi_name, "settings": repr(kw), "quantity": key,
                             "expected_voxels": exp["counts"][key], "float_voxels": rf.counts[key],
                             "dicom_voxels": rd.counts[key], "expected_cm3": exp["counts"][key] * voxel_cm3,
                             "dicom_cm3": rd.volumes_cm3[key]})
        ok = all(exp["counts"][k] == rd.counts[k] for k in exp["counts"])
        edge = on_threshold_mask(ph, roi, **kw)
        lbl = rd.ser_map if rd.params["use_ser_map"] else rd.washout_map
        off = ~edge
        ok_off = (np.array_equal(exp["label_map"][off], lbl[off].astype(int)) and
                  all(np.array_equal(exp["masks"][k][off], m[off]) for k, m in reference_masks(rd).items()))
        rows_md.append((name, roi_name,
                        f"{exp['counts']['FTV0']} / {rf.counts['FTV0']} / {rd.counts['FTV0']}",
                        f"{exp['counts']['FTVthresh']} / {rf.counts['FTVthresh']} / {rd.counts['FTVthresh']}",
                        f"{exp['counts']['ETV']} / {rd.counts['ETV']}",
                        f"{rd.volumes_cm3['FTV0']:.3f}", "yes" if ok else "no",
                        int(edge.sum()), "yes" if ok_off else "**NO**"))
        breakdowns[name] = region_breakdown(ph, roi, rd, rf, exp, edge)
    write_csv(args.out / "scenarios.csv", rows_csv)

    R("\n## Scenario results (voxel counts)\n")
    R("Columns give *expected / float / DICOM* voxel counts. *Expected* comes from the encoded PE/SER only; "
      "*float* and *DICOM* come from the reference implementation of the module applied to the noiseless "
      "floating-point signal and to the uint16 pixel values stored in DICOM. Unless stated otherwise: "
      "frames 0/1/7, background 60 %, PE ≥ 70 %, pre-defined SER range.\n")
    R("*On-threshold voxels* belong to regions whose encoded PE or SER equals a decision threshold of that "
      "scenario (e.g. SER = 0.9 with the 0.9 cut-off); their classification is decided by round-off and has no "
      "defined ground truth. The last column compares every other voxel individually (label map, FTV, FTVth "
      "and ETV membership).\n")
    R(md_table(["Scenario", "ROI", "FTV (SER>0)", "FTVth", "ETV exp / DICOM", "FTV cm³ (DICOM)",
                "All counts equal", "On-threshold voxels", "Voxel-wise agreement off-threshold"], rows_md))

    # ------------------------------------------------------------------------- region-level breakdown
    for name in (scenarios[0][0], scenarios[1][0], scenarios[6][0]):
        bd = breakdowns[name]
        write_csv(args.out / f"regions_{name.split()[0]}.csv", bd)
        R(f"\n## Region breakdown – scenario {name}\n")
        R("Each region is a (tissue, target PE, target SER) combination inside the ROI. Measured ranges are "
          "from the DICOM (uint16) data.\n")
        R(md_table(["Tissue", "PE", "SER", "Voxels", "S0 float→DICOM", "PE measured", "SER measured",
                    "SER float", "Expected", "DICOM", "On threshold", "Agree"],
                   [(r["tissue"], r["target_PE"], r["target_SER"], r["voxels"],
                     f"{r['S0_float']:.2f}→{r['S0_dicom']:.0f}",
                     f"{r['PE_dicom_min']:.2f}–{r['PE_dicom_max']:.2f}",
                     f"{r['SER_dicom_min']:.4f}–{r['SER_dicom_max']:.4f}", f"{r['SER_float']:.7f}",
                     r["expected"], r["dicom"], "yes" if r["on_threshold"] else "",
                     "yes" if r["agrees"] else "**no**") for r in bd]))

    # ------------------------------------------------------------------------------------- sweeps
    sweeps = [
        ("PE threshold (%)", "pe_threshold", list(np.arange(0, 141, 10.0)), "whole_phantom", {}),
        ("Background threshold (%)", "background_threshold", list(np.arange(0, 101, 10.0)), "whole_phantom",
         {"pe_threshold": 0.0}),
        ("Single SER threshold", "ser_threshold", [round(x, 2) for x in np.arange(0.0, 2.01, 0.1)], "breasts",
         {"display_ser_range": False, "pe_threshold": 50.0}),
    ]
    for title, key, values, roi_name, base in sweeps:
        roi = rois[roi_name]
        rows = []
        for v in values:
            kw = dict(base, **{key: float(v)})
            exp = expected_outcome(ph, roi, **kw)
            rd = run_reference(vol_dicom, ph, roi, **kw)
            rows.append({key: v, "expected_FTV0": exp["counts"]["FTV0"], "dicom_FTV0": rd.counts["FTV0"],
                         "expected_FTVth": exp["counts"]["FTVthresh"], "dicom_FTVth": rd.counts["FTVthresh"],
                         "dicom_FTV0_cm3": rd.volumes_cm3["FTV0"], "dicom_FTVth_cm3": rd.volumes_cm3["FTVthresh"],
                         "background_value_dicom": rd.background_threshold_value})
        write_csv(args.out / f"sweep_{key}.csv", rows)
        R(f"\n## Sweep: {title} (ROI {roi_name})\n")
        R(md_table(list(rows[0].keys()), [tuple(r.values()) for r in rows]))

    # ----------------------------------------------------------------- time-point selection sensitivity
    roi = rois["breasts"]
    rows = []
    for early in range(1, ph.n_frames - 1):
        rd = run_reference(vol_dicom, ph, roi, early=early, pe_threshold=50.0)
        m = rd.base_mask
        rows.append((early, rd.counts["FTV0"], rd.counts["FTVthresh"],
                     float(rd.pe_raw[m].mean()) if m.any() else float("nan"),
                     float(rd.ser_raw[m].mean()) if m.any() else float("nan")))
    R("\n## Sensitivity to the early post-contrast frame (ROI breasts, PE ≥ 50 %)\n")
    R("The phantom encodes PE/SER between frames 0, 1 and 7. The module's default early frame is 3 "
      "([q68]); any other choice measures different PE and SER.\n")
    R(md_table(["Early frame", "FTV voxels", "FTVth voxels", "mean PE in FTV", "mean SER in FTV"], rows))

    # ------------------------------------------------------------------------ full module tables for A
    rd = run_reference(vol_dicom, ph, rois["breasts"])
    R("\n## Module tables reproduced for scenario A\n")
    R("### Summary Table\n")
    R(md_table(["Parameter", "Value", "Units"], rd.summary_table))
    R("\n### SER Table\n")
    R(md_table(["Ranges", "Volume (cm3)", "Distribution (%)"], rd.bin_table))
    R(f"\nROI statistics: segment volume {rd.roi_stats['segment_volume_cm3']:.3f} cm³; "
      f"value reported as 'ROI Volume' (× π/6) {rd.roi_stats['reported_roi_volume_cm3']:.3f} cm³.\n")

    (args.out / "report.md").write_text("\n".join(report) + "\n")
    print(f"Report written to {args.out / 'report.md'}")
    for r in rows_md:
        print(" | ".join(map(str, r)))


if __name__ == "__main__":
    main()
