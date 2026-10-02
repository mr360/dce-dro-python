"""
Compare the output of the SeQ-DCEMRI 3D Slicer module (exported with export_plugin_results.py) with the
Slicer-free reference implementation (validation/reference_ftv.py) run on exactly the same input.

    python compare_with_plugin.py EXPORT_DIR [--phantom]

Checks, voxel by voxel inside the RefBox: PE map, SER label map, Washout label map; and every value of the
SER/FTV table, the Summary table and the TIC table. Exit status 0 means everything matched.

--phantom additionally checks that the exported input equals the regenerated FTV phantom inside the RefBox.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import reference_ftv as ref  # noqa: E402

# Tables are rounded by the module (np.round(x, 3) / np.round(x, 2)), hence the absolute tolerances.
TOL_MAP = 1e-9
TOL_TABLE = 1e-6


def read_table(path: Path) -> list[list[str]]:
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.reader(f))


def as_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


class Checker:
    def __init__(self):
        self.rows = []

    def check(self, what, ok, detail=""):
        self.rows.append((what, bool(ok), detail))

    def close(self, what, a, b, tol):
        a, b = as_float(a), as_float(b)
        ok = (np.isnan(a) and np.isnan(b)) or abs(a - b) <= tol * max(1.0, abs(b))
        self.check(what, ok, f"module={a!r} reference={b!r}")

    def report(self) -> bool:
        width = max(len(r[0]) for r in self.rows)
        for what, ok, detail in self.rows:
            print(f"{'PASS' if ok else 'FAIL'}  {what:<{width}}  {detail}")
        n_fail = sum(not r[1] for r in self.rows)
        print(f"\n{len(self.rows) - n_fail} passed, {n_fail} failed")
        return n_fail == 0


def run(export_dir: Path, check_phantom: bool = False) -> bool:
    params = json.loads((export_dir / "parameters.json").read_text())
    vol = np.load(export_dir / "input_box.npy")
    label = np.load(export_dir / "label_box.npy")
    maps = dict(np.load(export_dir / "maps_box.npz"))
    c = Checker()

    c.check("segment labelmap has the volume geometry", params.get("label_shape_matches_volume", True))
    nt, nk, nj, ni = vol.shape
    roi = {"IJKmin": (0, 0, 0), "IJKmax": (ni, nj, nk)}  # the export is already cropped to the RefBox
    res = ref.process(vol, roi, segment_mask=label, pre=params["pre"], early=params["early"],
                      late=params["late"], pe_threshold=params["pe_threshold"],
                      background_threshold=params["background_threshold"],
                      display_ser_range=params["display_ser_range"], ser_threshold=params["ser_threshold"],
                      use_ser_map=params["use_ser_map"], spacing_mm=tuple(params["spacing_mm"]),
                      timepoints_ms=np.array(params["timepoints_ms"]) if params["timepoints_ms"] else None,
                      injection_time_ms=params["injection_time_ms"])

    c.check("SER intervals", np.allclose(res.bins.intervals, params["ser_bins"] + (
        [] if params["ser_bins"][-1] == ref.SER_UPPER_THRESHOLD else [ref.SER_UPPER_THRESHOLD])),
        f"module={params['ser_bins']} reference={res.bins.intervals}")

    # ---------------------------------------------------------------------------------------------- maps
    if "PE" in maps:
        d = np.abs(maps["PE"] - res.pe)
        c.check("PE map (voxel-wise)", d.max() <= TOL_MAP * max(1.0, np.abs(res.pe).max()),
                f"max |diff| = {d.max():.3g}")
    for name, mine in (("SER", res.ser_map), ("Washout", res.washout_map)):
        if name in maps:
            n_diff = int(np.count_nonzero(maps[name].astype(int) != mine.astype(int)))
            c.check(f"{name} label map (voxel-wise)", n_diff == 0, f"{n_diff} voxels differ")

    # ---------------------------------------------------------------------------------------- SER table
    ser_rows = read_table(export_dir / "tables" / "ser_table.csv")[1:]
    module_ser = {r[0]: r[1:] for r in ser_rows if r and r[0]}
    for name, vol_cm3, dist in res.bin_table:
        if name not in module_ser:
            c.check(f"SER table row '{name}' present", False, "missing in module output")
            continue
        c.close(f"SER table '{name}' volume (cm3)", module_ser[name][0], vol_cm3, TOL_TABLE)
        c.close(f"SER table '{name}' distribution (%)", module_ser[name][1], dist, TOL_TABLE)
    extra = set(module_ser) - {r[0] for r in res.bin_table}
    c.check("SER table has no unexpected rows", not extra, ", ".join(sorted(extra)))

    # ------------------------------------------------------------------------------------ Summary table
    summary_rows = read_table(export_dir / "tables" / "summary_table.csv")[1:]
    module_summary = {r[0]: r[1] for r in summary_rows if r and r[0]}
    for name, value, _units in res.summary_table:
        if name == "ROI longest axis":
            continue  # oriented bounding box computed by Slicer; reference value is approximate
        if name == "ROI Volume":
            total = params.get("segment_voxels_total")
            if not params.get("label_shape_matches_volume", True) or (
                    total is not None and total != int(np.count_nonzero(label))):
                c.check("Summary 'ROI Volume'", True, "skipped: segment extends beyond the RefBox")
                continue
        c.close(f"Summary '{name}'", module_summary.get(name), value, TOL_TABLE)

    # ---------------------------------------------------------------------------------------- TIC table
    tic = read_table(export_dir / "tables" / "tic_table.csv")[1:]
    module_tic = np.array([[as_float(x) for x in r] for r in tic if r])
    if module_tic.shape == res.tic_table.shape:
        for j, col in enumerate(("time (min)", "mean PE (%)", "linear fit")):
            a, b = module_tic[:, j], res.tic_table[:, j]
            ok = np.allclose(a, b, rtol=TOL_TABLE, atol=TOL_TABLE, equal_nan=True)
            c.check(f"TIC table column '{col}'", ok, f"max |diff| = {np.nanmax(np.abs(a - b)):.3g}")
    else:
        c.check("TIC table shape", False, f"module {module_tic.shape} vs reference {res.tic_table.shape}")

    # ------------------------------------------------------------------------------ optional: phantom
    if check_phantom:
        from phantom import build_phantom, to_slicer_order
        ph = build_phantom()
        full = to_slicer_order(ph.saved_uint16)
        lo, hi = params["roi_ijk_min"], params["roi_ijk_max"]
        expected = full[:, lo[2]:hi[2], lo[1]:hi[1], lo[0]:hi[0]]
        same = expected.shape == vol.shape and np.array_equal(expected.astype(float), vol.astype(float))
        if not same and expected.shape == vol.shape:
            # Slicer may present rows/slices in the opposite order depending on the patient orientation
            for axes in ((2,), (1,), (1, 2)):
                if np.array_equal(np.flip(expected, axes).astype(float), vol.astype(float)):
                    same = True
                    break
        c.check("exported input equals the regenerated phantom", same)

    return c.report()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("export_dir", type=Path)
    ap.add_argument("--phantom", action="store_true")
    args = ap.parse_args()
    sys.exit(0 if run(args.export_dir, args.phantom) else 1)


if __name__ == "__main__":
    main()
