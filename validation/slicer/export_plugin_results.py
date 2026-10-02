"""
Export what the SeQ-DCEMRI module computed, so it can be checked against the reference implementation.

HOW TO USE (inside 3D Slicer, after pressing "Click to Process" in the Semi-Quantitative DCE-MRI module):
    1. Open the Python console (View > Python Console).
    2. Set the output folder and run this file:
           EXPORT_DIR = r"/path/to/export_folder"
           exec(open(r"/path/to/dce-dro-python/validation/slicer/export_plugin_results.py").read())
    3. On any machine with numpy/scipy, run:
           python validation/slicer/compare_with_plugin.py /path/to/export_folder

Everything is exported inside the RefBox only (the module crops to it before computing), so files stay small:
    input_box.npy       [t, k, j, i] input sequence inside the RefBox (values as held by Slicer)
    label_box.npy       [k, j, i] selected segment inside the RefBox (as used by the module)
    maps_box.npz        output maps PE, SER (label), Washout (label), MIP inside the RefBox
    tables/*.csv        TIC Table, Summary Table, SER Table (last run only)
    parameters.json     every setting passed to quantificationLogic.process, RefBox IJK, spacing, timings

Run the export immediately after processing and before changing any setting.
Tip: restart Slicer (or close the scene) before each validation run; the module appends table columns on
every run and the export only keeps the last three columns.
"""
import csv
import json
import os

import numpy as np
import slicer

try:
    EXPORT_DIR
except NameError:
    EXPORT_DIR = os.path.join(os.path.expanduser("~"), "seq_dcemri_export")

os.makedirs(os.path.join(EXPORT_DIR, "tables"), exist_ok=True)

widget = slicer.modules.quantification.widgetRepresentation().self()
logic = widget.logic
pn = widget._parameterNode
seq = pn.input4DVolume
if seq is None:
    raise RuntimeError("No input sequence selected in the module")

ref_node = seq.GetNthDataNode(0)
roi = logic.getBoxROIIJKCoordinates(widget.roiNode, ref_node)
imin = [int(v) for v in roi["IJKmin"]]
imax = [int(v) for v in roi["IJKmax"]]
box = (slice(imin[2], imax[2]), slice(imin[1], imax[1]), slice(imin[0], imax[0]))

# Input sequence inside the RefBox, in the dtype Slicer holds (the module converts to float64)
frames = [slicer.util.arrayFromVolume(seq.GetNthDataNode(t))[box] for t in range(seq.GetNumberOfDataNodes())]
np.save(os.path.join(EXPORT_DIR, "input_box.npy"), np.stack(frames, axis=0))

# Selected segment, read exactly as the module does (no reference volume) [q1663]
label = slicer.util.arrayFromSegmentBinaryLabelmap(pn.inputMaskVolume, widget.segmentID)
full_shape = slicer.util.arrayFromVolume(ref_node).shape
label_shape_matches = tuple(label.shape) == tuple(full_shape)
np.save(os.path.join(EXPORT_DIR, "label_box.npy"), label[box] if label_shape_matches else label)

# Output maps
maps = {}
out_seq = pn.outputSequenceMaps
for name in ("MIP", "PE", "SER", "Washout"):
    node = out_seq.GetDataNodeAtValue(name) if out_seq else None
    if node is not None:
        maps[name] = slicer.util.arrayFromVolume(node)[box]
np.savez_compressed(os.path.join(EXPORT_DIR, "maps_box.npz"), **maps)


def export_table(table_node, filename):
    table = table_node.GetTable()
    ncols = table.GetNumberOfColumns()
    cols = [table.GetColumn(c) for c in range(max(0, ncols - 3), ncols)]  # last run = last three columns
    with open(os.path.join(EXPORT_DIR, "tables", filename), "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([col.GetName() for col in cols])
        for r in range(table.GetNumberOfRows()):
            row = []
            for col in cols:
                v = col.GetValue(r) if r < col.GetNumberOfValues() else ""
                row.append(v if isinstance(v, str) else repr(float(v)))  # full double precision
            writer.writerow(row)


export_table(widget.TICTableNode, "tic_table.csv")
export_table(widget.SummaryTableNode, "summary_table.csv")
export_table(widget.SERDistributionTableNode, "ser_table.csv")

spacing = ref_node.GetSpacing()  # (i, j, k) in mm
timings = widget.timings or {}
params = {
    "pre": int(pn.indicesDCE.preContrast),
    "early": int(pn.indicesDCE.earlyPostContrast),
    "late": int(pn.indicesDCE.latePostContrast),
    "pe_threshold": float(pn.peakEnhancementThreshold),
    "background_threshold": float(pn.backgroundThreshold),
    "display_ser_range": bool(pn.displaySERrange),
    "ser_threshold": float(pn.signalEnhancementRatioThreshold),
    "use_ser_map": bool(widget.useSERmap),
    "ser_upper_threshold": float(widget.SER_UPPER_THRESHOLD),
    "ser_delta_factor": float(widget.SER_DELTA_FACTOR),
    "ser_bins": [float(x) for x in widget.serMapInterval],
    "roi_ijk_min": imin,
    "roi_ijk_max": imax,
    "omit_regions": [{"IJKmin": [int(v) for v in o["IJKmin"]], "IJKmax": [int(v) for v in o["IJKmax"]]}
                     for o in (logic.getBoxROIIJKCoordinates(n, ref_node) for n in widget.omitRoiList)],
    "spacing_mm": list(spacing),
    "volume_shape_kji": list(full_shape),
    "label_shape_matches_volume": label_shape_matches,
    "segment_voxels_total": int(np.count_nonzero(label)),  # ROI statistics use the whole segment [q1877]
    "timepoints_ms": [float(t) for t in np.asarray(timings.get("timepoints", []))],
    "injection_time_ms": float(timings.get("injectionTime", 0.0)),
    "segment_name": pn.inputMaskVolume.GetSegmentation().GetSegment(widget.segmentID).GetName(),
    "slicer_version": slicer.app.applicationVersion,
    "sequence_name": seq.GetName(),
}
with open(os.path.join(EXPORT_DIR, "parameters.json"), "w") as f:
    json.dump(params, f, indent=2)

print(f"SeQ-DCEMRI results exported to {EXPORT_DIR}")
if not label_shape_matches:
    print("WARNING: the segment labelmap does not have the volume's shape; see parameters.json")
