"""
Rebuild the FTV digital phantom (FTV_digitalPhantom.ipynb) without Jupyter, DICOM or plotting.

The phantom is NOT re-implemented here: the relevant code cells of the notebook are located by
a marker string and executed verbatim, so the data analysed by the validation is, by construction,
the data produced by the notebook. The only change is the number of slices, which can be reduced
for quick tests (all slices of the phantom are identical, so voxel counts scale linearly with it).

The notebook writes the 4D signal to DICOM as ``np.uint16(signal)`` (truncation towards zero,
RescaleSlope = 1, RescaleIntercept = 0; see the DICOM cell of the notebook). ``Phantom.saved_uint16``
reproduces exactly that step, so it holds the pixel values that 3D Slicer reads from the DICOM files.
"""
from __future__ import annotations

import hashlib
import io
import json
from contextlib import redirect_stdout
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

NOTEBOOK_PATH = Path(__file__).resolve().parents[1] / "FTV_digitalPhantom.ipynb"

# Each marker must identify exactly one code cell of the notebook. Cells are executed in notebook order.
CELL_MARKERS = (
    "def relaxationTerm(",            # E1 / E2 relaxation terms
    "def spgr_signal(",               # SPGR signal equation
    "def s2_ser(",                    # S(early) from PE, S(late) from SER
    "def get_coordinates_from_size(", # rectangle helper
    "tissue_pars = {",                # T1 / T2* / M0 per tissue
    "image_pars = {",                 # matrix size, resolution, TR / TE / FA, Tacq
    "baseline_images = {",            # drawing of the regions (T1, T2*, M0, S0, PE, SER maps)
    "StFull = np.zeros(",             # time axis: S0, S(early), S(late) and linear in-between frames
)


@dataclass
class Phantom:
    signal: np.ndarray            # float64 [row, col, slice, time], as computed by the notebook (StFull)
    target_pe: np.ndarray         # float64 [row, col], PE (%) the notebook encoded (baseline_images['PE'])
    target_ser: np.ndarray        # float64 [row, col], SER the notebook encoded (baseline_images['SER'])
    t1: np.ndarray                # float64 [row, col], baseline T1 (ms), used to name tissues
    s0: np.ndarray                # float64 [row, col], pre-contrast SPGR signal (baseline_images['St'])
    resolution_mm: tuple          # (row, col, slice) voxel size in mm
    pre_index: int
    early_index: int
    late_index: int
    frame_interval_s: float       # time between frames written to DICOM (AcquisitionTime / TriggerTime)
    tissue_t1: dict               # tissue name -> T1 (ms)
    sample_tissues: dict          # region definitions from the notebook
    notebook_sha256: str
    namespace: dict = field(repr=False, default_factory=dict)

    @property
    def n_frames(self) -> int:
        return self.signal.shape[-1]

    @property
    def voxel_volume_mm3(self) -> float:
        return float(np.prod(self.resolution_mm))

    @property
    def saved_uint16(self) -> np.ndarray:
        """Pixel values as written to DICOM by the notebook (np.uint16 cast, no rescaling)."""
        return np.uint16(self.signal)

    def tissue_name_map(self) -> np.ndarray:
        names = np.full(self.t1.shape, "background", dtype=object)
        names[self.t1 == 1] = "phantom_fill"  # interior of the outer frame (T1 = 1 ms in the notebook)
        for name, t1 in self.tissue_t1.items():
            names[self.t1 == t1] = name
        return names


def notebook_code_cells(notebook_path: Path = NOTEBOOK_PATH) -> list[str]:
    nb = json.loads(Path(notebook_path).read_text())
    return ["".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code"]


def build_phantom(n_slices: int | None = None, notebook_path: Path = NOTEBOOK_PATH) -> Phantom:
    """Execute the phantom-building cells of the notebook and return the result.

    n_slices: number of slices to generate (the notebook uses 120). Only the size of the 4D array
    changes; geometry in-plane, signal model and timing are untouched.
    """
    notebook_path = Path(notebook_path)
    cells = notebook_code_cells(notebook_path)
    selected = []
    for marker in CELL_MARKERS:
        hits = [i for i, src in enumerate(cells) if marker in src]
        if len(hits) != 1:
            raise RuntimeError(f"Marker {marker!r} found in {len(hits)} notebook cells (expected 1); "
                               "the notebook has changed and CELL_MARKERS must be reviewed.")
        selected.append(hits[0])
    if selected != sorted(selected):
        raise RuntimeError("Notebook cells are no longer in the expected order; review CELL_MARKERS.")

    import cv2  # the notebook draws the regions with OpenCV
    ns: dict = {"np": np, "cv2": cv2}
    with redirect_stdout(io.StringIO()):
        for idx in selected:
            src = cells[idx]
            if "image_pars = {" in src:
                exec(src, ns)
                if n_slices is not None:
                    # image_dim is a list used by the later cells; Tacq has already been computed from
                    # the original 120 slices, so the timing written to DICOM is unchanged.
                    ns["image_dim"][-1] = int(n_slices)
            else:
                exec(src, ns)

    images = ns["baseline_images"]
    return Phantom(
        signal=ns["StFull"],
        target_pe=np.array(images["PE"], dtype=float),
        target_ser=np.array(images["SER"], dtype=float),
        t1=np.array(images["T1"], dtype=float),
        s0=np.array(images["St"], dtype=float),
        resolution_mm=tuple(float(r) for r in ns["image_res"]),
        pre_index=int(ns["baseline_pre_contrast_index"]),
        early_index=int(ns["early_post_contrast_index"]),
        late_index=int(ns["late_post_contrast_index"]),
        # The DICOM cell writes AcquisitionTime = start + tind * Tacq seconds, TriggerTime = tind * Tacq * 1000 ms
        frame_interval_s=float(ns["image_pars"]["Tacq"]),
        tissue_t1={name: pars["T1"] for name, pars in ns["tissue_pars"].items()},
        sample_tissues=ns["sample_tissues"],
        notebook_sha256=hashlib.sha256(notebook_path.read_bytes()).hexdigest(),
        namespace=ns,
    )


def to_slicer_order(array_rcst: np.ndarray) -> np.ndarray:
    """[row, col, slice, time] -> [time, slice, row, col], the order of slicer.util.arrayFromVolume per frame
    (k, j, i) stacked by quantificationLogic.getVolumeDataFromSequence."""
    return np.transpose(array_rcst, (3, 2, 0, 1))


def load_dicom_series(dicom_dir: str | Path) -> np.ndarray:
    """Read the DICOM files written by the notebook into [time, slice, row, col] (stored pixel values).

    Used to prove that the dataset loaded into 3D Slicer is identical to the regenerated phantom.
    """
    import pydicom

    frames: dict[tuple[int, float], np.ndarray] = {}
    for path in sorted(Path(dicom_dir).rglob("*.dcm")):
        ds = pydicom.dcmread(path)
        key = (int(ds.TemporalPositionIdentifier), float(ds.ImagePositionPatient[2]))
        frames[key] = ds.pixel_array.astype(np.uint16)
    if not frames:
        raise FileNotFoundError(f"No .dcm files found under {dicom_dir}")
    times = sorted({k[0] for k in frames})
    zs = sorted({k[1] for k in frames})
    out = np.zeros((len(times), len(zs)) + next(iter(frames.values())).shape, dtype=np.uint16)
    for (t, z), img in frames.items():
        out[times.index(t), zs.index(z)] = img
    return out
