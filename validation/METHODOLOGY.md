# SeQ-DCEMRI: computation method, audit trail and validation

This document describes how the SeQ-DCEMRI 3D Slicer module (`parametricDCEMRI`) turns a DCE-MRI series
into PE, SER and FTV values. It ties each step to the line of code that performs it, and records how the
computation was validated with the FTV digital reference object (DRO) in this repository. It is written to
support the methods section of a paper. Section 7 gives draft text.

* **Software audited:** `mr360/parametricDCEMRI`, file `quantification/quantification.py`, commit
  [`7490646`](https://github.com/mr360/parametricDCEMRI/blob/7490646dfb019dbdfb9d8789054e5333d17837f7/quantification/quantification.py).
  Line references such as **L1696** point to that file at that commit.
* **Phantom:** `FTV_digitalPhantom.ipynb` (SHA-256 `33d815a1…229e`, commit `27ba0ae` of this repository).
* **Validation code:** `validation/` (this folder). Results: `validation/results/report.md`.

---

## 1. Processing pipeline

The module implements the three-time-point (3TP) method (Degani 1997; Furman-Haran 2002) as used for FTV
in ACRIN 6657 / I-SPY (Hylton 2012). All computations are in `quantificationLogic.process` (L1573–2033).
Inputs come from the GUI (`onApplyButton`, L639–692).

Notation: for voxel *v*, S₀ = signal in the pre-contrast frame, S₁ = early post-contrast frame,
S₂ = late post-contrast frame (frame indices chosen by the user). ε = 10⁻⁶ (L1234).

| Step | Operation | Equation / rule | Code |
|---|---|---|---|
| 1 | Load the 4D series as float64, `[t, k, j, i]` | — | L1424–1440, L1618 |
| 2 | Analysis mask *M*: the selected segment. If it is empty: the RefBox minus omit boxes | — | L1663–1677 |
| 3 | Crop the series and *M* to the RefBox | — | L1682–1688 |
| 4 | Background threshold over the **whole RefBox** | T_bkg = (B/100) · P₉₅(S₀ in RefBox) | L1696 |
| 5 | Keep voxels above background inside *M* | base = (S₀ ≥ T_bkg) ∧ M | L1698 |
| 6 | Percentage enhancement | PE = 100 · (S₁ − S₀) / (S₀ + ε) | L1705 |
| 7 | PE threshold (inclusive) | base ← base ∧ (PE ≥ PE_th) | L1707 |
| 8 | Signal enhancement ratio | SER = (S₁ − S₀) / (S₂ − S₀ + ε) | L1718 |
| 9 | Non-SER voxels | SER < 0 or SER > 3.0 → SER = 0 | L1719–1720 |
| 10 | SER label map, intervals **lb < SER ≤ ub** | preset: (0, 0.9], (0.9, 1.3], (1.3, 3.0]; single threshold T: (0, T], (T, 1.44T], (1.44T, 3.0] | L1739–1751, L1152–1212 |
| 11 | Kinetic-curve ("Washout") label map | Δ = 100 · (S₂ − S₁)/(S₁ − S₀ + ε): Persistent Δ > 10, Plateau \|Δ\| ≤ 10, Washout Δ < −10 | L1769–1784, L1128–1149 |
| 12 | Peak PE / peak SER | max over the means of non-overlapping 3×3×3 blocks of the masked PE / SER maps | L1812–1824 |
| 13 | FTV, FTVth, ETV voxel sets (below) | — | L1826–1846 |
| 14 | Volumes | voxel count × voxel volume (Slicer *Segment Statistics*) | L1849–1861 |
| 15 | Time–intensity curve | mean of 100 · (S(t) − S₀)/(S₀ + ε) over voxels in *base*; straight-line fit over frames 1…N−1 | L1866–1874 |

### Reported volumes

| Output | Voxel set (SER mode) | Voxel set (Washout mode) | Code |
|---|---|---|---|
| **FTV (SER>0)** | base ∧ 0 < SER ≤ 3.0 | base ∧ label ∈ {Persistent, Plateau, Washout} | L1828, L1838 |
| **FTVth (SER>SER_th)** | base ∧ SER_th < SER ≤ 3.0, with SER_th = 0.9 (preset) or T (single) | base ∧ label ∈ {Plateau, Washout} | L1829, L1839 |
| **ETV** | every RefBox voxel with PE > 0. No segment, background or PE-threshold restriction | same | L1832, L1845 |
| Bin volume / distribution | voxels in each label; distribution = count / FTV count × 100 | same | L1974–1986 |

### Parameters

| Parameter | GUI default | Range | Code |
|---|---|---|---|
| Pre / early / late frame | 0 / **3** / last | 0 … N−1 | L65–69 |
| Background threshold B | 60 % | 0–100 % | L133 |
| PE threshold | 70 % | 0–250 % | L132 |
| SER mode | pre-defined range | preset or single threshold | L140 |
| Single SER threshold T | 1.4 | 0 – 2.077 (= 3/(1 + 4/9)) | L134, L420 |
| SER upper limit | 3.0 (fixed) | — | L267 |
| Washout cut-off | 10 % (fixed) | — | L1146–1149 |

---

## 2. Validation design

### 2.1 Digital reference object
The DRO (`FTV_digitalPhantom.ipynb`) is a 528 × 528 × 120 × 8 series (0.68 × 0.68 × 1.5 mm). Rectangular
regions have known T1, T2* and M0. The baseline signal is the spoiled gradient-echo equation
(TR 5.2 ms, TE 2.5 ms, flip angle 10°):

S₀ = M₀ sin α · (1 − E₁)/(1 − cos α · E₁) · E₂,  E₁ = e^(−TR/T1), E₂ = e^(−TE/T2*)

Each region is assigned a target PE and SER, and the post-contrast signals are obtained by inverting the
module's equations:

S₁ = S₀ (1 + PE/100),  S₂ = S₀ + (S₁ − S₀)/SER   (`s1_pe`, `s2_ser`)

Frame 0 = S₀, frame 1 = S₁, frame 7 = S₂, and frames 2–6 lie on the straight line between S₁ and S₂. All
slices are identical. Each breast is divided into four quadrants by combining two PE values (60, 80 %) with
two SER values (left breast 0.7/1.2, right breast 0.9/1.3). The quadrants have glandular interiors and
fat contour rings. The liver (PE 50, SER 0.5) and heart (PE 120, SER 3.1) provide additional cases.
The series is written to DICOM as `uint16` with no rescaling, so the stored values are the signals
truncated to integers.

### 2.2 Three-tier comparison
1. **Expected (ground truth)** — `ground_truth.py`: the module's decision rules (Section 1) applied to the
   *encoded* PE and SER of every voxel. It needs no signal intensities except the noiseless S₀ for the
   background rule.
2. **Reference implementation** — `reference_ftv.py`: a line-by-line, Slicer-free transcription of
   `quantificationLogic.process`. It is applied (a) to the floating-point signal and (b) to the uint16
   values stored in DICOM.
3. **Module output** — `slicer/export_plugin_results.py` exports the maps and tables the module computed
   inside 3D Slicer. `slicer/compare_with_plugin.py` checks them voxel by voxel against tier 2 run on the
   same exported input (PE map, SER and Washout label maps, every table value).

Tiers 1 and 2 show that the method recovers the encoded values. Tiers 2 and 3 show that the software
performs the method as described.

### 2.3 Threshold ties
Some encoded values lie exactly on a decision threshold: SER 0.9 and 1.3 are bin edges, and PE 50 % is a
tie when the PE threshold is 50 %. For these voxels the outcome depends on round-off (ε in the module,
the 10⁻⁶ added in the notebook's `s2_ser`, and uint16 truncation), so ground truth is undefined. They are
flagged *on-threshold* and reported separately. Every other voxel is compared individually.

### 2.4 Scenarios
Ten parameter scenarios (defaults; PE thresholds 0 and 50 %; background 0 %; single SER thresholds 1.0 and
1.25; Washout mode; single-breast and default RefBoxes). Sweeps of the PE threshold (0–140 %), background
threshold (0–100 %) and single SER threshold (0–2.0). A test of the early-frame choice.
`python validation/validate_phantom.py` regenerates all of them.

---

## 3. Validation results (120-slice phantom, `results/report.md`)

* **Off-threshold voxels: 100 % agreement** between ground truth and the reference implementation on the
  DICOM values, in all 10 scenarios. This covers SER-bin label, FTV, FTVth and ETV membership of every voxel.
* **Recovered values** (DICOM data, all enhancing regions of the phantom). PE is within 1.32 percentage
  points of target (glandular tissue ≤ 0.43, e.g. 80 → 80.43). SER is within 0.018 of target (glandular
  tissue ≤ 0.0093, e.g. 0.7 → 0.7014, 1.3 → 1.2907). The largest errors occur in the low-signal fat and
  muscle contour rings, where integer truncation of S₀ ≈ 84–137 matters most. On the floating-point signal,
  PE and SER equal the targets to 10⁻⁶ relative (unit test `test_float_phantom_recovers_encoded_values`).
* **Default scenario** (frames 0/1/7, background 60 %, PE ≥ 70 %, preset SER, RefBox around both breasts):
  FTV = 778 320 voxels = **539.84 cm³**, identical to ground truth. FTVth = 593 520 voxels (411.67 cm³) against
  389 520 expected. The 204 000-voxel difference is exactly the glandular PE 80 / SER 0.9 quadrant. Its
  measured SER is 0.9024 after uint16 truncation (0.9000021 even without truncation), so it falls above the
  0.9 cut-off. This is a threshold tie (Section 2.3), not a computational error.
* Fat contour rings (S₀ = 84 against T_bkg ≈ 110) are always removed by the background rule, as the
  rule prescribes.
* **Module output (tier 3): pending.** Run `slicer/export_plugin_results.py` in 3D Slicer and
  `slicer/compare_with_plugin.py` (Section 6) to complete the software-equivalence check.

---

## 4. Behaviours to disclose (audit findings)

These follow from the code as written. They do not invalidate the FTV computation, but some affect how
results should be reported. Each is reproduced by the reference implementation and/or a unit test.

| # | Finding | Effect | Code |
|---|---|---|---|
| F1 | GUI default early frame is **3**, not 1 | With the DRO, any early frame ≠ 1 gives different PE/SER (report: *Sensitivity to the early post-contrast frame*). State the frames used. | L68 |
| F2 | Ties at thresholds: ε lowers PE slightly, so PE exactly equal to the threshold is **excluded**. SER exactly on a bin edge goes to the lower bin in exact arithmetic, but round-off decides in practice | Relevant for DRO design. Negligible for noisy clinical data | L1705, L1742 |
| F3 | Background threshold uses the 95th percentile of S₀ over the **RefBox**, not the segment | FTV depends on RefBox placement and on bright structures inside it (e.g. heart: T_bkg 139 vs 110) | L1696 |
| F4 | **ETV** = all RefBox voxels with PE > 0. Ignores the segment, the background threshold and the PE threshold | ETV is not a subset of the segment. Describe it accordingly | L1706, L1832 |
| F5 | Voxels with SER > 3 or SER < 0 are removed from FTV, but `base &= (SER ≥ 0)` is a no-op | These voxels still contribute to **Peak PE**, Maximum/Delta/First-pass enhancement and the time–intensity curve. In the DRO, Peak PE = 121.3 % comes from heart voxels excluded from FTV. In **Washout mode** they are counted in FTV (heart SER 3.1 → "Washout") | L1721, L1812–1874 |
| F6 | "ROI Volume" = Segment-Statistics volume × π/6, on the assumption that the statistic is a bounding-box volume | Slicer's `volume_cm3` is the voxel-counted segment volume, so the reported value is 52 % of it (DRO: 1178 vs 2251 cm³). Confirm with the tier-3 comparison | L1884–1890 |
| F7 | Summary row "SER Upper Threshold" shows the FTVth SER cut-off (0.9 or T), not 3.0 | Labelling only | L1904, L1918 |
| F8 | Peak PE/SER blocks start at index 1 of the cropped RefBox, and masked-out voxels count as 0 in the mean | Peak values depend on RefBox alignment and are lower at region edges. Not the maximum over all neighbourhoods | L1812–1824 |
| F9 | Without DICOM frame labels the time axis is `linspace(0, N min, N)`, so the spacing is N/(N−1) min | Affects "Enhancement Slope" and phase times only | L847 |
| F10 | RefBox thinner than 3 voxels along any axis → error | Usability | L1817–1824 |
| F11 | Washout comment says Δ is relative to S_early. The code divides by (S₁ − S₀), so Δ = 100 (1/SER − 1) | The Washout classes are SER classes with cut-offs 1/1.1 ≈ 0.909 and 1/0.9 ≈ 1.111 | L1768–1769 |
| F12 | FTV/ETV label volumes are written into a full-size node from a cropped array before import | Voxel counts (hence volumes) are unaffected; temporary segments are mispositioned. Confirmed by tier 3 if FTV matches | L1852 |
| F13 | Summary and SER tables gain three columns on every re-run | Read the last three columns, or restart the scene between analyses | L2011–2016 |

DRO-specific notes:

* The variable `Tacq` is labelled seconds but is computed in minutes. The DICOM files store a frame
  interval of 64.58 s. This affects only time labels.
* The outer frame uses T2* = 0.001 ms, which gives zero signal (the comment expects a bright frame).
* Because the SER targets 0.9 and 1.3 coincide with the preset bin edges, the right breast has no voxel
  with defined ground truth in preset mode. Moving those targets off the edges (e.g. 0.8 and 1.5), or
  scaling M₀ to reduce truncation, would make every quadrant a decisive test.

---

## 5. Reproducing the validation

```bash
pip install -r validation/requirements.txt
python -m pytest validation/tests              # 17 unit tests, ~5 s
python validation/validate_phantom.py           # full phantom, writes validation/results/ (~10 min)
python validation/validate_phantom.py --slices 6 --out /tmp/quick   # quick run (counts scale with slices)
python validation/validate_phantom.py --dicom-dir /path/to/FTV_DRO/800   # also verify the DICOM files
```

## 6. Software-equivalence check in 3D Slicer (tier 3)

1. Load the DRO DICOM series in 3D Slicer (the version you report) and open *Semi-Quantitative DCE-MRI*.
2. Set frames **0 / 1 / 7**, place the RefBox, choose the thresholds, press **Click to Process**.
3. In the Python console:
   ```python
   EXPORT_DIR = r"/path/to/export"
   exec(open(r"/path/to/dce-dro-python/validation/slicer/export_plugin_results.py").read())
   ```
4. `python validation/slicer/compare_with_plugin.py /path/to/export --phantom`
   Every line should read PASS. `--phantom` also checks that the data Slicer loaded is the regenerated DRO.
5. Repeat for each scenario you report, and keep the export folders and outputs as supplementary material.

---

## 7. Draft methods text

> *Semi-quantitative analysis.* PE, SER and FTV were computed with the SeQ-DCEMRI extension
> (commit 7490646) for 3D Slicer (version X), which implements the three-time-point method. For each voxel,
> PE = 100·(S₁ − S₀)/S₀ and SER = (S₁ − S₀)/(S₂ − S₀), where S₀, S₁ and S₂ are the signals in the
> pre-contrast, early (frame *e*, *t* = … min) and late (frame *l*, *t* = … min) post-contrast images.
> Voxels were analysed inside a user-drawn segmentation within a rectangular region of interest. Voxels were
> excluded if S₀ was below B % of the 95th percentile of S₀ within the region of interest, or if PE was
> below the PE threshold. Voxels with SER < 0 or SER > 3 were treated as non-enhancing for FTV. FTV was the
> volume of the remaining voxels with 0 < SER ≤ 3. FTVth additionally required SER > SER_th. Voxels were
> classified into SER ranges (0, 0.9], (0.9, 1.3] and (1.3, 3.0]. Volumes were computed as voxel count ×
> voxel volume.
>
> *Software validation.* The implementation was verified with a digital reference object in which every
> region has known PE and SER. Signals were generated from the spoiled gradient-echo equation, and the
> post-contrast signals were obtained by inverting the PE and SER definitions. Results were compared at
> three levels: (i) ground truth derived from the encoded values, (ii) an independent re-implementation of
> the algorithm applied to the stored DICOM values, and (iii) the output of the extension itself, voxel
> by voxel. Across ten parameter configurations and sweeps of the PE, background and SER thresholds, all
> voxels whose encoded values did not lie exactly on a decision threshold were classified identically to
> ground truth. After integer quantisation, encoded PE and SER were recovered to within 1.3 percentage
> points and 0.018 respectively (0.43 and 0.0093 in glandular-tissue regions).
