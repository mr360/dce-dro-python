# FTV validation tools

Tools to validate and audit the FTV, PE and SER values computed by the SeQ-DCEMRI 3D Slicer extension
([parametricDCEMRI](https://github.com/mr360/parametricDCEMRI)) using the FTV digital phantom
(`../FTV_digitalPhantom.ipynb`).

**Start with [METHODOLOGY.md](METHODOLOGY.md)**. It describes the method step by step with links to the
exact lines of code, lists the audit findings and gives draft methods text. The latest results are in
[results/report.md](results/report.md).

| File | Role |
|---|---|
| `phantom.py` | Rebuilds the phantom by running the notebook's own code cells (no copy of the phantom code). Gives the float signal, the uint16 values written to DICOM, and the encoded PE/SER maps |
| `ground_truth.py` | Expected outcome derived only from the encoded PE/SER values and the module's decision rules |
| `reference_ftv.py` | Slicer-free, line-by-line transcription of `quantificationLogic.process`, annotated with source line numbers |
| `validate_phantom.py` | Runs all scenarios and sweeps and writes `results/` (report and CSV files) |
| `slicer/export_plugin_results.py` | Run inside 3D Slicer after processing; exports the module's maps, tables and settings |
| `slicer/compare_with_plugin.py` | Checks an export voxel by voxel against the reference implementation |
| `tests/` | Unit tests for each decision rule, the phantom recovery and the comparison script |

```bash
pip install -r validation/requirements.txt
python -m pytest validation/tests
python validation/validate_phantom.py            # ~10 min for the full 120-slice phantom
```

The tools do not modify the extension. Where the extension behaves unexpectedly, the reference
implementation reproduces that behaviour so outputs can be compared exactly, and METHODOLOGY.md §4
records it.
