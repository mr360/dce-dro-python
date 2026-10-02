# Validation of SeQ-DCEMRI FTV computations with the FTV digital phantom

Generated 2026-10-02T13:47:03 by `validation/validate_phantom.py`.

## Provenance

| Item | Value |
|---|---|
| Phantom notebook SHA-256 | 33d815a1d1860612650f50fb793eb5eb4ae2d644d1a0ea74a195ee2c3436229e |
| dce-dro-python commit | 27ba0ae0e319d6f36893f716b0720a5364fc3eea |
| parametricDCEMRI commit audited | 7490646dfb019dbdfb9d8789054e5333d17837f7 |
| Python / numpy / scipy | 3.11.15 / 2.4.6 / 1.17.1 |
| Phantom size [row, col, slice, time] | 528×528×120×8 |
| Voxel size (mm) | 0.68 × 0.68 × 1.5 = 0.6936 mm³ |
| Pre / early / late frame | 0 / 1 / 7 |
| Frame interval written to DICOM | 64.577 s |

## ROI boxes used (IJK, exclusive upper bound)

| ROI | IJK min | IJK max |
|---|---|---|
| breasts | (134, 130, 0) | (394, 234, 120) |
| left_breast | (134, 130, 0) | (218, 234, 120) |
| right_breast | (310, 130, 0) | (394, 234, 120) |
| whole_phantom | (113, 130, 0) | (415, 432, 120) |
| module_default_box | (198, 198, 44) | (331, 331, 75) |

## Scenario results (voxel counts)

Columns give *expected / float / DICOM* voxel counts. *Expected* comes from the encoded PE/SER only; *float* and *DICOM* come from the reference implementation of the module applied to the noiseless floating-point signal and to the uint16 pixel values stored in DICOM. Unless stated otherwise: frames 0/1/7, background 60 %, PE ≥ 70 %, pre-defined SER range.

*On-threshold voxels* belong to regions whose encoded PE or SER equals a decision threshold of that scenario (e.g. SER = 0.9 with the 0.9 cut-off); their classification is decided by round-off and has no defined ground truth. The last column compares every other voxel individually (label map, FTV, FTVth and ETV membership).

| Scenario | ROI | FTV (SER>0) | FTVth | ETV exp / DICOM | FTV cm³ (DICOM) | All counts equal | On-threshold voxels | Voxel-wise agreement off-threshold |
|---|---|---|---|---|---|---|---|---|
| A  defaults, preset SER range | breasts | 778320 / 778320 / 778320 | 389520 / 593520 / 593520 | 2150400 / 2150400 | 539.843 | no | 552960 | yes |
| B  PE threshold 0 | breasts | 1563600 / 1563600 / 1563600 | 778320 / 1161840 / 1161840 | 2150400 / 2150400 | 1084.513 | no | 1046880 | yes |
| C  PE threshold 50 | whole_phantom | 1556640 / 1556640 / 1556640 | 778320 / 1161840 / 1161840 | 9349680 / 9349680 | 1079.686 | no | 7457520 | yes |
| D  background 0, PE 0 | whole_phantom | 8504400 / 8504400 / 8504400 | 1046880 / 1594560 / 1515720 | 9349680 / 9349680 | 5898.652 | no | 1046880 | yes |
| E  single SER threshold 1.0 | breasts | 778320 / 778320 / 778320 | 389520 / 389520 / 389520 | 2150400 / 2150400 | 539.843 | yes | 0 | yes |
| F  single SER threshold 1.25 | breasts | 1563600 / 1556640 / 1563600 | 394800 / 394800 / 394800 | 2150400 / 2150400 | 1084.513 | yes | 6960 | yes |
| G  Washout (CAD) map | breasts | 1613280 / 1606320 / 1613280 | 828000 / 828000 / 828000 | 2150400 / 2150400 | 1118.971 | yes | 6960 | yes |
| H  left breast only | left_breast | 779040 / 778320 / 779040 | 383520 / 383520 / 383520 | 1047600 / 1047600 | 540.342 | yes | 720 | yes |
| I  right breast only | right_breast | 779040 / 778320 / 779040 | 394800 / 778320 / 778320 | 1047600 / 1047600 | 540.342 | no | 1047600 | yes |
| J  module default RefBox | module_default_box | 15376 / 15376 / 15376 | 15376 / 15376 / 15376 | 459947 / 459947 | 10.665 | yes | 23343 | yes |

## Region breakdown – scenario A  defaults, preset SER range

Each region is a (tissue, target PE, target SER) combination inside the ROI. Measured ranges are from the DICOM (uint16) data.

| Tissue | PE | SER | Voxels | S0 float→DICOM | PE measured | SER measured | SER float | Expected | DICOM | On threshold | Agree |
|---|---|---|---|---|---|---|---|---|---|---|---|
| background | 0.000 | 0.000 | 1094400 | 14.25→14 | 0.00–0.00 | -0.0000–-0.0000 | -0.0000000 | excluded: 1094400 | excluded: 1094400 |  | yes |
| breast_fat | 60.000 | 0.700 | 53640 | 84.16→84 | 59.52–59.52 | 0.6944–0.6944 | 0.7000022 | excluded: 53640 | excluded: 53640 |  | yes |
| breast_fat | 60.000 | 0.900 | 78840 | 84.16→84 | 59.52–59.52 | 0.8929–0.8929 | 0.9000025 | excluded: 78840 | excluded: 78840 |  | yes |
| breast_fat | 60.000 | 1.200 | 85320 | 84.16→84 | 59.52–59.52 | 1.1905–1.1905 | 1.2000030 | excluded: 85320 | excluded: 85320 |  | yes |
| breast_fat | 60.000 | 1.300 | 50760 | 84.16→84 | 59.52–59.52 | 1.2821–1.2821 | 1.3000031 | excluded: 50760 | excluded: 50760 |  | yes |
| breast_fat | 80.000 | 0.700 | 50760 | 84.16→84 | 79.76–79.76 | 0.6979–0.6979 | 0.7000019 | excluded: 50760 | excluded: 50760 |  | yes |
| breast_fat | 80.000 | 0.900 | 85320 | 84.16→84 | 79.76–79.76 | 0.9054–0.9054 | 0.9000021 | excluded: 85320 | excluded: 85320 | yes | yes |
| breast_fat | 80.000 | 1.200 | 78840 | 84.16→84 | 79.76–79.76 | 1.1964–1.1964 | 1.2000025 | excluded: 78840 | excluded: 78840 |  | yes |
| breast_fat | 80.000 | 1.300 | 53640 | 84.16→84 | 79.76–79.76 | 1.3137–1.3137 | 1.3000026 | excluded: 53640 | excluded: 53640 | yes | yes |
| breast_glandular | 60.000 | 0.700 | 210000 | 184.75→184 | 60.33–60.33 | 0.6981–0.6981 | 0.7000022 | excluded: 210000 | excluded: 210000 |  | yes |
| breast_glandular | 60.000 | 0.900 | 179520 | 184.75→184 | 60.33–60.33 | 0.9024–0.9024 | 0.9000025 | excluded: 179520 | excluded: 179520 |  | yes |
| breast_glandular | 60.000 | 1.200 | 204000 | 184.75→184 | 60.33–60.33 | 1.1935–1.1935 | 1.2000030 | excluded: 204000 | excluded: 204000 |  | yes |
| breast_glandular | 60.000 | 1.300 | 184800 | 184.75→184 | 60.33–60.33 | 1.2907–1.2907 | 1.3000032 | excluded: 184800 | excluded: 184800 |  | yes |
| breast_glandular | 80.000 | 0.700 | 184800 | 184.75→184 | 80.43–80.43 | 0.7014–0.7014 | 0.7000019 | 0.00 < SER ≤ 0.90: 184800 | 0.00 < SER ≤ 0.90: 184800 |  | yes |
| breast_glandular | 80.000 | 0.900 | 204000 | 184.75→184 | 80.43–80.43 | 0.9024–0.9024 | 0.9000021 | 0.00 < SER ≤ 0.90: 204000 | 0.90 < SER ≤ 1.30: 204000 | yes | **no** |
| breast_glandular | 80.000 | 1.200 | 179520 | 184.75→184 | 80.43–80.43 | 1.2033–1.2033 | 1.2000025 | 0.90 < SER ≤ 1.30: 179520 | 0.90 < SER ≤ 1.30: 179520 |  | yes |
| breast_glandular | 80.000 | 1.300 | 210000 | 184.75→184 | 80.43–80.43 | 1.2982–1.2982 | 1.3000026 | 0.90 < SER ≤ 1.30: 210000 | 0.90 < SER ≤ 1.30: 210000 | yes | yes |
| major_muscle | 50.000 | 0.500 | 6960 | 136.82→136 | 50.74–50.74 | 0.5036–0.5036 | 0.5000020 | excluded: 6960 | excluded: 6960 |  | yes |
| major_muscle | 120.000 | 3.100 | 49680 | 136.82→136 | 121.32–121.32 | 3.1132–3.1132 | 3.1000035 | excluded: 49680 | excluded: 49680 |  | yes |

## Region breakdown – scenario B  PE threshold 0

Each region is a (tissue, target PE, target SER) combination inside the ROI. Measured ranges are from the DICOM (uint16) data.

| Tissue | PE | SER | Voxels | S0 float→DICOM | PE measured | SER measured | SER float | Expected | DICOM | On threshold | Agree |
|---|---|---|---|---|---|---|---|---|---|---|---|
| background | 0.000 | 0.000 | 1094400 | 14.25→14 | 0.00–0.00 | -0.0000–-0.0000 | -0.0000000 | excluded: 1094400 | excluded: 1094400 |  | yes |
| breast_fat | 60.000 | 0.700 | 53640 | 84.16→84 | 59.52–59.52 | 0.6944–0.6944 | 0.7000022 | excluded: 53640 | excluded: 53640 |  | yes |
| breast_fat | 60.000 | 0.900 | 78840 | 84.16→84 | 59.52–59.52 | 0.8929–0.8929 | 0.9000025 | excluded: 78840 | excluded: 78840 | yes | yes |
| breast_fat | 60.000 | 1.200 | 85320 | 84.16→84 | 59.52–59.52 | 1.1905–1.1905 | 1.2000030 | excluded: 85320 | excluded: 85320 |  | yes |
| breast_fat | 60.000 | 1.300 | 50760 | 84.16→84 | 59.52–59.52 | 1.2821–1.2821 | 1.3000031 | excluded: 50760 | excluded: 50760 | yes | yes |
| breast_fat | 80.000 | 0.700 | 50760 | 84.16→84 | 79.76–79.76 | 0.6979–0.6979 | 0.7000019 | excluded: 50760 | excluded: 50760 |  | yes |
| breast_fat | 80.000 | 0.900 | 85320 | 84.16→84 | 79.76–79.76 | 0.9054–0.9054 | 0.9000021 | excluded: 85320 | excluded: 85320 | yes | yes |
| breast_fat | 80.000 | 1.200 | 78840 | 84.16→84 | 79.76–79.76 | 1.1964–1.1964 | 1.2000025 | excluded: 78840 | excluded: 78840 |  | yes |
| breast_fat | 80.000 | 1.300 | 53640 | 84.16→84 | 79.76–79.76 | 1.3137–1.3137 | 1.3000026 | excluded: 53640 | excluded: 53640 | yes | yes |
| breast_glandular | 60.000 | 0.700 | 210000 | 184.75→184 | 60.33–60.33 | 0.6981–0.6981 | 0.7000022 | 0.00 < SER ≤ 0.90: 210000 | 0.00 < SER ≤ 0.90: 210000 |  | yes |
| breast_glandular | 60.000 | 0.900 | 179520 | 184.75→184 | 60.33–60.33 | 0.9024–0.9024 | 0.9000025 | 0.00 < SER ≤ 0.90: 179520 | 0.90 < SER ≤ 1.30: 179520 | yes | **no** |
| breast_glandular | 60.000 | 1.200 | 204000 | 184.75→184 | 60.33–60.33 | 1.1935–1.1935 | 1.2000030 | 0.90 < SER ≤ 1.30: 204000 | 0.90 < SER ≤ 1.30: 204000 |  | yes |
| breast_glandular | 60.000 | 1.300 | 184800 | 184.75→184 | 60.33–60.33 | 1.2907–1.2907 | 1.3000032 | 0.90 < SER ≤ 1.30: 184800 | 0.90 < SER ≤ 1.30: 184800 | yes | yes |
| breast_glandular | 80.000 | 0.700 | 184800 | 184.75→184 | 80.43–80.43 | 0.7014–0.7014 | 0.7000019 | 0.00 < SER ≤ 0.90: 184800 | 0.00 < SER ≤ 0.90: 184800 |  | yes |
| breast_glandular | 80.000 | 0.900 | 204000 | 184.75→184 | 80.43–80.43 | 0.9024–0.9024 | 0.9000021 | 0.00 < SER ≤ 0.90: 204000 | 0.90 < SER ≤ 1.30: 204000 | yes | **no** |
| breast_glandular | 80.000 | 1.200 | 179520 | 184.75→184 | 80.43–80.43 | 1.2033–1.2033 | 1.2000025 | 0.90 < SER ≤ 1.30: 179520 | 0.90 < SER ≤ 1.30: 179520 |  | yes |
| breast_glandular | 80.000 | 1.300 | 210000 | 184.75→184 | 80.43–80.43 | 1.2982–1.2982 | 1.3000026 | 0.90 < SER ≤ 1.30: 210000 | 0.90 < SER ≤ 1.30: 210000 | yes | yes |
| major_muscle | 50.000 | 0.500 | 6960 | 136.82→136 | 50.74–50.74 | 0.5036–0.5036 | 0.5000020 | 0.00 < SER ≤ 0.90: 6960 | 0.00 < SER ≤ 0.90: 6960 |  | yes |
| major_muscle | 120.000 | 3.100 | 49680 | 136.82→136 | 121.32–121.32 | 3.1132–3.1132 | 3.1000035 | excluded: 49680 | excluded: 49680 |  | yes |

## Region breakdown – scenario G  Washout (CAD) map

Each region is a (tissue, target PE, target SER) combination inside the ROI. Measured ranges are from the DICOM (uint16) data.

| Tissue | PE | SER | Voxels | S0 float→DICOM | PE measured | SER measured | SER float | Expected | DICOM | On threshold | Agree |
|---|---|---|---|---|---|---|---|---|---|---|---|
| background | 0.000 | 0.000 | 1094400 | 14.25→14 | 0.00–0.00 | -0.0000–-0.0000 | -0.0000000 | excluded: 1094400 | excluded: 1094400 |  | yes |
| breast_fat | 60.000 | 0.700 | 53640 | 84.16→84 | 59.52–59.52 | 0.6944–0.6944 | 0.7000022 | excluded: 53640 | excluded: 53640 |  | yes |
| breast_fat | 60.000 | 0.900 | 78840 | 84.16→84 | 59.52–59.52 | 0.8929–0.8929 | 0.9000025 | excluded: 78840 | excluded: 78840 |  | yes |
| breast_fat | 60.000 | 1.200 | 85320 | 84.16→84 | 59.52–59.52 | 1.1905–1.1905 | 1.2000030 | excluded: 85320 | excluded: 85320 |  | yes |
| breast_fat | 60.000 | 1.300 | 50760 | 84.16→84 | 59.52–59.52 | 1.2821–1.2821 | 1.3000031 | excluded: 50760 | excluded: 50760 |  | yes |
| breast_fat | 80.000 | 0.700 | 50760 | 84.16→84 | 79.76–79.76 | 0.6979–0.6979 | 0.7000019 | excluded: 50760 | excluded: 50760 |  | yes |
| breast_fat | 80.000 | 0.900 | 85320 | 84.16→84 | 79.76–79.76 | 0.9054–0.9054 | 0.9000021 | excluded: 85320 | excluded: 85320 |  | yes |
| breast_fat | 80.000 | 1.200 | 78840 | 84.16→84 | 79.76–79.76 | 1.1964–1.1964 | 1.2000025 | excluded: 78840 | excluded: 78840 |  | yes |
| breast_fat | 80.000 | 1.300 | 53640 | 84.16→84 | 79.76–79.76 | 1.3137–1.3137 | 1.3000026 | excluded: 53640 | excluded: 53640 |  | yes |
| breast_glandular | 60.000 | 0.700 | 210000 | 184.75→184 | 60.33–60.33 | 0.6981–0.6981 | 0.7000022 | Persistent: 210000 | Persistent: 210000 |  | yes |
| breast_glandular | 60.000 | 0.900 | 179520 | 184.75→184 | 60.33–60.33 | 0.9024–0.9024 | 0.9000025 | Persistent: 179520 | Persistent: 179520 |  | yes |
| breast_glandular | 60.000 | 1.200 | 204000 | 184.75→184 | 60.33–60.33 | 1.1935–1.1935 | 1.2000030 | Washout: 204000 | Washout: 204000 |  | yes |
| breast_glandular | 60.000 | 1.300 | 184800 | 184.75→184 | 60.33–60.33 | 1.2907–1.2907 | 1.3000032 | Washout: 184800 | Washout: 184800 |  | yes |
| breast_glandular | 80.000 | 0.700 | 184800 | 184.75→184 | 80.43–80.43 | 0.7014–0.7014 | 0.7000019 | Persistent: 184800 | Persistent: 184800 |  | yes |
| breast_glandular | 80.000 | 0.900 | 204000 | 184.75→184 | 80.43–80.43 | 0.9024–0.9024 | 0.9000021 | Persistent: 204000 | Persistent: 204000 |  | yes |
| breast_glandular | 80.000 | 1.200 | 179520 | 184.75→184 | 80.43–80.43 | 1.2033–1.2033 | 1.2000025 | Washout: 179520 | Washout: 179520 |  | yes |
| breast_glandular | 80.000 | 1.300 | 210000 | 184.75→184 | 80.43–80.43 | 1.2982–1.2982 | 1.3000026 | Washout: 210000 | Washout: 210000 |  | yes |
| major_muscle | 50.000 | 0.500 | 6960 | 136.82→136 | 50.74–50.74 | 0.5036–0.5036 | 0.5000020 | Persistent: 6960 | Persistent: 6960 | yes | yes |
| major_muscle | 120.000 | 3.100 | 49680 | 136.82→136 | 121.32–121.32 | 3.1132–3.1132 | 3.1000035 | Washout: 49680 | Washout: 49680 |  | yes |

## Sweep: PE threshold (%) (ROI whole_phantom)

| pe_threshold | expected_FTV0 | dicom_FTV0 | expected_FTVth | dicom_FTVth | dicom_FTV0_cm3 | dicom_FTVth_cm3 | background_value_dicom |
|---|---|---|---|---|---|---|---|
| 0.000 | 1556640 | 1556640 | 778320 | 1161840 | 1079.686 | 805.852 | 139.200 |
| 10.000 | 1556640 | 1556640 | 778320 | 1161840 | 1079.686 | 805.852 | 139.200 |
| 20.000 | 1556640 | 1556640 | 778320 | 1161840 | 1079.686 | 805.852 | 139.200 |
| 30.000 | 1556640 | 1556640 | 778320 | 1161840 | 1079.686 | 805.852 | 139.200 |
| 40.000 | 1556640 | 1556640 | 778320 | 1161840 | 1079.686 | 805.852 | 139.200 |
| 50.000 | 1556640 | 1556640 | 778320 | 1161840 | 1079.686 | 805.852 | 139.200 |
| 60.000 | 1556640 | 1556640 | 778320 | 1161840 | 1079.686 | 805.852 | 139.200 |
| 70.000 | 778320 | 778320 | 389520 | 593520 | 539.843 | 411.665 | 139.200 |
| 80.000 | 778320 | 778320 | 389520 | 593520 | 539.843 | 411.665 | 139.200 |
| 90.000 | 0 | 0 | 0 | 0 | 0.000 | 0.000 | 139.200 |
| 100.000 | 0 | 0 | 0 | 0 | 0.000 | 0.000 | 139.200 |
| 110.000 | 0 | 0 | 0 | 0 | 0.000 | 0.000 | 139.200 |
| 120.000 | 0 | 0 | 0 | 0 | 0.000 | 0.000 | 139.200 |
| 130.000 | 0 | 0 | 0 | 0 | 0.000 | 0.000 | 139.200 |
| 140.000 | 0 | 0 | 0 | 0 | 0.000 | 0.000 | 139.200 |

## Sweep: Background threshold (%) (ROI whole_phantom)

| background_threshold | expected_FTV0 | dicom_FTV0 | expected_FTVth | dicom_FTVth | dicom_FTV0_cm3 | dicom_FTVth_cm3 | background_value_dicom |
|---|---|---|---|---|---|---|---|
| 0.000 | 8504400 | 8504400 | 1046880 | 1515720 | 5898.652 | 1051.303 | 0.000 |
| 10.000 | 8504400 | 8504400 | 1046880 | 1515720 | 5898.652 | 1051.303 | 23.200 |
| 20.000 | 8504400 | 8504400 | 1046880 | 1515720 | 5898.652 | 1051.303 | 46.400 |
| 30.000 | 2537280 | 2537280 | 1046880 | 1515720 | 1759.857 | 1051.303 | 69.600 |
| 40.000 | 2000160 | 2000160 | 778320 | 1161840 | 1387.311 | 805.852 | 92.800 |
| 50.000 | 2000160 | 2000160 | 778320 | 1161840 | 1387.311 | 805.852 | 116.000 |
| 60.000 | 1556640 | 1556640 | 778320 | 1161840 | 1079.686 | 805.852 | 139.200 |
| 70.000 | 1556640 | 1556640 | 778320 | 1161840 | 1079.686 | 805.852 | 162.400 |
| 80.000 | 0 | 0 | 0 | 0 | 0.000 | 0.000 | 185.600 |
| 90.000 | 0 | 0 | 0 | 0 | 0.000 | 0.000 | 208.800 |
| 100.000 | 0 | 0 | 0 | 0 | 0.000 | 0.000 | 232.000 |

## Sweep: Single SER threshold (ROI breasts)

| ser_threshold | expected_FTV0 | dicom_FTV0 | expected_FTVth | dicom_FTVth | dicom_FTV0_cm3 | dicom_FTVth_cm3 | background_value_dicom |
|---|---|---|---|---|---|---|---|
| 0.000 | 1563600 | 1563600 | 1563600 | 1563600 | 1084.513 | 1084.513 | 110.400 |
| 0.100 | 1563600 | 1563600 | 1563600 | 1563600 | 1084.513 | 1084.513 | 110.400 |
| 0.200 | 1563600 | 1563600 | 1563600 | 1563600 | 1084.513 | 1084.513 | 110.400 |
| 0.300 | 1563600 | 1563600 | 1563600 | 1563600 | 1084.513 | 1084.513 | 110.400 |
| 0.400 | 1563600 | 1563600 | 1563600 | 1563600 | 1084.513 | 1084.513 | 110.400 |
| 0.500 | 1563600 | 1563600 | 1556640 | 1563600 | 1084.513 | 1084.513 | 110.400 |
| 0.600 | 1563600 | 1563600 | 1556640 | 1556640 | 1084.513 | 1079.686 | 110.400 |
| 0.700 | 1563600 | 1563600 | 1161840 | 1346640 | 1084.513 | 934.030 | 110.400 |
| 0.800 | 1563600 | 1563600 | 1161840 | 1161840 | 1084.513 | 805.852 | 110.400 |
| 0.900 | 1563600 | 1563600 | 778320 | 1161840 | 1084.513 | 805.852 | 110.400 |
| 1.000 | 1563600 | 1563600 | 778320 | 778320 | 1084.513 | 539.843 | 110.400 |
| 1.100 | 1563600 | 1563600 | 778320 | 778320 | 1084.513 | 539.843 | 110.400 |
| 1.200 | 1563600 | 1563600 | 394800 | 574320 | 1084.513 | 398.348 | 110.400 |
| 1.300 | 1563600 | 1563600 | 0 | 0 | 1084.513 | 0.000 | 110.400 |
| 1.400 | 1563600 | 1563600 | 0 | 0 | 1084.513 | 0.000 | 110.400 |
| 1.500 | 1563600 | 1563600 | 0 | 0 | 1084.513 | 0.000 | 110.400 |
| 1.600 | 1563600 | 1563600 | 0 | 0 | 1084.513 | 0.000 | 110.400 |
| 1.700 | 1563600 | 1563600 | 0 | 0 | 1084.513 | 0.000 | 110.400 |
| 1.800 | 1563600 | 1563600 | 0 | 0 | 1084.513 | 0.000 | 110.400 |
| 1.900 | 1563600 | 1563600 | 0 | 0 | 1084.513 | 0.000 | 110.400 |
| 2.000 | 1563600 | 1563600 | 0 | 0 | 1084.513 | 0.000 | 110.400 |

## Sensitivity to the early post-contrast frame (ROI breasts, PE ≥ 50 %)

The phantom encodes PE/SER between frames 0, 1 and 7. The module's default early frame is 3 ([q68]); any other choice measures different PE and SER.

| Early frame | FTV voxels | FTVth voxels | mean PE in FTV | mean SER in FTV |
|---|---|---|---|---|
| 1 | 1563600 | 1161840 | 71.864 | 1.085 |
| 2 | 1613280 | 1211520 | 71.925 | 1.072 |
| 3 | 1613280 | 1211520 | 71.998 | 1.058 |
| 4 | 1613280 | 1211520 | 71.951 | 1.043 |
| 5 | 1613280 | 1396320 | 71.950 | 1.028 |
| 6 | 1428480 | 1428480 | 75.108 | 1.012 |

## Module tables reproduced for scenario A

### Summary Table

| Parameter | Value | Units |
|---|---|---|
| Peak SER | 1.298 | [] |
| Peak PE | 121.324 | % |
| PE Threshold | 70.000 | % |
| SER Upper Threshold | 0.900 | [] |
| ROI longest axis | 180.000 | mm |
| ROI Volume | 1178.408 | cm3 |
| Bolus injection time | 0.000 | min |
| Early Phase Time | 1.076 | min |
| Late Phase Time | 7.534 | min |
| Maximum Enhancement | 92.672 | % |
| Delta Enhancement | -2.789 | % |
| First Pass Enhancement | 82.888 | % |
| Enhancement Slope | -0.418 | [] |

### SER Table

| Ranges | Volume (cm3) | Distribution (%) |
|---|---|---|
| 0.00 < SER ≤ 0.90 | 128.177 | 23.740 |
| 0.90 < SER ≤ 1.30 | 411.665 | 76.260 |
| FTVth (SER>SER_th) | 411.665 | 76.260 |
| FTV (SER>0) | 539.843 | 100.000 |
| ETV (Enhanced Tumour Volume) | 1491.517 | 100.000 |

ROI statistics: segment volume 2250.593 cm³; value reported as 'ROI Volume' (× π/6) 1178.408 cm³.

