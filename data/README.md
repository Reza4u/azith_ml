# Data (not included)

Participant data are **not** in this repository. De-identified data are available on
reasonable request to researchers who meet icddr,b's criteria for access to confidential
data, through the Research Administration of icddr,b (Ms Shiblee Sayeed,
shiblee_s@icddrb.org).

To rerun the analysis, place the two approved files here:

| File | Used by | Contents |
|---|---|---|
| `azithromycin_pediatric_clinical_data.rds` | all analyses | one row per child; outcome `azithromycin_failure` ("Yes"/"No"; missing = outcome not determined) and the baseline variables listed below |
| `modeldata.xlsx` (sheet `Sheet1`) | Supplementary Tables S6a–S6b only | trial dataset with the 48-hour outcome components `blood_48h`, `watery_stool_48h`, `fever_48h`, `clin_deterior_48h` ("Yes"/"No") |

Baseline variables read from the `.rds` file: `age` (months), `sex`, `resident_area`
("Urban Area"/"Rural Area"), `diarrhea_dur` (days), `fever`, `cough`, `resp_rate`,
`straining`, `haz`, `waz`, `whz`, `antibiotic_home` (recorded regimen, 11 levels),
`color_RE` (free text), `consist_RE`, `blood_RE`, `mucus_RE`, and the microscopy bands
`rbc_RE`, `pcell_RE`, `wbc_RE` ("1-5" … "31+") and `macro_RE` ("1-2", "2-5").
The coding applied to them is in `R/01_data.R` and Supplementary Section S1.1 of the paper.

A different folder can be used by setting the environment variable `AZITH_DATA_DIR`.
