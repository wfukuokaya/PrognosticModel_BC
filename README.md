# Prognostic model for patients with bladder cancer undergoing radical cystectomy

This repository contains the R code used in the study "Development and Validation of a Model to Predict Outcomes of Radical Cystectomy in Patients with Bladder Cancer." The study aimed to develop a model to predict cancer-specific survival in patients with operable bladder cancer undergoing radical cystectomy. The model was developed using data from patients with bladder cancer who had undergone radical cystectomy without receiving adjuvant chemotherapy across 36 institutions in the Japan Urological Oncology Group. The total cohort was split in a 2:1 ratio into the development cohort and the validation cohort.

## Repository contents

| File | Description |
|---|---|
| `main_script.R` | Analysis script (R) |

## Requirements

R (version 4.x) with the following packages: `broom`, `dcurves`, `ggfortify`, `ggsurvfit`, `glmnet`, `gtsummary`, `lubridate`, `missForest`, `patchwork`, `pec`, `rms`, `survival`, `survminer`, `tidymodels`, `tidyverse`, `timeROC`.

## Data availability

This repository contains analysis code only. The patient-level data analysed in the study are not distributed here.

## License

The code in this repository is released under the [MIT License](LICENSE).
