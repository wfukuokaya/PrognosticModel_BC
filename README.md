# Prognostic model for radical cystectomy in bladder cancer

R code for the analysis reported in:

> Development and Validation of a Model to Predict Outcomes of Radical Cystectomy in Patients with Bladder Cancer.

Development and validation of a model for cancer-specific survival after radical cystectomy without adjuvant chemotherapy, using Japan Urological Oncology Group data from 36 institutions split 2:1 into development and validation cohorts.

## Contents

`main_script.R`: analysis.

## Requirements

R 4.x with `broom`, `dcurves`, `ggfortify`, `ggsurvfit`, `glmnet`, `gtsummary`, `lubridate`, `missForest`, `patchwork`, `pec`, `rms`, `survival`, `survminer`, `tidymodels`, `tidyverse`, `timeROC`.

## Data

Patient-level data are not included in this repository.

## License

[MIT](LICENSE)
