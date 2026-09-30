# TMD sleep, pain, solar radiation, and vitamin D analyses

This repository reproduces the main and supplementary tables and figures for the study of sleep quality, psychological burden, pain, ambient solar radiation, and longitudinal serum 25-hydroxyvitamin D in patients with temporomandibular disorders.

Patient-level data are not included. The scripts expect the deidentified analysis workbook described in `data/README.md`.

## Repository structure

```text
R/                 Figure scripts
python/            Statistical analysis and table scripts
metadata/          Variable definitions used for Table S1
data/              Instructions for placing the analysis workbook
outputs/            Generated tables, figures, and supporting files
reference_outputs/  Aggregate tables used to verify reproduction
```

## Outputs

| Manuscript item | Script | Main output |
|---|---|---|
| Table 1 | `python/Table_1.py` | `outputs/tables/Table_1.csv` |
| Table 2 | `python/Table_2.py` | `outputs/tables/Table_2.csv` |
| Table 3 | `python/Table_3.py` | `outputs/tables/Table_3.csv` |
| Table 4 | `python/Table_4.py` | `outputs/tables/Table_4.csv` |
| Table S1 | `python/Table_S1.py` | `outputs/tables/Table_S1.csv` |
| Table S2 | `python/Table_S2.py` | `outputs/tables/Table_S2.csv` |
| Table S3 | `python/Table_S3.py` | `outputs/tables/Table_S3.csv` |
| Table S4 | `python/Table_S4.py` | `outputs/tables/Table_S4.csv` |
| Figure 1 | `R/Figure_1.R` | `outputs/figures/Figure_1.*` |
| Figure 2 | `R/Figure_2.R` | `outputs/figures/Figure_2.*` |
| Figure 3 | `R/Figure_3.R` | `outputs/figures/Figure_3.*` |
| Figure 4 | `R/Figure_4.R` | `outputs/figures/Figure_4.*` |
| Figure 5 | `R/Figure_5.R` | `outputs/figures/Figure_5.*` |
| Figure S1 | `R/Figure_S1.R` | `outputs/figures/Figure_S1.*` |
| Figure S2 | `R/Figure_S2.R` | `outputs/figures/Figure_S2.*` |

`Table_4.py` also generates Table S3 and the repeated cross-validation files required by Figure 5. Table S2 includes a separate file for the three-category VAS trajectory summary.

## Data setup

Place the workbook at:

```text
data/20260925_TMD_Sleep_VitaminD_Analysis_Master.xlsx
```

Alternatively, set the `TMD_SLEEP_DATA` environment variable to the workbook location.

PowerShell:

```powershell
$env:TMD_SLEEP_DATA = "C:\path\to\20260925_TMD_Sleep_VitaminD_Analysis_Master.xlsx"
```

macOS or Linux:

```bash
export TMD_SLEEP_DATA="/path/to/20260925_TMD_Sleep_VitaminD_Analysis_Master.xlsx"
```

All analyses read the `Analysis_Master` worksheet. The expected primary analysis cohort contains 120 participants: 62 good sleepers and 58 poor sleepers. SCL-90-R GSI and STOP-Bang scores are available for 92 and 101 participants, respectively. Missing questionnaire values are not imputed.

## Software

The table scripts were checked with Python 3.12. Install the required packages in an isolated environment:

```bash
python -m venv .venv
```

Windows:

```powershell
.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

macOS or Linux:

```bash
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Install the R packages with:

```bash
Rscript R/install_packages.R
```

The figure scripts use `readxl`, `readr`, `dplyr`, `tidyr`, `purrr`, `tibble`, `ggplot2`, `patchwork`, `scales`, `stringr`, `igraph`, `ggraph`, `ggforce`, `sandwich`, and `lmtest`. SVG export additionally uses `svglite`.

## Reproduction

Generate all tables:

```bash
python python/run_all_tables.py
```

Table 4 uses 100 repeated five-fold partitions and therefore takes longer than the other tables. After the table scripts finish, generate all figures:

```bash
Rscript R/run_all_figures.R
```

Figure S2 uses 2,000 stratified paired bootstrap resamples and is the slowest R analysis. Each script may also be run separately.

## Analysis conventions

Good sleepers were defined by PSQI global score ≤5 and poor sleepers by PSQI global score >5. Table 1 applies the Benjamini-Hochberg procedure across the 20 nondefinitional group comparisons; PSQI global score and sleep duration are shown descriptively because they define or contribute to sleep-group classification. Table 2 controls the three PSQI-score models and the three sleep-group models as separate analysis families. Table 3 controls separate and integrated follow-up 25(OH)D models as separate families.

Linear-model inference uses HC3 heteroscedasticity-consistent standard errors. Count and binary outcomes use robust Poisson and logistic models, respectively. Table 4 uses ridge regression with standardization and penalty selection confined to the training data. Model comparisons use paired Nadeau-Bengio corrected repeated-cross-validation tests, followed by Benjamini-Hochberg correction within each outcome.

Figure 2 excludes structural or part-whole correlations from FDR-based substantive interpretation. Figure S2 reports order-independent Shapley/dominance contributions with 95% intervals from 2,000 bootstrap resamples; it describes explanatory contribution rather than out-of-sample prediction.

## Reproduction checks

Successful execution should recover these denominators:

- PSQI analysis cohort: 120
- Good sleepers: 62
- Poor sleepers: 58
- SCL-90-R GSI available: 92
- STOP-Bang available: 101
- Fully integrated complete-case cohort for Figure S2: 81

Aggregate reference tables are provided in `reference_outputs/tables`. Small differences beyond the displayed precision usually indicate a different software version, data revision, missing-value rule, or cross-validation seed.

After generating the tables, compare them with the archived aggregate results:

```bash
python python/validate_outputs.py
```

## License

Code is released under the MIT License. Access to the study data remains subject to the approvals and data-sharing conditions described in the manuscript.
