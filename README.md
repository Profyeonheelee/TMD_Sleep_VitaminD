# TMD_Sleep_VitaminD

Analysis code accompanying the manuscript:

**Sleep Quality and Subsequent Vitamin D Status in Temporomandibular Disorders: Psychological, Pain, and Environmental Correlates**

This retrospective cohort study examines sleep quality, psychological burden, pain, ambient solar radiation, and subsequent serum 25-hydroxyvitamin D [25(OH)D] status in patients with temporomandibular disorders (TMD). The scripts generate the main and supplementary tables and figures.

Participant-level study data are not included in this public repository. Numerical reproduction requires access to the original analysis workbook, subject to the ethical approvals and data-sharing conditions described in the manuscript.

## Repository contents

All scripts are stored directly in the repository root.

| Files | Purpose |
|---|---|
| `common.py` | Shared Python data-loading, statistical, and output functions |
| `tables.py` | Shared functions used by the non-machine-learning table scripts |
| `Table_1.py`, `Table_2.py`, `Table_3.py` | Main Tables 1–3 |
| `Table_4.py` | Ridge-regression evaluation, Table 4, Table S3, and supporting results for Figure 5 |
| `Table_S1.py`, `Table_S2.py`, `Table_S3.py`, `Table_S4.py` | Supplementary Tables S1–S4 |
| `run_all_tables.py` | Runs the table analyses |
| `common.R` | Shared R input and output functions |
| `Figure_1.R` through `Figure_5.R` | Main Figures 1–5 |
| `Figure_S1.R`, `Figure_S2.R` | Supplementary Figures S1 and S2 |
| `install_packages.R` | Installs the R packages used by the figure scripts |

The shared table module must be named **`tables.py`**, in lowercase, and placed alongside the other Python scripts.

## Software and installation

The manuscript reports Python 3.12 and R 4.5.1. Package versions are not currently pinned in this repository.

Run the following commands from the repository root.

### Python

Create a virtual environment:

```bash
python -m venv .venv
```

Activate it on macOS or Linux:

```bash
source .venv/bin/activate
```

Or activate it in Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

Install the Python packages used by the available analysis scripts:

```bash
python -m pip install numpy pandas scipy scikit-learn openpyxl
```

### R

```bash
Rscript install_packages.R
```

The installer includes `readxl`, `readr`, `dplyr`, `tidyr`, `purrr`, `tibble`, `ggplot2`, `patchwork`, `scales`, `stringr`, `igraph`, `ggraph`, `ggforce`, `sandwich`, `lmtest`, and `svglite`.

## Project paths

For the flat repository layout, use the following project-root settings so that input data and generated outputs are located within the repository.

In `common.py`:

```python
PROJECT_ROOT = Path(__file__).resolve().parent
```

In `common.R`:

```r
project_dir <- normalizePath(script_dir, mustWork = FALSE)
```

The paths described below assume these settings. Helper versions using `parents[1]` in Python or `file.path(script_dir, "..")` in R instead resolve the project root one directory above the repository.

## Data setup

Place the analysis workbook at:

```text
data/20260925_TMD_Sleep_VitaminD_Analysis_Master.xlsx
```

Create the local `data/` directory if needed. All scripts read the `Analysis_Master` worksheet. The workbook must retain the variable names and coding expected by the scripts, including `Study_ID`, `Eligible_PSQI_analysis`, and `PSQI_global_score`.

Alternatively, set `TMD_SLEEP_DATA` to the workbook's absolute path. This overrides the input location without changing the output location.

macOS or Linux:

```bash
export TMD_SLEEP_DATA="/path/to/20260925_TMD_Sleep_VitaminD_Analysis_Master.xlsx"
```

Windows PowerShell:

```powershell
$env:TMD_SLEEP_DATA = "C:\path\to\20260925_TMD_Sleep_VitaminD_Analysis_Master.xlsx"
```

The expected PSQI analysis cohort contains 120 participants, comprising 62 good sleepers and 58 poor sleepers. SCL-90-R GSI is available for 92 participants and STOP-Bang for 101 participants. The fully integrated complete-case cohort contains 81 participants. Missing questionnaire values are not imputed, and analysis-specific sample sizes depend on the required variables.

## Running the analyses

### Tables

After adding `tables.py` and configuring the data path, generate the tables with:

```bash
python run_all_tables.py
```

Individual table scripts can also be run directly, for example:

```bash
python Table_1.py
python Table_4.py
```

`Table_4.py` also generates Table S3 and the supporting files required by Figure 5. It uses 100 repetitions of five-fold cross-validation and may take longer than the other table analyses.

### Figures

Run each figure script with `Rscript`:

```bash
Rscript Figure_1.R
Rscript Figure_2.R
Rscript Figure_3.R
Rscript Figure_4.R
Rscript Figure_5.R
Rscript Figure_S1.R
Rscript Figure_S2.R
```

Run `Table_4.py` before `Figure_5.R`, because Figure 5 reads its saved cross-validation results and participant-level predictions.

`Figure_S2.R` performs the dominance analysis and 2,000 stratified paired bootstrap resamples in R. This analysis may take additional time.

## Generated outputs

With the project-root settings above, the scripts create the following directories as needed:

| Directory | Contents |
|---|---|
| `outputs/tables/` | Main and supplementary tables in CSV format |
| `outputs/figures/` | Figure files and figure-specific statistical summaries |
| `outputs/supporting/` | Cross-validation metrics, model summaries, and predictions used by Figure 5 |

Figure export formats vary by script and include PDF, PNG, TIFF, and SVG. Supporting outputs can contain participant-level predictions and are intended for local use.

## Analysis notes

- Good sleepers have a PSQI global score of ≤5; poor sleepers have a score of >5.
- Correlation analyses use available pairs. Regression and prediction analyses use complete cases for the variables required by each analysis.
- Linear-model inference uses HC3 heteroscedasticity-consistent standard errors where specified.
- Ridge-regression standardization and penalty selection are performed within the training data. Table 4 uses fixed seeds and the same outer cross-validation partitions for models compared within each outcome.
- Repeated-cross-validation model comparisons use the Nadeau–Bengio correction, with Benjamini–Hochberg adjustment within each outcome.
- Figure S2 reports explanatory contributions to model R² and bootstrap uncertainty. These contributions do not represent out-of-sample predictive performance.

Refer to the manuscript Methods and table legends for the full model specifications, adjustment variables, and multiple-comparison families. Exact numerical reproduction depends on the original workbook and software environment.

