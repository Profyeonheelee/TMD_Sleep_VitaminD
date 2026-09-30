from __future__ import annotations

from pathlib import Path
from decimal import Decimal, ROUND_HALF_UP

import numpy as np
import pandas as pd
from scipy import stats

from common import (
    PROJECT_ROOT,
    SUPPORT_DIR,
    analysis_cohort,
    bh_adjust,
    format_effect,
    format_p,
    glm_hc3,
    ols_hc3,
    standardized_beta,
    write_table,
)


def _group_values(data: pd.DataFrame, variable: str):
    good = data.loc[data["PSQI_sleep_group"].eq("Good sleeper"), variable].dropna()
    poor = data.loc[data["PSQI_sleep_group"].eq("Poor sleeper"), variable].dropna()
    return good, poor


def _mean_sd(values: pd.Series, show_n: bool = False) -> str:
    def round_half_up(value):
        return Decimal(str(float(value))).quantize(Decimal("0.01"), rounding=ROUND_HALF_UP)

    text = f"{round_half_up(values.mean())} ± {round_half_up(values.std(ddof=1))}"
    return f"{text} (n={len(values)})" if show_n else text


def _count_percent(values: pd.Series, positive=1, denominator=None) -> str:
    values = values.dropna()
    denominator = len(values) if denominator is None else denominator
    count = int((values == positive).sum())
    return f"{count} ({100 * count / denominator:.1f})"


def table_1() -> pd.DataFrame:
    data = analysis_cohort()
    data["Female_sex"] = data["Sex"].eq("Female").astype(int)

    specification = [
        ("Demographic and clinical characteristics", "Age, years", "Age_years", "welch", "continuous", False),
        ("Demographic and clinical characteristics", "Female sex, n (%)", "Female_sex", "fisher", "binary", False),
        ("Demographic and clinical characteristics", "Symptom duration, months", "Symptom_duration_months", "welch", "continuous", False),
        ("Demographic and clinical characteristics", "Chronic TMD, n (%)", "Chronic_TMD_ge3m", "fisher", "binary", False),
        ("Demographic and clinical characteristics", "Clinical follow-up interval, months", "Clinical_followup_months", "welch", "continuous", False),
        ("Sleep and psychological characteristics", "PSQI global score", "PSQI_global_score", None, "continuous", False),
        ("Sleep and psychological characteristics", "Sleep duration, hours", "PSQI_sleep_duration_hours_std", None, "continuous", False),
        ("Sleep and psychological characteristics", "SCL-90-R GSI T-score", "SCL_GSI_T_harmonized", "welch", "continuous", True),
        ("Sleep and psychological characteristics", "STOP-Bang score", "STOPBANG_total_score", "welch", "continuous", True),
        ("Sleep and psychological characteristics", "Intermediate-to-high OSA risk (STOP-Bang ≥3), n (%)", "STOPBANG_score_ge3", "fisher", "binary_fraction", False),
        ("Baseline pain profile", "Baseline pain intensity, VAS", "Baseline_pain_intensity_VAS", "welch", "continuous", False),
        ("Baseline pain profile", "Baseline clinical pain burden, sites", "Baseline_pain_location_count", "wilcoxon", "continuous", False),
        ("Longitudinal pain outcomes", "Follow-up pain intensity, VAS", "Followup_pain_intensity_VAS", "welch", "continuous", False),
        ("Longitudinal pain outcomes", "Pain intensity improvement, VAS", "Pain_intensity_reduction_VAS", "welch", "continuous", False),
        ("Vitamin D characteristics", "Baseline serum 25(OH)D, ng/mL", "Baseline_25OHD_ng_mL", "welch", "continuous", False),
        ("Vitamin D characteristics", "Follow-up serum 25(OH)D, ng/mL", "Followup_25OHD_ng_mL", "welch", "continuous", False),
        ("Vitamin D characteristics", "Change in serum 25(OH)D, ng/mL", "Delta_25OHD_ng_mL", "welch", "continuous", False),
        ("Vitamin D characteristics", "Vitamin D prescription, n (%)", "Vitamin_D_prescription", "fisher", "binary", False),
        ("Environmental exposure", "90-day prebaseline sunshine, h/day", "Sunshine_preBL_90d_mean_hr", "welch", "continuous", False),
        ("Environmental exposure", "90-day prebaseline solar radiation, MJ/m²/day", "SolarRad_preBL_90d_mean_MJm2", "welch", "continuous", False),
        ("Environmental exposure", "Interval sunshine duration, h/day", "Sunshine_BL_to_FU_mean_hr", "welch", "continuous", False),
        ("Environmental exposure", "Interval solar radiation, MJ/m²/day", "SolarRad_BL_to_FU_mean_MJm2", "welch", "continuous", False),
    ]

    rows = []
    tested_p = []
    tested_row = []
    for section, label, variable, test, kind, show_n in specification:
        good, poor = _group_values(data, variable)
        if kind == "continuous":
            good_text = _mean_sd(good, show_n)
            poor_text = _mean_sd(poor, show_n)
        elif kind == "binary_fraction":
            good_text = f"{int(good.sum())}/{len(good)} ({100 * good.mean():.1f})"
            poor_text = f"{int(poor.sum())}/{len(poor)} ({100 * poor.mean():.1f})"
        else:
            good_text = _count_percent(good)
            poor_text = _count_percent(poor)

        p_value = np.nan
        if test == "welch":
            p_value = stats.ttest_ind(good, poor, equal_var=False).pvalue
        elif test == "wilcoxon":
            p_value = stats.mannwhitneyu(good, poor, alternative="two-sided").pvalue
        elif test == "fisher":
            contingency = pd.crosstab(
                data.loc[data["PSQI_sleep_group"].isin(["Good sleeper", "Poor sleeper"]), "PSQI_sleep_group"],
                data.loc[data["PSQI_sleep_group"].isin(["Good sleeper", "Poor sleeper"]), variable],
            )
            p_value = stats.fisher_exact(contingency.to_numpy()).pvalue

        row = {
            "Section": section,
            "Characteristic": label,
            "Good sleepers (n=62)": good_text,
            "Poor sleepers (n=58)": poor_text,
            "p-value": "—" if np.isnan(p_value) else format_p(p_value),
            "FDR q-value": "—",
        }
        rows.append(row)
        if not np.isnan(p_value):
            tested_p.append(p_value)
            tested_row.append(len(rows) - 1)

    adjusted = bh_adjust(tested_p)
    for row_index, q_value in zip(tested_row, adjusted):
        rows[row_index]["FDR q-value"] = format_p(q_value)

    result = pd.DataFrame(rows)
    write_table(result, "Table_1.csv")
    raw = pd.DataFrame({"row": tested_row, "p_value": tested_p, "FDR_q": adjusted})
    raw.to_csv(SUPPORT_DIR / "Table_1_statistics.csv", index=False)
    return result


def _table_2_row(data, outcome, predictor, model, outcome_label, predictor_label):
    covariates = [predictor, "Age_years", "Male", "Chronic_TMD_ge3m"]
    if model == "linear":
        fit = ols_hc3(data, outcome, covariates)
        effect = format_effect(
            "B", fit["coef"][predictor], fit["lower"][predictor], fit["upper"][predictor]
        )
        std_beta = standardized_beta(fit, predictor, outcome)
        std_text = f"{std_beta:.3f}"
    else:
        fit = glm_hc3(data, outcome, covariates, family="poisson")
        effect = format_effect(
            "IRR",
            np.exp(fit["coef"][predictor]),
            np.exp(fit["lower"][predictor]),
            np.exp(fit["upper"][predictor]),
        )
        std_text = "—"
    return {
        "Outcome": outcome_label,
        "Sleep measure": predictor_label,
        "n": fit["n"],
        "Adjusted effect (robust 95% CI)": effect,
        "Std. β": std_text,
        "p_raw": fit["p"][predictor],
    }


def table_2() -> pd.DataFrame:
    data = analysis_cohort()
    outcomes = [
        ("SCL_GSI_T_harmonized", "linear", "SCL-90-R GSI T-score"),
        ("Baseline_pain_intensity_VAS", "linear", "Baseline pain intensity, VAS"),
        ("Baseline_pain_location_count", "poisson", "Baseline clinical pain burden, sites"),
    ]
    predictor_sets = [
        ("PSQI_global_score", "PSQI global score"),
        ("PSQI_poor_sleeper_gt5", "Poor sleeper"),
    ]
    rows = []
    for predictor, predictor_label in predictor_sets:
        family_rows = [
            _table_2_row(data, outcome, predictor, model, label, predictor_label)
            for outcome, model, label in outcomes
        ]
        q_values = bh_adjust([row["p_raw"] for row in family_rows])
        for row, q_value in zip(family_rows, q_values):
            row["p-value"] = format_p(row.pop("p_raw"))
            row["FDR q-value"] = format_p(q_value)
            rows.append(row)
    result = pd.DataFrame(rows)
    write_table(result, "Table_2.csv")
    return result


def _vitamin_d_model(data, focal_predictors, seasonal_predictors=None):
    seasonal_predictors = seasonal_predictors or []
    adjustment = [
        "Baseline_25OHD_ng_mL",
        "Vitamin_D_prescription",
        "Age_years",
        "Male",
        "Clinical_followup_months",
    ]
    predictors = adjustment + focal_predictors + seasonal_predictors
    return ols_hc3(data, "Followup_25OHD_ng_mL", predictors)


def table_3() -> pd.DataFrame:
    data = analysis_cohort()
    focal = [
        ("PSQI_global_score", "PSQI global score"),
        ("PSQI_sleep_duration_hours_std", "Sleep duration"),
        ("SolarRad_BL_to_FU_mean_MJm2", "Interval solar radiation"),
    ]
    rows = []
    separate_rows = []
    for index, (variable, label) in enumerate(focal, start=1):
        fit = _vitamin_d_model(data, [variable])
        separate_rows.append((f"Separate model {index}", variable, label, fit))
    integrated = _vitamin_d_model(data, [item[0] for item in focal])
    integrated_rows = [("Integrated model", variable, label, integrated) for variable, label in focal]

    for family in (separate_rows, integrated_rows):
        q_values = bh_adjust([fit["p"][variable] for _, variable, _, fit in family])
        for (model_label, variable, predictor_label, fit), q_value in zip(family, q_values):
            rows.append({
                "Model": model_label,
                "Predictor": predictor_label,
                "n": fit["n"],
                "B (HC3 robust 95% CI), ng/mL": format_effect(
                    "B", fit["coef"][variable], fit["lower"][variable], fit["upper"][variable]
                ).replace("B=", ""),
                "Standardized β": f"{standardized_beta(fit, variable, 'Followup_25OHD_ng_mL'):.3f}",
                "p-value": format_p(fit["p"][variable]),
                "FDR q-value": format_p(q_value),
                "R² / adjusted R²": f"{fit['r2']:.3f} / {fit['adjusted_r2']:.3f}",
            })
    result = pd.DataFrame(rows)
    write_table(result, "Table_3.csv")
    return result


def table_s1() -> pd.DataFrame:
    data = analysis_cohort()
    definitions_path = PROJECT_ROOT / "metadata" / "variable_definitions.csv"
    definitions = pd.read_csv(definitions_path)
    availability = []
    for _, row in definitions.iterrows():
        variable = row["Analysis variable"]
        if variable == "Derived in analysis":
            available = len(data)
        elif variable == "VAS_reduction_ge40pct" or variable == "VAS_reduction_ge60pct":
            available = int(data["Baseline_pain_intensity_VAS"].gt(0).sum())
        else:
            available = int(data[variable].notna().sum())
        missing = len(data) - available
        availability.append((available, f"{missing} ({100 * missing / len(data):.1f})"))
    definitions["Available n"] = [x[0] for x in availability]
    definitions["Missing n (%)"] = [x[1] for x in availability]
    write_table(definitions, "Table_S1.csv")
    return definitions


def _observed_by_group(data, outcome, binary=False, positive=1):
    pieces = []
    for group, short in [("Good sleeper", "Good"), ("Poor sleeper", "Poor")]:
        values = data.loc[data["PSQI_sleep_group"].eq(group), outcome].dropna()
        if binary:
            pieces.append(f"{short}: {int((values == positive).sum())}/{len(values)} ({100 * (values == positive).mean():.1f}%)")
        else:
            pieces.append(f"{short}: {values.mean():.2f} ± {values.std(ddof=1):.2f}")
    return "; ".join(pieces)


def _pain_outcome_model(data, outcome, predictor, family):
    predictors = [
        predictor,
        "Baseline_pain_intensity_VAS",
        "Baseline_pain_location_count",
        "Age_years",
        "Male",
        "Chronic_TMD_ge3m",
        "Clinical_followup_months",
    ]
    if outcome == "Followup_pain_intensity_VAS":
        predictors = [p for p in predictors if p != outcome]
    if family == "linear":
        return ols_hc3(data, outcome, predictors)
    return glm_hc3(data, outcome, predictors, family="binomial")


def table_s2() -> pd.DataFrame:
    data = analysis_cohort()
    rows = []
    specifications = [
        ("Continuous pain outcome", "Follow-up VAS", "Followup_pain_intensity_VAS", "linear", False),
        ("VAS responder outcomes", "≥40% VAS improvement", "VAS_reduction_ge40pct", "binomial", True),
        ("VAS responder outcomes", "≥60% VAS improvement", "VAS_reduction_ge60pct", "binomial", True),
        ("Pain worsening", "Follow-up VAS greater than baseline", "Pain_worsening", "binomial", True),
    ]
    predictors = [
        ("PSQI_global_score", "PSQI global score"),
        ("PSQI_poor_sleeper_gt5", "Poor-sleeper status"),
    ]
    for section, outcome_label, outcome, family, binary in specifications:
        model_data = data.copy()
        if outcome.startswith("VAS_reduction"):
            model_data = model_data.loc[model_data["Baseline_pain_intensity_VAS"].gt(0)].copy()
        for predictor, predictor_label in predictors:
            fit = _pain_outcome_model(model_data, outcome, predictor, family)
            if family == "linear":
                effect = format_effect("B", fit["coef"][predictor], fit["lower"][predictor], fit["upper"][predictor])
            else:
                effect = format_effect(
                    "OR",
                    np.exp(fit["coef"][predictor]),
                    np.exp(fit["lower"][predictor]),
                    np.exp(fit["upper"][predictor]),
                )
            observed = "—" if predictor == "PSQI_global_score" else _observed_by_group(model_data, outcome, binary=binary)
            rows.append({
                "Section": section,
                "Outcome": outcome_label,
                "Sleep measure": predictor_label,
                "n": fit["n"],
                "Observed outcome by sleep group": observed,
                "Adjusted effect (robust 95% CI)": effect,
                "p_raw": fit["p"][predictor],
            })
    q_values = bh_adjust([row["p_raw"] for row in rows])
    for row, q_value in zip(rows, q_values):
        row["p-value"] = format_p(row.pop("p_raw"))
        row["FDR q-value"] = format_p(q_value)
    result = pd.DataFrame(rows)
    write_table(result, "Table_S2.csv")

    trajectory = []
    for label, relation in [
        ("Decreased", lambda x: x < 0),
        ("Unchanged", lambda x: x == 0),
        ("Increased", lambda x: x > 0),
    ]:
        row = {"VAS trajectory": label}
        change = data["Followup_pain_intensity_VAS"] - data["Baseline_pain_intensity_VAS"]
        for group, column in [("Good sleeper", "Good sleepers (n=62)"), ("Poor sleeper", "Poor sleepers (n=58)")]:
            mask = data["PSQI_sleep_group"].eq(group)
            count = int(relation(change[mask]).sum())
            total = int(mask.sum())
            row[column] = f"{count} ({100 * count / total:.1f}%)"
        trajectory.append(row)
    write_table(pd.DataFrame(trajectory), "Table_S2_VAS_trajectory.csv")
    trajectory_counts = pd.crosstab(
        data["PSQI_sleep_group"],
        pd.cut(
            data["Followup_pain_intensity_VAS"] - data["Baseline_pain_intensity_VAS"],
            bins=[-np.inf, -1e-12, 1e-12, np.inf],
            labels=["Decreased", "Unchanged", "Increased"],
        ),
    )
    chi_square, trajectory_p, degrees_freedom, _ = stats.chi2_contingency(trajectory_counts)
    write_table(
        pd.DataFrame([{
            "Comparison": "VAS trajectory by sleep group",
            "Chi-square": f"{chi_square:.2f}",
            "df": int(degrees_freedom),
            "p-value": format_p(trajectory_p),
        }]),
        "Table_S2_VAS_trajectory_test.csv",
    )
    return result


def table_s4() -> pd.DataFrame:
    data = analysis_cohort()
    focal = ["PSQI_global_score", "PSQI_sleep_duration_hours_std", "SolarRad_BL_to_FU_mean_MJm2"]
    models = [
        ("Primary integrated model", []),
        ("+ baseline-month sine/cosine", ["Season_BL_sine", "Season_BL_cosine"]),
        ("+ follow-up-month sine/cosine", ["Season_FU_sine", "Season_FU_cosine"]),
    ]
    rows = []
    variable = "SolarRad_BL_to_FU_mean_MJm2"
    for label, seasonal in models:
        fit = _vitamin_d_model(data, focal, seasonal_predictors=seasonal)
        rows.append({
            "Model": label,
            "n": fit["n"],
            "Solar-radiation B (95% CI)": f"{fit['coef'][variable]:.3f} ({fit['lower'][variable]:.3f} to {fit['upper'][variable]:.3f})",
            "β": f"{standardized_beta(fit, variable, 'Followup_25OHD_ng_mL'):.3f}",
            "p": format_p(fit["p"][variable]),
        })
    result = pd.DataFrame(rows)
    write_table(result, "Table_S4.csv")
    return result


def run_non_ml_tables() -> None:
    table_1()
    table_2()
    table_3()
    table_s1()
    table_s2()
    table_s4()


if __name__ == "__main__":
    run_non_ml_tables()
