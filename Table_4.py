from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import RidgeCV
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import KFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from common import SUPPORT_DIR, analysis_cohort, bh_adjust, format_p, write_table


warnings.filterwarnings("ignore", category=RuntimeWarning)

N_REPEATS = 100
N_FOLDS = 5
SEED = 20260926
ALPHAS = np.logspace(-3, 3, 13)


def corrected_repeated_cv_test(differences, n_folds=N_FOLDS):
    values = np.asarray(differences, dtype=float)
    values = values[np.isfinite(values)]
    mean_difference = values.mean()
    variance = values.var(ddof=1)
    correction = 1 / len(values) + 1 / (n_folds - 1)
    standard_error = np.sqrt(correction * variance)
    if standard_error == 0:
        return mean_difference, np.nan, 1.0
    t_statistic = mean_difference / standard_error
    p_value = 2 * stats.t.sf(abs(t_statistic), df=len(values) - 1)
    return mean_difference, t_statistic, p_value


def main():
    data = analysis_cohort()

    clinical = ["Age_years", "Male", "Chronic_TMD_ge3m"]
    sleep = ["PSQI_global_score", "PSQI_sleep_duration_hours_std"]
    pain = ["Baseline_pain_intensity_VAS", "Baseline_pain_location_count"]
    vitamin_d_core = [
        "Baseline_25OHD_ng_mL",
        "Vitamin_D_prescription",
        "Age_years",
        "Male",
        "Clinical_followup_months",
    ]

    tasks = {
        "Psychological burden": {
            "outcome": "SCL_GSI_T_harmonized",
            "models": {
                "Clinical benchmark": clinical,
                "+ Sleep": clinical + sleep,
                "+ Pain": clinical + pain,
                "+ Sleep and pain": clinical + sleep + pain,
            },
            "reference": "Clinical benchmark",
        },
        "Baseline pain intensity": {
            "outcome": "Baseline_pain_intensity_VAS",
            "models": {
                "Clinical benchmark": clinical,
                "+ Sleep": clinical + sleep,
            },
            "reference": "Clinical benchmark",
        },
        "Follow-up 25(OH)D": {
            "outcome": "Followup_25OHD_ng_mL",
            "models": {
                "Conditional benchmark": vitamin_d_core,
                "+ Sleep": vitamin_d_core + sleep,
                "+ Prebaseline environment": vitamin_d_core + ["SolarRad_preBL_60d_mean_MJm2"],
                "+ Sleep and prebaseline environment": vitamin_d_core + sleep + ["SolarRad_preBL_60d_mean_MJm2"],
                "Exposure-updated environment": vitamin_d_core + sleep + [
                    "SolarRad_preBL_60d_mean_MJm2",
                    "SolarRad_BL_to_FU_mean_MJm2",
                ],
            },
            "reference": "Conditional benchmark",
        },
    }

    fold_rows = []
    repeated_rows = []
    prediction_rows = []

    for task_index, (task_name, specification) in enumerate(tasks.items()):
        outcome = specification["outcome"]
        models = specification["models"]
        required = list(dict.fromkeys([outcome] + [v for variables in models.values() for v in variables]))
        task_data = data.loc[:, required].dropna().reset_index(drop=True)
        task_data.insert(0, "analysis_id", np.arange(1, len(task_data) + 1))
        y = task_data[outcome].to_numpy(dtype=float)

        for repeat_index in range(N_REPEATS):
            outer = KFold(
                n_splits=N_FOLDS,
                shuffle=True,
                random_state=SEED + task_index * 10000 + repeat_index,
            )
            fold_assignments = list(outer.split(task_data))

            for model_name, variables in models.items():
                out_of_fold = np.full(len(task_data), np.nan)
                for fold_index, (train_index, test_index) in enumerate(fold_assignments, start=1):
                    pipeline = Pipeline([
                        ("scale", StandardScaler()),
                        ("ridge", RidgeCV(alphas=ALPHAS, scoring="neg_root_mean_squared_error")),
                    ])
                    pipeline.fit(task_data.loc[train_index, variables], y[train_index])
                    predicted = pipeline.predict(task_data.loc[test_index, variables])
                    out_of_fold[test_index] = predicted
                    fold_rows.append({
                        "task": task_name,
                        "model": model_name,
                        "repeat_index": repeat_index + 1,
                        "fold_index": fold_index,
                        "n_test": len(test_index),
                        "RMSE": mean_squared_error(y[test_index], predicted) ** 0.5,
                        "MAE": mean_absolute_error(y[test_index], predicted),
                        "best_alpha": pipeline.named_steps["ridge"].alpha_,
                    })

                repeated_rows.append({
                    "task": task_name,
                    "model": model_name,
                    "repeat_index": repeat_index + 1,
                    "n": len(task_data),
                    "predictors": len(variables),
                    "R2": r2_score(y, out_of_fold),
                    "RMSE": mean_squared_error(y, out_of_fold) ** 0.5,
                    "MAE": mean_absolute_error(y, out_of_fold),
                })
                for participant_index, prediction in enumerate(out_of_fold):
                    prediction_rows.append({
                        "task": task_name,
                        "model": model_name,
                        "repeat_index": repeat_index + 1,
                        "analysis_id": task_data.loc[participant_index, "analysis_id"],
                        "observed": y[participant_index],
                        "predicted": prediction,
                    })

            if (repeat_index + 1) % 10 == 0:
                print(f"{task_name}: repeat {repeat_index + 1}/{N_REPEATS}", flush=True)

    folds = pd.DataFrame(fold_rows)
    repeated = pd.DataFrame(repeated_rows)
    predictions = pd.DataFrame(prediction_rows)

    summary_rows = []
    comparison_rows = []
    for task_name, specification in tasks.items():
        reference = specification["reference"]
        task_repeated = repeated.loc[repeated["task"].eq(task_name)]
        task_folds = folds.loc[folds["task"].eq(task_name)]
        for model_name in specification["models"]:
            model_repeated = task_repeated.loc[task_repeated["model"].eq(model_name)]
            row = {
                "task": task_name,
                "model": model_name,
                "n": int(model_repeated["n"].iloc[0]),
                "predictors": int(model_repeated["predictors"].iloc[0]),
            }
            for metric in ["R2", "RMSE", "MAE"]:
                values = model_repeated[metric].to_numpy()
                row[f"CV_{metric}"] = values.mean()
                row[f"CV_{metric}_lower"] = np.quantile(values, 0.025)
                row[f"CV_{metric}_upper"] = np.quantile(values, 0.975)

            if model_name == reference:
                row.update({
                    "Delta_RMSE": np.nan,
                    "Delta_RMSE_lower": np.nan,
                    "Delta_RMSE_upper": np.nan,
                    "corrected_p": np.nan,
                })
            else:
                paired_repeats = model_repeated[["repeat_index", "RMSE"]].merge(
                    task_repeated.loc[task_repeated["model"].eq(reference), ["repeat_index", "RMSE"]],
                    on="repeat_index",
                    suffixes=("_model", "_reference"),
                )
                repeat_difference = paired_repeats["RMSE_model"] - paired_repeats["RMSE_reference"]
                paired_folds = task_folds.loc[
                    task_folds["model"].eq(model_name), ["repeat_index", "fold_index", "RMSE"]
                ].merge(
                    task_folds.loc[
                        task_folds["model"].eq(reference), ["repeat_index", "fold_index", "RMSE"]
                    ],
                    on=["repeat_index", "fold_index"],
                    suffixes=("_model", "_reference"),
                )
                fold_difference = paired_folds["RMSE_model"] - paired_folds["RMSE_reference"]
                mean_difference, t_statistic, p_value = corrected_repeated_cv_test(fold_difference)
                row.update({
                    "Delta_RMSE": repeat_difference.mean(),
                    "Delta_RMSE_lower": np.quantile(repeat_difference, 0.025),
                    "Delta_RMSE_upper": np.quantile(repeat_difference, 0.975),
                    "corrected_p": p_value,
                })
                comparison_rows.append({
                    "Outcome": task_name,
                    "Reference model": reference,
                    "Comparison model": model_name,
                    "Mean fold ΔRMSE": mean_difference,
                    "Corrected t": t_statistic,
                    "Corrected p": p_value,
                })
            summary_rows.append(row)

    summary = pd.DataFrame(summary_rows)
    comparisons = pd.DataFrame(comparison_rows)
    summary["FDR_q"] = np.nan
    comparisons["FDR q"] = np.nan
    for task_name in comparisons["Outcome"].unique():
        mask = comparisons["Outcome"].eq(task_name)
        q_values = bh_adjust(comparisons.loc[mask, "Corrected p"])
        comparisons.loc[mask, "FDR q"] = q_values
        for comparison, q_value in zip(comparisons.loc[mask, "Comparison model"], q_values):
            summary.loc[
                summary["task"].eq(task_name) & summary["model"].eq(comparison), "FDR_q"
            ] = q_value

    formatted_rows = []
    for _, row in summary.iterrows():
        formatted_rows.append({
            "Outcome": row["task"],
            "Model": row["model"],
            "n": int(row["n"]),
            "Predictors, n": int(row["predictors"]),
            "CV R² (2.5th–97.5th)": f"{row['CV_R2']:.3f} ({row['CV_R2_lower']:.3f} to {row['CV_R2_upper']:.3f})",
            "CV RMSE (2.5th–97.5th)": f"{row['CV_RMSE']:.3f} ({row['CV_RMSE_lower']:.3f} to {row['CV_RMSE_upper']:.3f})",
            "ΔRMSE (2.5th–97.5th)": "—" if pd.isna(row["Delta_RMSE"]) else f"{row['Delta_RMSE']:.3f} ({row['Delta_RMSE_lower']:.3f} to {row['Delta_RMSE_upper']:.3f})",
            "Corrected p": "—" if pd.isna(row["corrected_p"]) else format_p(row["corrected_p"]),
            "FDR q": "—" if pd.isna(row["FDR_q"]) else format_p(row["FDR_q"]),
        })
    write_table(pd.DataFrame(formatted_rows), "Table_4.csv")

    table_s3 = comparisons.copy()
    for column in ["Mean fold ΔRMSE", "Corrected t"]:
        table_s3[column] = table_s3[column].map(lambda x: f"{x:.3f}")
    table_s3["Corrected p"] = table_s3["Corrected p"].map(format_p)
    table_s3["FDR q"] = table_s3["FDR q"].map(format_p)
    write_table(table_s3, "Table_S3.csv")

    SUPPORT_DIR.mkdir(parents=True, exist_ok=True)
    folds.to_csv(SUPPORT_DIR / "Table_4_fold_level_metrics.csv", index=False)
    repeated.to_csv(SUPPORT_DIR / "Table_4_repeated_cv_metrics.csv", index=False)
    summary.to_csv(SUPPORT_DIR / "Table_4_summary_raw.csv", index=False)
    comparisons.to_csv(SUPPORT_DIR / "Table_S3_comparisons_raw.csv", index=False)
    participant_predictions = (
        predictions.groupby(["task", "model", "analysis_id", "observed"], as_index=False)["predicted"]
        .mean()
    )
    participant_predictions.to_csv(
        SUPPORT_DIR / "Figure_5_participant_predictions.csv", index=False
    )
    print(f"Saved Table 4 supporting results to {SUPPORT_DIR}")


if __name__ == "__main__":
    main()
