from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA = PROJECT_ROOT / "data" / "20260925_TMD_Sleep_VitaminD_Analysis_Master.xlsx"
DATA_PATH = Path(os.environ.get("TMD_SLEEP_DATA", DEFAULT_DATA)).expanduser()
TABLE_DIR = PROJECT_ROOT / "outputs" / "tables"
SUPPORT_DIR = PROJECT_ROOT / "outputs" / "supporting"


def ensure_output_dirs() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    SUPPORT_DIR.mkdir(parents=True, exist_ok=True)


def load_master() -> pd.DataFrame:
    if not DATA_PATH.exists():
        raise FileNotFoundError(
            f"Data file not found: {DATA_PATH}\n"
            "Copy the analysis workbook to data/ or set TMD_SLEEP_DATA."
        )
    data = pd.read_excel(DATA_PATH, sheet_name="Analysis_Master")
    required = {"Study_ID", "Eligible_PSQI_analysis", "PSQI_global_score"}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Analysis_Master is missing: {', '.join(sorted(missing))}")
    return data


def analysis_cohort() -> pd.DataFrame:
    data = load_master()
    data = data.loc[data["Eligible_PSQI_analysis"].eq(1)].copy()
    data["Male"] = data["Sex"].eq("Male").astype(int)
    data["Female"] = data["Sex"].eq("Female").astype(int)
    data["Pain_worsening"] = (
        data["Followup_pain_intensity_VAS"] > data["Baseline_pain_intensity_VAS"]
    ).astype(int)
    data["Season_BL_sine"] = np.sin(2 * np.pi * data["VitD_draw_month_BL"] / 12)
    data["Season_BL_cosine"] = np.cos(2 * np.pi * data["VitD_draw_month_BL"] / 12)
    data["Season_FU_sine"] = np.sin(2 * np.pi * data["VitD_draw_month_FU"] / 12)
    data["Season_FU_cosine"] = np.cos(2 * np.pi * data["VitD_draw_month_FU"] / 12)
    if len(data) != 120:
        raise ValueError(f"Expected 120 PSQI-analysis participants; found {len(data)}")
    return data


def bh_adjust(p_values) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    result = np.empty_like(adjusted)
    result[order] = adjusted
    return result


def complete_frame(data: pd.DataFrame, outcome: str, predictors: list[str]):
    columns = [outcome] + predictors
    frame = data.loc[:, columns].dropna().copy()
    y = frame[outcome].to_numpy(dtype=float)
    x = np.column_stack([np.ones(len(frame)), frame[predictors].to_numpy(dtype=float)])
    return frame, y, x


def ols_hc3(data: pd.DataFrame, outcome: str, predictors: list[str]) -> dict:
    frame, y, x = complete_frame(data, outcome, predictors)
    xtx_inv = np.linalg.pinv(x.T @ x)
    beta = xtx_inv @ x.T @ y
    fitted = x @ beta
    residual = y - fitted
    hat = np.sum((x @ xtx_inv) * x, axis=1)
    adjusted_residual = residual / np.clip(1 - hat, 1e-8, None)
    meat = x.T @ ((adjusted_residual ** 2)[:, None] * x)
    covariance = xtx_inv @ meat @ xtx_inv
    se = np.sqrt(np.diag(covariance))
    df = len(y) - x.shape[1]
    t_value = beta / se
    p_value = 2 * stats.t.sf(np.abs(t_value), df=df)
    critical = stats.t.ppf(0.975, df=df)
    lower = beta - critical * se
    upper = beta + critical * se
    sse = np.sum(residual ** 2)
    sst = np.sum((y - y.mean()) ** 2)
    r2 = 1 - sse / sst
    adjusted_r2 = 1 - (1 - r2) * (len(y) - 1) / df
    names = ["Intercept"] + predictors
    return {
        "n": len(y),
        "names": names,
        "coef": dict(zip(names, beta)),
        "se": dict(zip(names, se)),
        "lower": dict(zip(names, lower)),
        "upper": dict(zip(names, upper)),
        "p": dict(zip(names, p_value)),
        "r2": r2,
        "adjusted_r2": adjusted_r2,
        "frame": frame,
        "outcome_sd": float(np.std(y, ddof=1)),
    }


def glm_hc3(
    data: pd.DataFrame,
    outcome: str,
    predictors: list[str],
    family: str,
    tolerance: float = 1e-10,
    max_iter: int = 200,
) -> dict:
    frame, y, x = complete_frame(data, outcome, predictors)
    beta = np.zeros(x.shape[1])
    beta[0] = np.log(np.clip(y.mean(), 1e-6, None)) if family == "poisson" else 0.0

    for _ in range(max_iter):
        eta = x @ beta
        if family == "poisson":
            mu = np.exp(np.clip(eta, -30, 30))
            weight = np.clip(mu, 1e-8, None)
            working = eta + (y - mu) / weight
        elif family == "binomial":
            mu = 1 / (1 + np.exp(-np.clip(eta, -30, 30)))
            weight = np.clip(mu * (1 - mu), 1e-8, None)
            working = eta + (y - mu) / weight
        else:
            raise ValueError("family must be 'poisson' or 'binomial'")

        xtwx = x.T @ (weight[:, None] * x)
        updated = np.linalg.pinv(xtwx) @ x.T @ (weight * working)
        if np.max(np.abs(updated - beta)) < tolerance:
            beta = updated
            break
        beta = updated

    eta = x @ beta
    if family == "poisson":
        mu = np.exp(np.clip(eta, -30, 30))
        weight = np.clip(mu, 1e-8, None)
    else:
        mu = 1 / (1 + np.exp(-np.clip(eta, -30, 30)))
        weight = np.clip(mu * (1 - mu), 1e-8, None)

    bread = np.linalg.pinv(x.T @ (weight[:, None] * x))
    root_w_x = np.sqrt(weight)[:, None] * x
    leverage = np.sum((root_w_x @ bread) * root_w_x, axis=1)
    score_residual = (y - mu) / np.clip(1 - leverage, 1e-8, None)
    meat = x.T @ ((score_residual ** 2)[:, None] * x)
    covariance = bread @ meat @ bread
    se = np.sqrt(np.diag(covariance))
    z_value = beta / se
    p_value = 2 * stats.norm.sf(np.abs(z_value))
    lower = beta - stats.norm.ppf(0.975) * se
    upper = beta + stats.norm.ppf(0.975) * se
    names = ["Intercept"] + predictors
    return {
        "n": len(y),
        "names": names,
        "coef": dict(zip(names, beta)),
        "se": dict(zip(names, se)),
        "lower": dict(zip(names, lower)),
        "upper": dict(zip(names, upper)),
        "p": dict(zip(names, p_value)),
        "frame": frame,
    }


def standardized_beta(fit: dict, predictor: str, outcome: str) -> float:
    frame = fit["frame"]
    return (
        fit["coef"][predictor]
        * frame[predictor].std(ddof=1)
        / frame[outcome].std(ddof=1)
    )


def format_p(value: float) -> str:
    return "<0.001" if value < 0.001 else f"{value:.3f}"


def format_effect(prefix: str, estimate: float, lower: float, upper: float) -> str:
    return f"{prefix}={estimate:.3f} ({lower:.3f} to {upper:.3f})"


def write_table(data: pd.DataFrame, filename: str) -> Path:
    ensure_output_dirs()
    path = TABLE_DIR / filename
    data.to_csv(path, index=False, encoding="utf-8-sig")
    print(f"Saved {path}")
    return path
