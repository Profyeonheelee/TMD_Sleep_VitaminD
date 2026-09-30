# Figure 4. Sleep, psychological, environmental, and vitamin D relationships

required_packages <- c(
  "readxl", "dplyr", "tidyr", "purrr", "tibble", "ggplot2",
  "patchwork", "sandwich", "lmtest"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(patchwork)
library(sandwich)
library(lmtest)

source(file.path(dirname(script_path <- sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)[1])), "common.R"))
output_dir <- figure_output_dir

dat <- read_excel(data_path, sheet = "Analysis_Master") |>
  filter(Eligible_PSQI_analysis == 1) |>
  mutate(Male = as.integer(Sex == "Male"))

stopifnot(nrow(dat) == 120)

format_p <- function(x) {
  ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
}

questionnaire_specification <- tribble(
  ~panel, ~outcome, ~panel_title, ~x_label, ~y_label, ~color,
  "A", "SCL_GSI_T_harmonized", "PSQI and psychological burden", "PSQI global score", "SCL-90-R GSI T-score", "#E64B35",
  "B", "STOPBANG_total_score", "PSQI and OSA risk", "PSQI global score", "STOP-Bang score", "#00A087",
  "C", "Baseline_pain_intensity_VAS", "PSQI and baseline pain", "PSQI global score", "Baseline VAS", "#6C4AA5"
)

# The FDR values for panels A-C are taken from the same complete correlation
# family used for Figure 2. This avoids recalculating FDR only among the three
# relationships selected for display.
correlation_variables <- c(
  "Baseline VAS" = "Baseline_pain_intensity_VAS",
  "Follow-up VAS" = "Followup_pain_intensity_VAS",
  "VAS improvement" = "Pain_intensity_reduction_VAS",
  "Baseline pain burden" = "Baseline_pain_location_count",
  "Follow-up pain burden" = "Followup_pain_location_count",
  "Pain-burden improvement" = "Pain_location_reduction",
  "Baseline 25(OH)D" = "Baseline_25OHD_ng_mL",
  "Follow-up 25(OH)D" = "Followup_25OHD_ng_mL",
  "Change in 25(OH)D" = "Delta_25OHD_ng_mL",
  "PSQI global score" = "PSQI_global_score",
  "Sleep duration" = "PSQI_sleep_duration_hours_std",
  "SCL-90-R GSI" = "SCL_GSI_T_harmonized",
  "STOP-Bang" = "STOPBANG_total_score",
  "Prebaseline sunshine" = "Sunshine_preBL_90d_mean_hr",
  "Prebaseline solar radiation" = "SolarRad_preBL_90d_mean_MJm2",
  "Interval sunshine" = "Sunshine_BL_to_FU_mean_hr",
  "Interval solar radiation" = "SolarRad_BL_to_FU_mean_MJm2",
  "Age" = "Age_years",
  "Symptom duration" = "Symptom_duration_months",
  "Follow-up interval" = "Clinical_followup_months"
)

pair_key <- function(a, b) {
  paste(sort(c(a, b)), collapse = " || ")
}

structural_pairs <- c(
  pair_key("Baseline VAS", "VAS improvement"),
  pair_key("Follow-up VAS", "VAS improvement"),
  pair_key("Baseline pain burden", "Pain-burden improvement"),
  pair_key("Follow-up pain burden", "Pain-burden improvement"),
  pair_key("Baseline 25(OH)D", "Change in 25(OH)D"),
  pair_key("Follow-up 25(OH)D", "Change in 25(OH)D"),
  pair_key("PSQI global score", "Sleep duration"),
  pair_key("Prebaseline sunshine", "Prebaseline solar radiation"),
  pair_key("Interval sunshine", "Interval solar radiation")
)

correlation_family <- map_dfr(
  combn(names(correlation_variables), 2, simplify = FALSE),
  function(pair) {
    x <- dat[[correlation_variables[[pair[1]]]]]
    y <- dat[[correlation_variables[[pair[2]]]]]
    keep <- complete.cases(x, y)
    result <- suppressWarnings(
      cor.test(x[keep], y[keep], method = "spearman", exact = FALSE)
    )

    tibble(
      variable_1 = pair[1],
      variable_2 = pair[2],
      pair = pair_key(pair[1], pair[2]),
      p_value = result$p.value,
      structural = pair_key(pair[1], pair[2]) %in% structural_pairs
    )
  }
)

tested_correlations <- which(!correlation_family$structural)
correlation_family$global_FDR_q <- NA_real_
correlation_family$global_FDR_q[tested_correlations] <- p.adjust(
  correlation_family$p_value[tested_correlations],
  method = "BH"
)

questionnaire_pair_names <- c(
  "SCL_GSI_T_harmonized" = pair_key("PSQI global score", "SCL-90-R GSI"),
  "STOPBANG_total_score" = pair_key("PSQI global score", "STOP-Bang"),
  "Baseline_pain_intensity_VAS" = pair_key("PSQI global score", "Baseline VAS")
)

run_spearman <- function(outcome) {
  analysis_data <- dat |>
    select(PSQI_global_score, all_of(outcome)) |>
    filter(if_all(everything(), ~ !is.na(.x)))

  test <- cor.test(
    analysis_data$PSQI_global_score,
    analysis_data[[outcome]],
    method = "spearman",
    exact = FALSE
  )

  tibble(
    n = nrow(analysis_data),
    rho = unname(test$estimate),
    p_value = test$p.value
  )
}

questionnaire_results <- questionnaire_specification |>
  mutate(result = map(outcome, run_spearman)) |>
  tidyr::unnest(result) |>
  mutate(
    pair = unname(questionnaire_pair_names[outcome]),
    FDR_q_value = map_dbl(
      pair,
      ~correlation_family$global_FDR_q[match(.x, correlation_family$pair)]
    )
  )

make_questionnaire_panel <- function(panel, outcome, panel_title, x_label,
                                     y_label, color, n, rho, p_value,
                                     FDR_q_value) {
  plot_data <- dat |>
    select(PSQI_global_score, value = all_of(outcome)) |>
    filter(if_all(everything(), ~ !is.na(.x)))

  annotation <- paste0(
    "n = ", n,
    "; Spearman ρ = ", sprintf("%.2f", rho),
    "\np = ", format_p(p_value),
    "; FDR q = ", format_p(FDR_q_value)
  )

  ggplot(plot_data, aes(x = PSQI_global_score, y = value)) +
    geom_jitter(
      width = 0.10,
      height = ifelse(
        outcome %in% c("STOPBANG_total_score", "Baseline_pain_intensity_VAS"),
        0.06,
        0
      ),
      size = 2.05,
      alpha = 0.64,
      color = color
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      linewidth = 0.95,
      color = color,
      fill = color,
      alpha = 0.16
    ) +
    annotate(
      "label",
      x = -Inf,
      y = Inf,
      label = annotation,
      hjust = -0.05,
      vjust = 1.15,
      size = 3.60,
      label.size = 0,
      fill = scales::alpha("white", 0.86),
      color = "#252A2E"
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 5),
      expand = expansion(mult = c(0.06, 0.16))
    ) +
    labs(
      title = paste(panel, panel_title),
      x = x_label,
      y = y_label
    ) +
    theme_classic(base_size = 12.5) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      axis.title = element_text(face = "bold", size = 11.5),
      axis.text = element_text(size = 10.5, color = "#202428"),
      axis.line = element_line(linewidth = 0.35, color = "#3C444B"),
      axis.ticks = element_line(linewidth = 0.35, color = "#3C444B"),
      plot.margin = margin(8, 10, 8, 8)
    )
}

questionnaire_panels <- pmap(
  questionnaire_results |>
    select(
      panel, outcome, panel_title, x_label, y_label, color,
      n, rho, p_value, FDR_q_value
    ),
  make_questionnaire_panel
)

vitaminD_specification <- tribble(
  ~panel, ~predictor, ~panel_title, ~x_label, ~color,
  "D", "PSQI_global_score", "Sleep quality and follow-up 25(OH)D", "Adjusted PSQI global score", "#E64B35",
  "E", "PSQI_sleep_duration_hours_std", "Sleep duration and follow-up 25(OH)D", "Adjusted sleep duration", "#3C5488",
  "F", "SolarRad_BL_to_FU_mean_MJm2", "Solar radiation and follow-up 25(OH)D", "Adjusted interval solar radiation", "#D98C2B"
)

vitaminD_covariates <- c(
  "Baseline_25OHD_ng_mL",
  "Vitamin_D_prescription",
  "Age_years",
  "Male",
  "Clinical_followup_months"
)

continuous_variables <- c(
  "Followup_25OHD_ng_mL", "Baseline_25OHD_ng_mL",
  "PSQI_global_score", "PSQI_sleep_duration_hours_std",
  "SolarRad_BL_to_FU_mean_MJm2", "Age_years",
  "Clinical_followup_months"
)

standardize <- function(x) {
  as.numeric(scale(x))
}

fit_vitaminD_model <- function(predictor) {
  model_variables <- unique(c(
    "Followup_25OHD_ng_mL",
    predictor,
    vitaminD_covariates
  ))

  analysis_data <- dat |>
    select(all_of(model_variables)) |>
    filter(if_all(everything(), ~ !is.na(.x)))

  variables_to_standardize <- intersect(
    model_variables,
    continuous_variables
  )

  standardized_data <- analysis_data |>
    mutate(across(all_of(variables_to_standardize), standardize))

  full_formula <- reformulate(
    c(predictor, vitaminD_covariates),
    response = "Followup_25OHD_ng_mL"
  )

  full_model <- lm(full_formula, data = standardized_data)

  robust_test <- lmtest::coeftest(
    full_model,
    vcov. = sandwich::vcovHC(full_model, type = "HC3")
  )

  beta <- robust_test[predictor, 1]
  robust_se <- robust_test[predictor, 2]
  degrees_freedom <- df.residual(full_model)
  critical_value <- qt(0.975, df = degrees_freedom)

  outcome_residual_model <- lm(
    reformulate(vitaminD_covariates, response = "Followup_25OHD_ng_mL"),
    data = standardized_data
  )

  predictor_residual_model <- lm(
    reformulate(vitaminD_covariates, response = predictor),
    data = standardized_data
  )

  plot_data <- tibble(
    x_residual = resid(predictor_residual_model),
    y_residual = resid(outcome_residual_model)
  )

  list(
    statistics = tibble(
      n = nrow(standardized_data),
      standardized_beta = beta,
      lower_95CI = beta - critical_value * robust_se,
      upper_95CI = beta + critical_value * robust_se,
      p_value = robust_test[predictor, 4]
    ),
    plot_data = plot_data
  )
}

vitaminD_fits <- map(vitaminD_specification$predictor, fit_vitaminD_model)

vitaminD_results <- vitaminD_specification |>
  mutate(
    n = map_int(vitaminD_fits, ~ .x$statistics$n),
    standardized_beta = map_dbl(
      vitaminD_fits,
      ~ .x$statistics$standardized_beta
    ),
    lower_95CI = map_dbl(vitaminD_fits, ~ .x$statistics$lower_95CI),
    upper_95CI = map_dbl(vitaminD_fits, ~ .x$statistics$upper_95CI),
    p_value = map_dbl(vitaminD_fits, ~ .x$statistics$p_value),
    FDR_q_value = p.adjust(p_value, method = "BH"),
    plot_data = map(vitaminD_fits, "plot_data")
  )

make_vitaminD_panel <- function(panel, predictor, panel_title, x_label,
                                color, n, standardized_beta, lower_95CI,
                                upper_95CI, p_value, FDR_q_value, plot_data) {
  annotation <- paste0(
    "n = ", n,
    "; β = ", sprintf("%.2f", standardized_beta),
    " (", sprintf("%.3f", lower_95CI),
    " to ", sprintf("%.3f", upper_95CI), ")",
    "\np = ", format_p(p_value),
    "; FDR q = ", format_p(FDR_q_value)
  )

  ggplot(plot_data, aes(x = x_residual, y = y_residual)) +
    geom_point(size = 2.05, alpha = 0.64, color = color) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      linewidth = 0.95,
      color = color,
      fill = color,
      alpha = 0.16
    ) +
    annotate(
      "label",
      x = -Inf,
      y = Inf,
      label = annotation,
      hjust = -0.05,
      vjust = 1.15,
      size = 3.55,
      label.size = 0,
      fill = scales::alpha("white", 0.86),
      color = "#252A2E"
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 5),
      expand = expansion(mult = c(0.06, 0.18))
    ) +
    labs(
      title = paste(panel, panel_title),
      x = x_label,
      y = "Adjusted follow-up 25(OH)D"
    ) +
    theme_classic(base_size = 12.5) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      axis.title = element_text(face = "bold", size = 11.3),
      axis.text = element_text(size = 10.5, color = "#202428"),
      axis.line = element_line(linewidth = 0.35, color = "#3C444B"),
      axis.ticks = element_line(linewidth = 0.35, color = "#3C444B"),
      plot.margin = margin(8, 10, 8, 8)
    )
}

vitaminD_panels <- pmap(
  vitaminD_results,
  make_vitaminD_panel
)

figure4 <- wrap_plots(
  c(questionnaire_panels, vitaminD_panels),
  ncol = 3
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure_4.png"
  ),
  plot = figure4,
  width = 13.5,
  height = 9.2,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure_4.tiff"
  ),
  plot = figure4,
  width = 13.5,
  height = 9.2,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

pdf_path <- file.path(
  output_dir,
  "Figure_4.pdf"
)

if (file.exists(pdf_path) && !file.remove(pdf_path)) {
  pdf_path <- file.path(
    output_dir,
      "Figure_4_new.pdf"
  )
  warning("The existing PDF appears to be open. Saving the new PDF as: ", pdf_path)
}

grDevices::pdf(
  file = pdf_path,
  width = 13.5,
  height = 9.2,
  useDingbats = FALSE,
  family = "Helvetica"
)
print(figure4)
grDevices::dev.off()

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    file.path(output_dir, "Figure_4.svg"),
    figure4,
    width = 13.5,
    height = 9.2,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}

write.csv(
  questionnaire_results |>
    select(
      panel, panel_title, n, rho, p_value, FDR_q_value
    ),
  file.path(output_dir, "Figure_4_questionnaire_statistics.csv"),
  row.names = FALSE
)

write.csv(
  vitaminD_results |>
    select(
      panel, panel_title, predictor, n, standardized_beta,
      lower_95CI, upper_95CI, p_value, FDR_q_value
    ),
  file.path(output_dir, "Figure_4_vitaminD_model_statistics.csv"),
  row.names = FALSE
)

cat("\nFigure 4: questionnaire relationships (panels A-C)\n")
for (i in seq_len(nrow(questionnaire_results))) {
  r <- questionnaire_results[i, ]
  cat(sprintf(
    "[%s] %s | n=%d | rho=%.3f | p=%s | FDR q=%s\n",
    r$panel, r$panel_title, r$n, r$rho,
    format_p(r$p_value), format_p(r$FDR_q_value)
  ))
}

cat("\nFigure 4: adjusted follow-up 25(OH)D models (panels D-F)\n")
cat(
  "Adjustment set: baseline 25(OH)D, vitamin D prescription, age, sex, ",
  "and follow-up interval.\n",
  sep = ""
)
for (i in seq_len(nrow(vitaminD_results))) {
  r <- vitaminD_results[i, ]
  cat(sprintf(
    paste0(
      "[%s] %s | n=%d\n",
      "    beta=%.3f | HC3 95%% CI %.3f to %.3f | p=%s | FDR q=%s\n"
    ),
    r$panel, r$panel_title, r$n,
    r$standardized_beta, r$lower_95CI, r$upper_95CI,
    format_p(r$p_value), format_p(r$FDR_q_value)
  ))
}

message("Figure 4 and its statistical results were saved to: ", output_dir)
message("PDF saved to: ", pdf_path)

# Display the completed figure in the RStudio Plots pane as well as saving it.
if (interactive()) {
  print(figure4)
}
