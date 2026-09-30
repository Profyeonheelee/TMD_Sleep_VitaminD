# Figure 3. Clinical, psychological, and biochemical distributions by sleep group

required_packages <- c(
  "readxl", "dplyr", "tidyr", "purrr", "ggplot2",
  "patchwork"
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
library(ggplot2)
library(patchwork)

source(file.path(dirname(script_path <- sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)[1])), "common.R"))
output_dir <- figure_output_dir

dat <- read_excel(data_path, sheet = "Analysis_Master") |>
  filter(Eligible_PSQI_analysis == 1) |>
  mutate(
    sleep_group = factor(
      PSQI_sleep_group,
      levels = c("Good sleeper", "Poor sleeper")
    ),
    Female_sex = as.integer(Sex == "Female")
  )

stopifnot(nrow(dat) == 120)
stopifnot(sum(dat$sleep_group == "Good sleeper", na.rm = TRUE) == 62)
stopifnot(sum(dat$sleep_group == "Poor sleeper", na.rm = TRUE) == 58)

variable_specification <- tribble(
  ~panel, ~variable, ~panel_title, ~y_label, ~test,
  "A", "Baseline_pain_intensity_VAS",  "Baseline VASᵃ",                    "VAS score",             "welch",
  "B", "Baseline_pain_location_count", "Baseline clinical pain burdenᵇ",   "Painful sites, n",       "wilcoxon",
  "C", "SCL_GSI_T_harmonized",         "SCL-90-R GSI T-scoreᵃ",            "T-score",               "welch",
  "D", "Followup_pain_intensity_VAS",  "Follow-up VASᵃ",                   "VAS score",             "welch",
  "E", "Pain_intensity_reduction_VAS", "VAS improvementᵃ",                "Baseline − follow-up",  "welch",
  "F", "Delta_25OHD_ng_mL",            "Change in serum 25(OH)Dᵃ",         "Change, ng/mL",         "welch"
)

format_p <- function(x) {
  ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
}

run_group_test <- function(variable, test) {
  analysis_data <- dat |>
    select(sleep_group, value = all_of(variable)) |>
    filter(!is.na(sleep_group), !is.na(value))

  good <- analysis_data$value[analysis_data$sleep_group == "Good sleeper"]
  poor <- analysis_data$value[analysis_data$sleep_group == "Poor sleeper"]

  if (test == "welch") {
    result <- t.test(poor, good, var.equal = FALSE)
    estimate <- mean(poor) - mean(good)
    estimate_label <- sprintf("Mean difference = %.2f", estimate)
  } else {
    result <- wilcox.test(poor, good, exact = FALSE, conf.int = FALSE)
    estimate <- median(poor) - median(good)
    estimate_label <- sprintf("Median difference = %.2f", estimate)
  }

  tibble(
    n_good = length(good),
    n_poor = length(poor),
    estimate = estimate,
    estimate_label = estimate_label,
    p_value = result$p.value
  )
}

# Table 1 contains 20 nondefinitional comparisons. The same FDR family is
# used here so identical comparisons have identical q-values in the table
# and figure.
table1_fdr_specification <- tribble(
  ~variable, ~test,
  "Age_years",                         "welch",
  "Female_sex",                       "fisher",
  "Symptom_duration_months",           "welch",
  "Chronic_TMD_ge3m",                  "fisher",
  "Clinical_followup_months",          "welch",
  "SCL_GSI_T_harmonized",              "welch",
  "STOPBANG_total_score",              "welch",
  "STOPBANG_score_ge3",                "fisher",
  "Baseline_pain_intensity_VAS",       "welch",
  "Baseline_pain_location_count",      "wilcoxon",
  "Followup_pain_intensity_VAS",       "welch",
  "Pain_intensity_reduction_VAS",      "welch",
  "Baseline_25OHD_ng_mL",              "welch",
  "Followup_25OHD_ng_mL",              "welch",
  "Delta_25OHD_ng_mL",                 "welch",
  "Vitamin_D_prescription",            "fisher",
  "Sunshine_preBL_90d_mean_hr",         "welch",
  "SolarRad_preBL_90d_mean_MJm2",       "welch",
  "Sunshine_BL_to_FU_mean_hr",          "welch",
  "SolarRad_BL_to_FU_mean_MJm2",        "welch"
)

run_table1_p_value <- function(variable, test) {
  analysis_data <- dat |>
    select(sleep_group, value = all_of(variable)) |>
    filter(!is.na(sleep_group), !is.na(value))

  if (test == "fisher") {
    return(fisher.test(table(analysis_data$sleep_group, analysis_data$value))$p.value)
  }

  good <- analysis_data$value[analysis_data$sleep_group == "Good sleeper"]
  poor <- analysis_data$value[analysis_data$sleep_group == "Poor sleeper"]

  if (test == "wilcoxon") {
    return(wilcox.test(poor, good, exact = FALSE)$p.value)
  }

  t.test(poor, good, var.equal = FALSE)$p.value
}

table1_fdr_results <- table1_fdr_specification |>
  mutate(
    table1_p_value = map2_dbl(variable, test, run_table1_p_value),
    FDR_q_value = p.adjust(table1_p_value, method = "BH")
  )

test_results <- variable_specification |>
  mutate(result = map2(variable, test, run_group_test)) |>
  unnest(result) |>
  left_join(
    table1_fdr_results |> select(variable, FDR_q_value),
    by = "variable"
  )

make_distribution_panel <- function(panel, variable, panel_title, y_label,
                                    test, n_good, n_poor, estimate,
                                    estimate_label, p_value, FDR_q_value) {
  plot_data <- dat |>
    transmute(
      sleep_group,
      value = .data[[variable]]
    ) |>
    filter(!is.na(sleep_group), !is.na(value)) |>
    mutate(
      group_label = factor(
        sleep_group,
        levels = c("Good sleeper", "Poor sleeper"),
        labels = c(
          sprintf("Good sleeper\n(n = %d)", n_good),
          sprintf("Poor sleeper\n(n = %d)", n_poor)
        )
      )
    )

  annotation <- paste0(
    "p = ", format_p(p_value),
    "; FDR q = ", format_p(FDR_q_value)
  )

  ggplot(plot_data, aes(x = group_label, y = value, fill = sleep_group)) +
    geom_violin(
      trim = FALSE,
      scale = "width",
      width = 0.88,
      alpha = 0.18,
      color = NA
    ) +
    geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      alpha = 0.58,
      linewidth = 0.48,
      color = "#202428"
    ) +
    geom_jitter(
      aes(color = sleep_group),
      width = 0.18,
      height = 0,
      size = 1.80,
      alpha = 0.62,
      show.legend = FALSE
    ) +
    annotate(
      "label",
      x = 1.5,
      y = Inf,
      label = annotation,
      vjust = 1.35,
      size = 3.65,
      label.size = 0,
      fill = scales::alpha("white", 0.82),
      color = "#252A2E"
    ) +
    scale_fill_manual(values = c(
      "Good sleeper" = "#3C5488",
      "Poor sleeper" = "#E64B35"
    )) +
    scale_color_manual(values = c(
      "Good sleeper" = "#3C5488",
      "Poor sleeper" = "#E64B35"
    )) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    labs(
      title = paste(panel, panel_title),
      x = NULL,
      y = y_label
    ) +
    theme_classic(base_size = 12.5) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      axis.title.y = element_text(face = "bold", size = 11.5, margin = margin(r = 6)),
      axis.text.x = element_text(size = 11, color = "#202428"),
      axis.text.y = element_text(size = 10.5, color = "#202428"),
      axis.line = element_line(linewidth = 0.35, color = "#3C444B"),
      axis.ticks = element_line(linewidth = 0.35, color = "#3C444B"),
      plot.margin = margin(8, 10, 8, 8)
    )
}

panels <- pmap(test_results, make_distribution_panel)

figure3 <- wrap_plots(panels, ncol = 3)

ggsave(
  file.path(output_dir, "Figure_3.png"),
  figure3,
  width = 13.5,
  height = 8.8,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(output_dir, "Figure_3.tiff"),
  figure3,
  width = 13.5,
  height = 8.8,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

pdf_path <- file.path(output_dir, "Figure_3.pdf")

if (file.exists(pdf_path) && !file.remove(pdf_path)) {
  pdf_path <- file.path(output_dir, "Figure_3_new.pdf")
  warning("The existing PDF appears to be open. Saving the new PDF as: ", pdf_path)
}

grDevices::pdf(
  file = pdf_path,
  width = 13.5,
  height = 8.8,
  useDingbats = FALSE,
  family = "Helvetica"
)

print(figure3)
grDevices::dev.off()

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    file.path(output_dir, "Figure_3.svg"),
    figure3,
    width = 13.5,
    height = 8.8,
    units = "in",
    device = svglite::svglite,
    bg = "white"
  )
}

write.csv(
  test_results |>
    select(
      panel, variable, panel_title, test, n_good, n_poor,
      estimate, p_value, FDR_q_value
    ),
  file.path(output_dir, "Figure_3_sleep_group_distributions_statistics.csv"),
  row.names = FALSE
)

message("Figure 3 and its statistical results were saved to: ", output_dir)
cat("\nFigure 3 statistical results\n")
print(test_results, n = Inf)
cat("\nTable 1 FDR family used for Figure 3 annotations\n")
print(table1_fdr_results, n = Inf)

if (interactive()) print(figure3)
