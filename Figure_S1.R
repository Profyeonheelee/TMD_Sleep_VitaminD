# Supplementary Figure S1. Correlation structures by sleep-quality group

required_packages <- c("readxl", "dplyr", "tidyr", "purrr", "ggplot2", "patchwork")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install the following packages: ", paste(missing_packages, collapse = ", "))
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
supplementary_output_dir <- figure_output_dir

dat <- read_excel(data_path, sheet = "Analysis_Master") |>
  filter(Eligible_PSQI_analysis == 1) |>
  mutate(
    sleep_group = factor(
      PSQI_sleep_group,
      levels = c("Good sleeper", "Poor sleeper")
    )
  )

stopifnot(
  nrow(dat) == 120,
  sum(dat$sleep_group == "Good sleeper") == 62,
  sum(dat$sleep_group == "Poor sleeper") == 58
)

# All variables included in the global FDR family.
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

# Fifteen core variables displayed in the heatmaps.
display_variables <- c(
  "Baseline VAS",
  "Follow-up VAS",
  "VAS improvement",
  "Baseline pain burden",
  "Follow-up pain burden",
  "Pain-burden improvement",
  "Baseline 25(OH)D",
  "Follow-up 25(OH)D",
  "Change in 25(OH)D",
  "PSQI global score",
  "Sleep duration",
  "SCL-90-R GSI",
  "STOP-Bang",
  "Prebaseline solar radiation",
  "Interval solar radiation"
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

calculate_correlations <- function(data, panel_name) {
  variable_pairs <- combn(names(correlation_variables), 2, simplify = FALSE)

  results <- map_dfr(variable_pairs, function(pair) {
    x <- data[[correlation_variables[[pair[1]]]]]
    y <- data[[correlation_variables[[pair[2]]]]]
    keep <- complete.cases(x, y)

    test <- suppressWarnings(
      cor.test(x[keep], y[keep], method = "spearman", exact = FALSE)
    )

    tibble(
      panel = panel_name,
      variable_1 = pair[1],
      variable_2 = pair[2],
      pairwise_n = sum(keep),
      rho = unname(test$estimate),
      p_value = test$p.value,
      structural = pair_key(pair[1], pair[2]) %in% structural_pairs
    )
  })

  tested <- which(!results$structural)
  results$global_FDR_q <- NA_real_
  results$global_FDR_q[tested] <- p.adjust(
    results$p_value[tested], method = "BH"
  )

  results
}

cor_good <- calculate_correlations(
  filter(dat, sleep_group == "Good sleeper"), "Good sleeper"
)
cor_poor <- calculate_correlations(
  filter(dat, sleep_group == "Poor sleeper"), "Poor sleeper"
)

all_correlations <- bind_rows(cor_good, cor_poor)

prepare_heatmap_data <- function(results) {
  symmetric <- bind_rows(
    results,
    results |>
      transmute(
        panel,
        variable_1_swapped = variable_2,
        variable_2_swapped = variable_1,
        pairwise_n,
        rho,
        p_value,
        structural,
        global_FDR_q
      ) |>
      rename(
        variable_1 = variable_1_swapped,
        variable_2 = variable_2_swapped
      )
  )

  diagonal <- tibble(
    panel = unique(results$panel),
    variable_1 = names(correlation_variables),
    variable_2 = names(correlation_variables),
    pairwise_n = map_int(
      correlation_variables,
      ~sum(!is.na(dat[[.x]]))
    ),
    rho = 1,
    p_value = NA_real_,
    structural = FALSE,
    global_FDR_q = NA_real_
  )

  bind_rows(symmetric, diagonal) |>
    filter(
      variable_1 %in% display_variables,
      variable_2 %in% display_variables
    ) |>
    mutate(
      row_position = match(variable_1, display_variables),
      column_position = match(variable_2, display_variables),
      variable_1 = factor(variable_1, levels = rev(display_variables)),
      variable_2 = factor(variable_2, levels = display_variables),
      show_cell = row_position >= column_position,
      significant = !is.na(global_FDR_q) & global_FDR_q < 0.05,
      rho_display = if_else(show_cell, rho, NA_real_),
      significance_symbol = case_when(
        !show_cell | row_position == column_position | structural ~ "",
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ ""
      ),
      cell_label = case_when(
        !show_cell ~ "",
        row_position == column_position ~ "1.00",
        TRUE ~ paste0(sprintf("%.2f", rho), significance_symbol)
      ),
      label_color = if_else(abs(rho_display) >= 0.55, "white", "#252A2E")
    )
}

make_heatmap <- function(results, panel_title, low_color, high_color,
                         text_size = 3.0) {
  plot_data <- prepare_heatmap_data(results)

  ggplot(plot_data, aes(x = variable_2, y = variable_1)) +
    geom_tile(
      aes(fill = rho_display),
      color = "white",
      linewidth = 0.35,
      na.rm = FALSE
    ) +
    geom_tile(
      data = filter(plot_data, show_cell, significant),
      aes(x = variable_2, y = variable_1),
      fill = NA,
      color = "#20262C",
      linewidth = 0.30,
      inherit.aes = FALSE
    ) +
    geom_text(
      aes(label = cell_label, color = label_color),
      size = text_size + 0.25,
      fontface = "bold",
      na.rm = TRUE
    ) +
    scale_color_identity() +
    scale_fill_gradient2(
      low = low_color,
      mid = "#FAF8F4",
      high = high_color,
      midpoint = 0,
      limits = c(-1, 1),
      breaks = c(-1, -0.5, 0, 0.5, 1),
      name = "Spearman ρ",
      na.value = "white"
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    coord_fixed() +
    labs(title = panel_title, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(
        face = "bold", size = 15, hjust = 0, margin = margin(b = 9)
      ),
      axis.text.x = element_text(
        angle = 48, hjust = 1, vjust = 1, size = 10.0, color = "#202428"
      ),
      axis.text.y = element_text(size = 10.0, color = "#202428"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      plot.margin = margin(8, 12, 8, 8)
    )
}

panel_b <- make_heatmap(
  cor_good, "A  Good sleepers (n = 62)",
  low_color = "#3C5488", high_color = "#E64B35", text_size = 4.05
)
panel_c <- make_heatmap(
  cor_poor, "B  Poor sleepers (n = 58)",
  low_color = "#3C5488", high_color = "#E64B35", text_size = 4.05
)

supplementary_figure_s1 <- (panel_b | panel_c) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

save_heatmap <- function(plot_object, stem, destination = supplementary_output_dir) {
  ggsave(
    file.path(destination, paste0(stem, ".png")),
    plot_object, width = 23.0, height = 11.5, units = "in", dpi = 600,
    bg = "white"
  )
  ggsave(
    file.path(destination, paste0(stem, ".tiff")),
    plot_object, width = 23.0, height = 11.5, units = "in", dpi = 600,
    compression = "lzw", bg = "white"
  )
  ggsave(
    file.path(destination, paste0(stem, ".pdf")),
    plot_object, width = 23.0, height = 11.5, units = "in",
    device = grDevices::pdf, bg = "white"
  )
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggsave(
      file.path(destination, paste0(stem, ".svg")),
      plot_object, width = 23.0, height = 11.5, units = "in",
      device = svglite::svglite, bg = "white"
    )
  }
  message(stem, " was saved to: ", destination)
}

save_heatmap(
  supplementary_figure_s1,
  "Figure_S1",
  supplementary_output_dir
)

write.csv(
  all_correlations,
  file.path(supplementary_output_dir, "Supplementary_Figure_S1_correlation_statistics.csv"),
  row.names = FALSE
)

if (interactive()) {
  print(supplementary_figure_s1)
}
