# Figure 5. Cross-validated prediction across the two principal research themes

required_packages <- c(
  "readr", "dplyr", "tidyr", "ggplot2", "patchwork", "scales", "stringr"
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

invisible(lapply(required_packages, library, character.only = TRUE))

source(file.path(dirname(script_path <- sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)[1])), "common.R"))

result_dir <- supporting_output_dir
figure_dir <- figure_output_dir
repeat_path <- file.path(result_dir, "Table_4_repeated_cv_metrics.csv")
summary_path <- file.path(result_dir, "Table_4_summary_raw.csv")
prediction_path <- file.path(result_dir, "Figure_5_participant_predictions.csv")

required_inputs <- c(repeat_path, summary_path, prediction_path)
if (any(!file.exists(required_inputs))) {
  stop(
    "Run python/Table_4.py before Figure_5.R. Missing: ",
    paste(basename(required_inputs[!file.exists(required_inputs)]), collapse = ", ")
  )
}

cv <- readr::read_csv(repeat_path, show_col_types = FALSE)
summary_data <- readr::read_csv(summary_path, show_col_types = FALSE)
predictions <- readr::read_csv(prediction_path, show_col_types = FALSE)

# -----------------------------------------------------------------------------
# Color palette
# -----------------------------------------------------------------------------
colors <- c(
  psychological = "#B23A48",  # muted, saturated red
  pain = "#6746A5",           # clear purple
  vitamin_d = "#A66A1F",      # ochre/brown
  environment = "#16705A",    # deep environmental green
  sleep = "#173F5F",          # navy
  neutral = "#596673"
)

light_fills <- c(
  psychological = "#F8EAEC",
  pain = "#F0ECF8",
  vitamin_d = "#F5EEE4",
  environment = "#E7F2ED",
  sleep = "#E8EFF4"
)

theme_publication <- theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(
      size = 17, face = "bold", color = "#15202B",
      margin = margin(b = 6)
    ),
    plot.subtitle = element_text(
      size = 12.5, color = "#53606C", margin = margin(b = 9)
    ),
    axis.title = element_text(size = 13.5, face = "bold", color = "#26313C"),
    axis.text = element_text(size = 11.3, color = "#26313C"),
    strip.background = element_rect(
      fill = "#F7F8F9", color = "#7D8790", linewidth = 0.55
    ),
    strip.text = element_text(
      size = 12.3, face = "bold", color = "#26313C",
      margin = margin(4, 4, 4, 4)
    ),
    panel.grid.major.y = element_line(color = "#E3E7EA", linewidth = 0.45),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.margin = margin(10, 13, 10, 10)
  )

format_p <- function(x) {
  ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
}

# -----------------------------------------------------------------------------
# Select the three models shown in the figure
# -----------------------------------------------------------------------------
selected_models <- tibble::tribble(
  ~task, ~model, ~target_label, ~domain, ~short_label,
  "Psychological burden", "+ Sleep and pain",
  "Psychological burden", "psychological", "GSI model",
  "Baseline pain intensity", "+ Sleep",
  "Baseline pain", "pain", "VAS model",
  "Follow-up 25(OH)D", "Exposure-updated environment",
  "Follow-up 25(OH)D", "environment", "25(OH)D model"
)

selected_summary <- summary_data |>
  inner_join(selected_models, by = c("task", "model"))

selected_cv <- cv |>
  inner_join(selected_models, by = c("task", "model")) |>
  mutate(
    target_label = factor(
      target_label,
      levels = c("Psychological burden", "Baseline pain", "Follow-up 25(OH)D")
    ),
    # ggplot draws discrete y-axis levels from bottom to top. Reverse the
    # intended reading order here so both performance panels display the
    # outcomes from top to bottom as psychological burden, baseline pain,
    # and follow-up 25(OH)D.
    short_label = factor(
      short_label,
      levels = c("25(OH)D model", "VAS model", "GSI model")
    )
  )

# Participant-level out-of-fold predictions were averaged across 100 repetitions
# before this plotting script was run.
selected_predictions <- predictions |>
  inner_join(selected_models, by = c("task", "model")) |>
  mutate(
    target_label = factor(
      target_label,
      levels = c("Psychological burden", "Baseline pain", "Follow-up 25(OH)D")
    )
  )

# -----------------------------------------------------------------------------
# Model summary cards
# -----------------------------------------------------------------------------
model_card <- function(panel_letter, heading, model_text, performance_text,
                       border_color, fill_color, heading_color = border_color) {
  ggplot() +
    annotation_custom(
      grid::roundrectGrob(
        x = 0.5, y = 0.5, width = 0.985, height = 0.92,
        r = grid::unit(0.045, "snpc"),
        gp = grid::gpar(
          fill = fill_color, col = border_color,
          lwd = 1.25, alpha = 0.92
        )
      ),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    ) +
    annotate(
      "text", x = 0.055, y = 0.86,
      label = paste0(panel_letter, "  ", heading),
      hjust = 0, vjust = 1, size = 6.2,
      fontface = "bold", color = heading_color
    ) +
    annotate(
      "text", x = 0.055, y = 0.67,
      label = model_text,
      hjust = 0, vjust = 1, size = 4.25,
      lineheight = 1.12, color = "#111820"
    ) +
    annotate(
      "text", x = 0.055, y = 0.16,
      label = performance_text,
      hjust = 0, vjust = 0, size = 4.5,
      lineheight = 1.12, fontface = "bold", color = heading_color
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(3, 5, 3, 5))
}

gsi_result <- selected_summary |>
  filter(task == "Psychological burden") |>
  slice(1)

vas_result <- selected_summary |>
  filter(task == "Baseline pain intensity") |>
  slice(1)

vitd_result <- selected_summary |>
  filter(task == "Follow-up 25(OH)D") |>
  slice(1)

card_a <- model_card(
  "A", "Psychological-burden model",
  paste0(
    "Outcome: SCL-90-R GSI T-score\n",
    "Clinical characteristics + PSQI + sleep duration\n",
    "+ baseline VAS + clinical pain burden"
  ),
  paste0(
    "CV R²  ", sprintf("%.3f", gsi_result$CV_R2),
    "     RMSE  ", sprintf("%.3f", gsi_result$CV_RMSE), "\n",
    "ΔRMSE  ", sprintf("%.3f", gsi_result$Delta_RMSE),
    "     p=", format_p(gsi_result$corrected_p),
    "     q=", format_p(gsi_result$FDR_q)
  ),
  colors[["psychological"]], light_fills[["psychological"]]
)

card_b <- model_card(
  "B", "Baseline-pain model",
  paste0(
    "Outcome: baseline VAS\n",
    "Age + sex + chronic TMD\n",
    "+ PSQI + sleep duration"
  ),
  paste0(
    "CV R²  ", sprintf("%.3f", vas_result$CV_R2),
    "     RMSE  ", sprintf("%.3f", vas_result$CV_RMSE), "\n",
    "ΔRMSE  ", sprintf("%.3f", vas_result$Delta_RMSE),
    "     p=", format_p(vas_result$corrected_p),
    "     q=", format_p(vas_result$FDR_q)
  ),
  colors[["pain"]], light_fills[["pain"]]
)

card_c <- model_card(
  "C", "Follow-up vitamin D model",
  paste0(
    "Outcome: follow-up serum 25(OH)D\n",
    "Clinical–biochemical benchmark + sleep\n",
    "+ prebaseline and interval solar radiation"
  ),
  paste0(
    "CV R²  ", sprintf("%.3f", vitd_result$CV_R2),
    "     RMSE  ", sprintf("%.3f", vitd_result$CV_RMSE), "\n",
    "ΔRMSE  ", sprintf("%.3f", vitd_result$Delta_RMSE),
    "     p=", format_p(vitd_result$corrected_p),
    "     q=", format_p(vitd_result$FDR_q)
  ),
  colors[["environment"]], light_fills[["environment"]]
)

# -----------------------------------------------------------------------------
# Panel D: observed versus cross-validated predicted values
# -----------------------------------------------------------------------------
domain_colors <- c(
  psychological = colors[["psychological"]],
  pain = colors[["pain"]],
  environment = colors[["environment"]]
)

panel_d <- ggplot(
  selected_predictions,
  aes(x = predicted, y = observed, color = domain, fill = domain)
) +
  geom_abline(
    slope = 1, intercept = 0,
    color = "#77818B", linewidth = 0.75, linetype = "22"
  ) +
  geom_smooth(
    method = "lm", formula = y ~ x,
    se = TRUE, alpha = 0.13, linewidth = 1.0
  ) +
  geom_point(
    size = 2.15, alpha = 0.58, stroke = 0.25,
    shape = 21, color = "white"
  ) +
  facet_wrap(~target_label, scales = "free", nrow = 1) +
  scale_color_manual(values = domain_colors) +
  scale_fill_manual(values = domain_colors) +
  labs(
    title = "D  Observed and cross-validated predicted outcomes",
    subtitle = "Participant-level predictions averaged across 100 repeated five-fold partitions",
    x = "Cross-validated predicted value",
    y = "Observed value"
  ) +
  theme_publication +
  theme(
    aspect.ratio = 0.92,
    panel.spacing.x = grid::unit(0.8, "lines")
  )

# -----------------------------------------------------------------------------
# Panel E: transparent distributions of CV R² and percentage RMSE change
# -----------------------------------------------------------------------------
reference_names <- c(
  "Psychological burden" = "Clinical benchmark",
  "Baseline pain intensity" = "Clinical benchmark",
  "Follow-up 25(OH)D" = "Conditional benchmark"
)

reference_rmse <- cv |>
  filter(
    (task == "Psychological burden" & model == "Clinical benchmark") |
      (task == "Baseline pain intensity" & model == "Clinical benchmark") |
      (task == "Follow-up 25(OH)D" & model == "Conditional benchmark")
  ) |>
  select(task, repeat_index, reference_RMSE = RMSE)

increment_data <- selected_cv |>
  left_join(reference_rmse, by = c("task", "repeat_index")) |>
  mutate(
    delta_pct = 100 * (RMSE - reference_RMSE) / reference_RMSE,
    short_label = factor(
      short_label,
      levels = c("25(OH)D model", "VAS model", "GSI model")
    )
  )

annotation_data <- selected_summary |>
  mutate(
    short_label = factor(
      short_label,
      levels = c("25(OH)D model", "VAS model", "GSI model")
    ),
    label = paste0(
      "p=", format_p(corrected_p),
      "; q=", format_p(FDR_q)
    )
  )

panel_e1 <- ggplot(
  selected_cv,
  aes(x = R2, y = short_label, fill = domain, color = domain)
) +
  geom_vline(xintercept = 0, color = "#505A64", linewidth = 0.65, linetype = "22") +
  geom_violin(
    orientation = "y", width = 0.77, trim = FALSE,
    alpha = 0.25, linewidth = 0.9, adjust = 1.05
  ) +
  geom_jitter(
    width = 0, height = 0.065, size = 0.72,
    alpha = 0.16, stroke = 0
  ) +
  geom_boxplot(
    orientation = "y", width = 0.15, outlier.shape = NA,
    alpha = 0.72, color = "#202830", linewidth = 0.5
  ) +
  stat_summary(
    fun = mean, geom = "point", shape = 23,
    size = 3.0, fill = "white", color = "#111820", stroke = 0.8
  ) +
  scale_color_manual(values = domain_colors) +
  scale_fill_manual(values = domain_colors) +
  labs(x = "Cross-validated R²", y = NULL) +
  theme_publication +
  theme(
    plot.title = element_blank(), plot.subtitle = element_blank(),
    axis.text.y = element_text(face = "bold")
  )

panel_e2 <- ggplot(
  increment_data,
  aes(x = delta_pct, y = short_label, fill = domain, color = domain)
) +
  geom_vline(xintercept = 0, color = "#505A64", linewidth = 0.65, linetype = "22") +
  geom_violin(
    orientation = "y", width = 0.77, trim = FALSE,
    alpha = 0.25, linewidth = 0.9, adjust = 1.05
  ) +
  geom_jitter(
    width = 0, height = 0.065, size = 0.72,
    alpha = 0.16, stroke = 0
  ) +
  geom_boxplot(
    orientation = "y", width = 0.15, outlier.shape = NA,
    alpha = 0.72, color = "#202830", linewidth = 0.5
  ) +
  stat_summary(
    fun = mean, geom = "point", shape = 23,
    size = 3.0, fill = "white", color = "#111820", stroke = 0.8
  ) +
  geom_text(
    data = annotation_data,
    aes(x = Inf, y = short_label, label = label),
    inherit.aes = FALSE, hjust = 1.03, vjust = -1.15,
    size = 3.25, color = "#202830"
  ) +
  scale_color_manual(values = domain_colors) +
  scale_fill_manual(values = domain_colors) +
  scale_x_continuous(
    labels = label_number(accuracy = 0.1, suffix = "%"),
    expand = expansion(mult = c(0.08, 0.21))
  ) +
  labs(
    x = "Change in CV RMSE versus benchmark\n(lower values indicate improvement)",
    y = NULL
  ) +
  theme_publication +
  theme(
    plot.title = element_blank(), plot.subtitle = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank()
  ) +
  coord_cartesian(clip = "off")

panel_e_header <- ggplot() +
  annotate(
    "text", x = 0, y = 0.82,
    label = "E  Repeated cross-validated performance",
    hjust = 0, size = 6.0, fontface = "bold", color = "#15202B"
  ) +
  annotate(
    "text", x = 0, y = 0.22,
    label = "Transparent distributions across 100 repeated five-fold partitions",
    hjust = 0, size = 4.05, color = "#53606C"
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void()

panel_e <- panel_e_header / (panel_e1 | panel_e2) +
  plot_layout(heights = c(0.12, 0.88), widths = c(1, 1.18))

# -----------------------------------------------------------------------------
# Final assembly
# -----------------------------------------------------------------------------
top_row <- (card_a | card_b | card_c) +
  plot_layout(widths = c(1, 1, 1.08))
bottom_row <- (panel_d | panel_e) +
  plot_layout(widths = c(1.26, 1))

figure5 <- top_row / bottom_row +
  plot_layout(heights = c(0.39, 0.61)) +
  plot_annotation(
    theme = theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(8, 8, 8, 8)
    )
  )

print(figure5)

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
png_path <- file.path(figure_dir, "Figure_5.png")
pdf_path <- file.path(figure_dir, "Figure_5.pdf")
svg_path <- file.path(figure_dir, "Figure_5.svg")

ggsave(
  png_path, figure5,
  width = 18.0, height = 11.7, units = "in", dpi = 450, bg = "white"
)

# Base PDF device avoids cairo_pdf output-stream errors on Windows.
grDevices::pdf(
  pdf_path, width = 18.0, height = 11.7,
  family = "sans", useDingbats = FALSE, onefile = TRUE
)
print(figure5)
grDevices::dev.off()

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    svg_path, figure5,
    width = 18.0, height = 11.7, units = "in",
    device = svglite::svglite, bg = "white"
  )
}

message("Figure 5 saved to: ", figure_dir)
