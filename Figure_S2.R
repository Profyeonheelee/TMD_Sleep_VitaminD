# Supplementary Figure S2. Domain contributions and explanatory performance

pkgs <- c("readxl", "dplyr", "purrr", "tibble", "tidyr", "ggplot2", "patchwork")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Install: ", paste(missing_pkgs, collapse = ", "))
invisible(lapply(pkgs, library, character.only = TRUE))

source(file.path(dirname(script_path <- sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)[1])), "common.R"))
output_dir <- figure_output_dir

B <- 2000
set.seed(20260909)

d <- read_excel(data_path, sheet = "Analysis_Master") |>
  filter(Eligible_PSQI_analysis == 1) |>
  mutate(
    Male = as.integer(Sex == "Male"),
    month_sin = sin(2 * pi * VitD_draw_month_BL / 12),
    month_cos = cos(2 * pi * VitD_draw_month_BL / 12)
  )

outcome <- "Followup_pain_intensity_VAS"
common <- list(
  Environment = c("SolarRad_BL_to_FU_mean_MJm2", "month_sin", "month_cos"),
  `Vitamin D` = c("Baseline_25OHD_ng_mL", "Delta_25OHD_ng_mL"),
  Sleep = c("PSQI_global_score", "PSQI_sleep_duration_hours_std"),
  `Demographic/clinical` = c("Age_years", "Male", "Chronic_TMD_ge3m", "Clinical_followup_months"),
  `Baseline pain burden` = "Baseline_pain_location_count",
  `Baseline VAS` = "Baseline_pain_intensity_VAS"
)
full <- list(
  Environment = common$Environment,
  `Vitamin D` = common$`Vitamin D`,
  Sleep = common$Sleep,
  Psychological = "SCL_GSI_T_harmonized",
  `STOP-Bang` = "STOPBANG_total_score",
  `Demographic/clinical` = common$`Demographic/clinical`,
  `Baseline pain burden` = common$`Baseline pain burden`,
  `Baseline VAS` = common$`Baseline VAS`
)

common_vars <- unique(c(outcome, unlist(common, use.names = FALSE)))
full_vars <- unique(c(outcome, unlist(full, use.names = FALSE)))
missing_vars <- setdiff(full_vars, names(d))
if (length(missing_vars)) stop("Missing variable(s): ", paste(missing_vars, collapse = ", "))

d120 <- d |>
  select(all_of(full_vars)) |>
  filter(if_all(all_of(common_vars), ~ !is.na(.x))) |>
  mutate(full_case = if_all(everything(), ~ !is.na(.x)))
d81 <- d120 |> filter(full_case) |> select(-full_case)
d39 <- d120 |> filter(!full_case)

if (nrow(d120) != 120 || nrow(d81) != 81) {
  warning("Observed cohort sizes: ", nrow(d120), " and ", nrow(d81))
}

r2_fit <- function(data, predictors) {
  x <- cbind(1, as.matrix(data[, predictors, drop = FALSE]))
  y <- data[[outcome]]
  fit <- lm.fit(x, y)
  1 - sum(fit$residuals^2) / sum((y - mean(y))^2)
}

model_metrics <- function(data, blocks) {
  fit <- lm(reformulate(unlist(blocks, use.names = FALSE), outcome), data = data)
  tibble(
    n = nobs(fit), predictors = length(coef(fit)) - 1,
    R2 = summary(fit)$r.squared,
    adjusted_R2 = summary(fit)$adj.r.squared,
    RMSE = sqrt(mean(residuals(fit)^2))
  )
}

dominance <- function(data, blocks) {
  nm <- names(blocks); m <- length(nm)
  key <- function(x) {
    if (length(x) == 0) "__EMPTY_MODEL__" else paste(sort(x), collapse = "|")
  }
  cache <- setNames(0, "__EMPTY_MODEL__")
  for (k in seq_len(m)) {
    for (s in combn(nm, k, simplify = FALSE)) {
      cache[key(s)] <- r2_fit(data, unlist(blocks[s], use.names = FALSE))
    }
  }
  weight <- map_dbl(nm, function(target) {
    others <- setdiff(nm, target); value <- 0
    for (k in 0:(m - 1)) {
      sets <- if (k == 0) list(character()) else combn(others, k, simplify = FALSE)
      w <- factorial(k) * factorial(m - k - 1) / factorial(m)
      for (s in sets) value <- value + w * (cache[key(c(s, target))] - cache[key(s)])
    }
    value
  })
  full_r2 <- cache[key(nm)]
  tibble(block = nm, dominance_R2 = weight,
         percent = 100 * weight / full_r2, full_R2 = full_r2)
}

scenario_names <- c(
  N120 = "n=120: common domains",
  N81_common = "n=81: common domains",
  N81_full = "n=81: full domains"
)

points <- bind_rows(
  dominance(d120, common) |> mutate(scenario = "N120"),
  dominance(d81, common) |> mutate(scenario = "N81_common"),
  dominance(d81, full) |> mutate(scenario = "N81_full")
)

models <- bind_rows(
  model_metrics(d120, common) |> mutate(scenario = "N120"),
  model_metrics(d81, common) |> mutate(scenario = "N81_common"),
  model_metrics(d81, full) |> mutate(scenario = "N81_full")
) |> mutate(label = unname(scenario_names[scenario]))

fit81_common <- lm(reformulate(unlist(common, use.names = FALSE), outcome), d81)
fit81_full <- lm(reformulate(unlist(full, use.names = FALSE), outcome), d81)
nested <- anova(fit81_common, fit81_full)
delta_R2 <- summary(fit81_full)$r.squared - summary(fit81_common)$r.squared
nested_F <- nested$F[2]
nested_p <- nested$`Pr(>F)`[2]

# Stratified paired bootstrap: fixed n=81 complete and n=39 incomplete cases.
boot <- vector("list", B)
for (b in seq_len(B)) {
  b81 <- d81[sample.int(nrow(d81), replace = TRUE), , drop = FALSE]
  b39 <- d39[sample.int(nrow(d39), replace = TRUE), , drop = FALSE]
  b120 <- bind_rows(b81, b39)
  x1 <- dominance(b120, common) |> select(block, N120 = dominance_R2)
  x2 <- dominance(b81, common) |> select(block, N81_common = dominance_R2)
  x3 <- dominance(b81, full) |> select(block, N81_full = dominance_R2)
  boot[[b]] <- full_join(full_join(x1, x2, by = "block"), x3, by = "block") |>
    mutate(iteration = b)
  if (b %% 100 == 0) message("Bootstrap ", b, "/", B)
}

boot_wide <- bind_rows(boot)
boot_long <- boot_wide |>
  pivot_longer(starts_with("N"), names_to = "scenario", values_to = "dominance_R2") |>
  filter(!is.na(dominance_R2))

ci <- boot_long |>
  group_by(scenario, block) |>
  summarise(
    lower = quantile(dominance_R2, .025),
    upper = quantile(dominance_R2, .975),
    .groups = "drop"
  )
points <- points |> left_join(ci, by = c("scenario", "block"))

point_diff <- points |>
  filter(block %in% names(common)) |>
  select(block, scenario, dominance_R2) |>
  pivot_wider(names_from = scenario, values_from = dominance_R2) |>
  transmute(
    block,
    Cohort = N81_common - N120,
    `Added information` = N81_full - N81_common
  ) |>
  pivot_longer(-block, names_to = "contrast", values_to = "difference")

boot_diff <- boot_wide |>
  filter(block %in% names(common)) |>
  transmute(
    iteration, block,
    Cohort = N81_common - N120,
    `Added information` = N81_full - N81_common
  ) |>
  pivot_longer(c(Cohort, `Added information`), names_to = "contrast", values_to = "boot_difference")

differences <- boot_diff |>
  group_by(contrast, block) |>
  summarise(
    lower = quantile(boot_difference, .025),
    upper = quantile(boot_difference, .975),
    SE = sd(boot_difference),
    .groups = "drop"
  ) |>
  left_join(point_diff, by = c("contrast", "block")) |>
  mutate(p_value = 2 * pnorm(-abs(difference / SE))) |>
  group_by(contrast) |>
  mutate(FDR_q = p.adjust(p_value, "BH")) |>
  ungroup()

block_levels <- rev(c(
  "Demographic/clinical", "Sleep", "Baseline VAS", "Psychological",
  "Vitamin D", "Environment", "Baseline pain burden", "STOP-Bang"
))
scenario_colors <- c(
  "n=120: common domains" = "#626B73",
  "n=81: common domains" = "#173F5F",
  "n=81: full domains" = "#7A3E1D"
)

points <- points |>
  mutate(
    block = factor(block, block_levels),
    label = factor(unname(scenario_names[scenario]), levels = scenario_names),
    estimate_label = sprintf("%.3f (%.1f%%)", dominance_R2, percent)
  )
differences <- differences |>
  mutate(block = factor(block, block_levels))

A_label_x <- max(points$upper, na.rm = TRUE) * 1.33
A_axis_max <- max(points$upper, na.rm = TRUE) * 1.38

A <- ggplot(points, aes(dominance_R2, block, color = label, shape = label)) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0,
                orientation = "y",
                position = position_dodge(.58), linewidth = .72) +
  geom_point(position = position_dodge(.58), size = 3.8) +
  geom_text(
    aes(x = A_label_x, label = estimate_label, group = label),
    position = position_dodge(.58),
    hjust = 1,
    size = 3.65,
    color = "#171A1D",
    show.legend = FALSE
  ) +
  scale_color_manual(values = scenario_colors) +
  scale_shape_manual(values = c(16, 17, 15)) +
  scale_x_continuous(
    limits = c(0, A_axis_max),
    breaks = scales::pretty_breaks(n = 5),
    labels = scales::number_format(accuracy = 0.01),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(title = "A  Order-independent domain contributions",
       subtitle = "Dominance R² with 95% CIs from 2,000 stratified paired bootstrap resamples; percentages refer to total model R²",
       x = "Order-independent contribution to R²", y = NULL,
       color = NULL, shape = NULL) +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(face = "bold", size = 17),
        plot.subtitle = element_text(size = 12, color = "#59636C"),
        axis.title.x = element_text(face = "bold", size = 13),
        axis.text = element_text(size = 12.5, color = "#202428"),
        panel.grid.major.y = element_line(color = "#ECEFF1", linewidth = .45),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        legend.text = element_text(size = 11.5),
        plot.margin = margin(7, 10, 7, 7))

models <- models |>
  mutate(label = factor(label, levels = rev(scenario_names)))
C <- ggplot(models, aes(y = label)) +
  geom_segment(aes(x = adjusted_R2, xend = R2, yend = label),
               linewidth = 1.5, color = "#CDD2D6") +
  geom_point(aes(x = R2, color = "R²"), size = 4.2) +
  geom_point(aes(x = adjusted_R2, color = "Adjusted R²"), size = 4.2, shape = 17) +
  geom_text(aes(x = R2, label = sprintf("%.3f", R2)),
            nudge_x = .025, hjust = 0, size = 3.85) +
  geom_text(aes(x = adjusted_R2, label = sprintf("%.3f", adjusted_R2)),
            nudge_x = -.025, hjust = 1, size = 3.85) +
  scale_color_manual(values = c("R²" = "#7A3E1D", "Adjusted R²" = "#173F5F")) +
  scale_x_continuous(limits = c(0, max(models$R2) + .10)) +
  labs(title = "B  Overall explanatory performance",
       x = "Model explanatory performance", y = NULL, color = NULL) +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(face = "bold", size = 17),
        axis.title.x = element_text(face = "bold", size = 13),
        axis.text = element_text(size = 12.5, color = "#202428"),
        panel.grid.major.y = element_line(color = "#ECEFF1", linewidth = .45),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        legend.text = element_text(size = 11.5))

supplementary_figure_s2 <- A + C + plot_layout(widths = c(1.65, 1))

stem <- file.path(output_dir, "Figure_S2")
ggsave(paste0(stem, ".png"), supplementary_figure_s2, width = 14.5, height = 7.5,
       units = "in", dpi = 600, bg = "white")
ggsave(paste0(stem, ".tiff"), supplementary_figure_s2, width = 14.5, height = 7.5,
       units = "in", dpi = 600,
       compression = "lzw", bg = "white")
grDevices::pdf(paste0(stem, ".pdf"), width = 14.5, height = 7.5,
               useDingbats = FALSE, family = "Helvetica")
print(supplementary_figure_s2); grDevices::dev.off()

if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    paste0(stem, ".svg"), supplementary_figure_s2,
    width = 14.5, height = 7.5, units = "in",
    device = svglite::svglite, bg = "white"
  )
}

write.csv(points, paste0(stem, "_domain_contributions.csv"), row.names = FALSE)
write.csv(differences, paste0(stem, "_domain_differences.csv"), row.names = FALSE)
write.csv(models, paste0(stem, "_model_performance.csv"), row.names = FALSE)

cat("\nModel performance\n"); print(models, n = Inf)
cat(sprintf("\nAdded SCL-90-R GSI and STOP-Bang: ΔR²=%.3f, F=%.3f, p=%.3f\n",
            delta_R2, nested_F, nested_p))
cat("\nDomain contributions\n"); print(points, n = Inf)
cat("\nBootstrap differences, p-values, and FDR q-values\n"); print(differences, n = Inf)
message("Supplementary Figure S2 and results saved to: ", output_dir)
if (interactive()) print(supplementary_figure_s2)
