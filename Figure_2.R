# Figure 2. Overall correlation matrix and association network

required_packages <- c(
  "readxl", "dplyr", "tidyr", "purrr", "tibble", "ggplot2",
  "igraph", "ggraph", "ggforce", "patchwork"
)
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
library(tibble)
library(ggplot2)
library(igraph)
library(ggraph)
library(ggforce)
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

cor_overall <- calculate_correlations(dat, "Overall")

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

panel_a <- make_heatmap(
  cor_overall, "A  Overall correlation matrix (n = 120)",
  low_color = "#148F77", high_color = "#6C4AA5", text_size = 3.05
)

# Domain colors follow the visual language used in Figure 4.
domain_colors <- c(
  "Sleep" = "#3C5488",
  "Psychological" = "#E64B35",
  "Pain" = "#6C4AA5",
  "Vitamin D" = "#8A5A2B",
  "Environment" = "#D98C2B",
  "OSA risk" = "#00A087"
)

network_nodes <- tibble(
  name = display_variables,
  domain = case_when(
    name %in% c(
      "Baseline VAS", "Follow-up VAS", "VAS improvement",
      "Baseline pain burden", "Follow-up pain burden",
      "Pain-burden improvement"
    ) ~ "Pain",
    name %in% c(
      "Baseline 25(OH)D", "Follow-up 25(OH)D", "Change in 25(OH)D"
    ) ~ "Vitamin D",
    name %in% c("PSQI global score", "Sleep duration") ~ "Sleep",
    name == "SCL-90-R GSI" ~ "Psychological",
    name == "STOP-Bang" ~ "OSA risk",
    name %in% c(
      "Prebaseline solar radiation", "Interval solar radiation"
    ) ~ "Environment",
    TRUE ~ NA_character_
  ),
  domain = factor(domain, levels = names(domain_colors))
)

# Only FDR-significant, nonstructural correlations are drawn as edges.
network_edges <- cor_overall |>
  filter(
    variable_1 %in% display_variables,
    variable_2 %in% display_variables,
    !structural,
    !is.na(global_FDR_q),
    global_FDR_q < 0.05
  ) |>
  transmute(
    from = variable_1,
    to = variable_2,
    rho,
    p_value,
    FDR_q = global_FDR_q,
    direction = factor(
      if_else(rho >= 0, "Positive", "Negative"),
      levels = c("Positive", "Negative")
    ),
    edge_strength = abs(rho),
    rho_label = sprintf("%.2f", rho)
  )

association_graph <- graph_from_data_frame(
  network_edges,
  directed = FALSE,
  vertices = network_nodes
)

set.seed(20260910)
network_layout <- create_layout(
  association_graph,
  layout = "fr",
  weights = E(association_graph)$edge_strength
)

node_positions <- as_tibble(network_layout) |>
  select(name, x, y)
component_membership <- components(association_graph)$membership
node_positions <- node_positions |>
  mutate(
    component = unname(component_membership[name]),
    component_size = ave(component, component, FUN = length)
  )
network_edge_labels <- network_edges |>
  left_join(
    transmute(node_positions, from = name, x_from = x, y_from = y),
    by = "from"
  ) |>
  left_join(
    transmute(node_positions, to = name, x_to = x, y_to = y),
    by = "to"
  ) |>
  mutate(
    x_mid = (x_from + x_to) / 2,
    y_mid = (y_from + y_to) / 2
  )

panel_network <- ggraph(network_layout) +
  ggforce::geom_mark_hull(
    data = filter(node_positions, component_size >= 2),
    aes(x = x, y = y, group = component),
    inherit.aes = FALSE,
    expand = grid::unit(4.5, "mm"),
    radius = grid::unit(2.5, "mm"),
    concavity = 3,
    fill = scales::alpha("#8C969F", 0.055),
    color = scales::alpha("#6F7A83", 0.42),
    linewidth = 0.65,
    show.legend = FALSE
  ) +
  geom_edge_link(
    aes(width = edge_strength, color = direction),
    alpha = 0.72,
    lineend = "round",
    show.legend = TRUE
  ) +
  geom_label(
    data = network_edge_labels,
    aes(x = x_mid, y = y_mid, label = rho_label),
    inherit.aes = FALSE,
    size = 3.15,
    linewidth = 0.15,
    label.padding = grid::unit(0.10, "lines"),
    fill = grDevices::adjustcolor("white", alpha.f = 0.88),
    color = "#252A2E",
    show.legend = FALSE
  ) +
  geom_node_point(
    aes(fill = domain),
    shape = 21,
    size = 9.0,
    stroke = 0.65,
    color = "#30363B"
  ) +
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 3.8,
    lineheight = 0.95,
    color = "#202428"
  ) +
  scale_fill_manual(values = domain_colors, drop = FALSE, name = "Domain") +
  scale_edge_color_manual(
    values = c("Positive" = "#B24A4A", "Negative" = "#3F6F9F"),
    name = "Association"
  ) +
  scale_edge_width(range = c(0.65, 2.50), guide = "none") +
  labs(
    title = "B  FDR-significant association network",
    subtitle = "Edges represent nonstructural correlations with FDR q<0.05"
  ) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0),
    plot.subtitle = element_text(size = 10.5, color = "#59636B"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(face = "bold", size = 10.5),
    legend.text = element_text(size = 9.5),
    plot.margin = margin(12, 24, 12, 24)
  )

figure_2 <- panel_a | panel_network
figure_2 <- figure_2 + plot_layout(widths = c(1.22, 0.90))

ggsave(
  file.path(output_dir, "Figure_2.png"),
  figure_2, width = 21.0, height = 11.5, units = "in", dpi = 600,
  bg = "white"
)
ggsave(
  file.path(output_dir, "Figure_2.tiff"),
  figure_2, width = 21.0, height = 11.5, units = "in", dpi = 600,
  compression = "lzw", bg = "white"
)
ggsave(
  file.path(output_dir, "Figure_2.pdf"),
  figure_2, width = 21.0, height = 11.5, units = "in",
  device = grDevices::pdf, bg = "white"
)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(
    file.path(output_dir, "Figure_2.svg"),
    figure_2, width = 21.0, height = 11.5, units = "in",
    device = svglite::svglite, bg = "white"
  )
}

write.csv(
  cor_overall,
  file.path(output_dir, "Figure_2_overall_correlation_statistics.csv"),
  row.names = FALSE
)

write.csv(
  network_edges,
  file.path(output_dir, "Figure_2_network_edges.csv"),
  row.names = FALSE
)

message("Figure 2 and supporting results were saved to: ", output_dir)

if (interactive()) {
  print(figure_2)
}
