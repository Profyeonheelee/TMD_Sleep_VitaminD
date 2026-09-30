# Figure 1. Cohort selection and longitudinal assessment framework
# -----------------------------------------------------------------------------
# This script uses only base R graphics packages (grid and grDevices).
# It draws the figure in RStudio and saves PDF/PNG copies. Files are first
# written to an ASCII-only temporary path and then copied to the destination,
# which avoids Cairo errors caused by long paths or non-ASCII folder names.

required_packages <- c("grid", "readxl")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install the following packages: ", paste(missing_packages, collapse = ", "))
}

library(grid)
source(file.path(dirname(script_path <- sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)[1])), "common.R"))
output_dir <- figure_output_dir

master <- readxl::read_excel(data_path, sheet = "Analysis_Master")
analysis_data <- master[master$Eligible_PSQI_analysis == 1, ]
n_recorded <- sum(master$PSQI_available == 1, na.rm = TRUE)
n_analysis <- nrow(analysis_data)
n_incomplete <- n_recorded - n_analysis
n_good <- sum(analysis_data$PSQI_sleep_group == "Good sleeper", na.rm = TRUE)
n_poor <- sum(analysis_data$PSQI_sleep_group == "Poor sleeper", na.rm = TRUE)
n_scl <- sum(!is.na(analysis_data$SCL_GSI_T_harmonized))
n_stopbang <- sum(!is.na(analysis_data$STOPBANG_total_score))

analysis_data$Male <- as.integer(analysis_data$Sex == "Male")
analysis_data$month_sin <- sin(2 * pi * analysis_data$VitD_draw_month_BL / 12)
analysis_data$month_cos <- cos(2 * pi * analysis_data$VitD_draw_month_BL / 12)
complete_case_variables <- c(
  "Followup_pain_intensity_VAS", "Baseline_pain_intensity_VAS",
  "Age_years", "Male", "Chronic_TMD_ge3m", "Clinical_followup_months",
  "Baseline_pain_location_count", "PSQI_global_score",
  "PSQI_sleep_duration_hours_std", "SCL_GSI_T_harmonized",
  "STOPBANG_total_score", "Baseline_25OHD_ng_mL", "Delta_25OHD_ng_mL",
  "SolarRad_BL_to_FU_mean_MJm2", "month_sin", "month_cos"
)
n_complete <- sum(complete.cases(analysis_data[, complete_case_variables]))

stopifnot(
  n_recorded == 136, n_analysis == 120, n_incomplete == 16,
  n_good == 62, n_poor == 58, n_scl == 92, n_stopbang == 101,
  n_complete == 81
)

# ---- Nature-inspired palette -------------------------------------------------
col_ink        <- "#17212B"
col_line       <- "#52616B"
col_light_line <- "#B8C2CC"
col_panel      <- "#F7F9FB"
col_neutral    <- "#EEF1F4"

col_sleep      <- "#244A73"  # navy
col_psych      <- "#B65353"  # muted red
col_pain       <- "#6A4C93"  # purple
col_vitd       <- "#9A6726"  # ochre/brown
col_env        <- "#3E7652"  # green
col_osa        <- "#6C857B"  # grey-green

# High-chroma colours used only for the circular domain markers.
dot_sleep      <- "#0067B9"
dot_psych      <- "#E64B35"
dot_pain       <- "#7A3EB1"
dot_vitd       <- "#E69F00"
dot_env        <- "#009E73"
dot_osa        <- "#00A087"

fill_sleep     <- "#E7EEF5"
fill_poor      <- "#F6E7E6"
fill_vitd      <- "#F5EDE1"
fill_env       <- "#E8F1EB"
fill_theme1    <- "#EEEAF5"
fill_theme2    <- "#EEF3EA"

# ---- Drawing helpers ---------------------------------------------------------
box <- function(x, y, w, h, label = NULL,
                fill = "white", border = col_line, lwd = 1.2,
                fontsize = 11.5, fontface = "plain", text_col = col_ink,
                radius = 0.014, lineheight = 1.08) {
  grid.roundrect(
    x = unit(x, "npc"), y = unit(y, "npc"),
    width = unit(w, "npc"), height = unit(h, "npc"),
    r = unit(radius, "snpc"),
    gp = gpar(fill = fill, col = border, lwd = lwd)
  )
  if (!is.null(label)) {
    grid.text(
      label, x = unit(x, "npc"), y = unit(y, "npc"),
      gp = gpar(col = text_col, fontsize = fontsize,
                fontface = fontface, lineheight = lineheight)
    )
  }
}

connector <- function(x1, y1, x2, y2, dashed = FALSE,
                      arrow_end = TRUE, col = col_line, lwd = 1.65) {
  line_gp <- gpar(
    col = col, lwd = lwd, lty = if (dashed) 2 else 1,
    lineend = "butt", linejoin = "mitre"
  )

  if (!arrow_end) {
    grid.lines(
      x = unit(c(x1, x2), "npc"), y = unit(c(y1, y2), "npc"),
      gp = line_gp
    )
    return(invisible(NULL))
  }

  # Custom device-independent arrowheads. The shaft stops at the base of the
  # triangle and only the triangle tip touches the target box boundary.
  # This prevents the small gaps/overlaps produced by grid::arrow().
  head_length <- unit(3.0, "mm")
  half_width  <- unit(1.75, "mm")

  if (abs(x2 - x1) < 1e-9 && y2 < y1) {             # downward arrow
    tip_y  <- unit(y2, "npc")
    base_y <- tip_y + head_length
    grid.lines(
      x = unit.c(unit(x1, "npc"), unit(x1, "npc")),
      y = unit.c(unit(y1, "npc"), base_y),
      gp = line_gp
    )
    grid.polygon(
      x = unit.c(unit(x2, "npc") - half_width,
                 unit(x2, "npc") + half_width,
                 unit(x2, "npc")),
      y = unit.c(base_y, base_y, tip_y),
      gp = gpar(fill = col, col = col, lwd = 0.4, linejoin = "mitre")
    )
  } else if (abs(y2 - y1) < 1e-9 && x2 > x1) {      # rightward arrow
    tip_x  <- unit(x2, "npc")
    base_x <- tip_x - head_length
    grid.lines(
      x = unit.c(unit(x1, "npc"), base_x),
      y = unit.c(unit(y1, "npc"), unit(y1, "npc")),
      gp = line_gp
    )
    grid.polygon(
      x = unit.c(base_x, base_x, tip_x),
      y = unit.c(unit(y2, "npc") - half_width,
                 unit(y2, "npc") + half_width,
                 unit(y2, "npc")),
      gp = gpar(fill = col, col = col, lwd = 0.4, linejoin = "mitre")
    )
  } else {
    stop("connector() currently supports downward or rightward arrows only.")
  }
}

down_arrowhead <- function(x, y_tip, col = col_line) {
  tip_y  <- unit(y_tip, "npc")
  base_y <- tip_y + unit(3.0, "mm")
  grid.polygon(
    x = unit.c(unit(x, "npc") - unit(1.75, "mm"),
               unit(x, "npc") + unit(1.75, "mm"),
               unit(x, "npc")),
    y = unit.c(base_y, base_y, tip_y),
    gp = gpar(fill = col, col = col, lwd = 0.4, linejoin = "mitre")
  )
}

right_arrowhead <- function(x_tip, y, col = col_line) {
  tip_x  <- unit(x_tip, "npc")
  base_x <- tip_x - unit(3.0, "mm")
  grid.polygon(
    x = unit.c(base_x, base_x, tip_x),
    y = unit.c(unit(y, "npc") - unit(1.75, "mm"),
               unit(y, "npc") + unit(1.75, "mm"),
               unit(y, "npc")),
    gp = gpar(fill = col, col = col, lwd = 0.4, linejoin = "mitre")
  )
}

time_axis_arrow <- function(x1, x2, y,
                            fill_col = col_sleep,
                            fill_alpha = 0.20,
                            border_alpha = 0.60) {
  y_mid      <- unit(y, "npc")
  x_start    <- unit(x1, "npc")
  x_tip      <- unit(x2, "npc")
  head_base  <- x_tip - unit(8.0, "mm")
  shaft_half <- unit(1.7, "mm")
  head_half  <- unit(4.3, "mm")

  grid.polygon(
    x = unit.c(x_start, head_base, head_base, x_tip,
               head_base, head_base, x_start),
    y = unit.c(y_mid - shaft_half, y_mid - shaft_half,
               y_mid - head_half, y_mid,
               y_mid + head_half, y_mid + shaft_half,
               y_mid + shaft_half),
    gp = gpar(
      fill = grDevices::adjustcolor(fill_col, alpha.f = fill_alpha),
      col  = grDevices::adjustcolor(fill_col, alpha.f = border_alpha),
      lwd = 1.0, linejoin = "mitre"
    )
  )
}

panel_title <- function(letter, title, y) {
  grid.text(
    paste0(letter, "  ", title),
    x = unit(0.035, "npc"), y = unit(y, "npc"), just = "left",
    gp = gpar(col = col_ink, fontsize = 16.5, fontface = "bold")
  )
}

domain_line <- function(x, y, domain, detail, colour, dot_colour,
                        fontsize = 11.0) {
  grid.circle(
    x = unit(x - 0.108, "npc"), y = unit(y, "npc"),
    r = unit(1.85, "mm"),
    gp = gpar(fill = dot_colour, col = "white", lwd = 0.9)
  )
  grid.text(
    domain, x = unit(x - 0.096, "npc"), y = unit(y, "npc"), just = "left",
    gp = gpar(col = col_ink, fontsize = fontsize, fontface = "bold")
  )
  grid.text(
    detail, x = unit(x - 0.010, "npc"), y = unit(y, "npc"), just = "left",
    gp = gpar(col = col_ink, fontsize = fontsize)
  )
}

draw_figure_1 <- function() {
  grid.newpage()
  pushViewport(viewport(gp = gpar(fontfamily = "sans")))

  # ---------------------------------------------------------------------------
  # A. Cohort selection and data availability
  # ---------------------------------------------------------------------------
  panel_title("A", "Cohort selection and data availability", 0.968)

  # CONNECTION LAYER. Drawn first as continuous paths; boxes are drawn later
  # and therefore mask the inner half of each boundary stroke cleanly.
  connector(0.48, 0.8505, 0.48, 0.806)
  connector(0.63, 0.765, 0.725, 0.765, dashed = TRUE,
            col = "#6F7A83", lwd = 1.45)

  # One continuous polyline for the complete upper split. Retracing the short
  # central segment avoids any raster gap at the T-junction.
  child_base_y <- unit(0.665, "npc") + unit(3.0, "mm")
  grid.lines(
    x = unit(c(0.335, 0.335, 0.48, 0.48, 0.48, 0.625, 0.625), "npc"),
    y = unit.c(child_base_y, unit(0.695, "npc"), unit(0.695, "npc"),
               unit(0.724, "npc"), unit(0.695, "npc"),
               unit(0.695, "npc"), child_base_y),
    gp = gpar(col = col_line, lwd = 1.55, lineend = "butt", linejoin = "mitre")
  )
  down_arrowhead(0.335, 0.665)
  down_arrowhead(0.625, 0.665)

  # One continuous polyline for the lower merge.
  availability_base_y <- unit(0.520, "npc") + unit(3.0, "mm")
  grid.lines(
    x = unit(c(0.335, 0.335, 0.48, 0.48, 0.48, 0.625, 0.625), "npc"),
    y = unit.c(unit(0.575, "npc"), unit(0.550, "npc"),
               unit(0.550, "npc"), availability_base_y,
               unit(0.550, "npc"), unit(0.550, "npc"),
               unit(0.575, "npc")),
    gp = gpar(col = col_line, lwd = 1.55, lineend = "butt", linejoin = "mitre")
  )
  down_arrowhead(0.48, 0.520)

  box(
    0.48, 0.888, 0.30, 0.075,
    sprintf("PSQI data recorded\n(n = %d)", n_recorded),
    fill = col_panel, border = col_sleep,
    fontsize = 12.2, fontface = "bold"
  )
  box(
    0.48, 0.765, 0.30, 0.082,
    sprintf("PSQI global score calculable\nPrimary analysis cohort (n = %d)", n_analysis),
    fill = "white", border = col_sleep,
    fontsize = 12.0, fontface = "bold"
  )

  box(
    0.825, 0.765, 0.20, 0.082,
    sprintf("Insufficient data to calculate\nPSQI global score (n = %d)", n_incomplete),
    fill = col_neutral, border = "#8B9298",
    fontsize = 11.0, text_col = col_ink
  )

  box(
    0.335, 0.620, 0.245, 0.090,
    sprintf("Good sleepers\nPSQI <=5  (n = %d)", n_good),
    fill = fill_sleep, border = col_sleep,
    fontsize = 12.0, fontface = "bold", text_col = col_ink
  )
  box(
    0.625, 0.620, 0.245, 0.090,
    sprintf("Poor sleepers\nPSQI >5  (n = %d)", n_poor),
    fill = fill_poor, border = col_psych,
    fontsize = 12.0, fontface = "bold", text_col = col_ink
  )

  box(0.48, 0.475, 0.64, 0.090,
      fill = "#F4F7F8", border = col_osa, fontsize = 10.0)
  grid.lines(
    x = unit(c(0.505, 0.505), "npc"), y = unit(c(0.442, 0.508), "npc"),
    gp = gpar(col = "#C8D0D6", lwd = 0.9)
  )
  grid.text(
    sprintf(
      "Core domains available in all %d participants\nSleep | pain | serum 25(OH)D | solar radiation",
      n_analysis
    ),
    x = unit(0.335, "npc"), y = unit(0.475, "npc"),
    gp = gpar(col = col_ink, fontsize = 11.2, fontface = "bold", lineheight = 1.08)
  )
  grid.text(
    sprintf(
      "Additional data availability\nSCL-90-R GSI: n = %d | STOP-Bang: n = %d\nFully integrated complete-case cohort: n = %d",
      n_scl, n_stopbang, n_complete
    ),
    x = unit(0.675, "npc"), y = unit(0.475, "npc"),
    gp = gpar(col = col_ink, fontsize = 10.5, lineheight = 1.08)
  )

  # ARROWHEAD LAYER. Redrawn last so every head remains fully visible above
  # rounded-box fills and borders, while the shafts stay behind the boxes.
  down_arrowhead(0.48, 0.806)
  right_arrowhead(0.725, 0.765, col = "#6F7A83")
  down_arrowhead(0.335, 0.665)
  down_arrowhead(0.625, 0.665)
  down_arrowhead(0.48, 0.520)

  # Fine separator between panels.
  grid.lines(
    x = unit(c(0.035, 0.965), "npc"), y = unit(c(0.422, 0.422), "npc"),
    gp = gpar(col = "#D6DDE3", lwd = 0.9)
  )

  # ---------------------------------------------------------------------------
  # B. Longitudinal assessment framework
  # ---------------------------------------------------------------------------
  panel_title("B", "Longitudinal assessment framework", 0.392)

  # Main time line and markers.
  time_axis_arrow(0.105, 0.910, 0.265)
  grid.lines(x = unit(c(0.255, 0.255), "npc"), y = unit(c(0.248, 0.282), "npc"),
             gp = gpar(col = col_ink, lwd = 1.5))
  grid.lines(x = unit(c(0.775, 0.775), "npc"), y = unit(c(0.248, 0.282), "npc"),
             gp = gpar(col = col_ink, lwd = 1.5))

  grid.text("Baseline assessment", x = unit(0.255, "npc"), y = unit(0.305, "npc"),
            gp = gpar(col = col_ink, fontsize = 12.0, fontface = "bold"))
  grid.text("same-day questionnaires, clinical examination, and blood sampling",
            x = unit(0.255, "npc"), y = unit(0.286, "npc"),
            gp = gpar(col = col_ink, fontsize = 9.8))
  grid.text("Follow-up assessment", x = unit(0.775, "npc"), y = unit(0.305, "npc"),
            gp = gpar(col = col_ink, fontsize = 12.0, fontface = "bold"))
  grid.text("clinical examination and blood sampling",
            x = unit(0.775, "npc"), y = unit(0.286, "npc"),
            gp = gpar(col = col_ink, fontsize = 9.8))

  # Environmental exposure ribbon.
  box(
    0.515, 0.350, 0.52, 0.036,
    "Date-linked interval solar radiation",
    fill = fill_env, border = col_env, lwd = 1.0,
    fontsize = 11.0, fontface = "bold", text_col = col_ink,
    radius = 0.010
  )
  # The ribbon begins and ends directly above the baseline and follow-up
  # markers, respectively; no crossing guide lines are needed.

  # Baseline and follow-up domain cards.
  connector(0.255, 0.248, 0.255, 0.2425,
            arrow_end = TRUE, col = col_line, lwd = 1.45)
  box(0.255, 0.165, 0.36, 0.155, fill = "white", border = "#9CA8B3", lwd = 1.0)
  domain_line(0.255, 0.215, "Sleep", "PSQI global score and sleep duration", col_sleep, dot_sleep)
  domain_line(0.255, 0.188, "Psychological", "SCL-90-R GSI", col_psych, dot_psych)
  domain_line(0.255, 0.161, "OSA risk", "STOP-Bang", col_osa, dot_osa)
  domain_line(0.255, 0.134, "Pain", "VAS and clinical pain burden", col_pain, dot_pain)
  domain_line(0.255, 0.107, "Vitamin D", "serum 25(OH)D", col_vitd, dot_vitd)

  connector(0.775, 0.248, 0.775, 0.225,
            arrow_end = TRUE, col = col_line, lwd = 1.45)
  box(0.775, 0.165, 0.28, 0.120, fill = "white", border = "#9CA8B3", lwd = 1.0)
  domain_line(0.775, 0.190, "Pain", "VAS and clinical pain burden", col_pain, dot_pain, 11.0)
  domain_line(0.775, 0.148, "Vitamin D", "serum 25(OH)D", col_vitd, dot_vitd, 11.0)

  # Keep temporal assessment arrowheads above the domain-card borders.
  down_arrowhead(0.255, 0.2425)
  down_arrowhead(0.775, 0.225)

  # Analytic themes: deliberately phrased as associations, not causal pathways.
  box(
    0.32125, 0.045, 0.4925, 0.058,
    "Theme 1 | Baseline sleep quality, psychological burden, and pain",
    fill = fill_theme1, border = col_pain, lwd = 1.0,
    fontsize = 11.2, fontface = "bold", text_col = col_ink
  )
  box(
    0.75375, 0.045, 0.3225, 0.058,
    "Theme 2 | Sleep/environment and follow-up 25(OH)D",
    fill = fill_theme2, border = col_env, lwd = 1.0,
    fontsize = 11.0, fontface = "bold", text_col = col_ink
  )

  popViewport()
}

# ---- Show in the RStudio Plots pane -----------------------------------------
draw_figure_1()

# ---- Safe export -------------------------------------------------------------
stem <- "Figure_1"

tmp_pdf <- file.path(tempdir(), paste0(stem, ".pdf"))
grDevices::pdf(tmp_pdf, width = 16, height = 10.2, useDingbats = FALSE)
draw_figure_1()
grDevices::dev.off()
pdf_path <- file.path(output_dir, paste0(stem, ".pdf"))
if (!file.copy(tmp_pdf, pdf_path, overwrite = TRUE)) {
  warning("PDF could not be copied to the requested output folder: ", pdf_path)
}

tmp_png <- file.path(tempdir(), paste0(stem, ".png"))
grDevices::png(tmp_png, width = 4800, height = 3060, res = 300, bg = "white")
draw_figure_1()
grDevices::dev.off()
png_path <- file.path(output_dir, paste0(stem, ".png"))
if (!file.copy(tmp_png, png_path, overwrite = TRUE)) {
  warning("PNG could not be copied to the requested output folder: ", png_path)
}

if (requireNamespace("svglite", quietly = TRUE)) {
  tmp_svg <- file.path(tempdir(), paste0(stem, ".svg"))
  svglite::svglite(tmp_svg, width = 16, height = 10.2, bg = "white")
  draw_figure_1()
  grDevices::dev.off()
  svg_path <- file.path(output_dir, paste0(stem, ".svg"))
  if (!file.copy(tmp_svg, svg_path, overwrite = TRUE)) {
    warning("SVG could not be copied to the requested output folder: ", svg_path)
  }
}

cat("\nFigure 1 saved to:\n", normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n")
cat("PDF: ", pdf_path, "\n", sep = "")
cat("PNG: ", png_path, "\n", sep = "")
