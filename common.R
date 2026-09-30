script_path <- function() {
  argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(argument) == 0) return(normalizePath(getwd()))
  normalizePath(sub("^--file=", "", argument[[1]]), mustWork = FALSE)
}

script_dir <- dirname(script_path())
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)

default_data_path <- file.path(
  project_dir, "data", "20260925_TMD_Sleep_VitaminD_Analysis_Master.xlsx"
)
data_path <- Sys.getenv("TMD_SLEEP_DATA", unset = default_data_path)

if (!file.exists(data_path)) {
  stop(
    "Data file not found: ", data_path, "\n",
    "Copy the analysis workbook to data/ or set TMD_SLEEP_DATA."
  )
}

figure_output_dir <- file.path(project_dir, "outputs", "figures")
supporting_output_dir <- file.path(project_dir, "outputs", "supporting")
dir.create(figure_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supporting_output_dir, recursive = TRUE, showWarnings = FALSE)

save_plot_formats <- function(plot_object, stem, width, height, dpi = 450) {
  png_path <- paste0(stem, ".png")
  pdf_path <- paste0(stem, ".pdf")
  svg_path <- paste0(stem, ".svg")

  ggplot2::ggsave(
    png_path, plot_object, width = width, height = height,
    units = "in", dpi = dpi, bg = "white"
  )
  grDevices::pdf(
    pdf_path, width = width, height = height,
    family = "sans", useDingbats = FALSE, onefile = TRUE
  )
  print(plot_object)
  grDevices::dev.off()

  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(
      svg_path, plot_object, width = width, height = height,
      units = "in", device = svglite::svglite, bg = "white"
    )
  }
  invisible(c(png = png_path, pdf = pdf_path, svg = svg_path))
}
