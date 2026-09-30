packages <- c(
  "readxl", "readr", "dplyr", "tidyr", "purrr", "tibble",
  "ggplot2", "patchwork", "scales", "stringr", "igraph",
  "ggraph", "ggforce", "sandwich", "lmtest", "svglite"
)

missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
} else {
  message("All required R packages are installed.")
}
