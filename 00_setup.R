if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

renv::init(bare = TRUE, restart = FALSE)

packages <- c(
  "pscl", "MASS",
  "dplyr", "tidyr", "purrr", "ggplot2",
  "httr2", "jsonlite", "digest",
  "knitr"
)

renv::install(packages, prompt = FALSE)
renv::snapshot(type = "all", prompt = FALSE)

cat("\nrenv.lock written to:", file.path(getwd(), "renv.lock"), "\n")
cat("Commit renv.lock, .Rprofile, and renv/ to the repository before submission.\n")