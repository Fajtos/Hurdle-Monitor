# To document the retrieved FAERS panel,
# README for the data folder, and SHA-256 checksums for reproducibility.
#Prerequisites: data/faers/{id}.csv for each id in data/panel.csv
# Usage (from the project root):
#   source("application/00_document_data.R")
# Outputs:
#   out/application/tables/data_overview.tex and .csv
#   data/faers/README.md
#   data/faers/checksums.txt

suppressPackageStartupMessages({
  library(digest)
})

source("R/io.R")
set_out_root("out/application")

if (!file.exists("data/panel.csv")) stop("data/panel.csv is missing.")
panel <- read.csv("data/panel.csv", stringsAsFactors = FALSE)
TRAIN_END <- as.Date("2017-12-31")

rows <- list()
for (i in seq_len(nrow(panel))) {
  id <- panel$id[i]
  f <- sprintf("data/faers/%s.csv", id)
  if (!file.exists(f)) {
    message("Skipping ", id, ": ", f, " not found")
    next
  }
  d <- read.csv(f, stringsAsFactors = FALSE)
  d$date <- as.Date(d$date)
  weekday <- format(d$date, "%u") %in% c("1", "2", "3", "4", "5")
  wd_train <- d$count[weekday & d$date <= TRAIN_END]
  pos_train <- wd_train[wd_train > 0]
  rows[[length(rows) + 1]] <- data.frame(
    id = id,
    drug = panel$drug[i],
    reaction = panel$reaction[i],
    total = sum(d$count),
    weekday_zero_train = round(mean(wd_train == 0), 3),
    positive_mean_train = round(mean(pos_train), 2),
    max_count = max(d$count),
    first_date = as.character(min(d$date)),
    last_date = as.character(max(d$date)),
    n_days = nrow(d)
  )
}

overview <- do.call(rbind, rows)
write_table(overview, "data_overview",
            caption = "Panel of FAERS drug-reaction series. Weekday zero share and positive-count mean are computed on the training window, 2012-01-01 to 2017-12-31.")

cat("\nData overview:\n\n")
print(overview, row.names = FALSE)

readme <- c(
  "# FAERS panel data",
  "",
  sprintf("Retrieved: %s", format(Sys.Date(), "%Y-%m-%d")),
  sprintf("Source: openFDA drug-event endpoint, https://api.fda.gov/drug/event.json"),
  "",
  "## Files",
  "",
  paste0("- ", list.files("data/faers", pattern = "\\.csv$")),
  "",
  "## Columns",
  "",
  "- `date`: calendar date, daily from 2012-01-01 to 2024-12-31 (4749 rows).",
  "- `count`: integer count of adverse-event reports for the series on that date.",
  "",
  "## Panel definition",
  "",
  "Each series is defined by the generic drug name and, where present, a",
  "MedDRA reaction term. See `data/panel.csv` for the exact definitions.",
  "",
  "## Retrieval protocol",
  "",
  "The series were retrieved by `application/01_retrieve_faers.R` using the",
  "openFDA drug-event endpoint with `count=receivedate`, one request per",
  "two-year window, and days without any reports filled as zero. Retrieval",
  "used an authenticated api.data.gov key stored outside the repository.",
  "",
  "## Verification",
  "",
  "SHA-256 checksums of the CSVs are in `checksums.txt`. If the retrieval is",
  "re-run and the checksums differ, the upstream data have changed or the",
  "panel definition has been modified."
)
writeLines(readme, "data/faers/README.md")

files <- list.files("data/faers", pattern = "\\.csv$", full.names = TRUE)
checks <- sapply(files, function(f) digest::digest(f, algo = "sha256", file = TRUE))
checksum_lines <- sprintf("%s  %s", checks, basename(files))
writeLines(checksum_lines, "data/faers/checksums.txt")
