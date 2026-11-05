#  To retrieve daily adverse event report counts from the public openFDA FAERS
# drug event endpoint for the five drug reaction pairs defined in
# data/panel.csv. Each pair produces one CSV in data/faers/{id}.csv with
# two columns: date (daily, 2012-01-01 through 2024-12-31, 4749 rows)
# and count (integer, zero on days with no reports).
# Prerequisites:
#   data/panel.csv exists with columns id, drug, reaction, role, note
#   OPENFDA_API_KEY environment variable is set
#
# Usage (from the project root):
#   source("application/01_retrieve_faers.R")
# Outputs:
#   data/faers/{id}.csv  for each row of data/panel.csv
#   console diagnostics per series and a summary table

suppressPackageStartupMessages({
  library(httr2)
  library(jsonlite)
})

api_key <- Sys.getenv("OPENFDA_API_KEY")
if (api_key == "") {
  stop("OPENFDA_API_KEY is not set. Set it in .Renviron or via Sys.setenv().")
}
if (!file.exists("data/panel.csv")) {
  stop("data/panel.csv is missing.")
}

panel <- read.csv("data/panel.csv", stringsAsFactors = FALSE)
required_cols <- c("id", "drug", "reaction")
missing_cols <- setdiff(required_cols, colnames(panel))
if (length(missing_cols) > 0) {
  stop(sprintf("data/panel.csv missing columns: %s",
               paste(missing_cols, collapse = ", ")))
}

dir.create("data/faers", recursive = TRUE, showWarnings = FALSE)

url_escape <- function(x) {
  x <- gsub(" ", "%20", x, fixed = TRUE)
  x <- gsub("\t", "%09", x, fixed = TRUE)
  x <- gsub("\n", "%0A", x, fixed = TRUE)
  x <- gsub("\r", "%0D", x, fixed = TRUE)
  x <- gsub("%", "%25", x, fixed = TRUE)
  x <- gsub("%2520", "%20", x, fixed = TRUE)
  x <- gsub("%2509", "%09", x, fixed = TRUE)
  x <- gsub("%250A", "%0A", x, fixed = TRUE)
  x <- gsub("%250D", "%0D", x, fixed = TRUE)
  x
}

fetch_one <- function(drug, reaction, start = 2012, end = 2024,
                      sleep = 0.4, max_retries = 3, debug = FALSE) {
  drug_enc <- url_escape(drug)
  reaction_enc <- if (!is.na(reaction) && reaction != "") url_escape(reaction) else ""
  rows <- list()
  for (y0 in seq(start, end, by = 2)) {
    y1 <- min(y0 + 1, end)
    parts <- sprintf('patient.drug.openfda.generic_name:"%s"', drug_enc)
    if (reaction_enc != "") {
      parts <- c(parts, sprintf('patient.reaction.reactionmeddrapt:"%s"', reaction_enc))
    }
    parts <- c(parts, sprintf("receivedate:[%d0101+TO+%d1231]", y0, y1))
    search <- paste(parts, collapse = "+AND+")
    url <- sprintf(
      "https://api.fda.gov/drug/event.json?search=%s&count=receivedate&limit=1000&api_key=%s",
      search, api_key
    )
    if (debug && y0 == start) {
      cat("\n    [debug] URL:\n    ", url, "\n\n", sep = "")
    }
    resp <- NULL
    for (attempt in seq_len(max_retries)) {
      resp <- tryCatch(
        request(url) |> req_perform(),
        error = function(e) {
          if (attempt == max_retries) {
            message("  HTTP failure on ", drug, " ", reaction, " ", y0,
                    ": ", conditionMessage(e))
            return(NULL)
          }
          Sys.sleep(2 * attempt)
          NULL
        }
      )
      if (!is.null(resp)) break
    }
    if (!is.null(resp)) {
      js <- tryCatch(resp_body_json(resp, simplifyVector = TRUE),
                     error = function(e) NULL)
      if (!is.null(js$results) && length(js$results) > 0 &&
          is.data.frame(js$results)) {
        rows[[length(rows) + 1]] <- js$results
      }
    }
    Sys.sleep(sleep)
  }
  
  full <- data.frame(
    date = seq(as.Date(sprintf("%d-01-01", start)),
               as.Date(sprintf("%d-12-31", end)), by = "day"),
    count = 0L
  )
  if (length(rows) > 0) {
    df <- do.call(rbind, rows)
    df$date_parsed <- as.Date(df$time, format = "%Y%m%d")
    idx <- match(df$date_parsed, full$date)
    valid <- !is.na(idx)
    full$count[idx[valid]] <- as.integer(df$count[valid])
  }
  full
}

diagnose <- function(d, train_end = as.Date("2017-12-31")) {
  weekday <- format(d$date, "%u") %in% c("1", "2", "3", "4", "5")
  wd_train <- d$count[weekday & d$date <= train_end]
  pos_train <- wd_train[wd_train > 0]
  list(
    n_rows = nrow(d),
    total_reports = sum(d$count),
    zero_share_train = if (length(wd_train) > 0) mean(wd_train == 0) else NA_real_,
    positive_mean_train = if (length(pos_train) > 0) mean(pos_train) else NA_real_,
    max_count = max(d$count)
  )
}

cat("Retrieving FAERS report counts for", nrow(panel), "series.\n\n")

summary_rows <- list()
for (i in seq_len(nrow(panel))) {
  id <- panel$id[i]
  cat(sprintf("[%d/%d] %s: %s / %s ... ",
              i, nrow(panel), id, panel$drug[i], panel$reaction[i]))
  df <- fetch_one(panel$drug[i], panel$reaction[i],
                  debug = (i == 1))
  write.csv(df, sprintf("data/faers/%s.csv", id), row.names = FALSE)
  d <- diagnose(df)
  cat(sprintf("rows=%d, weekday_zero=%.3f, positive_mean=%.2f\n",
              d$n_rows, d$zero_share_train, d$positive_mean_train))
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    id = id,
    drug = panel$drug[i],
    reaction = panel$reaction[i],
    n_rows = d$n_rows,
    total_reports = d$total_reports,
    weekday_zero_share_train = round(d$zero_share_train, 3),
    positive_mean_train = round(d$positive_mean_train, 2),
    max_count = d$max_count
  )
}

summary_df <- do.call(rbind, summary_rows)
cat("\nSummary:\n\n")
print(summary_df, row.names = FALSE)

cat("\nFiles in data/faers/:\n")
for (f in list.files("data/faers", pattern = "\\.csv$")) {
  cat("  ", f, "\n", sep = "")
}

keep <- summary_df$weekday_zero_share_train >= 0.15 &
  summary_df$weekday_zero_share_train <= 0.85
cat("\nScreening rule (weekday zero share in [0.15, 0.85]):\n")
for (i in seq_len(nrow(summary_df))) {
  cat(sprintf("  %s %s (%.3f)\n",
              ifelse(keep[i], "RETAIN ", "EXCLUDE"),
              summary_df$id[i],
              summary_df$weekday_zero_share_train[i]))
}
cat(sprintf("\n%d of %d series retained.\n", sum(keep), length(keep)))