# Download daily FAERS counts for the verified events in events/events.csv.
# Needs: OPENFDA_API_KEY in .Renviron (never commit the key), internet.
# Usage from the project root:  source("application/08_event_fetch.R")
# Writes: data/faers_events/{event_id}.csv, checksums.txt, retrieved.txt
source("R/event_study.R")
suppressPackageStartupMessages({ library(httr2); library(digest) })

api_key <- Sys.getenv("OPENFDA_API_KEY")
if (api_key == "") stop("OPENFDA_API_KEY is not set. Add OPENFDA_API_KEY=... to ~/.Renviron and restart R.")
ev <- read_events("events/events.csv")
ev <- ev[ev$type == "event", ]
todo <- ev[tolower(ev$verified) == "yes", ]
if (nrow(ev) > nrow(todo)) message("Skipping ", nrow(ev) - nrow(todo), " event(s) with verified != yes: ", paste(ev$event_id[tolower(ev$verified) != "yes"], collapse = ", "))
if (nrow(todo) == 0) stop("No verified events. Set verified = yes in events/events.csv after checking date, URL and term.")

dir.create("data/faers_events", recursive = TRUE, showWarnings = FALSE)
log_lines <- character()
for (i in seq_len(nrow(todo))) {
  r <- todo[i, ]; out <- file.path("data/faers_events", paste0(r$event_id, ".csv"))
  if (file.exists(out)) { cat("exists, skipping:", out, "\n"); next }
  cat("fetching", r$event_id, "(", r$drug, "/", r$reaction, ")\n")
  drugs <- trimws(strsplit(r$drug, "|", fixed = TRUE)[[1]])
  ser <- fetch_faers_series(drugs, r$reaction, 2012, 2024, api_key)
  write.csv(ser, out, row.names = FALSE)
  line <- sprintf("%s: %d reports over %d days, %.1f%% zero days", r$event_id, sum(ser$count), nrow(ser), 100 * mean(ser$count == 0))
  cat("  ", line, "\n"); log_lines <- c(log_lines, line)
}
files <- list.files("data/faers_events", pattern = "\\.csv$", full.names = TRUE)
sums <- vapply(files, function(f) digest(f, algo = "sha256", file = TRUE), character(1))
writeLines(sprintf("%s  %s", sums, basename(files)), "data/faers_events/checksums.txt")
writeLines(c(paste("Retrieved:", format(Sys.time(), "%Y-%m-%d %H:%M %Z")), paste("R:", R.version.string),
             "Source: openFDA drug-event endpoint, count=receivedate, 2012-01-01 to 2024-12-31",
             "Multiple drugs in one event are summed per day (a report listing two of them counts twice).", log_lines),
           "data/faers_events/retrieved.txt", useBytes = TRUE)
cat("\nDone. Next: source('application/09_event_screen.R')\n")
