# Event-window study: head-wise monitor and comparators on every event and on matched placebo windows.
# Usage from the project root:  source("application/10_event_study.R")
# Writes: out/events/tables/event_runs.csv, event_summary.{csv,tex}, head_attribution.csv,
#         out/events/raw/event_traces.rds, out/events/session_info.txt, out/events/constants.csv
source("R/io.R"); source("R/eprocess.R"); source("R/helpers.R"); source("R/forecast.R"); source("R/event_study.R")
suppressPackageStartupMessages(library(pscl))
set_out_root("out/events"); dir.create("out/events/tables", recursive = TRUE, showWarnings = FALSE)

ev <- read_events("events/events.csv")
ev <- ev[ev$type == "event" & tolower(ev$verified) == "yes", ]
if (nrow(ev) == 0) stop("No verified events in events/events.csv")
panel <- read.csv("data/panel.csv", stringsAsFactors = FALSE)

runs <- list(); traces <- list()
do_run <- function(run_id, type, event_id, series, d, event_date, expected_head) {
  cat(sprintf("  %-40s", run_id))
  res <- run_one(d, event_date, run_seed(run_id))
  cat(res$row$status, "\n")
  runs[[length(runs) + 1]] <<- as.data.frame(c(list(run_id = run_id, type = type, event_id = event_id, series = series,
                                                    expected_head = expected_head), res$row), stringsAsFactors = FALSE)
  if (!is.null(res$trace)) traces[[run_id]] <<- res$trace
  res$row$status
}

cat("Events\n")
ran_events <- character()
for (i in seq_len(nrow(ev))) {
  r <- ev[i, ]; f <- file.path("data/faers_events", paste0(r$event_id, ".csv"))
  if (!file.exists(f)) stop("Missing ", f, ". Run application/08_event_fetch.R first.")
  st <- do_run(r$event_id, "event", r$event_id, paste(r$drug, r$reaction, sep = " / "), read_series(f), r$event_date, r$expected_head)
  if (st == "run") ran_events <- c(ran_events, r$event_id)
}

cat("Placebo windows (same calendar dates, panel series assumed to have no event in the window)\n")
for (eid in ran_events) {
  r <- ev[ev$event_id == eid, ]
  for (j in seq_len(nrow(panel))) {
    same <- grepl(tolower(panel$drug[j]), tolower(r$drug), fixed = TRUE) && tolower(panel$reaction[j]) == tolower(r$reaction)
    if (same) next
    f <- file.path("data/faers", paste0(panel$id[j], ".csv"))
    do_run(paste0(eid, "@", panel$id[j]), "placebo", eid, paste(panel$drug[j], panel$reaction[j], sep = " / "), read_series(f), r$event_date, NA_character_)
  }
}

runs <- do.call(rbind, runs)
write.csv(runs, "out/events/tables/event_runs.csv", row.names = FALSE)
summ <- summarise_runs(runs)
write_table(summ, "event_summary", caption = "First-crossing outcome by method. For events: false alarm = first crossing in the quiet period, buffer = in the 30 days before the event, detected = within 180 days of the event, missed = none. For placebo windows any alarm is a false alarm. Comparator thresholds are calibrated by simulation under the frozen forecast.")
cat("\nSummary\n"); print(summ, row.names = FALSE)

det <- runs[runs$type == "event" & runs$status == "run" & runs$any_head_class == "detected" & runs$expected_head %in% c("occ", "mag"), ]
att <- data.frame(n_detected_with_expected_head = nrow(det), n_correct_head = sum(det$head_first == det$expected_head, na.rm = TRUE))
write.csv(att, "out/events/tables/head_attribution.csv", row.names = FALSE)
cat("\nHead attribution (first head to cross equals expected head):\n"); print(att, row.names = FALSE)

save_raw(traces, "event_traces")
writeLines(capture.output(sessionInfo()), "out/events/session_info.txt")
write.csv(data.frame(name = c(names(EV), "alpha", "tol_occ", "tol_mag", "calib_B", "pit_seed_pos", "pit_seed_whole"),
                     value = c(unlist(EV), EV_ALPHA, EV_TOL, EV_CALIB_B, EV_PIT_SEEDS)), "out/events/constants.csv", row.names = FALSE)
cat("\nDone. Results are in out/events/. Commit them (see README_EVENTS.md).\n")
