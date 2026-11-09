# Apply the series screen (protocol, Part B and D2) to each verified event, using that event's own training window.
# Usage from the project root:  source("application/09_event_screen.R")
# Writes: out/events/tables/event_screen.{csv,tex}
source("R/io.R"); source("R/event_study.R"); set_out_root("out/events")
ev <- read_events("events/events.csv")
ev <- ev[ev$type == "event" & tolower(ev$verified) == "yes", ]
if (nrow(ev) == 0) stop("No verified events in events/events.csv")
rows <- lapply(seq_len(nrow(ev)), function(i) {
  r <- ev[i, ]; f <- file.path("data/faers_events", paste0(r$event_id, ".csv"))
  if (!file.exists(f)) stop("Missing ", f, ". Run application/08_event_fetch.R first.")
  d <- read_series(f); w <- event_windows(r$event_date); sc <- screen_event(d, w)
  data.frame(event_id = r$event_id, event_date = r$event_date, train_end = as.character(w$train_end), n_train = sc$n_train,
             weekday_zero_share = round(sc$zero_share_wd, 3), mean_positive = round(sc$mean_pos_wd, 2),
             result = sc$reason, stringsAsFactors = FALSE)
})
scr <- do.call(rbind, rows)
write_table(scr, "event_screen", caption = "Series screen for each event, on that event's own training window (weekday zero share in [0.15, 0.85], mean positive count at least 1.3, at least 1500 training days).")
print(scr, row.names = FALSE)
npass <- sum(scr$result == "pass")
cat(sprintf("\n%d of %d events pass the screen.\n", npass, nrow(scr)))
if (npass < 10) cat("Fewer than 10 pass: the study is reported as a small case series with every screened event listed (no change of design).\n")
cat("Next: source('application/10_event_study.R')\n")
