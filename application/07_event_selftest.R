# Self-test of the event-study code on simulated data. No internet, writes nothing to data/ or out/.
# Run before the event study. From the project root:  source("application/07_event_selftest.R")
source("R/io.R"); source("R/eprocess.R"); source("R/helpers.R"); source("R/forecast.R"); source("R/event_study.R")
suppressPackageStartupMessages(library(pscl))

dates <- seq(as.Date("2012-01-01"), as.Date("2024-12-31"), by = "day")
n <- length(dates); doy <- as.integer(format(dates, "%j")); wk <- as.integer(format(dates, "%u")) <= 5
sim_series <- function(shift_date = NULL, shift = 1.5, seed = 1) {
  set.seed(seed)
  eta <- ifelse(wk, 0, -3.5)
  if (!is.null(shift_date)) eta <- eta + ifelse(dates >= shift_date & wk, shift, 0)   # occurrence shift only
  p_occ <- plogis(eta); mu <- exp(0.8 + 0.3 * sin(2 * base::pi * doy / 365.25)); theta <- 3
  q0 <- (theta / (theta + mu))^theta
  y_pos <- qnbinom(q0 + runif(n) * (1 - q0), size = theta, mu = mu)
  data.frame(date = dates, count = as.integer(ifelse(runif(n) < p_occ, y_pos, 0L)))
}
chk <- function(label, ok) cat(sprintf("%-70s %s\n", label, if (isTRUE(ok)) "OK" else "CHECK"))

ev_date <- as.Date("2019-06-01")
w <- event_windows(ev_date)
chk("quiet period is 365 days", as.integer(w$quiet[2] - w$quiet[1]) + 1L == 365L)
chk("buffer is 30 days", as.integer(w$buffer[2] - w$buffer[1]) + 1L == 30L)
chk("detection window is 180 days", as.integer(w$detect[2] - w$detect[1]) + 1L == 180L)
chk("training ends the day before the quiet period", w$train_end == w$quiet[1] - 1)

d0 <- sim_series(NULL, seed = 11); d1 <- sim_series(ev_date, seed = 12)
fc <- make_forecast(d1, w)
prev <- d1$count[match(fc$te$date[1] - 1, d1$date)]
chk("no look-ahead: first lag feature equals log(1 + previous day's count)", isTRUE(all.equal(fc$te$y_lag[1], log1p(prev))))

Z1 <- as.integer(fc$te$count > 0)
chk("restart-free mixture matches R/eprocess.R (occurrence)", isTRUE(all.equal(mix_log(occ_steps(fc$pi, Z1, EV_TOL[1])), eprocess_occurrence(fc$pi, Z1, EV_TOL[1]))))

cat("\nRunning the no-change series and the shifted series (about a minute)...\n")
r0 <- run_one(d0, ev_date, seed = 1)$row
r1 <- run_one(d1, ev_date, seed = 2)$row
cat("\nNo-change series: first crossing by method\n")
print(unlist(r0[paste0(EV_METHODS, "_class")]))
cat("\nSeries with an occurrence shift at the event date: first crossing by method\n")
print(unlist(r1[paste0(EV_METHODS, "_class")])); cat("Delay in days (occ head):", r1$occ_delay, "\n")

chk("no-change series: occurrence and magnitude heads stay quiet", all(c(r0$occ_class, r0$mag_class) == "none"))
chk("shifted series: occurrence head detects inside the window", identical(r1$occ_class, "detected"))
chk("shifted series: the occurrence head fires first, not the magnitude head", identical(r1$head_first, "occ"))
cat("\nIf every line says OK, the pipeline works as intended. A CHECK on the last three lines can be chance in one\n",
    "random series; a CHECK on the first six lines indicates a code problem.\n", sep = "")
