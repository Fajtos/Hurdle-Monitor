# Head-specific refit, done properly (replaces 05_refit.R as the paper's refit experiment).


source("R/io.R"); set_out_root("out/application")
source("R/eprocess.R"); source("R/helpers.R")
suppressPackageStartupMessages({ library(dplyr); library(pscl) })

REFIT_IDS <- c("r2", "r3", "r6")
TRAIN_END <- as.Date("2017-12-31"); ALPHA <- 0.05; TOL <- c(0.02, 0.03)
THR_HEAD <- log(2 / ALPHA); REFIT_L <- 500; PIT_SEED <- 2025
form <- count ~ sin_year + cos_year + weekend + y_lag + y_roll7

forecast_parts <- function(fit, df) {
  X <- model.matrix(form, df)
  list(pi = as.numeric(plogis(X %*% fit$coefficients$zero)),
       mu = as.numeric(exp(X %*% fit$coefficients$count)),
       theta = as.numeric(fit$theta))
}
monitor <- function(fit, df, seed = PIT_SEED) {
  fp <- forecast_parts(fit, df); y <- df$count; Z <- as.integer(y > 0)
  U <- rep(0.5, length(y)); pos <- which(y > 0)
  if (length(pos)) {
    set.seed(seed)
    q0 <- (fp$theta / (fp$theta + fp$mu[pos]))^fp$theta
    Fy   <- (pnbinom(y[pos],     size = fp$theta, mu = fp$mu[pos]) - q0) / (1 - q0)
    Fym1 <- (pnbinom(y[pos] - 1, size = fp$theta, mu = fp$mu[pos]) - q0) / (1 - q0)
    U[pos] <- Fym1 + runif(length(pos)) * (Fy - Fym1)
  }
  lo <- eprocess_occurrence(fp$pi, Z, TOL[1]); lm <- eprocess_pit(U, Z, TOL[2])
  list(occ = lo, mag = lm,
       score = -mean(hurdle_logpmf(y, fp$pi, fp$mu, fp$theta)))
}
first_str <- function(logE, dates) { i <- first_cross(logE, THR_HEAD); if (is.na(i)) "none" else as.character(dates[i]) }
mix_fit <- function(base, win, occ_from_win, mag_from_win) {
  f <- base
  if (occ_from_win) f$coefficients$zero  <- win$coefficients$zero
  if (mag_from_win) { f$coefficients$count <- win$coefficients$count; f$theta <- win$theta }
  f
}

rows <- list()
for (id in REFIT_IDS) {
  d <- read.csv(sprintf("data/faers/%s.csv", id)); d$date <- as.Date(d$date)
  d <- build_features(arrange(filter(d, !is.na(count)), date))
  df_tr <- d[d$date <= TRAIN_END, ]; df_te <- d[d$date > TRAIN_END, ]
  fit0 <- pscl::hurdle(form, data = df_tr, dist = "negbin")
  m0 <- monitor(fit0, df_te)
  io <- first_cross(m0$occ, THR_HEAD); im <- first_cross(m0$mag, THR_HEAD)
  if (is.na(io) && is.na(im)) { message(id, ": no baseline alarm, skipped"); next }
  flagged <- if (is.na(io)) "magnitude" else if (is.na(im)) "occurrence" else if (im <= io) "magnitude" else "occurrence"
  other_alarmed <- if (flagged == "magnitude") !is.na(io) else !is.na(im)
  adate <- df_te$date[min(c(io, im), na.rm = TRUE)]
  win <- tryCatch(pscl::hurdle(form, data = d[d$date > adate - REFIT_L & d$date <= adate, ], dist = "negbin"),
                  error = function(e) NULL)
  if (is.null(win)) { message(id, ": refit window fit failed, skipped"); next }
  post <- d[d$date > adate, ]
  fm <- flagged == "magnitude"
  arms <- list(
    "A0 no refit"            = fit0,
    "B refit flagged head"   = mix_fit(fit0, win, occ_from_win = !fm, mag_from_win = fm),
    "C refit other head"     = mix_fit(fit0, win, occ_from_win = fm,  mag_from_win = !fm),
    "D refit both heads"     = mix_fit(fit0, win, TRUE, TRUE))
  for (a in names(arms)) {
    m <- monitor(arms[[a]], post)
    rows[[length(rows) + 1]] <- data.frame(
      id = id, flagged = flagged, other_also_alarmed = other_alarmed, alarm_date = as.character(adate),
      arm = a, post_days = nrow(post),
      occ_alarm = first_str(m$occ, post$date), mag_alarm = first_str(m$mag, post$date),
      occ_terminal = round(tail(m$occ, 1), 2), mag_terminal = round(tail(m$mag, 1), 2),
      post_logscore = round(m$score, 3))
  }
}
tab <- do.call(rbind, rows)
write_table(tab, "refit_panel",
  caption = "Head-specific refit on the series with a baseline alarm. From the first-alarm date, four arms are monitored on the same post-alarm window with fresh e-processes and the same PIT seed. Refit window: last 500 days ending at the alarm. post\\_logscore is the mean negative log score on the post-alarm window (lower is better).")
print(tab, row.names = FALSE)
