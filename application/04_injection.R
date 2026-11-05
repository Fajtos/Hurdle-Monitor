# Semi synthetic validation to inject a known failure into each panel series
# and confirm that the monitor localizes it to the correct head.
#
# For each series:
#   1. Fit the hurdle model on the original training window.
#   2. Freeze parameters and the covariate trajectory at pre-injection values.
#   3. Apply the occurrence injection (or magnitude injection) from 2019-06-01.
#   4. Run the two-head monitor and record first-alarm dates.
#
# Reads:    data/faers/*.csv, data/panel.csv, R/io.R, R/forecast.R, R/eprocess.R
# Writes:   out/application/tables/injection.{csv,tex}
#           out/application/raw/injection_raw.rds
#
# Usage (from the project root):
#   source("application/04_injection.R")


source("R/io.R")
set_out_root("out/application")
source("R/eprocess.R")
source("R/forecast.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(pscl)
})

TRAIN_END    <- as.Date("2017-12-31")
INJECT_FROM  <- as.Date("2019-06-01")
ALPHA        <- 0.05
TOL          <- c(0.02, 0.03)
THR_HEAD     <- log(2 / ALPHA)
PIT_SEED     <- 2025

inject_occ <- function(y, dates, inject_from, seed, p_drop = 0.4) {
  set.seed(seed)
  y_inj <- y
  m <- dates >= inject_from & y > 0
  drop <- m & (runif(length(y)) < p_drop)
  y_inj[drop] <- 0L
  y_inj
}

inject_mag <- function(y, dates, inject_from, seed, rate = 0.5) {
  set.seed(seed)
  y_inj <- y
  m <- dates >= inject_from & y > 0
  y_inj[m] <- y[m] + rpois(sum(m), rate * y[m])
  y_inj
}

fit_and_forecast <- function(df_tr, df_te, form) {
  fit_hur <- tryCatch(pscl::hurdle(form, data = df_tr, dist = "negbin"),
                      error = function(e) {
                        message("  Hurdle fit failed: ", conditionMessage(e))
                        NULL
                      })
  if (is.null(fit_hur)) return(NULL)
  X_te <- model.matrix(form, df_te)
  coef_zero  <- fit_hur$coefficients$zero
  coef_count <- fit_hur$coefficients$count
  theta_hur  <- as.numeric(fit_hur$theta)
  pi_pos_te  <- as.numeric(plogis(X_te %*% coef_zero))
  mu_ut_te   <- as.numeric(exp(X_te %*% coef_count))
  list(pi_pos = pi_pos_te, mu_ut = mu_ut_te, theta = theta_hur, fit = fit_hur)
}

run_monitor <- function(y_eval, pi_pos_te, mu_ut_te, theta, dates_te, seed) {
  Z_te <- as.integer(y_eval > 0)
  U_pos <- rep(0.5, length(y_eval))
  pos   <- which(y_eval > 0)
  if (length(pos) > 0) {
    set.seed(seed)
    q0    <- (theta / (theta + mu_ut_te[pos]))^theta
    F_y   <- pnbinom(y_eval[pos],     size = theta, mu = mu_ut_te[pos])
    F_ym1 <- pnbinom(y_eval[pos] - 1, size = theta, mu = mu_ut_te[pos])
    Fy_tr   <- (F_y   - q0) / (1 - q0)
    Fym1_tr <- (F_ym1 - q0) / (1 - q0)
    V <- runif(length(pos))
    U_pos[pos] <- Fym1_tr + V * (Fy_tr - Fym1_tr)
  }
  logE_occ <- eprocess_occurrence(pi_pos_te, Z_te, TOL[1])
  logE_mag <- eprocess_pit(U_pos, Z_te, TOL[2])
  i_occ <- first_cross(logE_occ, THR_HEAD)
  i_mag <- first_cross(logE_mag, THR_HEAD)
  list(
    occ_alarm    = if (is.na(i_occ)) "none" else as.character(dates_te[i_occ]),
    mag_alarm    = if (is.na(i_mag)) "none" else as.character(dates_te[i_mag]),
    occ_terminal = round(logE_occ[length(logE_occ)], 2),
    mag_terminal = round(logE_mag[length(logE_mag)], 2)
  )
}

panel <- read.csv("data/panel.csv", stringsAsFactors = FALSE)
rows  <- list()

for (i in seq_len(nrow(panel))) {
  id <- panel$id[i]
  f  <- sprintf("data/faers/%s.csv", id)
  if (!file.exists(f)) {
    message("Skipping ", id, ": file missing")
    next
  }
  
  cat(sprintf("\n=== %s: %s / %s ===\n", id, panel$drug[i], panel$reaction[i]))
  
  d <- read.csv(f, stringsAsFactors = FALSE)
  d$date <- as.Date(d$date)
  d <- d |> filter(!is.na(count)) |> arrange(date)
  n <- nrow(d)
  
  doy <- as.integer(format(d$date, "%j"))
  d$sin_year <- sin(2 * pi * doy / 365.25)
  d$cos_year <- cos(2 * pi * doy / 365.25)
  d$weekend  <- as.integer(format(d$date, "%u") %in% c("6", "7"))
  d$y_lag    <- c(0, log1p(d$count[-n]))
  roll7 <- as.numeric(stats::filter(d$count, rep(1/7, 7), sides = 1))
  roll7[is.na(roll7)] <- 0
  d$y_roll7 <- log1p(c(0, roll7[-n]))
  
  n0    <- sum(d$date <= TRAIN_END)
  df_tr <- d[1:n0, ]
  df_te <- d[(n0 + 1):n, ]
  
  form <- count ~ sin_year + cos_year + weekend + y_lag + y_roll7
  fc   <- fit_and_forecast(df_tr, df_te, form)
  if (is.null(fc)) { message("  Skipping ", id); next }
  
  y_orig   <- df_te$count
  dates_te <- df_te$date
  seed_base <- 5000 + i * 10
  
  base <- run_monitor(y_orig, fc$pi_pos, fc$mu_ut, fc$theta,
                      dates_te, PIT_SEED)
  
  y_occ <- inject_occ(y_orig, dates_te, INJECT_FROM, seed_base + 1)
  occ   <- run_monitor(y_occ, fc$pi_pos, fc$mu_ut, fc$theta,
                       dates_te, PIT_SEED)
  
  y_mag <- inject_mag(y_orig, dates_te, INJECT_FROM, seed_base + 2)
  mag   <- run_monitor(y_mag, fc$pi_pos, fc$mu_ut, fc$theta,
                       dates_te, PIT_SEED)
  
  cat(sprintf("  Baseline:  occ=%s mag=%s\n", base$occ_alarm, base$mag_alarm))
  cat(sprintf("  Occ.inj.:  occ=%s mag=%s\n", occ$occ_alarm, occ$mag_alarm))
  cat(sprintf("  Mag.inj.:  occ=%s mag=%s\n", mag$occ_alarm, mag$mag_alarm))
  
  rows[[length(rows) + 1]] <- data.frame(
    id                 = id,
    series             = paste(panel$drug[i], panel$reaction[i], sep = " / "),
    baseline_occ       = base$occ_alarm,
    baseline_mag       = base$mag_alarm,
    occinj_occ         = occ$occ_alarm,
    occinj_mag         = occ$mag_alarm,
    occinj_occ_term    = occ$occ_terminal,
    occinj_mag_term    = occ$mag_terminal,
    maginj_occ         = mag$occ_alarm,
    maginj_mag         = mag$mag_alarm,
    maginj_occ_term    = mag$occ_terminal,
    maginj_mag_term    = mag$mag_terminal
  )
}

inj_tbl <- do.call(rbind, rows)

write_table(inj_tbl, "injection",
            caption = "Semi-synthetic validation. Occurrence injection sets each positive day to zero with probability $0.4$ from 2019-06-01; magnitude injection adds a $\\mathrm{Poisson}(0.5\\,y_t)$ count to each positive day from 2019-06-01. The forecast is frozen at its pre-injection parameters and covariate trajectory. A well-calibrated monitor flags the injected head.")

cat("\n=== Injection summary ===\n\n")
print(inj_tbl, row.names = FALSE)
cat("\nTable written to out/application/tables/injection.{csv,tex}\n")

save_raw(list(inj_tbl = inj_tbl, panel = panel), "injection_raw")