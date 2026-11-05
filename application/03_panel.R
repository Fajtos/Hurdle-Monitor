#Panel analysis, we run the two head monitor on every retained series.
# Reads:    data/faers/*.csv, data/panel.csv, R/io.R, R/forecast.R, R/eprocess.R
# Writes:   out/application/tables/panel.{csv,tex}
#           out/application/raw/panel_raw.rds

source("R/io.R")
set_out_root("out/application")
source("R/eprocess.R")
source("R/helpers.R")
source("R/forecast.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(MASS)
  library(pscl)
})

TRAIN_END <- as.Date("2017-12-31")
ALPHA     <- 0.05
TOL       <- c(0.02, 0.03)
THR_HEAD  <- log(2 / ALPHA)
THR_WHOLE <- log(1 / ALPHA)

rqr_hurdle_vec <- function(y, pi_pos, mu_ut, theta, rng_seed = NULL) {
  if (!is.null(rng_seed)) set.seed(rng_seed)
  n  <- length(y)
  q0 <- (theta / (theta + mu_ut))^theta
  U  <- numeric(n)
  pos <- which(y > 0)
  if (length(pos) > 0) {
    F_nb_y   <- pnbinom(y[pos],     size = theta, mu = mu_ut[pos])
    F_nb_ym1 <- pnbinom(y[pos] - 1, size = theta, mu = mu_ut[pos])
    F_tr_y   <- (F_nb_y   - q0[pos]) / (1 - q0[pos])
    F_tr_ym1 <- (F_nb_ym1 - q0[pos]) / (1 - q0[pos])
    F_full_y   <- (1 - pi_pos[pos]) + pi_pos[pos] * F_tr_y
    F_full_ym1 <- (1 - pi_pos[pos]) + pi_pos[pos] * F_tr_ym1
    V <- runif(length(pos))
    U[pos] <- F_full_ym1 + V * (F_full_y - F_full_ym1)
  }
  zero <- which(y == 0)
  if (length(zero) > 0) U[zero] <- runif(length(zero)) * (1 - pi_pos[zero])
  qnorm(U)
}

predict_hurdle <- function(fit_hur, X_te) {
  coef_zero  <- fit_hur$coefficients$zero
  coef_count <- fit_hur$coefficients$count
  theta_hur  <- as.numeric(fit_hur$theta)
  pi_pos_te <- as.numeric(plogis(X_te %*% coef_zero))
  mu_ut_te  <- as.numeric(exp(X_te %*% coef_count))
  list(pi_pos = pi_pos_te, mu_ut = mu_ut_te, theta = theta_hur)
}

first_str <- function(logE, thr, dates) {
  i <- first_cross(logE, thr)
  if (is.na(i)) "none" else as.character(dates[i])
}

panel <- read.csv("data/panel.csv", stringsAsFactors = FALSE)
rows  <- list()

for (i in seq_len(nrow(panel))) {
  id <- panel$id[i]
  cat(sprintf("\n=== %s: %s / %s ===\n", id, panel$drug[i], panel$reaction[i]))
  
  f <- sprintf("data/faers/%s.csv", id)
  if (!file.exists(f)) {
    message("  Skipping ", id, ": file missing")
    next
  }
  
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
  
  weekday_tr    <- format(df_tr$date, "%u") %in% c("1","2","3","4","5")
  zero_share_wd <- mean(df_tr$count[weekday_tr] == 0)
  
  if (zero_share_wd < 0.15 || zero_share_wd > 0.85) {
    message(sprintf("  Skipping %s: weekday zero share %.3f outside [0.15, 0.85]",
                    id, zero_share_wd))
    next
  }
  
  form  <- count ~ sin_year + cos_year + weekend + y_lag + y_roll7
  X_te  <- model.matrix(form, df_te)
  
  fit_ols  <- lm(form, data = df_tr)
  fit_pois <- glm(form, data = df_tr, family = poisson())
  fit_nb   <- tryCatch(MASS::glm.nb(form, data = df_tr),
                       error = function(e) NULL)
  fit_hur  <- tryCatch(pscl::hurdle(form, data = df_tr, dist = "negbin"),
                       error = function(e) {
                         message("  Hurdle fit failed for ", id, ": ",
                                 conditionMessage(e))
                         NULL
                       })
  
  if (is.null(fit_hur)) {
    message("  Skipping ", id, ": hurdle fit failed")
    next
  }
  
  sigma_ols <- summary(fit_ols)$sigma
  theta_nb  <- if (!is.null(fit_nb)) fit_nb$theta else NA_real_
  
  pred_ols  <- predict(fit_ols, newdata = df_te)
  pred_pois <- predict(fit_pois, newdata = df_te, type = "response")
  pred_nb   <- if (!is.null(fit_nb)) {
    predict(fit_nb, newdata = df_te, type = "response")
  } else rep(NA_real_, nrow(df_te))
  
  pred_h    <- predict_hurdle(fit_hur, X_te)
  pi_pos_te <- pred_h$pi_pos
  mu_ut_te  <- pred_h$mu_ut
  theta_hur <- pred_h$theta
  q0_te     <- (theta_hur / (theta_hur + mu_ut_te))^theta_hur
  
  # OLS log score on the same discretized basis as 02_case_study.R
  oos_ols <- ols_disc_logscore(df_te$count, pred_ols, sigma_ols)
  oos_pois <- -mean(dpois(df_te$count, pred_pois, log = TRUE))
  oos_nb   <- if (!is.null(fit_nb)) {
    -mean(dnbinom(df_te$count, size = theta_nb, mu = pred_nb, log = TRUE))
  } else NA_real_
  oos_hur  <- -mean(ifelse(
    df_te$count == 0,
    log(1 - pi_pos_te),
    log(pi_pos_te) +
      log(dnbinom(df_te$count, size = theta_hur, mu = mu_ut_te) / (1 - q0_te))
  ))
  
  oos_vec <- c(OLS = oos_ols, Poisson = oos_pois, NB = oos_nb, Hurdle = oos_hur)
  winner  <- names(oos_vec)[which.min(oos_vec)]
  if (winner == "Hurdle") winner <- "Hurdle-NB"
  
  cat(sprintf("  OOS: OLS=%.3f Poisson=%.3f NB=%s Hurdle=%.3f -> %s\n",
              oos_ols, oos_pois,
              ifelse(is.na(oos_nb), "NA", sprintf("%.3f", oos_nb)),
              oos_hur, winner))
  
  Z_te      <- as.integer(df_te$count > 0)
  U_pos_te  <- randomized_pit_positive(df_te$count, mu_ut_te, theta_hur, seed = 2025)
  rqr_te    <- rqr_hurdle_vec(df_te$count, pi_pos_te, mu_ut_te, theta_hur,
                              rng_seed = 2024)
  U_full_te <- pnorm(rqr_te)
  
  logE_occ_tol    <- eprocess_occurrence(pi_pos_te, Z_te, TOL[1])
  logE_mag_tol    <- eprocess_pit(U_pos_te,   Z_te, TOL[2])
  logE_whole_tol  <- eprocess_pit(U_full_te,  rep(1, length(U_full_te)), TOL[2])
  logE_occ_strict <- eprocess_occurrence(pi_pos_te, Z_te, 0)
  logE_mag_strict <- eprocess_pit(U_pos_te,   Z_te, 0)
  
  occ_tol_date    <- first_str(logE_occ_tol,    THR_HEAD,  df_te$date)
  mag_tol_date    <- first_str(logE_mag_tol,    THR_HEAD,  df_te$date)
  whole_tol_date  <- first_str(logE_whole_tol,  THR_WHOLE, df_te$date)
  occ_strict_date <- first_str(logE_occ_strict, THR_HEAD,  df_te$date)
  mag_strict_date <- first_str(logE_mag_strict, THR_HEAD,  df_te$date)
  
  cat(sprintf("  Alarms (tolerant): occ=%s mag=%s whole=%s\n",
              occ_tol_date, mag_tol_date, whole_tol_date))
  cat(sprintf("  Alarms (strict):   occ=%s mag=%s\n",
              occ_strict_date, mag_strict_date))
  
  # Diagnostic columns for the dilution discussion
  T_monitor <- length(df_te$count)
  rate_mag  <- logE_mag_tol[T_monitor] / T_monitor
  rate_whole <- logE_whole_tol[T_monitor] / T_monitor
  pi_mean   <- mean(pi_pos_te)
  pi_inv_cubed <- pi_mean^(-3)
  rate_ratio <- if (abs(rate_whole) > 1e-9) rate_mag / rate_whole else NA_real_
  
  rows[[length(rows) + 1]] <- data.frame(
    id = id,
    series = paste(panel$drug[i], panel$reaction[i], sep = " / "),
    zero_share = round(zero_share_wd, 3),
    oos_winner = winner,
    occ_tol    = occ_tol_date,
    mag_tol    = mag_tol_date,
    whole_tol  = whole_tol_date,
    occ_strict = occ_strict_date,
    mag_strict = mag_strict_date,
    occ_terminal   = round(logE_occ_tol[length(logE_occ_tol)], 2),
    mag_terminal   = round(logE_mag_tol[length(logE_mag_tol)], 2),
    whole_terminal = round(logE_whole_tol[length(logE_whole_tol)], 2),
    pi_mean        = round(pi_mean, 3),
    rate_mag       = round(rate_mag, 5),
    rate_whole     = round(rate_whole, 5),
    rate_ratio     = round(rate_ratio, 2),
    pi_inv_cubed   = round(pi_inv_cubed, 2)
  )
}

panel_tbl <- do.call(rbind, rows)

write_table(panel_tbl, "panel",
            caption = "Panel of first-alarm dates and terminal log e-process for the five retained series under tolerant and strict margins. The OOS winner column reports the model with the lowest out-of-sample log score among OLS, Poisson, NB, and Hurdle-NB. The rate columns give the mean per-day log e-process growth (terminal log E divided by the monitoring-window length) for the magnitude head and the whole-PIT test, and rate\\_ratio is their ratio. $\\pi$\\_mean is the mean of the fitted occurrence probability over the monitoring window; $\\pi^{-3}$ is the dilution factor Proposition 1 predicts for the growth-rate ratio under a magnitude-only alternative with constant $\\pi$. The observed ratio only matches the theoretical prediction when both rates are positive; where the whole-PIT rate is negative, the theoretical comparison does not apply and the whole-PIT has accumulated evidence \\emph{against} a drift.")

cat("\n=== Panel summary ===\n\n")
print(panel_tbl, row.names = FALSE)
cat("\nPanel table written to out/application/tables/panel.{csv,tex}\n")

save_raw(list(panel_tbl = panel_tbl, panel = panel), "panel_raw")