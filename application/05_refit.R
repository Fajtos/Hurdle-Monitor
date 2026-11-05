# Head-specific refit experiment.


source("R/io.R")
set_out_root("out/application")
source("R/eprocess.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(pscl)
})

CASE_ID    <- "r2"
TRAIN_END  <- as.Date("2017-12-31")
ALPHA      <- 0.05
TOL        <- c(0.02, 0.03)
THR_HEAD   <- log(2 / ALPHA)
REFIT_L    <- 500
PIT_SEED   <- 2025

d <- read.csv(sprintf("data/faers/%s.csv", CASE_ID), stringsAsFactors = FALSE)
d$date <- as.Date(d$date)
d <- d |> filter(!is.na(count)) |> arrange(date)
n <- nrow(d)

build_features <- function(d) {
  n <- nrow(d)
  doy <- as.integer(format(d$date, "%j"))
  d$sin_year <- sin(2 * pi * doy / 365.25)
  d$cos_year <- cos(2 * pi * doy / 365.25)
  d$weekend  <- as.integer(format(d$date, "%u") %in% c("6", "7"))
  d$y_lag    <- c(0, log1p(d$count[-n]))
  roll7 <- as.numeric(stats::filter(d$count, rep(1/7, 7), sides = 1))
  roll7[is.na(roll7)] <- 0
  d$y_roll7 <- log1p(c(0, roll7[-n]))
  d
}
d <- build_features(d)

n0    <- sum(d$date <= TRAIN_END)
df_tr <- d[1:n0, ]
form  <- count ~ sin_year + cos_year + weekend + y_lag + y_roll7

fit_orig <- pscl::hurdle(form, data = df_tr, dist = "negbin")

monitor_full <- function(fit, df_eval, pit_seed = PIT_SEED) {
  X_eval <- model.matrix(form, df_eval)
  pi_pos <- as.numeric(plogis(X_eval %*% fit$coefficients$zero))
  mu_ut  <- as.numeric(exp(X_eval %*% fit$coefficients$count))
  theta  <- as.numeric(fit$theta)
  Z      <- as.integer(df_eval$count > 0)
  
  U_pos <- rep(0.5, nrow(df_eval))
  pos <- which(df_eval$count > 0)
  if (length(pos) > 0) {
    set.seed(pit_seed)
    q0  <- (theta / (theta + mu_ut[pos]))^theta
    F_y   <- pnbinom(df_eval$count[pos],     size = theta, mu = mu_ut[pos])
    F_ym1 <- pnbinom(df_eval$count[pos] - 1, size = theta, mu = mu_ut[pos])
    Fy_tr   <- (F_y   - q0) / (1 - q0)
    Fym1_tr <- (F_ym1 - q0) / (1 - q0)
    V <- runif(length(pos))
    U_pos[pos] <- Fym1_tr + V * (Fy_tr - Fym1_tr)
  }
  logE_occ <- eprocess_occurrence(pi_pos, Z, TOL[1])
  logE_mag <- eprocess_pit(U_pos, Z, TOL[2])
  list(logE_occ = logE_occ, logE_mag = logE_mag, dates = df_eval$date)
}

# Original baseline: monitor post-2017 window
df_te   <- d |> filter(date > TRAIN_END)
m_orig  <- monitor_full(fit_orig, df_te)

occ_alarm_idx <- first_cross(m_orig$logE_occ, THR_HEAD)
mag_alarm_idx <- first_cross(m_orig$logE_mag, THR_HEAD)

if (is.na(occ_alarm_idx) && is.na(mag_alarm_idx)) {
  stop("No alarm in the baseline; cannot run the refit experiment.")
}

first_alarm_idx  <- min(c(occ_alarm_idx, mag_alarm_idx), na.rm = TRUE)
first_alarm_date <- df_te$date[first_alarm_idx]
flagged_head <- if (first_alarm_idx == mag_alarm_idx) "magnitude" else "occurrence"

cat(sprintf("First alarm: %s (head: %s)\n", first_alarm_date, flagged_head))
cat(sprintf("Refit window: last %d days ending %s, i.e. %s to %s\n",
            REFIT_L, first_alarm_date,
            first_alarm_date - REFIT_L, first_alarm_date))

# Build the sliding refit window
df_refit <- d |> filter(date > (first_alarm_date - REFIT_L), date <= first_alarm_date)
cat(sprintf("Refit window length: %d observations\n", nrow(df_refit)))

fit_refit_win <- pscl::hurdle(form, data = df_refit, dist = "negbin")

# Arm B: refit only the flagged head (magnitude in this case)
fit_B <- fit_orig
if (flagged_head == "magnitude") {
  fit_B$coefficients$count <- fit_refit_win$coefficients$count
  fit_B$theta              <- fit_refit_win$theta
} else {
  fit_B$coefficients$zero <- fit_refit_win$coefficients$zero
}

# Arm C: refit only the wrong head (control)
fit_C <- fit_orig
if (flagged_head == "magnitude") {
  fit_C$coefficients$zero <- fit_refit_win$coefficients$zero
} else {
  fit_C$coefficients$count <- fit_refit_win$coefficients$count
  fit_C$theta              <- fit_refit_win$theta
}

# Post-alarm window
df_post <- d |> filter(date > first_alarm_date)
cat(sprintf("Post-alarm window: %d observations\n", nrow(df_post)))

m_B <- monitor_full(fit_B, df_post, pit_seed = PIT_SEED + 1L)
m_C <- monitor_full(fit_C, df_post, pit_seed = PIT_SEED + 2L)

occ_term_orig <- round(m_orig$logE_occ[length(m_orig$logE_occ)], 2)
mag_term_orig <- round(m_orig$logE_mag[length(m_orig$logE_mag)], 2)
occ_term_B    <- round(m_B$logE_occ[length(m_B$logE_occ)], 2)
mag_term_B    <- round(m_B$logE_mag[length(m_B$logE_mag)], 2)
occ_term_C    <- round(m_C$logE_occ[length(m_C$logE_occ)], 2)
mag_term_C    <- round(m_C$logE_mag[length(m_C$logE_mag)], 2)

first_str <- function(logE, thr, dates) {
  i <- first_cross(logE, thr)
  if (is.na(i)) "none" else as.character(dates[i])
}

refit_tbl <- data.frame(
  Arm = c("A. Original (no refit)",
          "B. Refit flagged head (magnitude) on sliding window",
          "C. Refit control head (occurrence) on sliding window"),
  Refit_window = c("none",
                   sprintf("%s to %s", first_alarm_date - REFIT_L, first_alarm_date),
                   sprintf("%s to %s", first_alarm_date - REFIT_L, first_alarm_date)),
  Occ_first_alarm = c(
    first_str(m_orig$logE_occ, THR_HEAD, m_orig$dates),
    first_str(m_B$logE_occ,    THR_HEAD, m_B$dates),
    first_str(m_C$logE_occ,    THR_HEAD, m_C$dates)
  ),
  Mag_first_alarm = c(
    first_str(m_orig$logE_mag, THR_HEAD, m_orig$dates),
    first_str(m_B$logE_mag,    THR_HEAD, m_B$dates),
    first_str(m_C$logE_mag,    THR_HEAD, m_C$dates)
  ),
  Occ_terminal = c(occ_term_orig, occ_term_B, occ_term_C),
  Mag_terminal = c(mag_term_orig, mag_term_B, mag_term_C)
)

write_table(refit_tbl, "refit",
            caption = sprintf("Head-specific refit experiment. The baseline first alarm is the %s head on %s. Arms B and C refit one head each on a sliding window of the last %d observations ending at the alarm date, keep the other head frozen at the original 2012--2017 fit, and monitor the post-alarm window. Arm B refits the flagged head; arm C refits the wrong head as a control. If attribution matters, arm B should quiet the flagged head substantially more than arm C does.",
                              flagged_head, first_alarm_date, REFIT_L))

cat("\n=== Refit summary ===\n\n")
print(refit_tbl, row.names = FALSE)

save_raw(list(fit_orig = fit_orig, fit_B = fit_B, fit_C = fit_C,
              fit_refit_win = fit_refit_win,
              m_orig = m_orig, m_B = m_B, m_C = m_C,
              first_alarm_date = first_alarm_date,
              flagged_head = flagged_head,
              refit_tbl = refit_tbl), "refit_raw")

cat("\nTable:   out/application/tables/refit.{csv,tex}\n")