# Case-study analysis for the metformin / LACTIC ACIDOSIS series (r2).
# Produces: four-model ladder (OLS, Poisson, NB, hurdle-NB), coefficients,
#           classical diagnostics (RQR, envelope, Cook's distance), and the
#           two-head monitor (occurrence, magnitude, whole-forecast PIT).

source("R/io.R")
set_out_root("out/application")
source("R/eprocess.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(MASS)
  library(pscl)
})

CASE_ID   <- "r2"
TRAIN_END <- as.Date("2017-12-31")
ALPHA     <- 0.05
TOL       <- c(0.02, 0.03)
THR_HEAD  <- log(2 / ALPHA)
THR_WHOLE <- log(1 / ALPHA)

d <- read.csv(sprintf("data/faers/%s.csv", CASE_ID), stringsAsFactors = FALSE)
d$date <- as.Date(d$date)
d <- d |> filter(!is.na(count)) |> arrange(date)

n <- nrow(d)
doy <- as.integer(format(d$date, "%j"))
d$sin_year <- sin(2 * pi * doy / 365.25)
d$cos_year <- cos(2 * pi * doy / 365.25)
d$weekend  <- as.integer(format(d$date, "%u") %in% c("6", "7"))
d$y_lag    <- c(0, log1p(d$count[-n]))
roll7 <- as.numeric(stats::filter(log1p(d$count), rep(1/7, 7), sides = 1))
roll7[is.na(roll7)] <- 0
d$y_roll7 <- roll7

n0    <- sum(d$date <= TRAIN_END)
df_tr <- d[1:n0, ]
df_te <- d[(n0 + 1):n, ]

cat(sprintf("Case-study series: %s\n", CASE_ID))
cat(sprintf("Training observations:   %d (2012-01-01 to 2017-12-31)\n", nrow(df_tr)))
cat(sprintf("Monitoring observations: %d (2018-01-01 to 2024-12-31)\n", nrow(df_te)))
cat(sprintf("Total report count:      %d\n\n", sum(d$count)))

form <- count ~ sin_year + cos_year + weekend + y_lag + y_roll7

fit_ols  <- lm(form, data = df_tr)
fit_pois <- glm(form, data = df_tr, family = poisson())
fit_nb   <- tryCatch(MASS::glm.nb(form, data = df_tr),
                     error = function(e) { message("NB failed: ", conditionMessage(e)); NULL })
fit_hur  <- tryCatch(pscl::hurdle(form, data = df_tr, dist = "negbin"),
                     error = function(e) { message("Hurdle failed: ", conditionMessage(e)); NULL })

if (is.null(fit_hur)) stop("Hurdle fit failed; cannot proceed.")

ll_ols  <- as.numeric(logLik(fit_ols))
ll_pois <- as.numeric(logLik(fit_pois))
ll_nb   <- if (!is.null(fit_nb)) as.numeric(logLik(fit_nb)) else NA_real_
ll_hur  <- as.numeric(logLik(fit_hur))

k_ols  <- length(coef(fit_ols))
k_pois <- length(coef(fit_pois))
k_nb   <- if (!is.null(fit_nb)) length(coef(fit_nb)) + 1 else NA_integer_
k_hur  <- length(fit_hur$coefficients$zero) + length(fit_hur$coefficients$count) + 1

aic_ols  <- -2 * ll_ols  + 2 * k_ols
aic_pois <- -2 * ll_pois + 2 * k_pois
aic_nb   <- if (!is.null(fit_nb)) -2 * ll_nb + 2 * k_nb else NA_real_
aic_hur  <- -2 * ll_hur  + 2 * k_hur

pred_ols_te  <- predict(fit_ols, newdata = df_te)
pred_pois_te <- predict(fit_pois, newdata = df_te, type = "response")
pred_nb_te   <- if (!is.null(fit_nb)) predict(fit_nb, newdata = df_te, type = "response") else rep(NA_real_, nrow(df_te))

X_te <- model.matrix(form, df_te)
X_tr <- model.matrix(form, df_tr)

coef_zero  <- fit_hur$coefficients$zero
coef_count <- fit_hur$coefficients$count
theta_hur  <- as.numeric(fit_hur$theta)

pi_raw_tr <- as.numeric(plogis(X_tr %*% coef_zero))
emp_pi_pos_tr <- mean(df_tr$count > 0)
if (abs(mean(pi_raw_tr) - emp_pi_pos_tr) < abs(mean(1 - pi_raw_tr) - emp_pi_pos_tr)) {
  pi_pos_te <- as.numeric(plogis(X_te %*% coef_zero))
} else {
  pi_pos_te <- 1 - as.numeric(plogis(X_te %*% coef_zero))
}
mu_ut_te <- as.numeric(exp(X_te %*% coef_count))
q0_te    <- (theta_hur / (theta_hur + mu_ut_te))^theta_hur
mpos_te  <- mu_ut_te / (1 - q0_te)

cat(sprintf("Model matrix columns: %d, occurrence coefficients: %d, count coefficients: %d\n\n",
            ncol(X_te), length(coef_zero), length(coef_count)))

sigma_ols <- summary(fit_ols)$sigma
theta_nb  <- if (!is.null(fit_nb)) fit_nb$theta else NA_real_

log_score_ols  <- -mean(dnorm(df_te$count, pred_ols_te, sigma_ols, log = TRUE))
log_score_pois <- -mean(dpois(df_te$count, pred_pois_te, log = TRUE))
log_score_nb   <- if (!is.null(fit_nb)) {
  -mean(dnbinom(df_te$count, size = theta_nb, mu = pred_nb_te, log = TRUE))
} else NA_real_

ll_hur_te <- ifelse(
  df_te$count == 0,
  log(1 - pi_pos_te),
  log(pi_pos_te) + log(dnbinom(df_te$count, size = theta_hur, mu = mu_ut_te) / (1 - q0_te))
)
log_score_hur <- -mean(ll_hur_te)

disp_pois <- sum((df_tr$count - fitted(fit_pois))^2 / fitted(fit_pois)) /
  (nrow(df_tr) - length(coef(fit_pois)))
disp_nb <- if (!is.null(fit_nb)) {
  mu_nb_tr <- fitted(fit_nb)
  sum((df_tr$count - mu_nb_tr)^2 / (mu_nb_tr + mu_nb_tr^2 / fit_nb$theta)) /
    (nrow(df_tr) - length(coef(fit_nb)) - 1)
} else NA_real_

ladder <- data.frame(
  Model = c("OLS (Gaussian)", "Poisson GLM", "NB GLM", "Hurdle-NB (proposed)"),
  Params = c(k_ols, k_pois, k_nb, k_hur),
  LogLik = round(c(ll_ols, ll_pois, ll_nb, ll_hur), 1),
  AIC    = round(c(aic_ols, aic_pois, aic_nb, aic_hur), 1),
  Pearson_dispersion = c(NA, round(disp_pois, 2), round(disp_nb, 2), NA),
  OOS_log_score = round(c(log_score_ols, log_score_pois, log_score_nb, log_score_hur), 3)
)
write_table(ladder, "ladder",
            caption = "Model ladder on the training window. In-sample log-likelihood, AIC and Pearson dispersion; out-of-sample mean negative log score on the monitoring window.")
cat("Ladder:\n"); print(ladder); cat("\n")

weekday_tr <- format(df_tr$date, "%u") %in% c("1","2","3","4","5")
zero_share_all   <- mean(df_tr$count == 0)
zero_share_wd    <- mean(df_tr$count[weekday_tr] == 0)
var_to_mean      <- var(df_tr$count) / mean(df_tr$count)
pos_tr           <- df_tr$count[df_tr$count > 0]
acf1             <- acf(df_tr$count, lag.max = 1, plot = FALSE)$acf[2]
acf7             <- acf(df_tr$count, lag.max = 7, plot = FALSE)$acf[8]
weekend_mean     <- mean(df_tr$count[df_tr$weekend == 1])
weekday_mean     <- mean(df_tr$count[df_tr$weekend == 0])
runs             <- rle(as.integer(df_tr$count == 0))
zero_runs        <- runs$lengths[runs$values == 1]
mean_zero_run    <- if (length(zero_runs) > 0) mean(zero_runs) else 0
max_zero_run     <- if (length(zero_runs) > 0) max(zero_runs) else 0

eda <- data.frame(
  Statistic = c(
    "Training observations",
    "Monitoring observations",
    "Overall zero share (training)",
    "Weekday zero share (training)",
    "Variance-to-mean ratio (training)",
    "Mean of positive counts (training)",
    "Maximum count (training)",
    "Autocorrelation, lag 1",
    "Autocorrelation, lag 7",
    "Weekend mean count",
    "Weekday mean count",
    "Mean zero-run length",
    "Maximum zero-run length"
  ),
  Value = c(
    nrow(df_tr),
    nrow(df_te),
    sprintf("%.3f", zero_share_all),
    sprintf("%.3f", zero_share_wd),
    sprintf("%.2f", var_to_mean),
    sprintf("%.2f", mean(pos_tr)),
    max(df_tr$count),
    sprintf("%.3f", acf1),
    sprintf("%.3f", acf7),
    sprintf("%.2f", weekend_mean),
    sprintf("%.2f", weekday_mean),
    sprintf("%.1f", mean_zero_run),
    max_zero_run
  )
)
write_table(eda, "eda",
            caption = "Exploratory statistics, case-study series.")
cat("EDA:\n"); print(eda); cat("\n")

sum_zero  <- summary(fit_hur)$coefficients$zero
sum_count_full <- summary(fit_hur)$coefficients$count
n_count_coef <- length(fit_hur$coefficients$count)
sum_count <- sum_count_full[seq_len(n_count_coef), , drop = FALSE]

coefs <- data.frame(
  Term      = rownames(sum_zero),
  Occ_est   = round(sum_zero[, 1], 3),
  Occ_se    = round(sum_zero[, 2], 3),
  Occ_p     = round(sum_zero[, 4], 4),
  Occ_or    = round(exp(sum_zero[, 1]), 3),
  Mag_est   = round(sum_count[, 1], 3),
  Mag_se    = round(sum_count[, 2], 3),
  Mag_p     = round(sum_count[, 4], 4),
  Mag_mr    = round(exp(sum_count[, 1]), 3)
)
write_table(coefs, "coefs",
            caption = "Hurdle-NB coefficients, both heads. Odds ratios use $e^{\\hat\\alpha_j}$; mean ratios use $e^{\\hat\\gamma_j}$.")
cat("Coefficients:\n"); print(coefs); cat("\n")

rqr_hurdle_vec <- function(y, pi_pos, mu_ut, theta, rng_seed = NULL) {
  if (!is.null(rng_seed)) set.seed(rng_seed)
  n <- length(y)
  q0 <- (theta / (theta + mu_ut))^theta
  U <- numeric(n)
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

simulate_hurdle_obs <- function(pi_pos, mu_ut, theta) {
  n <- length(pi_pos)
  Z <- rbinom(n, 1, pi_pos)
  y <- integer(n)
  for (i in which(Z == 1)) {
    v <- 0L
    while (v == 0L) v <- rnbinom(1, size = theta, mu = mu_ut[i])
    y[i] <- v
  }
  y
}

set.seed(2024)
rqr_obs <- rqr_hurdle_vec(df_te$count, pi_pos_te, mu_ut_te, theta_hur)

B <- 250
set.seed(12345)
rqr_sim_mat <- matrix(NA_real_, nrow = nrow(df_te), ncol = B)
for (b in seq_len(B)) {
  y_sim <- simulate_hurdle_obs(pi_pos_te, mu_ut_te, theta_hur)
  rqr_sim_mat[, b] <- rqr_hurdle_vec(y_sim, pi_pos_te, mu_ut_te, theta_hur)
}

fitted_mean_te <- pi_pos_te * mpos_te
ord      <- order(fitted_mean_te)
env_lo   <- apply(rqr_sim_mat[ord, ], 1, quantile, probs = 0.025, na.rm = TRUE)
env_hi   <- apply(rqr_sim_mat[ord, ], 1, quantile, probs = 0.975, na.rm = TRUE)
inside   <- mean(rqr_obs[ord] >= env_lo & rqr_obs[ord] <= env_hi, na.rm = TRUE)

diag_tbl <- data.frame(
  Statistic = c(
    "RQR mean",
    "RQR variance",
    "RQR lag-1 autocorrelation",
    "Fraction inside 95% envelope",
    "Observed zeros in monitoring window",
    "Expected zeros in monitoring window"
  ),
  Value = c(
    sprintf("%.3f", mean(rqr_obs, na.rm = TRUE)),
    sprintf("%.3f", var(rqr_obs, na.rm = TRUE)),
    sprintf("%.3f", acf(rqr_obs[!is.na(rqr_obs)], lag.max = 1, plot = FALSE)$acf[2]),
    sprintf("%.3f", inside),
    as.character(sum(df_te$count == 0)),
    sprintf("%.1f", sum(1 - pi_pos_te))
  )
)
write_table(diag_tbl, "diag",
            caption = "Randomized quantile residual diagnostics, monitoring window. Under a correct fit, RQR has mean $0$ and variance $1$, and approximately $95\\%$ of points fall inside the envelope.")
cat("Diagnostics:\n"); print(diag_tbl); cat("\n")

Z_te <- as.integer(df_te$count > 0)

U_pos_te <- rep(0.5, nrow(df_te))
pos_idx <- which(df_te$count > 0)
if (length(pos_idx) > 0) {
  set.seed(2025)
  U_pos_te[pos_idx] <- {
    F_y   <- pnbinom(df_te$count[pos_idx],     size = theta_hur, mu = mu_ut_te[pos_idx])
    F_ym1 <- pnbinom(df_te$count[pos_idx] - 1, size = theta_hur, mu = mu_ut_te[pos_idx])
    F_0   <- pnbinom(0,                        size = theta_hur, mu = mu_ut_te[pos_idx])
    Fy_tr    <- (F_y   - F_0) / (1 - F_0)
    Fym1_tr  <- (F_ym1 - F_0) / (1 - F_0)
    V <- runif(length(pos_idx))
    Fym1_tr + V * (Fy_tr - Fym1_tr)
  }
}

U_full_te <- pnorm(rqr_obs)

logE_occ   <- eprocess_occurrence(pi_pos_te, Z_te, TOL[1])
logE_mag   <- eprocess_pit(U_pos_te,   Z_te, TOL[2])
logE_whole <- eprocess_pit(U_full_te,  rep(1, length(U_full_te)), TOL[2])

first_str <- function(logE, thr, dates) {
  i <- first_cross(logE, thr)
  if (is.na(i)) "none" else as.character(dates[i])
}

monitor_tbl <- data.frame(
  Head = c("Occurrence", "Magnitude", "Whole-forecast PIT"),
  First_alarm = c(
    first_str(logE_occ,   THR_HEAD,  df_te$date),
    first_str(logE_mag,   THR_HEAD,  df_te$date),
    first_str(logE_whole, THR_WHOLE, df_te$date)
  ),
  Terminal_log_E = round(c(
    logE_occ[length(logE_occ)],
    logE_mag[length(logE_mag)],
    logE_whole[length(logE_whole)]
  ), 2)
)
write_table(monitor_tbl, "monitor",
            caption = "First-alarm dates and terminal log e-process, case-study series, tolerant margin.")
cat("Monitor:\n"); print(monitor_tbl); cat("\n")

monthly <- d |>
  mutate(ym = as.Date(paste0(format(date, "%Y-%m"), "-01"))) |>
  group_by(ym) |>
  summarise(mean_count = mean(count), .groups = "drop")
p_series <- ggplot(monthly, aes(x = ym, y = mean_count)) +
  geom_line(colour = "steelblue", linewidth = 0.7) +
  geom_vline(xintercept = TRAIN_END, linetype = "dashed", colour = "grey40") +
  labs(x = NULL, y = "Monthly mean daily count") +
  theme_minimal()
save_figure(p_series, "series", width = 9, height = 4)

rqr_df <- data.frame(fitted = fitted_mean_te[ord], rqr = rqr_obs[ord],
                     env_lo = env_lo, env_hi = env_hi)
p_rqr <- ggplot(rqr_df, aes(x = fitted)) +
  geom_ribbon(aes(ymin = env_lo, ymax = env_hi), fill = "grey80", alpha = 0.5) +
  geom_point(aes(y = rqr), size = 0.6, alpha = 0.5, colour = "steelblue") +
  geom_hline(yintercept = c(-1.96, 0, 1.96),
             linetype = c("dashed", "solid", "dashed")) +
  labs(x = "Fitted mean", y = "Randomized quantile residual") +
  theme_minimal()
save_figure(p_rqr, "rqr_envelope", width = 8, height = 4.5)

Z_tr <- as.integer(df_tr$count > 0)
df_tr_z <- df_tr
df_tr_z$Z <- Z_tr
fit_occ_glm <- glm(update(form, Z ~ .), data = df_tr_z, family = binomial())
inf_df <- data.frame(
  index = seq_len(nrow(df_tr)),
  cook  = cooks.distance(fit_occ_glm)
)
p_inf <- ggplot(inf_df, aes(x = index, y = cook)) +
  geom_point(size = 0.5, colour = "steelblue", alpha = 0.7) +
  labs(x = "Training observation index", y = "Cook's distance (occurrence head)") +
  theme_minimal()
save_figure(p_inf, "influence", width = 9, height = 4)

monitor_df <- data.frame(
  date      = df_te$date,
  logE_occ  = logE_occ,
  logE_mag  = logE_mag,
  logE_whole = logE_whole
)
p_mon <- ggplot(monitor_df, aes(x = date)) +
  geom_line(aes(y = logE_occ,   colour = "Occurrence"), linewidth = 0.7) +
  geom_line(aes(y = logE_mag,   colour = "Magnitude"),  linewidth = 0.7) +
  geom_line(aes(y = logE_whole, colour = "Whole PIT"),  linewidth = 0.7, linetype = "dashed") +
  geom_hline(yintercept = THR_HEAD, linetype = "dotted") +
  scale_colour_manual(values = c("Occurrence" = "red",
                                 "Magnitude"  = "blue",
                                 "Whole PIT"  = "grey40")) +
  labs(x = NULL, y = expression(log~E[t]), colour = NULL) +
  theme_minimal() + theme(legend.position = "top")
save_figure(p_mon, "monitor", width = 9, height = 4)

save_raw(list(
  d = d, df_tr = df_tr, df_te = df_te,
  fit_ols = fit_ols, fit_pois = fit_pois, fit_nb = fit_nb, fit_hur = fit_hur,
  pi_pos_te = pi_pos_te, mu_ut_te = mu_ut_te, theta_hur = theta_hur,
  logE_occ = logE_occ, logE_mag = logE_mag, logE_whole = logE_whole
), sprintf("case_%s", CASE_ID))
