# Extra checks requested in review.
#  (1) Vuong test, hurdle-NB vs NB GLM (non-nested), AIC-corrected, per series, on the training window.
#  (2) Whole-forecast PIT test with three margins: strict (0), the paper's delta (0.03), and a
#      pi^2-scaled margin delta*pi_t^2 (the whole-forecast location signal is diluted by pi^2, so the
#      same delta handicaps it). pi_t is predictable, so validity is unchanged.
# Usage: source("application/07_extra_checks.R")
# Writes: out/application/tables/vuong.{csv,tex}, whole_margin.{csv,tex}

source("R/io.R"); set_out_root("out/application")
source("R/eprocess.R"); source("R/helpers.R")
suppressPackageStartupMessages({ library(dplyr); library(MASS); library(pscl) })
TRAIN_END <- as.Date("2017-12-31"); ALPHA <- 0.05; TOL <- c(0.02, 0.03)
THR_HEAD <- log(2 / ALPHA); THR_WHOLE <- log(1 / ALPHA)
form <- count ~ sin_year + cos_year + weekend + y_lag + y_roll7
panel <- read.csv("data/panel.csv", stringsAsFactors = FALSE)
first_str <- function(logE, thr, dates) { i <- first_cross(logE, thr); if (is.na(i)) "none" else as.character(dates[i]) }

vu <- list(); wm <- list()
for (i in seq_len(nrow(panel))) {
  id <- panel$id[i]
  d <- read.csv(sprintf("data/faers/%s.csv", id)); d$date <- as.Date(d$date)
  d <- build_features(arrange(filter(d, !is.na(count)), date))
  df_tr <- d[d$date <= TRAIN_END, ]; df_te <- d[d$date > TRAIN_END, ]
  fit_h  <- pscl::hurdle(form, data = df_tr, dist = "negbin")
  fit_nb <- MASS::glm.nb(form, data = df_tr)
  y <- df_tr$count; n <- length(y)
  P <- predict(fit_h, type = "prob")                       # n x (max+1) matrix
  ll_h  <- log(pmax(P[cbind(seq_len(n), y + 1)], 1e-300))
  ll_nb <- dnbinom(y, size = fit_nb$theta, mu = fitted(fit_nb), log = TRUE)
  if (abs(sum(ll_h) - as.numeric(logLik(fit_h))) > 1e-3) warning(id, ": hurdle loglik mismatch")
  k_h <- length(fit_h$coefficients$zero) + length(fit_h$coefficients$count) + 1
  k_n <- length(coef(fit_nb)) + 1
  m <- ll_h - ll_nb
  z <- (sum(m) - (k_h - k_n)) / (sqrt(n) * sd(m))           # AIC-corrected Vuong statistic
  vu[[length(vu) + 1]] <- data.frame(id = id, series = paste(panel$drug[i], panel$reaction[i], sep = " / "),
    vuong_z = round(z, 2), p_two_sided = round(2 * pnorm(-abs(z)), 4),
    favours = ifelse(abs(z) < 1.96, "neither", ifelse(z > 0, "hurdle-NB", "NB")))

  # whole-PIT with three margins
  X_tr <- model.matrix(form, df_tr); X_te <- model.matrix(form, df_te)
  pi_te <- hurdle_pi(fit_h, X_tr, X_te, df_tr$count)
  mu_te <- as.numeric(exp(X_te %*% fit_h$coefficients$count)); th <- as.numeric(fit_h$theta)
  U <- pnorm(rqr_hurdle_vec(df_te$count, pi_te, mu_te, th, rng_seed = 2024))
  one <- rep(1, length(U)); T_ <- length(U)
  e0 <- eprocess_pit(U, one, 0); e1 <- eprocess_pit(U, one, TOL[2]); e2 <- eprocess_pit(U, one, TOL[2] * pi_te^2)
  wm[[length(wm) + 1]] <- data.frame(id = id, pi_mean = round(mean(pi_te), 3),
    whole_strict = first_str(e0, THR_WHOLE, df_te$date), whole_delta = first_str(e1, THR_WHOLE, df_te$date),
    whole_pi2delta = first_str(e2, THR_WHOLE, df_te$date),
    rate_strict = round(tail(e0, 1) / T_, 5), rate_delta = round(tail(e1, 1) / T_, 5), rate_pi2delta = round(tail(e2, 1) / T_, 5))
}
vu <- do.call(rbind, vu); wm <- do.call(rbind, wm)
write_table(vu, "vuong", caption = "AIC-corrected Vuong test of hurdle-NB against the NB GLM on the training window (two-sided; |z|<1.96 means neither model is preferred).")
write_table(wm, "whole_margin", caption = "Whole-forecast PIT test under three margins: strict, delta=0.03, and delta scaled by the predictable pi\\_t squared. rate = terminal log E divided by monitoring days.")
print(vu, row.names = FALSE); print(wm, row.names = FALSE)
