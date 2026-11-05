# Shared helpers (added in patch 2026-11-05)

# Log score of the discretised Gaussian: P(Y=k) = Phi((k+.5-mu)/s) - Phi((k-.5-mu)/s).
# Computed in the upper tail when lo > 0 to avoid the cancellation that
# made pnorm(hi) - pnorm(lo) underflow to 0 for large counts.
ols_disc_logscore <- function(y, mu, s) {
  lo <- (y - 0.5 - mu) / s
  hi <- (y + 0.5 - mu) / s
  lsf_lo <- pnorm(lo, lower.tail = FALSE, log.p = TRUE)
  lsf_hi <- pnorm(hi, lower.tail = FALSE, log.p = TRUE)
  lp_up  <- lsf_lo + log1p(-exp(lsf_hi - lsf_lo))
  lp_lo  <- log(pmax(pnorm(hi) - pnorm(lo), 1e-300))
  lp <- ifelse(lo > 0, lp_up, lp_lo)
  -mean(lp)
}

# Per-observation hurdle-NB log probabilities given pi = P(Y>0), untruncated mean mu, theta
hurdle_logpmf <- function(y, pi_pos, mu_ut, theta) {
  q0 <- (theta / (theta + mu_ut))^theta
  ifelse(y == 0, log1p(-pi_pos),
         log(pi_pos) + dnbinom(y, size = theta, mu = mu_ut, log = TRUE) - log1p(-q0))
}

build_features <- function(d) {
  n <- nrow(d)
  doy <- as.integer(format(d$date, "%j"))
  d$sin_year <- sin(2 * pi * doy / 365.25)
  d$cos_year <- cos(2 * pi * doy / 365.25)
  d$weekend  <- as.integer(format(d$date, "%u") %in% c("6", "7"))
  d$y_lag    <- c(0, log1p(d$count[-n]))
  roll7 <- as.numeric(stats::filter(d$count, rep(1/7, 7), sides = 1))
  roll7[is.na(roll7)] <- 0
  d$y_roll7 <- log1p(c(0, roll7[-n]))   # uses days t-7..t-1 only (no leak)
  d
}

rqr_hurdle_vec <- function(y, pi_pos, mu_ut, theta, rng_seed = NULL) {
  if (!is.null(rng_seed)) set.seed(rng_seed)
  n <- length(y); q0 <- (theta / (theta + mu_ut))^theta; U <- numeric(n)
  pos <- which(y > 0)
  if (length(pos) > 0) {
    Fy   <- (pnbinom(y[pos],     size = theta, mu = mu_ut[pos]) - q0[pos]) / (1 - q0[pos])
    Fym1 <- (pnbinom(y[pos] - 1, size = theta, mu = mu_ut[pos]) - q0[pos]) / (1 - q0[pos])
    Ffy   <- (1 - pi_pos[pos]) + pi_pos[pos] * Fy
    Ffym1 <- (1 - pi_pos[pos]) + pi_pos[pos] * Fym1
    U[pos] <- Ffym1 + runif(length(pos)) * (Ffy - Ffym1)
  }
  zero <- which(y == 0)
  if (length(zero) > 0) U[zero] <- runif(length(zero)) * (1 - pi_pos[zero])
  qnorm(U)
}

# Same sign convention check used in 02/03: pscl's zero part models P(Y>0)
hurdle_pi <- function(fit_hur, X_tr, X_te, y_tr) {
  cz <- fit_hur$coefficients$zero
  pr <- as.numeric(plogis(X_tr %*% cz)); emp <- mean(y_tr > 0)
  if (abs(mean(pr) - emp) < abs(mean(1 - pr) - emp)) as.numeric(plogis(X_te %*% cz))
  else 1 - as.numeric(plogis(X_te %*% cz))
}
