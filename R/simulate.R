simulate_hurdle <- function(T_total, n_train, shift_time = NULL,
                            d_occ = 0, d_mag = 0, theta = 5,
                            alpha = c(-1.0, 0.8, 0.5),
                            gamma = c(1.0, 0.4, 0.2),
                            rho_x = 0.7, seed = 1) {
  set.seed(seed)
  T <- T_total
  x <- numeric(T)
  eps <- rnorm(T)
  x[1] <- eps[1] / sqrt(1 - rho_x^2)
  if (T > 1) {
    for (t in 2:T) {
      x[t] <- rho_x * x[t - 1] + sqrt(1 - rho_x^2) * eps[t]
    }
  }
  y <- integer(T)
  shift_idx <- if (is.null(shift_time)) T + 1 else n_train + shift_time
  
  for (t in seq_len(T)) {
    y_lag <- if (t == 1) 0 else log1p(y[t - 1])
    d_occ_t <- if (t >= shift_idx) d_occ else 0
    d_mag_t <- if (t >= shift_idx) d_mag else 0
    pi_t <- plogis(alpha[1] + alpha[2] * x[t] + alpha[3] * y_lag + d_occ_t)
    if (rbinom(1, 1, pi_t) == 1) {
      mu_t <- exp(gamma[1] + gamma[2] * x[t] + gamma[3] * y_lag + d_mag_t)
      v <- 0L
      while (v == 0L) v <- rnbinom(1, size = theta, mu = mu_t)
      y[t] <- v
    }
  }
  data.frame(t = seq_len(T), x = x, y = y, y_lag = c(0, log1p(y[-T])))
}