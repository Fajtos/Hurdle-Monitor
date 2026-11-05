BETS <- c(-0.4, -0.2, -0.1, -0.05, 0.05, 0.1, 0.2, 0.4)

eprocess_occurrence <- function(pi_hat, Z, delta = 0) {
  T <- length(Z)
  features <- list(
    rep(1, T),
    2 * pi_hat - 1,
    c(0, 2 * Z[-T] - 1)
  )
  logE_matrix <- matrix(0, T, length(features) * length(BETS))
  k <- 0
  for (f in features) {
    g <- f * (Z - pi_hat)
    for (b in BETS) {
      k <- k + 1
      logE_matrix[, k] <- cumsum(log1p(b * (g - sign(b) * delta)))
    }
  }
  apply(logE_matrix, 1, log_sum_exp) - log(ncol(logE_matrix))
}

eprocess_pit <- function(U, active, delta = 0) {
  T <- length(U)
  U_last <- c(0.5, U[-T])
  features <- list(
    2 * U - 1,
    (12 * (U - 0.5)^2 - 1) / 2,
    (2 * U - 1) * (2 * U_last - 1)
  )
  logE_matrix <- matrix(0, T, length(features) * length(BETS))
  k <- 0
  for (f in features) {
    g <- active * f
    for (b in BETS) {
      k <- k + 1
      logE_matrix[, k] <- cumsum(log1p(b * (g - sign(b) * delta * active)))
    }
  }
  apply(logE_matrix, 1, log_sum_exp) - log(ncol(logE_matrix))
}

first_cross <- function(logE, threshold) {
  idx <- which(logE >= threshold)
  if (length(idx) == 0) NA_integer_ else as.integer(idx[1])
}