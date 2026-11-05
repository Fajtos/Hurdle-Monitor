fit_and_predict <- function(df_train, df_test,
                            oracle = FALSE,
                            alpha_true = c(-1.0, 0.8, 0.5),
                            gamma_true = c(1.0, 0.4, 0.2),
                            theta_true = 5,
                            formula_occ = y ~ x + y_lag,
                            formula_mag = y ~ x + y_lag) {
  df_train_z <- df_train
  df_train_z$Z <- as.integer(df_train$y > 0)
  df_test_z <- df_test
  df_test_z$Z <- as.integer(df_test$y > 0)
  
  if (oracle) {
    pi_hat <- plogis(alpha_true[1] + alpha_true[2] * df_test$x +
                       alpha_true[3] * df_test$y_lag)
    mu_ut <- exp(gamma_true[1] + gamma_true[2] * df_test$x +
                   gamma_true[3] * df_test$y_lag)
    theta_hat <- theta_true
  } else {
    occ_formula <- update(formula_occ, Z ~ .)
    fit_occ <- tryCatch(
      glm(occ_formula, data = df_train_z, family = binomial()),
      error = function(e) {
        message("glm occurrence failed: ", conditionMessage(e))
        NULL
      }
    )
    fit_mag <- tryCatch(
      pscl::hurdle(formula_mag, data = df_train, dist = "negbin"),
      error = function(e) {
        message("pscl::hurdle failed: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(fit_occ) || is.null(fit_mag)) return(NULL)
    pi_hat <- as.numeric(predict(fit_occ, newdata = df_test_z, type = "response"))
    X_mag <- model.matrix(formula_mag, data = df_test)
    coef_mag <- fit_mag$coefficients$count
    mu_ut <- as.numeric(exp(X_mag %*% coef_mag))
    theta_hat <- as.numeric(fit_mag$theta)
  }
  q0 <- (theta_hat / (theta_hat + mu_ut)) ^ theta_hat
  m_pos <- mu_ut / (1 - q0)
  list(pi = pi_hat, mu_ut = mu_ut, m_pos = m_pos, theta = theta_hat)
}

randomized_pit_positive <- function(y, mu_ut, theta, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(y)
  U <- rep(0.5, n)
  pos <- which(y > 0)
  if (length(pos) > 0) {
    F_y <- pnbinom(y[pos], size = theta, mu = mu_ut[pos])
    F_ym1 <- pnbinom(y[pos] - 1, size = theta, mu = mu_ut[pos])
    F_0 <- pnbinom(0, size = theta, mu = mu_ut[pos])
    Fy_tr <- (F_y - F_0) / (1 - F_0)
    Fym1_tr <- (F_ym1 - F_0) / (1 - F_0)
    V <- runif(length(pos))
    U[pos] <- Fym1_tr + V * (Fy_tr - Fym1_tr)
  }
  U
}

randomized_pit_full <- function(y, pi_hat, mu_ut, theta, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(y)
  U <- numeric(n)
  for (i in seq_len(n)) {
    if (y[i] == 0) {
      U[i] <- runif(1, 0, 1 - pi_hat[i])
    } else {
      F_0 <- 1 - pi_hat[i]
      F_y <- F_0 + pi_hat[i] * pnbinom(y[i], size = theta, mu = mu_ut[i])
      F_ym1 <- F_0 + pi_hat[i] * pnbinom(y[i] - 1, size = theta, mu = mu_ut[i])
      U[i] <- F_ym1 + runif(1) * (F_y - F_ym1)
    }
  }
  U
}