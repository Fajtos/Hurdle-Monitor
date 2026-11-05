dispersion <- function(y, mu, n_params) {
  sum((y - mu)^2 / mu) / (length(y) - n_params)
}

wald_and_lr <- function(fit_occ, fit_mag) {
  list(
    occ_wald = summary(fit_occ)$coefficients,
    mag_wald = summary(fit_mag)$coefficients,
    occ_null = lmtest::lrtest(fit_occ),
    mag_null = lmtest::lrtest(fit_mag)
  )
}

influence_heads <- function(fit_occ, fit_mag) {
  list(
    occ = data.frame(
      hat = hatvalues(fit_occ),
      cook = cooks.distance(fit_occ),
      resid = rstandard(fit_occ)
    ),
    mag = data.frame(
      hat = hatvalues(fit_mag),
      cook = cooks.distance(fit_mag),
      resid = rstandard(fit_mag)
    )
  )
}

simulate_rqr <- function(hurdle_fit, n_sim = 250) {
  DHARMa::simulateResiduals(hurdle_fit, n = n_sim)
}