run_one_replication <- function(scenario, seed, oracle, config) {
  df_full <- simulate_hurdle(
    T_total = config$n_train + config$T_monitor,
    n_train = config$n_train,
    shift_time = config$shift_time,
    d_occ = scenario$d_occ,
    d_mag = scenario$d_mag,
    theta = scenario$theta,
    seed = seed
  )
  df_train <- df_full[1:config$n_train, ]
  df_test <- df_full[(config$n_train + 1):(config$n_train + config$T_monitor), ]
  
  fc <- fit_and_predict(df_train, df_test, oracle = oracle)
  if (is.null(fc)) return(NULL)
  
  Z <- as.integer(df_test$y > 0)
  U_pos <- randomized_pit_positive(df_test$y, fc$mu_ut, fc$theta, seed = seed + 1L)
  U_full <- randomized_pit_full(df_test$y, fc$pi, fc$mu_ut, fc$theta, seed = seed + 2L)
  
  out <- list(scenario = scenario$name, seed = seed, oracle = oracle)
  threshold_head <- log(2 / config$alpha)
  threshold_whole <- log(1 / config$alpha)
  
  for (m in names(config$margins)) {
    margin <- config$margins[[m]]
    delta_occ <- if (length(margin) == 1) margin else margin[1]
    delta_mag <- if (length(margin) == 1) margin else margin[2]
    
    logE_occ <- eprocess_occurrence(fc$pi, Z, delta_occ)
    logE_mag <- eprocess_pit(U_pos, Z, delta_mag)
    logE_whole <- eprocess_pit(U_full, rep(1, length(U_full)), delta_mag)
    
    t_occ <- first_cross(logE_occ, threshold_head)
    t_mag <- first_cross(logE_mag, threshold_head)
    t_whole <- first_cross(logE_whole, threshold_whole)
    
    out[[m]] <- list(
      t_occ = t_occ,
      t_mag = t_mag,
      t_whole = t_whole,
      occ_flagged = !is.na(t_occ),
      mag_flagged = !is.na(t_mag),
      occ_pre_shift = !is.na(t_occ) && t_occ < config$shift_time,
      mag_pre_shift = !is.na(t_mag) && t_mag < config$shift_time,
      rate_occ = logE_occ[config$T_monitor] / config$T_monitor,
      rate_mag = logE_mag[config$T_monitor] / config$T_monitor,
      rate_whole = logE_whole[config$T_monitor] / config$T_monitor,
      pi_mean = mean(fc$pi)
    )
  }
  
  rows <- list()
  for (m in names(config$margins)) {
    r <- out[[m]]
    rows[[m]] <- data.frame(
      scenario = scenario$name,
      seed = seed,
      oracle = oracle,
      margin = m,
      occ_flagged = r$occ_flagged,
      mag_flagged = r$mag_flagged,
      occ_pre_shift = r$occ_pre_shift,
      mag_pre_shift = r$mag_pre_shift,
      t_occ = r$t_occ,
      t_mag = r$t_mag,
      t_whole = r$t_whole,
      rate_occ = r$rate_occ,
      rate_mag = r$rate_mag,
      rate_whole = r$rate_whole,
      pi_mean = r$pi_mean
    )
  }
  do.call(rbind, rows)
}