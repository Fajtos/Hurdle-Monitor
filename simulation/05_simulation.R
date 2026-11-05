suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(ggplot2)
  library(tidyr)
})

source("R/io.R")
set_out_root("out/simulation")
source("R/simulate.R")
source("R/forecast.R")
source("R/eprocess.R")
source("R/replicate.R")

N_REP <- 400
ALPHA <- 0.05
PROGRESS_EVERY <- 200

CONFIG <- list(
  n_train = 2000,
  T_monitor = 3000,
  shift_time = 1000,
  alpha = ALPHA,
  margins = list(strict = 0, tolerant = c(0.02, 0.03))
)

SCENARIOS <- list(
  list(name = "null", d_occ = 0, d_mag = 0, theta = 5),
  list(name = "occ_pos", d_occ = 0.6, d_mag = 0, theta = 5),
  list(name = "occ_neg", d_occ = -0.6, d_mag = 0, theta = 5),
  list(name = "mag_moderate", d_occ = 0, d_mag = 0.22, theta = 5),
  list(name = "mag_large", d_occ = 0, d_mag = 0.50, theta = 5),
  list(name = "dispersion", d_occ = 0, d_mag = 0, theta = 2.5),
  list(name = "both", d_occ = 0.6, d_mag = 0.22, theta = 5)
)

counter <- new.env()
counter$n <- 0
start_time <- Sys.time()

tick <- function(total) {
  counter$n <- counter$n + 1
  if (counter$n %% PROGRESS_EVERY == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    rate <- counter$n / max(elapsed, 0.01)
    eta <- (total - counter$n) / rate
    cat(sprintf("Progress: %d / %d (%.1f%%)  elapsed %.1f min  ETA %.1f min\n",
                counter$n, total, 100 * counter$n / total, elapsed, eta))
  }
}

main_rds <- "out/simulation/raw/simulation_raw.rds"

if (file.exists(main_rds)) {
  cat("Found", main_rds, "- loading main simulation, skipping replications.\n")
  results <- readRDS(main_rds)
  cat(sprintf("Loaded %d main simulation rows.\n\n", nrow(results)))
} else {
  cat("=== Main simulation ===\n")
  n_main <- length(SCENARIOS) * N_REP * 2
  cat(sprintf("Running %d main simulation jobs.\n", n_main))
  jobs <- expand.grid(
    scenario_idx = seq_along(SCENARIOS),
    seed = seq_len(N_REP) * 1000 + 7,
    oracle = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  results <- purrr::pmap_dfr(
    list(jobs$scenario_idx, jobs$seed, jobs$oracle),
    function(idx, sd, orc) {
      tick(n_main)
      tryCatch(
        run_one_replication(SCENARIOS[[idx]], sd, orc, CONFIG),
        error = function(e) {
          message("Replication failed: scenario=", SCENARIOS[[idx]]$name,
                  " seed=", sd, " oracle=", orc, " -> ", conditionMessage(e))
          NULL
        }
      )
    }
  )
  if (is.null(results) || nrow(results) == 0) {
    stop("All main replications failed.")
  }
  save_raw(results, "simulation_raw")
  save_raw(CONFIG, "simulation_config")
  save_raw(SCENARIOS, "simulation_scenarios")
  cat("\n")
}

results$forecaster <- ifelse(results$oracle, "Oracle", "Estimated")

null_rows <- results %>% filter(scenario == "null")

sim_null <- null_rows %>%
  group_by(forecaster, margin) %>%
  summarise(
    occ_false_alarm_rate = round(mean(occ_flagged), 4),
    mag_false_alarm_rate = round(mean(mag_flagged), 4),
    occ_early_rate = round(mean(occ_pre_shift), 4),
    mag_early_rate = round(mean(mag_pre_shift), 4),
    n = n(),
    .groups = "drop"
  )
write_table(sim_null, "sim_null",
            caption = "Null calibration. Per-head false-alarm rates on a null stream.")

alt_rows <- results %>%
  filter(scenario != "null") %>%
  mutate(
    scenario_label = recode(scenario,
                            occ_pos = "Occurrence positive",
                            occ_neg = "Occurrence negative",
                            mag_moderate = "Magnitude moderate",
                            mag_large = "Magnitude large",
                            dispersion = "Dispersion",
                            both = "Both shifted"
    )
  )

sim_detection <- alt_rows %>%
  filter(forecaster == "Estimated", margin == "tolerant") %>%
  group_by(scenario_label, scenario) %>%
  summarise(
    occ_detect_rate = round(mean(occ_flagged), 3),
    mag_detect_rate = round(mean(mag_flagged), 3),
    occ_delay_median = safe_median_delay(t_occ, occ_flagged, CONFIG$shift_time),
    mag_delay_median = safe_median_delay(t_mag, mag_flagged, CONFIG$shift_time),
    localization_correct = round(mean(
      (scenario %in% c("occ_pos", "occ_neg") & occ_flagged & !mag_flagged) |
        (scenario %in% c("mag_moderate", "mag_large") & mag_flagged & !occ_flagged) |
        (scenario == "both" & occ_flagged & mag_flagged) |
        (scenario == "dispersion" & !occ_flagged & !mag_flagged)
    ), 3),
    n = n(),
    .groups = "drop"
  ) %>%
  dplyr::select(-scenario)
write_table(sim_detection, "sim_detection",
            caption = "Detection and localization under the tolerant margin, estimated forecaster.")

cat("\n=== Tolerance grid ===\n")
counter$n <- 0
start_time <- Sys.time()
n_tolerance <- 5 * 3 * 100
cat(sprintf("Running %d tolerance grid jobs.\n", n_tolerance))

delta_grid <- c(0, 0.01, 0.02, 0.03, 0.05)
tolerance_rows <- list()

for (d in delta_grid) {
  config_d <- CONFIG
  config_d$margins <- list(single = c(d, d))
  seeds <- seq_len(100) * 1000 + 11
  for (scen in c("null", "occ_pos", "mag_large")) {
    scenario <- SCENARIOS[[which(sapply(SCENARIOS, `[[`, "name") == scen)]]
    reps <- purrr::map_dfr(seeds, function(sd) {
      tick(n_tolerance)
      tryCatch(run_one_replication(scenario, sd, FALSE, config_d),
               error = function(e) NULL)
    })
    if (!is.null(reps) && nrow(reps) > 0) {
      is_null_scen <- scen == "null"
      tolerance_rows[[length(tolerance_rows) + 1]] <- data.frame(
        delta = d,
        scenario = scen,
        occ_rate = round(mean(reps$occ_flagged), 3),
        mag_rate = round(mean(reps$mag_flagged), 3),
        occ_delay_median = if (is_null_scen) NA_integer_ else
          safe_median_delay(reps$t_occ, reps$occ_flagged, CONFIG$shift_time),
        mag_delay_median = if (is_null_scen) NA_integer_ else
          safe_median_delay(reps$t_mag, reps$mag_flagged, CONFIG$shift_time)
      )
    }
  }
}

sim_tolerance <- if (length(tolerance_rows) > 0) {
  do.call(rbind, tolerance_rows)
} else {
  data.frame(
    delta = numeric(0), scenario = character(0),
    occ_rate = numeric(0), mag_rate = numeric(0),
    occ_delay_median = integer(0), mag_delay_median = integer(0)
  )
}
write_table(sim_tolerance, "sim_tolerance",
            caption = "Sensitivity to the tolerance margin.")

cat("\n=== Dependence grid ===\n")
counter$n <- 0
start_time <- Sys.time()
n_dependence <- 2 * 2 * 200
cat(sprintf("Running %d dependence grid jobs.\n", n_dependence))

dep_rows <- list()

for (label in c("with_lag", "without_lag")) {
  formula_use <- if (label == "with_lag") y ~ x + y_lag else y ~ x
  for (scen in c("null", "occ_pos")) {
    scenario <- SCENARIOS[[which(sapply(SCENARIOS, `[[`, "name") == scen)]]
    config_d <- CONFIG
    config_d$margins <- list(tolerant = c(0.02, 0.03))
    reps <- purrr::map_dfr(seq_len(200) * 1000 + 17, function(sd) {
      tick(n_dependence)
      tryCatch({
        df_full <- simulate_hurdle(
          config_d$n_train + config_d$T_monitor,
          config_d$n_train, config_d$shift_time,
          d_occ = scenario$d_occ, d_mag = scenario$d_mag,
          theta = scenario$theta, seed = sd)
        df_train <- df_full[1:config_d$n_train, ]
        df_test <- df_full[(config_d$n_train + 1):(config_d$n_train + config_d$T_monitor), ]
        
        df_train_z <- df_train
        df_train_z$Z <- as.integer(df_train$y > 0)
        df_test_z <- df_test
        df_test_z$Z <- as.integer(df_test$y > 0)
        
        fit_occ <- glm(update(formula_use, Z ~ .),
                       data = df_train_z, family = binomial())
        fit_mag <- pscl::hurdle(formula_use, data = df_train, dist = "negbin")
        
        pi_hat <- as.numeric(predict(fit_occ, newdata = df_test_z, type = "response"))
        X_mag <- model.matrix(formula_use, data = df_test)
        mu_ut <- as.numeric(exp(X_mag %*% fit_mag$coefficients$count))
        theta_hat <- fit_mag$theta
        Z <- as.integer(df_test$y > 0)
        U_pos <- randomized_pit_positive(df_test$y, mu_ut, theta_hat, seed = sd + 1L)
        logE_occ <- eprocess_occurrence(pi_hat, Z, 0.02)
        logE_mag <- eprocess_pit(U_pos, Z, 0.03)
        thr <- log(2 / ALPHA)
        data.frame(
          forecaster = label,
          scenario = scen,
          occ_flagged = !is.na(first_cross(logE_occ, thr)),
          mag_flagged = !is.na(first_cross(logE_mag, thr))
        )
      }, error = function(e) NULL)
    })
    if (!is.null(reps) && nrow(reps) > 0) dep_rows[[length(dep_rows) + 1]] <- reps
  }
}

sim_dependence <- if (length(dep_rows) > 0) {
  do.call(rbind, dep_rows) %>%
    group_by(forecaster, scenario) %>%
    summarise(
      occ_rate = round(mean(occ_flagged), 3),
      mag_rate = round(mean(mag_flagged), 3),
      n = n(),
      .groups = "drop"
    )
} else {
  data.frame(
    forecaster = character(0), scenario = character(0),
    occ_rate = numeric(0), mag_rate = numeric(0), n = integer(0)
  )
}
write_table(sim_dependence, "sim_dependence",
            caption = "Dependence calibration. Forecaster with the serial covariate included versus omitted.")

cat("\n=== Figures ===\n")

plot_null <- results %>%
  filter(scenario == "null", margin == "tolerant") %>%
  transmute(forecaster = forecaster, head = "Occurrence", rate = rate_occ) %>%
  bind_rows(
    results %>%
      filter(scenario == "null", margin == "tolerant") %>%
      transmute(forecaster = forecaster, head = "Magnitude", rate = rate_mag)
  )

fig_null <- ggplot(plot_null, aes(x = rate, fill = head)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ forecaster) +
  labs(x = "log E_T / T under the null", y = "Count", fill = NULL,
       title = "Growth rate under the null") +
  theme_minimal() + theme(legend.position = "top")
save_figure(fig_null, "sim_null_histogram", width = 8, height = 4)

plot_det <- alt_rows %>%
  filter(forecaster == "Estimated", margin == "tolerant") %>%
  dplyr::select(scenario_label, t_occ, t_mag) %>%
  pivot_longer(cols = c("t_occ", "t_mag"),
               names_to = "head", values_to = "cross_time") %>%
  filter(!is.na(cross_time)) %>%
  mutate(head = recode(head, t_occ = "Occurrence", t_mag = "Magnitude"))

fig_det <- ggplot(plot_det, aes(x = cross_time, fill = head)) +
  geom_histogram(bins = 40, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = CONFIG$shift_time, linetype = "dashed") +
  facet_wrap(~ scenario_label) +
  labs(x = "First crossing time", y = "Count", fill = NULL) +
  theme_minimal() + theme(legend.position = "top")
save_figure(fig_det, "sim_detection_histogram", width = 10, height = 5)
