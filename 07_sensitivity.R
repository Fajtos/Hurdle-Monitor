suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(ggplot2)
})

source("R/io.R")
source("R/simulate.R")
source("R/forecast.R")
source("R/eprocess.R")
source("R/replicate.R")

ALPHA <- 0.05
N_REP <- 400
PROGRESS_EVERY <- 100

SCENARIO_NULL <- list(name = "null", d_occ = 0, d_mag = 0, theta = 5)

BASE_CONFIG <- list(
  n_train = 2000,
  T_monitor = 3000,
  shift_time = 1000,
  alpha = ALPHA,
  margins = list(strict = 0, tolerant = c(0.02, 0.03))
)

n0_grid <- c(500, 1000, 2000, 4000)
T_grid  <- c(1000, 2000, 3000, 5000)

SEED_STREAM <- seq_len(N_REP) * 1000 + 7

cat("=== Part A: training-window sensitivity ===\n")
cat(sprintf("Running %d jobs (4 n0 values x %d replications, matching main seeds).\n\n",
            4 * N_REP, N_REP))
start_time <- Sys.time()
counter <- new.env(); counter$n <- 0

results_A <- purrr::pmap_dfr(
  expand.grid(n0 = n0_grid, seed = SEED_STREAM),
  function(n0, seed) {
    counter$n <- counter$n + 1
    if (counter$n %% PROGRESS_EVERY == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      rate <- counter$n / max(elapsed, 0.01)
      eta <- (4 * N_REP - counter$n) / rate
      cat(sprintf("  Progress %d / %d (%.1f%%), elapsed %.1f min, ETA %.1f min\n",
                  counter$n, 4 * N_REP, 100 * counter$n / (4 * N_REP),
                  elapsed, eta))
    }
    cfg <- BASE_CONFIG; cfg$n_train <- n0
    out <- tryCatch(run_one_replication(SCENARIO_NULL, seed, FALSE, cfg),
                    error = function(e) {
                      message("  Failed: n0=", n0, " seed=", seed,
                              " -> ", conditionMessage(e))
                      NULL
                    })
    if (!is.null(out)) out$n0 <- n0
    out
  }
)

n_succ_A <- nrow(results_A) / 2
cat(sprintf("\nPart A complete: %d / %d jobs succeeded.\n\n",
            n_succ_A, 4 * N_REP))

sens_n0 <- results_A %>%
  group_by(n0, margin) %>%
  summarise(
    occ_false_alarm_rate = round(mean(occ_flagged), 4),
    mag_false_alarm_rate = round(mean(mag_flagged), 4),
    n = n(),
    .groups = "drop"
  )
write_table(sens_n0, "sens_n0",
            caption = "Sensitivity of the null false-alarm rate to the training-window size. Null scenario, estimated forecaster, monitoring window fixed at 3000 days, 400 replications per cell, same seed stream as the main simulation.")

cat("=== Part B: monitoring-window sensitivity ===\n")
cat(sprintf("Running %d jobs (4 T values x %d replications, matching main seeds).\n\n",
            4 * N_REP, N_REP))
start_time <- Sys.time()
counter$n <- 0

results_B <- purrr::pmap_dfr(
  expand.grid(T_monitor = T_grid, seed = SEED_STREAM),
  function(T_monitor, seed) {
    counter$n <- counter$n + 1
    if (counter$n %% PROGRESS_EVERY == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      rate <- counter$n / max(elapsed, 0.01)
      eta <- (4 * N_REP - counter$n) / rate
      cat(sprintf("  Progress %d / %d (%.1f%%), elapsed %.1f min, ETA %.1f min\n",
                  counter$n, 4 * N_REP, 100 * counter$n / (4 * N_REP),
                  elapsed, eta))
    }
    cfg <- BASE_CONFIG; cfg$T_monitor <- T_monitor
    out <- tryCatch(run_one_replication(SCENARIO_NULL, seed, FALSE, cfg),
                    error = function(e) {
                      message("  Failed: T=", T_monitor, " seed=", seed,
                              " -> ", conditionMessage(e))
                      NULL
                    })
    if (!is.null(out)) out$T_monitor <- T_monitor
    out
  }
)

n_succ_B <- nrow(results_B) / 2
cat(sprintf("\nPart B complete: %d / %d jobs succeeded.\n\n",
            n_succ_B, 4 * N_REP))

sens_T <- results_B %>%
  group_by(T_monitor, margin) %>%
  summarise(
    occ_false_alarm_rate = round(mean(occ_flagged), 4),
    mag_false_alarm_rate = round(mean(mag_flagged), 4),
    n = n(),
    .groups = "drop"
  )
write_table(sens_T, "sens_T",
            caption = "Sensitivity of the null false-alarm rate to the monitoring-window length. Null scenario, estimated forecaster, training window fixed at 2000 days, 400 replications per cell, same seed stream as the main simulation.")

cat("=== Consistency check: (n0=2000, T=3000) vs main simulation ===\n")
main_rds <- "out/raw/simulation_raw.rds"
if (file.exists(main_rds)) {
  main_null <- readRDS(main_rds) %>%
    filter(scenario == "null", oracle == FALSE) %>%
    group_by(margin) %>%
    summarise(
      occ_main = round(mean(occ_flagged), 4),
      mag_main = round(mean(mag_flagged), 4),
      .groups = "drop"
    )
  overlap <- sens_n0 %>%
    filter(n0 == 2000) %>%
    dplyr::select(margin, occ_sens = occ_false_alarm_rate,
                  mag_sens = mag_false_alarm_rate)
  comparison <- main_null %>% inner_join(overlap, by = "margin")
  print(comparison)
  all_match <- all(comparison$occ_main == comparison$occ_sens) &&
    all(comparison$mag_main == comparison$mag_sens)
  if (all_match) {
    cat("\nPASS: overlapping cell matches the main simulation exactly.\n")
    cat("The sensitivity tables and Table 5.3 are arithmetically identical\n")
    cat("at the (n0 = 2000, T = 3000) point.\n\n")
  } else {
    cat("\nFAIL: overlapping cell differs from the main simulation.\n")
    cat("Check that BASE_CONFIG$n_train, BASE_CONFIG$T_monitor, BASE_CONFIG$margins\n")
    cat("and the seed stream all match 05_simulation.R exactly.\n\n")
  }
} else {
  cat("Main simulation raw file not found. Skipping consistency check.\n\n")
}

cat("=== Figures ===\n")

plot_A <- sens_n0 %>%
  mutate(head = "Occurrence", rate = occ_false_alarm_rate) %>%
  bind_rows(sens_n0 %>% mutate(head = "Magnitude", rate = mag_false_alarm_rate)) %>%
  mutate(margin_label = ifelse(margin == "strict", "Strict", "Tolerant"))

fig_A <- ggplot(plot_A, aes(x = n0, y = rate, colour = head, linetype = margin_label)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  geom_hline(yintercept = ALPHA / 2, linetype = "dashed", colour = "grey40") +
  scale_x_log10(breaks = n0_grid) +
  scale_colour_manual(values = c(Occurrence = "steelblue", Magnitude = "indianred")) +
  labs(x = expression(paste("Training-window size ", n[0])),
       y = "False-alarm rate", colour = "Head", linetype = "Margin") +
  theme_minimal()
save_figure(fig_A, "sens_n0", width = 8, height = 4.5)

plot_B <- sens_T %>%
  mutate(head = "Occurrence", rate = occ_false_alarm_rate) %>%
  bind_rows(sens_T %>% mutate(head = "Magnitude", rate = mag_false_alarm_rate)) %>%
  mutate(margin_label = ifelse(margin == "strict", "Strict", "Tolerant"))

fig_B <- ggplot(plot_B, aes(x = T_monitor, y = rate, colour = head, linetype = margin_label)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  geom_hline(yintercept = ALPHA / 2, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c(Occurrence = "steelblue", Magnitude = "indianred")) +
  labs(x = expression(paste("Monitoring-window length ", T)),
       y = "False-alarm rate", colour = "Head", linetype = "Margin") +
  theme_minimal()
save_figure(fig_B, "sens_T", width = 8, height = 4.5)

cat("Outputs written:\n")
cat("  out/tables/sens_n0.csv and .tex\n")
cat("  out/tables/sens_T.csv and .tex\n")
cat("  out/figures/sens_n0.pdf\n")
cat("  out/figures/sens_T.pdf\n")