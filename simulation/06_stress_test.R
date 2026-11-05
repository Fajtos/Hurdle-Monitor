suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(ggplot2)
})

source("R/io.R")
set_out_root("out/simulation")
source("R/simulate.R")
source("R/forecast.R")
source("R/eprocess.R")
source("R/replicate.R")

N_REP <- 400
ALPHA <- 0.05
PROGRESS_EVERY <- 100

CONFIG <- list(
  n_train = 2000,
  T_monitor = 3000,
  shift_time = 1000,
  alpha = ALPHA,
  margins = list(tolerant = c(0.02, 0.03))
)

SHIFTS <- list(
  list(name = "moderate", d_occ = 0, d_mag = 0.22, theta = 5),
  list(name = "small",    d_occ = 0, d_mag = 0.12, theta = 5)
)

run_shift <- function(shift, seed, config) {
  df_full <- simulate_hurdle(
    T_total = config$n_train + config$T_monitor,
    n_train = config$n_train,
    shift_time = config$shift_time,
    d_occ = shift$d_occ, d_mag = shift$d_mag, theta = shift$theta,
    seed = seed
  )
  df_train <- df_full[1:config$n_train, ]
  df_test <- df_full[(config$n_train + 1):(config$n_train + config$T_monitor), ]
  fc <- fit_and_predict(df_train, df_test, oracle = FALSE)
  if (is.null(fc)) return(NULL)
  Z <- as.integer(df_test$y > 0)
  U_pos <- randomized_pit_positive(df_test$y, fc$mu_ut, fc$theta, seed = seed + 1L)
  U_full <- randomized_pit_full(df_test$y, fc$pi, fc$mu_ut, fc$theta, seed = seed + 2L)
  logE_mag <- eprocess_pit(U_pos, Z, 0.03)
  logE_whole <- eprocess_pit(U_full, rep(1, length(U_full)), 0.03)
  data.frame(
    shift = shift$name,
    seed = seed,
    pi_mean = mean(fc$pi),
    rate_mag = logE_mag[config$T_monitor] / config$T_monitor,
    rate_whole = logE_whole[config$T_monitor] / config$T_monitor,
    t_mag = first_cross(logE_mag, log(2 / config$alpha)),
    t_whole = first_cross(logE_whole, log(1 / config$alpha))
  )
}

n_total <- length(SHIFTS) * N_REP
cat(sprintf("Stress test: %d jobs.\n\n", n_total))

counter <- new.env()
counter$n <- 0
start_time <- Sys.time()

jobs <- expand.grid(
  shift_idx = seq_along(SHIFTS),
  seed = seq_len(N_REP) * 1000 + 23,
  stringsAsFactors = FALSE
)

results <- purrr::pmap_dfr(
  list(jobs$shift_idx, jobs$seed),
  function(idx, sd) {
    counter$n <- counter$n + 1
    if (counter$n %% PROGRESS_EVERY == 0) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      rate <- counter$n / max(elapsed, 0.01)
      eta <- (n_total - counter$n) / rate
      cat(sprintf("Progress: %d / %d (%.1f%%)  elapsed %.1f min  ETA %.1f min\n",
                  counter$n, n_total, 100 * counter$n / n_total, elapsed, eta))
    }
    tryCatch(
      run_shift(SHIFTS[[idx]], sd, CONFIG),
      error = function(e) {
        message("Stress test failed: shift=", SHIFTS[[idx]]$name,
                " seed=", sd, " -> ", conditionMessage(e))
        NULL
      }
    )
  }
)

if (is.null(results) || nrow(results) == 0) {
  stop("All stress test replications failed.")
}

save_raw(results, "dilution_raw")

dilution <- results %>%
  group_by(shift) %>%
  summarise(
    mag_detect_rate = round(mean(!is.na(t_mag)), 3),
    whole_detect_rate = round(mean(!is.na(t_whole)), 3),
    detection_ratio = round(mean(!is.na(t_mag)) / max(mean(!is.na(t_whole)), 1e-6), 2),
    mag_delay_median = safe_median_delay(t_mag, !is.na(t_mag), CONFIG$shift_time),
    whole_delay_median = safe_median_delay(t_whole, !is.na(t_whole), CONFIG$shift_time),
    mean_logE_mag = round(mean(rate_mag) * CONFIG$T_monitor, 2),
    mean_logE_whole = round(mean(rate_whole) * CONFIG$T_monitor, 2),
    theoretical_ratio = round(mean(1 / pi_mean^3), 2),
    n = n(),
    .groups = "drop"
  )
write_table(dilution, "dilution",
            caption = "Dilution stress test. Detection rates, delays, and mean terminal log e-process under magnitude-only shifts.")

print(dilution)

sample_seed <- 1023
shift_for_figure <- SHIFTS[[1]]
df_full <- simulate_hurdle(
  T_total = CONFIG$n_train + CONFIG$T_monitor,
  n_train = CONFIG$n_train,
  shift_time = CONFIG$shift_time,
  d_occ = shift_for_figure$d_occ, d_mag = shift_for_figure$d_mag,
  theta = shift_for_figure$theta, seed = sample_seed
)
df_train <- df_full[1:CONFIG$n_train, ]
df_test <- df_full[(CONFIG$n_train + 1):(CONFIG$n_train + CONFIG$T_monitor), ]
fc <- fit_and_predict(df_train, df_test, oracle = FALSE)

Z <- as.integer(df_test$y > 0)
U_pos <- randomized_pit_positive(df_test$y, fc$mu_ut, fc$theta, seed = sample_seed + 1L)
U_full <- randomized_pit_full(df_test$y, fc$pi, fc$mu_ut, fc$theta, seed = sample_seed + 2L)
logE_mag <- eprocess_pit(U_pos, Z, 0.03)
logE_whole <- eprocess_pit(U_full, rep(1, length(U_full)), 0.03)

figure_data <- data.frame(
  t = 1:CONFIG$T_monitor,
  logE_mag = logE_mag,
  logE_whole = logE_whole
)

fig_dilution <- ggplot(figure_data, aes(x = t)) +
  geom_line(aes(y = logE_mag, colour = "Magnitude head"), linewidth = 0.8) +
  geom_line(aes(y = logE_whole, colour = "Whole-forecast PIT"),
            linewidth = 0.8, linetype = "dashed") +
  geom_hline(yintercept = log(2 / CONFIG$alpha), linetype = "dotted") +
  geom_vline(xintercept = CONFIG$shift_time, linetype = "dotdash") +
  scale_colour_manual(values = c("Magnitude head" = "steelblue",
                                 "Whole-forecast PIT" = "darkgray")) +
  labs(x = "Monitoring time (days)", y = expression(log~E[t]), colour = NULL) +
  theme_minimal() + theme(legend.position = "top")
save_figure(fig_dilution, "dilution", width = 8, height = 4)

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))