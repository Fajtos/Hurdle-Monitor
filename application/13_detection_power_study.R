# Detection power of the frozen forecaster and the rolling refit forecaster under injected shifts (data/protocol_F.md).
# From the project root:  source("application/13_detection_power_study.R")
# Smoke test on two windows:  Sys.setenv(FA_SMOKE = "1"); source("application/13_detection_power_study.R")
# The run can be interrupted and resumed: finished jobs are kept in jobs.csv.
# Reads: out/events/falsealarm/windows.csv.  Writes: out/events/power/jobs.csv, power_by_scenario.csv, power_pooled.csv, decision_F.txt, run_log.txt, session_info.txt
source("R/io.R"); source("R/eprocess.R"); source("R/helpers.R"); source("R/forecast.R"); source("R/event_study.R")
suppressPackageStartupMessages({ library(pscl); library(digest) })

SMOKE <- Sys.getenv("FA_SMOKE") == "1"
OUT <- if (SMOKE) "out/events/power_smoke" else "out/events/power"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
N_MON <- 365L; ONSET <- 181L; REFIT <- 30L; ROLL <- 730L
BOOT_B <- 2000L; BOOT_SEED <- 20261110L; INJ_SEED <- 20261110L
FC_NAMES <- c("F0", "F3")
SCEN <- c("none", "mag_up15", "mag_up2", "mag_down", "occ_drop", "occ_up")
EFFECT <- setdiff(SCEN, "none")
RETAIN_RATIO <- 0.8

# counter for fit warnings (a job whose fit fails is excluded and reported on the console)
CNT <- new.env(); CNT$warn <- 0L

fit_one <- function(df, tr, te) {
  fit <- tryCatch(withCallingHandlers(pscl::hurdle(EV_FORM, data = df[tr, ], dist = "negbin"),
                                      warning = function(w) { CNT$warn <- CNT$warn + 1L; invokeRestart("muffleWarning") }),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  X <- model.matrix(EV_FORM, df[te, ])
  list(pi = as.numeric(plogis(X %*% fit$coefficients$zero)),
       mu = as.numeric(exp(X %*% fit$coefficients$count)),
       theta = rep(as.numeric(fit$theta), length(te)))
}
make_adaptive <- function(df, t0, window = NULL, refit = NULL) {
  te_all <- which(df$date > t0 & df$date <= t0 + N_MON)
  off <- as.integer(df$date[te_all] - t0) - 1L
  block <- if (is.null(refit)) rep(0L, length(te_all)) else off %/% refit
  pi <- mu <- theta <- numeric(length(te_all))
  for (b in unique(block)) {
    sel <- which(block == b)
    fit_date <- t0 + if (is.null(refit)) 0L else b * refit
    lo <- if (is.null(window)) as.Date("1900-01-01") else fit_date - window
    tr <- which(df$date <= fit_date & df$date > lo)
    r <- fit_one(df, tr, te_all[sel])
    if (is.null(r)) return(NULL)
    pi[sel] <- r$pi; mu[sel] <- r$mu; theta[sel] <- r$theta
  }
  list(te = df[te_all, ], pi = pi, mu = mu, theta = theta)
}
pit_pos_vec <- function(y, mu, theta, seed) {
  set.seed(seed); n <- length(y); U <- rep(0.5, n); pos <- which(y > 0)
  if (length(pos) > 0) {
    th <- theta[pos]; m <- mu[pos]
    F_y <- pnbinom(y[pos], size = th, mu = m); F_ym1 <- pnbinom(y[pos] - 1, size = th, mu = m); F_0 <- pnbinom(0, size = th, mu = m)
    Fy <- (F_y - F_0) / (1 - F_0); Fym1 <- (F_ym1 - F_0) / (1 - F_0)
    V <- runif(length(pos)); U[pos] <- Fym1 + V * (Fy - Fym1)
  }
  U
}
rqr_vec <- function(y, pi_pos, mu, theta, seed) {
  set.seed(seed); n <- length(y); q0 <- (theta / (theta + mu))^theta; U <- numeric(n); pos <- which(y > 0)
  if (length(pos) > 0) {
    th <- theta[pos]; m <- mu[pos]
    Fy <- (pnbinom(y[pos], size = th, mu = m) - q0[pos]) / (1 - q0[pos]); Fym1 <- (pnbinom(y[pos] - 1, size = th, mu = m) - q0[pos]) / (1 - q0[pos])
    Ffy <- (1 - pi_pos[pos]) + pi_pos[pos] * Fy; Ffym1 <- (1 - pi_pos[pos]) + pi_pos[pos] * Fym1
    U[pos] <- Ffym1 + runif(length(pos)) * (Ffy - Ffym1)
  }
  zero <- which(y == 0)
  if (length(zero) > 0) U[zero] <- runif(length(zero)) * (1 - pi_pos[zero])
  qnorm(U)
}
monitor_heads_vec <- function(fc) {
  y <- fc$te$count; n <- length(y); Z <- as.integer(y > 0)
  U_pos  <- pit_pos_vec(y, fc$mu, fc$theta, EV_PIT_SEEDS[["pos"]])
  U_full <- pnorm(rqr_vec(y, fc$pi, fc$mu, fc$theta, EV_PIT_SEEDS[["whole"]]))
  list(occ   = edetector_log(occ_steps(fc$pi, Z, EV_TOL[1])),
       mag   = edetector_log(pit_steps(U_pos, Z, EV_TOL[2])),
       whole = edetector_log(pit_steps(U_full, rep(1, n), EV_TOL[2])))
}
# first alarm day (index in the monitored period) for either head and for the whole forecast
alarm_idx <- function(fc) {
  hd <- monitor_heads_vec(fc); n <- length(fc$te$count)
  io <- first_idx(hd$occ, EV_THR_HEAD(n)); im <- first_idx(hd$mag, EV_THR_HEAD(n)); iw <- first_idx(hd$whole, EV_THR_WHOLE(n))
  ia <- if (all(is.na(c(io, im)))) NA_integer_ else min(c(io, im), na.rm = TRUE)
  c(any_head = ia, whole = iw)
}

# shifts are applied to the raw counts from the onset day on; the forecasters then see the shifted series
inject <- function(d, t0, kind, seed) {
  if (kind == "none") return(d)
  idx <- which(d$date >= t0 + ONSET & d$date <= t0 + N_MON)
  y <- d$count[idx]; n <- length(y); set.seed(seed)
  weekday <- as.integer(format(d$date[idx], "%u")) < 6
  d$count[idx] <- switch(kind,
                         mag_up15 = y + rbinom(n, y, 0.5),
                         mag_up2  = 2L * y,
                         mag_down = rbinom(n, y, 0.5),
                         occ_drop = ifelse(y > 0 & runif(n) < 0.3, 0L, y),
                         occ_up   = ifelse(y == 0 & weekday & runif(n) < 0.3, 1L, y))
  d$count <- as.integer(d$count)
  d
}
job_seed <- function(sid, t0, kind) strtoi(substr(digest(paste(INJ_SEED, sid, t0, kind), algo = "md5"), 1, 6), 16L)

pe <- read.csv("out/events/falsealarm/windows.csv", stringsAsFactors = FALSE)
wins <- unique(pe[, c("series", "t0")])
wins <- wins[order(wins$series, wins$t0), ]
if (SMOKE) wins <- wins[c(1, nrow(wins)), ]
files <- c(list.files("data/faers", pattern = "^r[0-9]+\\.csv$", full.names = TRUE),
           list.files("data/faers_events", pattern = "^E[0-9]+_.*\\.csv$", full.names = TRUE))
path_of <- function(sid) files[sub("\\.csv$", "", basename(files)) == sid]

store <- file.path(OUT, "jobs.csv")
done <- if (file.exists(store)) read.csv(store, stringsAsFactors = FALSE) else NULL
key_done <- if (is.null(done)) character() else unique(paste(done$series, done$t0, done$scenario))
cat(nrow(wins), "windows,", length(SCEN), "scenarios\n")

for (i in seq_len(nrow(wins))) {
  sid <- wins$series[i]; t0 <- wins$t0[i]; t0d <- as.Date(t0)
  d0 <- read_series(path_of(sid))
  for (sc in SCEN) {
    if (paste(sid, t0, sc) %in% key_done) next
    d <- inject(d0, t0d, sc, job_seed(sid, t0, sc))
    df <- build_features(d)
    w0 <- CNT$warn
    fcs <- list(F0 = make_adaptive(df, t0d), F3 = make_adaptive(df, t0d, window = ROLL, refit = REFIT))
    if (any(vapply(fcs, is.null, logical(1)))) { cat("  ", sid, t0, sc, "fit failed, job excluded\n"); next }
    rows <- do.call(rbind, lapply(FC_NAMES, function(nm) {
      a <- alarm_idx(fcs[[nm]])
      data.frame(series = sid, t0 = t0, scenario = sc, forecaster = nm, idx_any = a[["any_head"]], idx_whole = a[["whole"]],
                 warn = CNT$warn - w0)
    }))
    write.table(rows, store, sep = ",", row.names = FALSE, col.names = !file.exists(store), append = file.exists(store), qmethod = "double")
  }
  cat("  ", sid, t0, "done\n")
}

j <- read.csv(store, stringsAsFactors = FALSE)
j$pre_any <- !is.na(j$idx_any) & j$idx_any < ONSET
j$hit_any <- !is.na(j$idx_any) & j$idx_any >= ONSET
j$pre_whole <- !is.na(j$idx_whole) & j$idx_whole < ONSET
j$hit_whole <- !is.na(j$idx_whole) & j$idx_whole >= ONSET
j$delay_any <- ifelse(j$hit_any, j$idx_any - ONSET, NA)

# consistency check: the no-injection jobs must reproduce the alarms of Part E
n0 <- j[j$scenario == "none", ]
pe_any <- pe[, c("series", "t0", "forecaster", "any_head")]
chk <- merge(n0, pe_any, by = c("series", "t0", "forecaster"))
chk_msg <- paste("No-injection alarms equal Part E in", sum((!is.na(chk$idx_any)) == chk$any_head), "of", nrow(chk), "series-window-forecaster rows")
cat(chk_msg, "\n")

ser <- unique(j$series)
by_scen <- do.call(rbind, lapply(SCEN, function(sc) do.call(rbind, lapply(FC_NAMES, function(fc) {
  a <- j[j$forecaster == fc & j$scenario == sc, ]
  data.frame(scenario = sc, forecaster = fc, windows = nrow(a), alarm_before_onset = mean(a$pre_any), correct_detection = mean(a$hit_any),
             conditional_power = if (sum(!a$pre_any) > 0) sum(a$hit_any) / sum(!a$pre_any) else NA_real_,
             median_delay_days = if (any(a$hit_any)) median(a$delay_any, na.rm = TRUE) else NA_real_,
             whole_correct = mean(a$hit_whole))
}))))
write.csv(by_scen, file.path(OUT, "power_by_scenario.csv"), row.names = FALSE)
print(by_scen, row.names = FALSE, digits = 3)

# pooled over the five shift scenarios (each scenario has equal weight); series are resampled for the intervals
set.seed(BOOT_SEED)
boot <- replicate(BOOT_B, sample(ser, length(ser), replace = TRUE), simplify = FALSE)
pool_stat <- function(s_idx) {
  a <- do.call(rbind, lapply(s_idx, function(s) j[j$series == s, ]))
  f <- function(fc, col) mean(tapply(a[[col]][a$forecaster == fc & a$scenario %in% EFFECT], a$scenario[a$forecaster == fc & a$scenario %in% EFFECT], mean))
  c(F0 = f("F0", "hit_any"), F3 = f("F3", "hit_any"))
}
pt <- pool_stat(ser)
bs <- sapply(boot, pool_stat)
pooled <- data.frame(quantity = c("F0 correct detection", "F3 correct detection", "F3 minus F0", "F3 over F0"),
                     estimate = c(pt[["F0"]], pt[["F3"]], pt[["F3"]] - pt[["F0"]], pt[["F3"]] / pt[["F0"]]),
                     lo = c(quantile(bs["F0", ], 0.025), quantile(bs["F3", ], 0.025), quantile(bs["F3", ] - bs["F0", ], 0.025), quantile(bs["F3", ] / bs["F0", ], 0.025, na.rm = TRUE)),
                     hi = c(quantile(bs["F0", ], 0.975), quantile(bs["F3", ], 0.975), quantile(bs["F3", ] - bs["F0", ], 0.975), quantile(bs["F3", ] / bs["F0", ], 0.975, na.rm = TRUE)))
write.csv(pooled, file.path(OUT, "power_pooled.csv"), row.names = FALSE)
print(pooled, row.names = FALSE, digits = 3)

ratio <- pt[["F3"]] / pt[["F0"]]
msg <- if (ratio >= RETAIN_RATIO) {
  paste0("Power retained: F3 correct detection is ", round(ratio, 3), " times that of F0 (rule: at least ", RETAIN_RATIO, ").")
} else {
  paste0("Power cost: F3 correct detection is ", round(ratio, 3), " times that of F0 (rule: at least ", RETAIN_RATIO, ").")
}
cat("\n", msg, "\n")
writeLines(c(chk_msg, capture.output(print(pooled, row.names = FALSE)), msg), file.path(OUT, "decision_F.txt"))
writeLines(c(paste("jobs:", nrow(j) / length(FC_NAMES)), paste("expected jobs:", nrow(wins) * length(SCEN)), paste("fit warnings:", sum(j$warn))), file.path(OUT, "run_log.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info.txt"))