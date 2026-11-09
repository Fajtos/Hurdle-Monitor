# Comparison with the rank based PIT monitor of Farran (2026) on the windows and shifts of Parts E and F (data/protocol_G.md).
# From the project root:  source("application/14_pit_monitor_comparison.R")
# Smoke test on two windows:  Sys.setenv(FA_SMOKE = "1"); source("application/14_pit_monitor_comparison.R")
# The run can be interrupted and resumed: finished jobs are kept in jobs.csv.
# Reads: out/events/falsealarm/windows.csv, out/events/power/jobs.csv.  Writes: out/events/pitmon/jobs.csv, by_scenario.csv, pooled.csv, attribution.csv, decision_G.txt, run_log.txt, session_info.txt
source("R/io.R"); source("R/eprocess.R"); source("R/helpers.R"); source("R/forecast.R"); source("R/event_study.R")
suppressPackageStartupMessages({ library(pscl); library(digest) })

SMOKE <- Sys.getenv("FA_SMOKE") == "1"
OUT <- if (SMOKE) "out/events/pitmon_smoke" else "out/events/pitmon"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
N_MON <- 365L; ONSET <- 181L; REFIT <- 30L; ROLL <- 730L
BOOT_B <- 2000L; BOOT_SEED <- 20261111L; INJ_SEED <- 20261110L
PM_B <- 100L; PM_ALPHA <- 0.05; PM_SEEDS <- c(whole = 31L, occ = 32L, mag = 33L)
FC_NAMES <- c("F0", "F3")
SCEN <- c("none", "mag_up15", "mag_up2", "mag_down", "occ_drop", "occ_up")
EFFECT <- setdiff(SCEN, "none")
RETAIN_RATIO <- 0.8
TARGET <- c(mag_up15 = "mag", mag_up2 = "mag", occ_drop = "occ")

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

# Rank based PIT monitor (Farran 2026): conformal p-values from the ranks of the PITs, histogram density bets (B bins, pseudocount 1),
# mixture over change points with weights 1 / (t (t + 1)); alarm when M_t reaches 1 / alpha.
pitmon_stat <- function(u, seed, B = PM_B) {
  n <- length(u); if (n == 0) return(numeric(0))
  set.seed(seed); V <- runif(n)
  cnt <- rep(1, B); tot <- B; M <- 0; out <- numeric(n)
  for (t in seq_len(n)) {
    R <- sum(u[seq_len(t)] <= u[t]); p <- (R - 1 + V[t]) / t
    b <- min(B, floor(p * B) + 1L)
    M <- B * cnt[b] / tot * (M + 1 / (t * (t + 1)))
    out[t] <- M; cnt[b] <- cnt[b] + 1; tot <- tot + 1
  }
  out
}
first_ge <- function(stat, thr) { i <- which(stat >= thr); if (length(i) > 0) i[1] else NA_integer_ }

# alarm days for both monitors on one forecast; indices refer to the monitored period
alarms_all <- function(fc) {
  y <- fc$te$count; n <- length(y); Z <- as.integer(y > 0)
  hd <- monitor_heads_vec(fc)
  ed <- c(occ = first_idx(hd$occ, EV_THR_HEAD(n)), mag = first_idx(hd$mag, EV_THR_HEAD(n)), whole = first_idx(hd$whole, EV_THR_WHOLE(n)))
  U_full <- pnorm(rqr_vec(y, fc$pi, fc$mu, fc$theta, EV_PIT_SEEDS[["whole"]]))
  pm_whole <- first_ge(pitmon_stat(U_full, PM_SEEDS[["whole"]]), 1 / PM_ALPHA)
  set.seed(PM_SEEDS[["occ"]] + 1000L); Vo <- runif(n)
  U_occ <- ifelse(Z == 1L, (1 - fc$pi) + Vo * fc$pi, Vo * (1 - fc$pi))
  pm_occ <- first_ge(pitmon_stat(U_occ, PM_SEEDS[["occ"]]), 2 / PM_ALPHA)
  pos <- which(Z == 1L)
  U_pos <- pit_pos_vec(y, fc$mu, fc$theta, EV_PIT_SEEDS[["pos"]])[pos]
  k <- first_ge(pitmon_stat(U_pos, PM_SEEDS[["mag"]]), 2 / PM_ALPHA)
  pm_mag <- if (is.na(k)) NA_integer_ else pos[k]
  c(ed_occ = ed[["occ"]], ed_mag = ed[["mag"]], ed_whole = ed[["whole"]], pm_whole = pm_whole, pm_occ = pm_occ, pm_mag = pm_mag)
}
min_na <- function(a, b) { v <- c(a, b); if (all(is.na(v))) NA_integer_ else as.integer(min(v, na.rm = TRUE)) }

pe <- read.csv("out/events/falsealarm/windows.csv", stringsAsFactors = FALSE)
wins <- unique(pe[, c("series", "t0")]); wins <- wins[order(wins$series, wins$t0), ]
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
      a <- alarms_all(fcs[[nm]])
      data.frame(series = sid, t0 = t0, scenario = sc, forecaster = nm, t(a), warn = CNT$warn - w0, row.names = NULL)
    }))
    write.table(rows, store, sep = ",", row.names = FALSE, col.names = !file.exists(store), append = file.exists(store), qmethod = "double")
  }
  cat("  ", sid, t0, "done\n")
}

j <- read.csv(store, stringsAsFactors = FALSE)
j$ed_any <- mapply(min_na, j$ed_occ, j$ed_mag)
j$pm_any <- mapply(min_na, j$pm_occ, j$pm_mag)
METHODS <- c(ed_any = "Headwise e-detector, either head", ed_whole = "E-detector, whole forecast", pm_whole = "PITMonitor, whole forecast", pm_any = "PIT rank monitor, either head")
for (m in names(METHODS)) {
  j[[paste0("pre_", m)]] <- !is.na(j[[m]]) & j[[m]] < ONSET
  j[[paste0("hit_", m)]] <- !is.na(j[[m]]) & j[[m]] >= ONSET
  j[[paste0("any_", m)]] <- !is.na(j[[m]])
}

# consistency check: the e-detector alarms must reproduce Part F
pf <- read.csv("out/events/power/jobs.csv", stringsAsFactors = FALSE)
chk <- merge(j, pf[, c("series", "t0", "scenario", "forecaster", "idx_any", "idx_whole")], by = c("series", "t0", "scenario", "forecaster"))
chk_msg <- paste("E-detector alarms equal Part F in", sum(mapply(identical, as.integer(chk$ed_any), as.integer(chk$idx_any)) & mapply(identical, as.integer(chk$ed_whole), as.integer(chk$idx_whole))), "of", nrow(chk), "rows")
cat(chk_msg, "\n")

by_scen <- do.call(rbind, lapply(names(METHODS), function(m) do.call(rbind, lapply(FC_NAMES, function(fc) do.call(rbind, lapply(SCEN, function(sc) {
  a <- j[j$forecaster == fc & j$scenario == sc, ]
  data.frame(method = METHODS[[m]], forecaster = fc, scenario = sc, windows = nrow(a), any_alarm = mean(a[[paste0("any_", m)]]),
             alarm_before_onset = mean(a[[paste0("pre_", m)]]), correct_detection = mean(a[[paste0("hit_", m)]]))
}))))))
write.csv(by_scen, file.path(OUT, "by_scenario.csv"), row.names = FALSE)
print(by_scen, row.names = FALSE, digits = 3)

ser <- unique(j$series)
set.seed(BOOT_SEED)
boot <- replicate(BOOT_B, sample(ser, length(ser), replace = TRUE), simplify = FALSE)
# false alarm rate: any alarm in the no shift scenario. Pooled detection: mean over the five shifts of the correct detection rate minus the no shift rate after day 181.
pool_stat <- function(s_idx, m, fc) {
  a <- do.call(rbind, lapply(s_idx, function(s) j[j$series == s & j$forecaster == fc, ]))
  hit <- tapply(a[[paste0("hit_", m)]], a$scenario, mean)
  c(false_alarm = mean(a[[paste0("any_", m)]][a$scenario == "none"]),
    correct = mean(hit[EFFECT]), net_correct = mean(hit[EFFECT]) - hit[["none"]])
}
pooled <- do.call(rbind, lapply(names(METHODS), function(m) do.call(rbind, lapply(FC_NAMES, function(fc) {
  pt <- pool_stat(ser, m, fc); bs <- sapply(boot, pool_stat, m = m, fc = fc)
  data.frame(method = METHODS[[m]], forecaster = fc, quantity = names(pt), estimate = unname(pt),
             lo = unname(apply(bs, 1, quantile, 0.025, na.rm = TRUE)), hi = unname(apply(bs, 1, quantile, 0.975, na.rm = TRUE)))
}))))
write.csv(pooled, file.path(OUT, "pooled.csv"), row.names = FALSE)
print(pooled, row.names = FALSE, digits = 3)

# attribution on the shifts that act on one head: counts times 1.5 and 2 (magnitude), positive days set to zero (occurrence)
attr_tab <- do.call(rbind, lapply(c("ed", "pm"), function(pfx) do.call(rbind, lapply(FC_NAMES, function(fc) do.call(rbind, lapply(names(TARGET), function(sc) {
  a <- j[j$forecaster == fc & j$scenario == sc, ]
  tgt <- a[[paste0(pfx, "_", TARGET[[sc]])]]; oth <- a[[paste0(pfx, "_", setdiff(c("occ", "mag"), TARGET[[sc]]))]]
  det <- a[[paste0("hit_", pfx, "_any")]]
  ok <- det & (!is.na(tgt)) & (is.na(oth) | tgt <= oth)
  data.frame(monitor = pfx, forecaster = fc, scenario = sc, windows_with_correct_time_alarm = sum(det), target_head_first_or_tied = sum(ok))
}))))))
write.csv(attr_tab, file.path(OUT, "attribution.csv"), row.names = FALSE)
print(attr_tab, row.names = FALSE)

# decision rule of the protocol: for each forecaster, the headwise e-detector is competitive with PITMonitor if its pooled net correct detection is at least RETAIN_RATIO times that of PITMonitor
lines <- chk_msg
for (fc in FC_NAMES) {
  a <- pooled$estimate[pooled$method == METHODS[["ed_any"]] & pooled$forecaster == fc & pooled$quantity == "net_correct"]
  b <- pooled$estimate[pooled$method == METHODS[["pm_whole"]] & pooled$forecaster == fc & pooled$quantity == "net_correct"]
  r <- a / b
  lines <- c(lines, paste0(fc, ": net correct detection of the headwise e-detector over PITMonitor is ", round(r, 3), if (is.finite(r) && r >= RETAIN_RATIO) " (competitive)" else " (not competitive)", " under the rule of at least ", RETAIN_RATIO, "."))
}
cat(lines, sep = "\n")
writeLines(lines, file.path(OUT, "decision_G.txt"))
writeLines(c(paste("jobs:", nrow(j) / length(FC_NAMES)), paste("expected jobs:", nrow(wins) * length(SCEN)), paste("fit warnings:", sum(j$warn))), file.path(OUT, "run_log.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info.txt"))
