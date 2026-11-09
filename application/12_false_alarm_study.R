# False alarm calibration of four forecasters on many windows (data/protocol_E.md).
# From the project root:  source("application/12_false_alarm_study.R")
# Smoke test on two series and two windows each:  Sys.setenv(FA_SMOKE = "1"); source("application/12_false_alarm_study.R")
# The run can be interrupted and resumed: finished windows are kept in windows.csv.
# Writes: out/events/falsealarm/windows.csv, summary.csv, decision.txt, session_info.txt
source("R/io.R"); source("R/eprocess.R"); source("R/helpers.R"); source("R/forecast.R"); source("R/event_study.R")
suppressPackageStartupMessages(library(pscl))

SMOKE <- Sys.getenv("FA_SMOKE") == "1"
OUT <- if (SMOKE) "out/events/falsealarm_smoke" else "out/events/falsealarm"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
N_MON <- 365L; STRIDE <- 182L; T0_FIRST <- as.Date("2015-12-31"); LAST <- as.Date("2024-12-31")
REFIT <- 30L; ROLL <- 730L; BOOT_B <- 2000L; BOOT_SEED <- 20261109L
FC_NAMES <- c("F0", "F1", "F2", "F3")

fit_one <- function(df, tr, te) {
  fit <- tryCatch(pscl::hurdle(EV_FORM, data = df[tr, ], dist = "negbin"), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  X <- model.matrix(EV_FORM, df[te, ])
  list(pi = as.numeric(plogis(X %*% fit$coefficients$zero)),
       mu = as.numeric(exp(X %*% fit$coefficients$count)),
       theta = rep(as.numeric(fit$theta), length(te)))
}
# window = NULL uses all data up to the refit date; refit = NULL never refits after t0
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

# PIT functions that accept a theta that changes between refits (same random draws as R/forecast.R and R/helpers.R when theta is constant)
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
alarm_row <- function(fc) {
  hd <- monitor_heads_vec(fc); n <- length(fc$te$count)
  i <- c(occ = first_idx(hd$occ, EV_THR_HEAD(n)), mag = first_idx(hd$mag, EV_THR_HEAD(n)), whole = first_idx(hd$whole, EV_THR_WHOLE(n)))
  any_head <- !(is.na(i[["occ"]]) && is.na(i[["mag"]]))
  c(occ = !is.na(i[["occ"]]), mag = !is.na(i[["mag"]]), any_head = any_head, whole = !is.na(i[["whole"]]))
}

ev <- read_events("events/events.csv")
files <- c(list.files("data/faers", pattern = "^r[0-9]+\\.csv$", full.names = TRUE),
           list.files("data/faers_events", pattern = "^E[0-9]+_.*\\.csv$", full.names = TRUE))
if (SMOKE) files <- files[c(1, length(files))]
store <- file.path(OUT, "windows.csv")
done <- if (file.exists(store)) read.csv(store, stringsAsFactors = FALSE) else NULL
key_done <- if (is.null(done)) character() else paste(done$series, done$t0)

for (f in files) {
  sid <- sub("\\.csv$", "", basename(f))
  d <- read_series(f); df <- build_features(d)
  ed <- ev$event_date[ev$event_id == sid]
  last_end <- if (length(ed) == 1) ed - 30 else LAST
  t0s <- seq(T0_FIRST, by = STRIDE, length.out = 40); t0s <- t0s[t0s + N_MON <= min(last_end, LAST)]
  if (SMOKE) t0s <- head(t0s, 2)
  cat(sid, ":", length(t0s), "candidate windows\n")
  for (t0 in as.character(t0s)) {
    t0d <- as.Date(t0)
    if (paste(sid, t0) %in% key_done) next
    sc <- screen_event(d, list(train_end = t0d, mon_end = t0d + N_MON))
    if (!isTRUE(sc$pass)) { cat("  ", t0, "screened out\n"); next }
    fcs <- list(F0 = make_adaptive(df, t0d), F1 = make_adaptive(df, t0d, window = ROLL),
                F2 = make_adaptive(df, t0d, refit = REFIT), F3 = make_adaptive(df, t0d, window = ROLL, refit = REFIT))
    if (any(vapply(fcs, is.null, logical(1)))) { cat("  ", t0, "fit failed, window excluded\n"); next }
    rows <- do.call(rbind, lapply(FC_NAMES, function(nm) {
      a <- alarm_row(fcs[[nm]])
      data.frame(series = sid, t0 = t0, forecaster = nm, occ = a[["occ"]], mag = a[["mag"]], any_head = a[["any_head"]], whole = a[["whole"]])
    }))
    write.table(rows, store, sep = ",", row.names = FALSE, col.names = !file.exists(store), append = file.exists(store), qmethod = "double")
    cat("  ", t0, "done\n")
  }
}

w <- read.csv(store, stringsAsFactors = FALSE)
metrics <- c("occ", "mag", "any_head", "whole")
agg <- aggregate(w[metrics], by = list(series = w$series, forecaster = w$forecaster), FUN = sum)
cnt <- aggregate(w["occ"], by = list(series = w$series, forecaster = w$forecaster), FUN = length); names(cnt)[3] <- "n"
agg <- merge(agg, cnt, by = c("series", "forecaster"))
ser <- unique(agg$series)
rate <- function(fc, m, s = ser) { a <- agg[agg$forecaster == fc & agg$series %in% s, ]; sum(a[[m]]) / sum(a$n) }
set.seed(BOOT_SEED)
boot <- replicate(BOOT_B, sample(ser, length(ser), replace = TRUE), simplify = FALSE)
brate <- function(fc, m) vapply(boot, function(s) { a <- agg[agg$forecaster == fc, ]; a <- a[match(s, a$series), ]; sum(a[[m]], na.rm = TRUE) / sum(a$n, na.rm = TRUE) }, numeric(1))
summ <- do.call(rbind, lapply(FC_NAMES, function(fc) do.call(rbind, lapply(metrics, function(m) {
  b <- brate(fc, m); data.frame(forecaster = fc, outcome = m, windows = sum(agg$n[agg$forecaster == fc]), rate = rate(fc, m),
                                lo = unname(quantile(b, 0.025, na.rm = TRUE)), hi = unname(quantile(b, 0.975, na.rm = TRUE)))
}))))
write.csv(summ, file.path(OUT, "summary.csv"), row.names = FALSE)
print(summ, row.names = FALSE, digits = 3)

r0 <- rate("F0", "any_head")
cand <- do.call(rbind, lapply(c("F1", "F2", "F3"), function(fc) data.frame(forecaster = fc, rate = rate(fc, "any_head"))))
cand$qualifies <- cand$rate <= 0.10 & cand$rate <= 0.5 * r0
msg <- if (!any(cand$qualifies)) "No forecaster qualifies under the decision rule of data/protocol_E.md: negative result." else {
  q <- cand[cand$qualifies, ]; best <- min(q$rate); q <- q[q$rate <= best + 0.01, ]
  paste("Carried forward:", q$forecaster[1], "(lowest rate; ties within 0.01 resolved in the order F1, F2, F3)")
}
cat("\nF0 either head rate:", round(r0, 3), "\n"); print(cand, row.names = FALSE, digits = 3); cat(msg, "\n")
writeLines(c(paste("F0 either-head rate:", round(r0, 4)), capture.output(print(cand, row.names = FALSE)), msg), file.path(OUT, "decision.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info.txt"))
