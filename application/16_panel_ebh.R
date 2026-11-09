# Panel monitoring with e-BH on the e-detector statistics (data/protocol_H.md).
# From the project root:  source("application/16_panel_ebh.R")
# Reads: out/events/panel/traj/*.rds (from application/15_panel_trajectories.R).  Writes: out/events/panel/panel_summary.csv, decision_H.txt, session_info.txt
OUT <- if (Sys.getenv("FA_SMOKE") == "1") "out/events/panel_smoke" else "out/events/panel"
N <- 365L; ONSET <- 181L; ALPHA <- 0.05; P_SHIFT <- 0.3; B <- 1000L; SEED <- 20261112L; BOOT_B <- 2000L; MIN_M <- 3L
FC_NAMES <- c("F0", "F3"); EFFECT <- c("mag_up15", "mag_up2", "mag_down", "occ_drop", "occ_up"); SCEN <- c("none", EFFECT)

fl <- list.files(file.path(OUT, "traj"), full.names = TRUE)
if (length(fl) == 0) stop("no trajectory files in ", file.path(OUT, "traj"))
lse2 <- function(a, b) { m <- pmax(a, b); m + log(exp(a - m) + exp(b - m)) }
LE <- list(); meta <- NULL
for (f in fl) {
  x <- readRDS(f)
  for (fc in FC_NAMES) LE[[paste(x$series, x$t0, x$scenario, fc, sep = "|")]] <- lse2(x$traj[[fc]]$occ, x$traj[[fc]]$mag) - log(2 * N)
  meta <- rbind(meta, data.frame(series = x$series, t0 = x$t0, scenario = x$scenario))
}
cnt <- table(meta$series, meta$t0, dnn = c("series", "t0"))
cmp <- aggregate(scenario ~ series + t0, meta, function(s) all(SCEN %in% s))
cmp <- cmp[cmp$scenario, ]
origins <- sort(unique(cmp$t0)); origins <- origins[vapply(origins, function(o) sum(cmp$t0 == o) >= MIN_M, logical(1))]
cat(length(origins), "origins with at least", MIN_M, "complete series\n")

ebh_reject <- function(Mat) {
  M <- ncol(Mat); thr <- log(M / (ALPHA * seq_len(M)))
  S <- t(apply(Mat, 1, sort, decreasing = TRUE))
  ok <- S >= matrix(rep(thr, each = nrow(Mat)), nrow(Mat))
  W <- ok * col(ok); khat <- W[cbind(seq_len(nrow(W)), max.col(W, "first"))]
  R <- matrix(FALSE, nrow(Mat), M)
  for (i in seq_len(M)) R[, i] <- (1 + rowSums(Mat > Mat[, i])) <= khat
  R
}
rules <- function(Mat) {
  M <- ncol(Mat)
  list(ebh = ebh_reject(Mat), bonferroni = Mat >= log(M / ALPHA), uncorrected = Mat >= log(1 / ALPHA))
}
score <- function(R, shifted) {
  truth <- outer(seq_len(N) >= ONSET, shifted, "&")
  nrej <- rowSums(R); any_t <- which(nrej > 0)
  fdp_end <- if (nrej[N] > 0) sum(R[N, ] & !shifted) / nrej[N] else 0
  pw_end <- if (any(shifted)) sum(R[N, ] & shifted) / sum(shifted) else NA_real_
  if (length(any_t) > 0) { tau <- any_t[1]; fdp_tau <- sum(R[tau, ] & !truth[tau, ]) / nrej[tau]; tp_tau <- sum(R[tau, ] & truth[tau, ]) } else { fdp_tau <- NA_real_; tp_tau <- NA_real_ }
  c(fdp_end = fdp_end, power_end = pw_end, reject_any = as.numeric(length(any_t) > 0), fdp_tau = fdp_tau, tp_tau = tp_tau)
}

set.seed(SEED)
res <- list()
for (o in origins) {
  ser <- cmp$series[cmp$t0 == o]; M <- length(ser)
  for (fc in FC_NAMES) {
    acc <- array(NA_real_, c(B, 3, 5), dimnames = list(NULL, c("ebh", "bonferroni", "uncorrected"), c("fdp_end", "power_end", "reject_any", "fdp_tau", "tp_tau")))
    for (b in seq_len(B)) {
      sc <- ifelse(runif(M) < P_SHIFT, sample(EFFECT, M, replace = TRUE), "none")
      Mat <- vapply(seq_len(M), function(i) LE[[paste(ser[i], o, sc[i], fc, sep = "|")]], numeric(N))
      rl <- rules(Mat); sh <- sc != "none"
      for (r in names(rl)) acc[b, r, ] <- score(rl[[r]], sh)
    }
    for (r in dimnames(acc)[[2]]) res[[length(res) + 1]] <- data.frame(origin = o, forecaster = fc, rule = r, M = M, t(colMeans(acc[, r, ], na.rm = TRUE)))
  }
}
res <- do.call(rbind, res)
set.seed(SEED + 1)
bt <- replicate(BOOT_B, sample(origins, length(origins), replace = TRUE), simplify = FALSE)
summ <- do.call(rbind, lapply(FC_NAMES, function(fc) do.call(rbind, lapply(c("ebh", "bonferroni", "uncorrected"), function(r) do.call(rbind, lapply(c("fdp_end", "power_end", "reject_any", "fdp_tau", "tp_tau"), function(m) {
  a <- res[res$forecaster == fc & res$rule == r, ]
  est <- mean(a[[m]], na.rm = TRUE)
  b <- vapply(bt, function(s) mean(a[[m]][match(s, a$origin)], na.rm = TRUE), numeric(1))
  data.frame(forecaster = fc, rule = r, quantity = m, estimate = est, lo = unname(quantile(b, 0.025, na.rm = TRUE)), hi = unname(quantile(b, 0.975, na.rm = TRUE)))
}))))))
write.csv(summ, file.path(OUT, "panel_summary.csv"), row.names = FALSE)
print(summ, row.names = FALSE, digits = 3)

# decision rule of the protocol: e-BH controls the false discovery proportion at the first rejection time if its estimate is at most ALPHA, and it is more powerful than Bonferroni if its power at day 365 is larger
lines <- character()
for (fc in FC_NAMES) {
  g <- function(r, q) summ$estimate[summ$forecaster == fc & summ$rule == r & summ$quantity == q]
  lines <- c(lines, paste0(fc, ": false discovery proportion at the first rejection, e-BH ", round(g("ebh", "fdp_tau"), 3), " (", if (g("ebh", "fdp_tau") <= ALPHA) "within" else "above", " ", ALPHA, "); uncorrected ", round(g("uncorrected", "fdp_tau"), 3),
                            ". Power at day ", N, ": e-BH ", round(g("ebh", "power_end"), 3), ", Bonferroni ", round(g("bonferroni", "power_end"), 3), ", uncorrected ", round(g("uncorrected", "power_end"), 3), "."))
}
cat(lines, sep = "\n"); writeLines(lines, file.path(OUT, "decision_H.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info.txt"))
