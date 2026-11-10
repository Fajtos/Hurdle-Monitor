# Post hoc summaries of the panel study of Part H. This script is exploratory and is not part of the protocol in data/protocol_H.md.
# From the project root:  source("application/17_panel_posthoc.R")
# Reads: out/events/panel/traj/*.rds.  Writes: out/events/panel/posthoc.txt
# It repeats the random draws of application/16_panel_ebh.R in the same order with the same seed, so the instances are the same.
OUT <- "out/events/panel"
N <- 365L; ONSET <- 181L; ALPHA <- 0.05; P_SHIFT <- 0.3; B <- 1000L; SEED <- 20261112L; MIN_M <- 3L
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
cmp <- aggregate(scenario ~ series + t0, meta, function(s) all(SCEN %in% s))
cmp <- cmp[cmp$scenario, ]
origins <- sort(unique(cmp$t0)); origins <- origins[vapply(origins, function(o) sum(cmp$t0 == o) >= MIN_M, logical(1))]

ebh_reject <- function(Mat) {
  M <- ncol(Mat); thr <- log(M / (ALPHA * seq_len(M)))
  S <- t(apply(Mat, 1, sort, decreasing = TRUE))
  ok <- S >= matrix(rep(thr, each = nrow(Mat)), nrow(Mat))
  W <- ok * col(ok); khat <- W[cbind(seq_len(nrow(W)), max.col(W, "first"))]
  R <- matrix(FALSE, nrow(Mat), M)
  for (i in seq_len(M)) R[, i] <- (1 + rowSums(Mat > Mat[, i])) <= khat
  R
}

set.seed(SEED)
rows <- list()
for (o in origins) {
  ser <- cmp$series[cmp$t0 == o]; M <- length(ser)
  for (fc in FC_NAMES) {
    v <- rep(NA_real_, B)
    for (b in seq_len(B)) {
      sc <- ifelse(runif(M) < P_SHIFT, sample(EFFECT, M, replace = TRUE), "none")
      Mat <- vapply(seq_len(M), function(i) LE[[paste(ser[i], o, sc[i], fc, sep = "|")]], numeric(N))
      R <- ebh_reject(Mat); sh <- sc != "none"
      truth <- outer(seq_len(N) >= ONSET, sh, "&")
      nrej <- rowSums(R); any_t <- which(nrej > 0)
      if (length(any_t) > 0) { tau <- any_t[1]; v[b] <- sum(R[tau, ] & !truth[tau, ]) / nrej[tau] }
    }
    rows[[length(rows) + 1]] <- data.frame(origin = as.character(o), forecaster = fc, M = M, share = mean(!is.na(v)),
                                           fdp_given_flag = if (any(!is.na(v))) mean(v, na.rm = TRUE) else NA_real_,
                                           fdp_unconditional = mean(ifelse(is.na(v), 0, v)))
  }
}
res <- do.call(rbind, rows)

mx <- do.call(rbind, lapply(origins, function(o) {
  ser <- cmp$series[cmp$t0 == o]
  data.frame(origin = as.character(o), series = ser,
             max_log_e_F3_no_shift = vapply(ser, function(s) max(LE[[paste(s, o, "none", "F3", sep = "|")]]), numeric(1)))
}))

lines <- c("Number of origins: ", length(origins), "Series per origin:", capture.output(print(table(res$M[res$forecaster == "F3"]))),
           "", "Average over origins of the false discovery proportion at the first flag, zero when nothing is flagged:",
           capture.output(print(tapply(res$fdp_unconditional, res$forecaster, mean))),
           "", "Average over origins of the proportion among instances that flag something:",
           capture.output(print(tapply(res$fdp_given_flag, res$forecaster, mean, na.rm = TRUE))),
           "", "By origin:", capture.output(print(res, row.names = FALSE, digits = 3)),
           "", "Largest log e-value of F3 without an injected shift, by origin and series:", capture.output(print(mx, row.names = FALSE, digits = 3)))
writeLines(lines, file.path(OUT, "posthoc.txt"))
cat(lines, sep = "\n")
