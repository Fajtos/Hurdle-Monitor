# Helper functions for the event-window study (data/protocol.md, Part D and D2).
# Needs R/io.R, R/eprocess.R, R/helpers.R and R/forecast.R to be sourced first.
# All design constants live here so that nothing is typed in two places.

EV         <- list(quiet = 365L, buffer = 30L, detect = 180L, min_train = 1500L)
EV_ALPHA   <- 0.05
EV_TOL     <- c(0.02, 0.03)                     # (occurrence, magnitude/whole) margins, as in the panel analysis
# e-detector thresholds depend on the monitored length n: log(2n/alpha) per head
EV_THR_HEAD  <- function(n) log(2 * n / EV_ALPHA)
EV_THR_WHOLE <- function(n) log(n / EV_ALPHA)            # whole forecast: log(n/alpha)
EV_FORM    <- count ~ sin_year + cos_year + weekend + y_lag + y_roll7
EV_PIT_SEEDS <- c(pos = 2025L, whole = 2024L)   # same PIT seeds as the panel analysis
EV_CALIB_B <- 300L                              # null simulations used to calibrate comparator thresholds
EV_SCREEN  <- list(zero_lo = 0.15, zero_hi = 0.85, min_pos_mean = 1.3, min_pos_days = 30L)
EV_METHODS <- c("occ", "mag", "whole", "any_head", "cusum", "ewma", "maxsprt")

# data
read_series <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$date <- as.Date(d$date)
  d <- d[order(d$date), ]
  d$count <- as.integer(d$count)
  if (anyNA(d$count) || anyNA(d$date)) stop("Missing values in ", path)
  if (any(diff(as.integer(d$date)) != 1L)) stop("Dates are not a complete daily sequence in ", path)
  d
}

# windows (Part D2)
# First monitored day is T_e - 395. The quiet period is the first 365 monitored days,
# the buffer is the 30 days before T_e, the detection window is T_e .. T_e + 179.
event_windows <- function(event_date, ev = EV) {
  te <- as.Date(event_date)
  mon_start <- te - (ev$quiet + ev$buffer)
  list(event = te,
       train_end = mon_start - 1,
       quiet  = c(mon_start, mon_start + ev$quiet - 1),
       buffer = c(te - ev$buffer, te - 1),
       detect = c(te, te + ev$detect - 1),
       mon_end = te + ev$detect - 1)
}

screen_event <- function(d, w) {
  out <- list(n_train = sum(d$date <= w$train_end), zero_share_wd = NA_real_, mean_pos_wd = NA_real_,
              pass = FALSE, reason = "")
  if (w$mon_end > max(d$date)) { out$reason <- "follow-up extends beyond the data"; return(out) }
  if (out$n_train < EV$min_train) { out$reason <- "training shorter than 1500 days"; return(out) }
  tr  <- d[d$date <= w$train_end, ]
  wd  <- tr$count[as.POSIXlt(tr$date)$wday %in% 1:5]
  pos <- wd[wd > 0]
  out$zero_share_wd <- mean(wd == 0)
  out$mean_pos_wd   <- if (length(pos) > 0) mean(pos) else NA_real_
  ok <- out$zero_share_wd >= EV_SCREEN$zero_lo && out$zero_share_wd <= EV_SCREEN$zero_hi &&
        length(pos) > EV_SCREEN$min_pos_days && !is.na(out$mean_pos_wd) && out$mean_pos_wd >= EV_SCREEN$min_pos_mean
  out$pass   <- isTRUE(ok)
  out$reason <- if (isTRUE(ok)) "pass" else "screen rule not met"
  out
}

# forecast (frozen at the end of training)
make_forecast <- function(d, w) {
  d  <- build_features(d)                       # lag features use the previous day only
  tr <- d[d$date <= w$train_end, ]
  te <- d[d$date >= w$quiet[1] & d$date <= w$mon_end, ]
  fit <- pscl::hurdle(EV_FORM, data = tr, dist = "negbin")
  X <- model.matrix(EV_FORM, te)
  list(te = te,
       pi = as.numeric(plogis(X %*% fit$coefficients$zero)),
       mu = as.numeric(exp(X %*% fit$coefficients$count)),
       theta = as.numeric(fit$theta))
}

# the head-wise monitor
# Per-day log-increments of every betting component (same bets and features as R/eprocess.R).
occ_steps <- function(pi_hat, Z, delta) {
  T <- length(Z); fl <- list(rep(1, T), 2 * pi_hat - 1, c(0, 2 * Z[-T] - 1))
  do.call(cbind, lapply(fl, function(f) { g <- f * (Z - pi_hat)
    sapply(BETS, function(b) log1p(b * (g - sign(b) * delta))) }))
}
pit_steps <- function(U, active, delta) {
  T <- length(U); U_last <- c(0.5, U[-T])
  fl <- list(2 * U - 1, (12 * (U - 0.5)^2 - 1) / 2, (2 * U - 1) * (2 * U_last - 1))
  do.call(cbind, lapply(fl, function(f) { g <- active * f
    sapply(BETS, function(b) log1p(b * (g - sign(b) * delta * active))) }))
}
# Mixture e-process without restarts (equals eprocess_occurrence / eprocess_pit) and the e-detector:
# log D_t = log sum_{s<=t} E_{s:t}, each E_{s:t} the same mixture started fresh at day s (Saha-Ramdas).
mix_log <- function(steps) { C <- apply(steps, 2, cumsum); apply(C, 1, log_sum_exp) - log(ncol(C)) }
edetector_log <- function(steps) {
  C <- apply(steps, 2, cumsum); C0 <- rbind(0, C); K <- ncol(C); T <- nrow(C); out <- numeric(T)
  for (t in seq_len(T)) {
    M <- sweep(C0[1:t, , drop = FALSE], 2, C[t, ], function(a, b) b - a)
    m <- apply(M, 1, max); out[t] <- log_sum_exp(m + log(rowSums(exp(M - m))) - log(K))
  }
  out
}
monitor_heads <- function(fc) {
  y <- fc$te$count; n <- length(y); Z <- as.integer(y > 0)
  U_pos  <- randomized_pit_positive(y, fc$mu, fc$theta, seed = EV_PIT_SEEDS[["pos"]])
  U_full <- pnorm(rqr_hurdle_vec(y, fc$pi, fc$mu, fc$theta, rng_seed = EV_PIT_SEEDS[["whole"]]))
  list(occ   = edetector_log(occ_steps(fc$pi, Z, EV_TOL[1])),
       mag   = edetector_log(pit_steps(U_pos, Z, EV_TOL[2])),
       whole = edetector_log(pit_steps(U_full, rep(1, n), EV_TOL[2])))
}

# comparators (calibrated by simulation under the frozen model)
pred_moments <- function(pi, mu, theta) {
  q0 <- (theta / (theta + mu))^theta
  m1 <- mu / (1 - q0)                           # E[Y | Y > 0]
  m2 <- ((mu + mu^2 / theta) + mu^2) / (1 - q0) # E[Y^2 | Y > 0]
  mean <- pi * m1
  list(mean = mean, var = pmax(pi * m2 - mean^2, 1e-8))
}
cusum_stat <- function(z, k = 0.5) {
  n <- length(z); out <- numeric(n); up <- 0; dn <- 0
  for (t in seq_len(n)) {
    up <- max(0, up + z[t] - k); dn <- max(0, dn - z[t] - k)
    out[t] <- max(up, dn)
  }
  out
}
ewma_stat <- function(z, lambda = 0.1) {
  n <- length(z); out <- numeric(n); s <- 0
  for (t in seq_len(n)) { s <- (1 - lambda) * s + lambda * z[t]; out[t] <- abs(s) }
  out / sqrt(lambda / (2 - lambda))
}
maxsprt_stat <- function(y, expected) {         # one-sided Poisson-type MaxSPRT on cumulative observed vs expected
  C <- cumsum(y); E <- cumsum(expected)
  ifelse(C > E, C * log(C / E) - (C - E), 0)
}
simulate_null <- function(pi, mu, theta) {      # draws from the frozen forecast, conditional on the observed covariates
  n <- length(pi); q0 <- (theta / (theta + mu))^theta
  y_pos <- qnbinom(q0 + runif(n) * (1 - q0), size = theta, mu = mu)
  as.integer(ifelse(runif(n) < pi, y_pos, 0L))
}
comparator_stats <- function(fc, seed, B = EV_CALIB_B) {
  mom <- pred_moments(fc$pi, fc$mu, fc$theta)
  sdv <- sqrt(mom$var)
  y <- fc$te$count
  obs <- list(cusum = cusum_stat((y - mom$mean) / sdv), ewma = ewma_stat((y - mom$mean) / sdv),
              maxsprt = maxsprt_stat(y, mom$mean))
  set.seed(seed)
  mx <- matrix(NA_real_, B, 3)
  for (b in seq_len(B)) {
    ys <- simulate_null(fc$pi, fc$mu, fc$theta)
    mx[b, 1] <- max(cusum_stat((ys - mom$mean) / sdv))
    mx[b, 2] <- max(ewma_stat((ys - mom$mean) / sdv))
    mx[b, 3] <- max(maxsprt_stat(ys, mom$mean))
  }
  thr <- apply(mx, 2, quantile, probs = 1 - EV_ALPHA, names = FALSE)
  list(stat = obs, thr = c(cusum = thr[1], ewma = thr[2], maxsprt = thr[3]))
}

# classification of the first crossing
first_idx <- function(stat, thr) { i <- which(stat >= thr); if (length(i) > 0) i[1] else NA_integer_ }
classify_idx <- function(i, dates, w) {
  if (is.na(i)) return(list(cls = "none", date = NA_character_, delay = NA_integer_))
  d <- dates[i]
  cls <- if (d <= w$quiet[2]) "false_alarm" else if (d < w$detect[1]) "buffer_alarm" else "detected"
  list(cls = cls, date = as.character(d), delay = if (cls == "detected") as.integer(d - w$detect[1]) else NA_integer_)
}

# one event (or placebo) window
run_one <- function(d, event_date, seed) {
  w  <- event_windows(event_date)
  sc <- screen_event(d, w)
  row <- list(event_date = as.character(w$event), n_train = sc$n_train,
              zero_share_wd = sc$zero_share_wd, mean_pos_wd = sc$mean_pos_wd,
              screen = sc$reason, status = if (sc$pass) "run" else "screened_out", head_first = NA_character_)
  for (m in EV_METHODS) { row[[paste0(m, "_class")]] <- NA_character_; row[[paste0(m, "_date")]] <- NA_character_; row[[paste0(m, "_delay")]] <- NA_integer_ }
  if (!sc$pass) return(list(row = row, trace = NULL))
  fc <- tryCatch(make_forecast(d, w), error = function(e) { message("  hurdle fit failed: ", conditionMessage(e)); NULL })
  if (is.null(fc)) { row$status <- "fit_failed"; return(list(row = row, trace = NULL)) }
  dates <- fc$te$date
  hd <- monitor_heads(fc)
  cp <- comparator_stats(fc, seed)
  idx <- list(occ = first_idx(hd$occ, EV_THR_HEAD(length(dates))), mag = first_idx(hd$mag, EV_THR_HEAD(length(dates))),
              whole = first_idx(hd$whole, EV_THR_WHOLE(length(dates))),
              cusum = first_idx(cp$stat$cusum, cp$thr[["cusum"]]),
              ewma = first_idx(cp$stat$ewma, cp$thr[["ewma"]]),
              maxsprt = first_idx(cp$stat$maxsprt, cp$thr[["maxsprt"]]))
  idx$any_head <- if (all(is.na(c(idx$occ, idx$mag)))) NA_integer_ else min(c(idx$occ, idx$mag), na.rm = TRUE)
  for (m in EV_METHODS) {
    r <- classify_idx(idx[[m]], dates, w)
    row[[paste0(m, "_class")]] <- r$cls; row[[paste0(m, "_date")]] <- r$date; row[[paste0(m, "_delay")]] <- r$delay
  }
  if (!is.na(idx$any_head)) {
    row$head_first <- if (is.na(idx$occ)) "mag" else if (is.na(idx$mag)) "occ" else if (idx$occ < idx$mag) "occ" else if (idx$mag < idx$occ) "mag" else "same_day"
  }
  trace <- data.frame(date = dates, occ = hd$occ, mag = hd$mag, whole = hd$whole,
                      cusum = cp$stat$cusum, ewma = cp$stat$ewma, maxsprt = cp$stat$maxsprt)
  attr(trace, "thresholds") <- c(head = EV_THR_HEAD(length(dates)), whole = EV_THR_WHOLE(length(dates)), cp$thr)
  list(row = row, trace = trace)
}

run_seed <- function(run_id) {
  v <- utf8ToInt(run_id)
  as.integer(sum(as.numeric(v) * seq_along(v)) %% 1000000)
}

# summary tables
summarise_runs <- function(runs) {
  runs <- runs[runs$status == "run", , drop = FALSE]
  out <- list()
  for (tp in unique(runs$type)) for (m in EV_METHODS) {
    sub <- runs[runs$type == tp, , drop = FALSE]
    cl <- sub[[paste0(m, "_class")]]; dl <- sub[[paste0(m, "_delay")]]
    out[[length(out) + 1]] <- data.frame(
      type = tp, method = m, n = length(cl),
      false_alarm = mean(cl == "false_alarm"), buffer_alarm = mean(cl == "buffer_alarm"),
      detected = mean(cl == "detected"), missed = mean(cl == "none"),
      any_alarm = mean(cl != "none"),
      median_delay = if (any(!is.na(dl))) median(dl, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

# openFDA download (same query pattern as application/01_retrieve_faers.R)
fetch_faers_series <- function(drugs, reaction, start = 2012, end = 2024, api_key, sleep = 0.4, max_retries = 3) {
  esc <- function(x) gsub("'", "%27", gsub(" ", "%20", x, fixed = TRUE), fixed = TRUE)
  full <- data.frame(date = seq(as.Date(sprintf("%d-01-01", start)), as.Date(sprintf("%d-12-31", end)), by = "day"), count = 0L)
  for (drug in drugs) {
    for (y0 in seq(start, end, by = 2)) {
      y1 <- min(y0 + 1, end)
      parts <- sprintf('patient.drug.openfda.generic_name:"%s"', esc(drug))
      if (!is.na(reaction) && nzchar(reaction)) parts <- c(parts, if (grepl("\\*$", reaction)) sprintf("patient.reaction.reactionmeddrapt:%s", esc(reaction))   # wildcard term, unquoted
                         else sprintf('patient.reaction.reactionmeddrapt:"%s"', esc(reaction)))
      parts <- c(parts, sprintf("receivedate:[%d0101+TO+%d1231]", y0, y1))
      url <- sprintf("https://api.fda.gov/drug/event.json?search=%s&count=receivedate&limit=1000&api_key=%s", paste(parts, collapse = "+AND+"), api_key)
      resp <- NULL
      for (attempt in seq_len(max_retries)) {
        resp <- tryCatch(httr2::req_perform(httr2::request(url)),
                         httr2_http_404 = function(e) "none",       # openFDA answers 404 when nothing matches
                         error = function(e) { if (attempt < max_retries) Sys.sleep(2 * attempt); NULL })
        if (!is.null(resp)) break
      }
      if (is.null(resp)) stop("Request failed after retries: ", drug, " / ", reaction, " / ", y0)
      if (!identical(resp, "none")) {
        js <- httr2::resp_body_json(resp, simplifyVector = TRUE)
        if (!is.null(js$results) && is.data.frame(js$results) && nrow(js$results) > 0) {
          idx <- match(as.Date(js$results$time, format = "%Y%m%d"), full$date)
          ok <- !is.na(idx)
          full$count[idx[ok]] <- full$count[idx[ok]] + as.integer(js$results$count[ok])
        }
      }
      Sys.sleep(sleep)
    }
  }
  full
}

# events file reader (accepts 2019-04-30 or 4/30/2019, e.g. after Excel)
read_events <- function(path = "events/events.csv") {
  ev <- read.csv(path, stringsAsFactors = FALSE, strip.white = TRUE)
  ev$event_date <- as.Date(ev$event_date, tryFormats = c("%Y-%m-%d", "%m/%d/%Y"))
  if (any(is.na(ev$event_date))) stop("Some event_date values could not be read in ", path)
  ev
}
