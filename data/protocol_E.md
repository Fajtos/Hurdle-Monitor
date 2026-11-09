# Protocol E. False alarm calibration of forecaster choices

Decided 2026-11-09, before any run of `application/12_false_alarm_study.R`. This file is separate from `data/protocol.md` and does not change it.

## Question

In the confirmatory event study (Part D2) the monitor raised false alarms on real series at rates above its nominal level. The monthly comparison of forecast and data (`out/events/exploratory/`, run after the confirmatory study) shows a persistent level offset in the occurrence head from the first monitored months in several series. This part asks whether forecasters that update over time reduce false alarms on real series, using many windows and no events.

## Series and windows

- Series: every file in `data/faers/` (panel series r2, r3, r4, r6, r8) and in `data/faers_events/` (event series E04 to E20 that were downloaded), each used as one series regardless of its event.
- Origins: for a series, the origin t0 runs over 2015-12-31 and every 182 days after it. The monitored period is t0 + 1 to t0 + 365 (n = 365 days).
- Event series: only windows whose monitored period ends before the event date minus 30 days are used, so that no window contains the event or its buffer. Panel series have no restriction, and their windows may contain unrecorded events, which makes false alarm rates conservative.
- Screen: a window is analysed only if, on the training data up to t0, the weekday zero share lies in [0.15, 0.85], there are more than 30 positive weekdays, the mean positive weekday count is at least 1.3, and the training data have at least 1500 days.
- A window is excluded, and counted, if any of the four forecasters fails to fit.

## Forecasters

All use the hurdle negative binomial model of Part B with the same five covariates, built from past days only. Each forecast for day t uses data up to a refit date before t, so every forecast is determined by information available before day t.

| Name | Training data | Refit |
|---|---|---|
| F0 frozen | all data up to t0 | never |
| F1 recent | the last 730 days up to t0 | never |
| F2 expanding | all data up to the refit date | every 30 days from t0 |
| F3 rolling | the last 730 days up to the refit date | every 30 days from t0 |

No other forecaster is tried in this part.

## Monitor

The e-detector of Part D2, with bets, features, tolerances (0.02 occurrence, 0.03 magnitude and whole) and PIT seeds (2025 positive part, 2024 whole forecast) unchanged. Thresholds use n = 365: log(2n/0.05) for each head and log(n/0.05) for the whole forecast. Any alarm in the monitored period is a false alarm.

## Outcomes and analysis

- For each forecaster: the proportion of windows with an alarm on the occurrence head, on the magnitude head, on either head, and on the whole forecast.
- Uncertainty: 95% percentile intervals from a bootstrap that resamples series with replacement (2000 resamples, seed 20261109), because windows of the same series overlap.
- Reference levels: 0.025 per head and 0.05 for either head and for the whole forecast.

## Decision rule

A forecaster F1, F2 or F3 is carried forward to a separate detection study (a later protocol, with new events chosen by a rule fixed beforehand) if its false alarm proportion for either head is at most 0.10 and at most half the proportion for F0. If several qualify, the one with the lowest proportion is carried forward, and among proportions within 0.01 the order of preference is F1, F2, F3. If none qualifies, the result is reported as a negative result, and no further forecaster is tried in this part.

## Not included

Detection power, the nine events of the confirmatory study, and any change to thresholds, windows or alarm rules. Any such change requires a new protocol.

## Deviations log

| Date | What changed | Why |
|---|---|---|
| (none yet) | | |
