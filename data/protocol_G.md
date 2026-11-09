# Protocol G. Comparison with the rank based PIT monitor of Farran (2026)

Decided 2026-11-09, after Parts E and F and before any run of `application/14_pit_monitor_comparison.R`. The protocol and the script are committed and tagged `protocol-G-frozen` before the first run, including the smoke test. This file does not change the earlier protocols.

## Question

The closest published monitor of forecast calibration is PITMonitor (Farran, arXiv 2603.13156). It converts the probability integral transform of each observation into a conformal p-value by ranking it against the earlier transforms of the same stream, bets with a histogram density, and mixes over change points. It monitors the whole forecast, and its null is that the transforms are exchangeable, so a forecast that is stably miscalibrated does not raise an alarm. The headwise e-detector of this project tests absolute calibration of each head. This part compares the two on the same windows, forecasts and injected shifts as Parts E and F, and reports false alarms, detection and attribution.

## Data

The 141 windows of Part E (`out/events/falsealarm/windows.csv`), the six scenarios of Part F (none, counts times 1.5, counts doubled, counts halved, 30% of positive days set to zero, 30% of zero weekdays set to one), the same injection seeds, and the forecasters F0 and F3. The e-detector alarms recomputed here must equal those of `out/events/power/jobs.csv`, and the number of agreeing rows is reported.

## Monitors

- E-detector, either head and whole forecast: exactly as in Parts D2, E and F.
- PITMonitor, whole forecast: the randomized probability integral transform of the whole hurdle forecast (the one used for the whole forecast e-detector, seed 2024) is ranked against the earlier values in the monitored period only, p_t = (R_t - 1 + V_t) / t with V_t uniform (seed 31). The e-value of day t is the histogram density of the p-values seen before day t, evaluated at p_t, with 100 bins and a pseudocount of 1 per bin. The statistic is M_t = e_t (M_{t-1} + 1 / (t (t + 1))) with M_0 = 0, and the alarm is the first day with M_t at least 1 / 0.05.
- PIT rank monitor, either head: the same construction applied separately to the randomized transform of the occurrence head (all days, seed 32) and of the positive part (positive days only, seed 33), each with threshold 2 / 0.05, so that the two heads share the level 0.05 by the union bound. This monitor is a variant built for this comparison and is not part of the published method.

## Outcomes

Per scenario, forecaster and monitor: the proportion of windows with an alarm on any day, with an alarm before day 181, and with a first alarm on or after day 181 (correct detection). The false alarm rate of a monitor is its proportion of windows with an alarm on any day in the scenario without a shift.

Attribution is assessed on the scenarios that act on one head: counts times 1.5 and counts doubled (magnitude) and positive days set to zero (occurrence). Among windows with a correct detection, it is the number in which the target head alarms first or on the same day as the other head. It applies to the two headwise monitors only, because the whole forecast monitors do not name a head.

## Analysis and decision rule

- Pooled correct detection is the mean over the five shift scenarios. Net correct detection subtracts the proportion of windows with an alarm on or after day 181 in the scenario without a shift. Intervals are 95% percentile intervals from a bootstrap that resamples series (2000 resamples, seed 20261111).
- For each forecaster, the headwise e-detector is competitive with PITMonitor if its pooled net correct detection is at least 0.8 times that of PITMonitor. Otherwise it is not competitive. The ratio is reported in both cases.
- The monitors test different null hypotheses, and this is stated in the report: PITMonitor tests that the transforms are exchangeable, and the e-detector tests that the forecast is calibrated. A monitor with a lower false alarm rate on real series is not thereby better, because they do not answer the same question.

## Not included

Other implementations or settings of PITMonitor (bins, weights), other forecasters, other shifts, new events. Each requires a new protocol.

## Deviations log

| Date | What changed | Why |
| :- | :- | :- |
| (none yet) | | |
