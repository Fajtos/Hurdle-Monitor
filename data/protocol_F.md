# Protocol F. Detection power of the frozen and the rolling refit forecaster

Decided 2026-11-09, after the results of Part E and before any run of `application/13_detection_power_study.R`. This file is separate from `data/protocol.md` and `data/protocol_E.md` and does not change them.

## Question

Part E found that the rolling refit forecaster F3 (last 730 days, refit every 30 days) raised false alarms in 3.5% of real windows against 13.5% for the frozen forecaster F0. A forecaster that adapts can also absorb a real shift and so detect less. This part measures detection on real series with shifts of known size.

## Windows

The 141 series-windows of Part E (`out/events/falsealarm/windows.csv`), unchanged: origin t0, training data up to t0, monitored period t0 + 1 to t0 + 365. A job is one window and one scenario.

## Scenarios

A shift starts on monitored day 181 (date t0 + 181) and lasts to the end of the monitored period. It is applied to the raw daily counts, and both forecasters then see the shifted series, including in their features and, for F3, in their refits. Random draws use a seed computed from the series, the origin and the scenario.

| Scenario | Change to the counts from day 181 |
| :- | :- |
| none | none, used as a consistency check against Part E |
| mag_up15 | y plus a Binomial(y, 0.5) draw, so about 1.5 times the count |
| mag_up2 | 2 times the count |
| mag_down | each count replaced by a Binomial(y, 0.5) draw, so about 0.5 times the count |
| occ_drop | each positive day set to zero with probability 0.3 |
| occ_up | each zero weekday set to 1 with probability 0.3 |

No other scenario is tried in this part.

## Forecasters and monitor

F0 and F3 exactly as in Part E, with the e-detector, bets, tolerances, seeds and thresholds of Part D2 and Part E unchanged. F1 and F2 are not run here.

## Outcomes

- Correct detection: the first alarm of either head falls on day 181 or later. A window whose first alarm falls before day 181 is a failure for that window, because the monitor would have stopped.
- Alarm before onset: the first alarm of either head falls before day 181.
- Conditional power: correct detections divided by windows without an alarm before onset.
- Delay: days from day 181 to the first alarm, median among correct detections.
- The same quantities for the whole forecast are reported as secondary.

## Analysis and decision rule

- Per scenario and forecaster, all of the outcomes above.
- Pooled correct detection: the mean of the five shift scenarios, each with equal weight, for F0 and for F3. Intervals are 95% percentile intervals from a bootstrap that resamples series with replacement (2000 resamples, seed 20261110).
- Statement made in the paper: F3 retains power if the pooled correct detection of F3 is at least 0.8 times that of F0, and otherwise F3 has a power cost. The ratio, the difference and their intervals are reported in both cases.
- The scenario none must reproduce the alarms of Part E. The number of agreeing rows is reported.

## Not included

New events, any forecaster other than F0 and F3, shifts other than those listed, any change to thresholds, windows or alarm rules. Each requires a new protocol.

## Deviations log

| Date | What changed | Why |
| :- | :- |
| (none yet) | | |
