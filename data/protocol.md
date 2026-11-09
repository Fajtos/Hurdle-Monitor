# Analysis protocol

This file records every analysis choice and says **when** each one was fixed. It has four parts, because the work was done in stages and the protocol must not claim more than is true.

- **Part A: simulation.** Fixed before the 400-replication runs. Matches the code in `simulation/05_simulation.R` and `simulation/06_stress_test.R`.
- **Part B: FAERS application.** Choices made in the application scripts before the results were read.
- **Part C: checks added after results were seen.** Labelled as such. They are not pre-specified.
- **Part D: Track A comparators.** To be frozen by a dated commit before they are run.

## Part A. Simulation (pre-specified)

| Item | Value | Where in code |
|---|---|---|
| Replications per scenario | 400 (main), 400 (stress test) | `N_REP` |
| Training window | 2000 days | `CONFIG$n_train` |
| Monitoring window | 3000 days | `CONFIG$T_monitor` |
| Shift time | day 1000 of monitoring | `CONFIG$shift_time` |
| Level alpha | 0.05 | `ALPHA` |
| Per-head threshold | 2/alpha | `R/eprocess.R` |
| Whole-forecast threshold | 1/alpha | `R/eprocess.R` |
| Margins | strict (0) and tolerant (0.02 occurrence, 0.03 magnitude) | `CONFIG$margins` |
| Bet grid | plus and minus 0.05, 0.1, 0.2, 0.4 | `BETS` in `R/eprocess.R` |
| Forecaster | estimated hurdle-NB; "oracle" (true parameters) run alongside | `run_one_replication` |

Main-simulation scenarios (all in `simulation/05_simulation.R`, `SCENARIOS`):

| Name | Occurrence shift | Magnitude shift | Dispersion theta |
|---|---|---|---|
| null | 0 | 0 | 5 |
| occ_pos | +0.6 | 0 | 5 |
| occ_neg | -0.6 | 0 | 5 |
| mag_moderate | 0 | 0.22 | 5 |
| mag_large | 0 | 0.50 | 5 |
| dispersion | 0 | 0 | 2.5 |
| both | +0.6 | 0.22 | 5 |

Seeds. Replication number r (1 to 400) uses seed r x 1000 + 7 in the main simulation and r x 1000 + 23 in the stress test. The randomised PIT for the positive part uses seed + 1; the full-forecast PIT in the stress test uses seed + 2. The seed stream r x 1000 + 7 is reused in the sensitivity study.

Status: the replication settings above were fixed before the 400-replication runs. The earlier version of this file mentioned a `TEST_MODE`; there is no such switch in the code, so that sentence has been removed.

## Part B. FAERS application

Data
- Source: openFDA drug-event endpoint (FAERS), daily report counts by `receivedate`, retrieved by `application/01_retrieve_faers.R`.
- Retrieved 2026-11-05. SHA-256 checksums in `data/faers/checksums.txt`; if a re-download changes them, the upstream data or the panel definition changed.
- Window: 2012-01-01 to 2024-12-31 (4749 days); days with no report are zero.
- Series (`data/panel.csv`): r2 metformin / LACTIC ACIDOSIS; r3 lisinopril / COUGH; r4 atorvastatin / MYALGIA; r6 isotretinoin / DEPRESSION; r8 levothyroxine / ALOPECIA.

Design
- Training window ends 2017-12-31. The forecast is fitted once on it and then frozen. Monitoring starts 2018-01-01.
- Forecaster: hurdle with negative binomial positive part (`pscl::hurdle`, dist = "negbin").
- Covariates known before day t: sin and cos of day of year, weekend indicator, log(1 + y[t-1]), log(1 + mean of days t-7 to t-1) (`build_features` in `R/helpers.R`).
- Ladder for comparison: OLS, Poisson, negative binomial, hurdle-NB.
- Monitor: alpha 0.05, head flag at 2/alpha, whole-forecast alarm at 1/alpha, tolerance (0.02, 0.03).
- Case study: r2 (`application/02_case_study.R`). Panel: all five series (`application/03_panel.R`).
- Injection check (`application/04_injection.R`): from 2019-06-01, drop each positive day to zero with probability 0.4 (occurrence failure), or add Poisson(0.5 x count) to each positive count (magnitude failure). Checks that the right head alarms.
- Seeds: PIT seed 2025 in the panel and refit scripts (2024 for the residual diagnostics); case-study set.seed values 2024, 12345, 2025 are in `application/02_case_study.R`.

Series selection (screening step, done before the panel was fetched)
- Eight candidate drug/reaction pairs (r1 to r8) were screened with `quick_check.py` (archived in `application/archive/quick_check.py`; a later R version, `application/00_quick_check.R`, was removed). Its raw output is saved in `data/faers/screening_output.txt`. Run date: 2026-11-04 (author's recollection; not independently timestamped), the day before the panel retrieval of 2026-11-05.
- Rule of thumb applied: keep a series if, in the **training** window, the weekday zero share is in [0.15, 0.85] and the mean positive count is at least 1.3 (so positive days are not all equal to 1).
- Result:

| Id | Series | Training weekday zero share | Mean positive count | Decision |
|---|---|---|---|---|
| r1 | valsartan (all reactions) | 0.00 | 35.88 | dropped: too dense, no zeros |
| r2 | metformin / LACTIC ACIDOSIS | 0.22 | 3.07 | kept |
| r3 | lisinopril / COUGH | 0.19 | 2.67 | kept |
| r4 | atorvastatin / MYALGIA | 0.22 | 2.27 | kept |
| r5 | semaglutide (all reactions) | 0.97 | 1.00 | dropped: too sparse, all positives equal 1 |
| r6 | isotretinoin / DEPRESSION | 0.56 | 2.95 | kept |
| r7 | ranitidine / PROSTATE CANCER | 0.98 | 1.00 | dropped: too sparse, all positives equal 1 |
| r8 | levothyroxine / ALOPECIA | 0.40 | 2.39 | kept |

- Limitation, stated openly: the screening output also printed yearly means for 2012 to 2024, so the monitoring-period behaviour of each series was visible when the panel was chosen, and the `note` column in `data/panel.csv` mentions later changes. The rule itself uses training-window statistics only, but the panel should be read as a screened, hand-assembled panel, not a random or blind sample of series.

Reference events
- No dated regulatory events have been assigned to the series yet. Until they are, every series is used only to report alarm dates and false-alarm behaviour, not detection delay against a known event. Dated events (with source links) will be added to the table below **before** any comparator in Part D is run.

| Series | Event | Date | Source (URL) | Expected head |
|---|---|---|---|---|
| (to be filled) | | | | |

## Part C. Added after results were seen (exploratory)

These were added in response to review, after the main results existed. They are reported as additional checks, not as pre-specified tests.
- `application/06_refit_panel.R`: head-specific refit on series r2, r3, r6. From the first alarm date, four arms (no refit, refit the flagged head, refit the other head, refit both) are monitored on the post-alarm window with fresh e-processes. Refit window is the last 500 days up to the alarm.
- `application/07_extra_checks.R`: Vuong test (hurdle-NB versus NB, training window) and a whole-forecast test with three margins (strict, 0.03, and 0.03 times the occurrence probability squared).
- `simulation/07_sensitivity.R`: null false-alarm rate against training length (four values) and monitoring length (four values), 400 replications each.

## Part D. Track A comparators

Frozen here so that the comparison is not tuned on results. The event study (Part D2) compares, on the same frozen forecast and windows:
1. The head-wise monitor (occurrence head, magnitude head, either head).
2. A one-part monitor (whole-forecast PIT only, no heads).
3. CUSUM on standardized residuals (k = 0.5).
4. EWMA on standardized residuals (lambda = 0.1).
5. One-sided Poisson-type MaxSPRT.

Outcomes: false-alarm rate in quiet periods and placebo windows, detection rate and median delay after dated events, and head-attribution accuracy.

Refit strategies (never refit, rolling refit, head-specific refit) and the sequential score-difference gate are not part of this run. They belong to separate analyses that will be pre-specified in their own protocol before any such run.

## Part D2. Event-window design (decided 2026-11-08, before any event-based run; final, no further edits planned)

The Part B results keep the fixed split (training to 2017-12-31). The event study uses one window per event so that every event has its own quiet period and none falls in a shared control year.

For an event with date T_e (all dates inclusive):
- **Training:** the hurdle-NB is fitted on all data from 2012-01-01 to T_e minus 396 days and never refitted. Minimum training length 1500 days.
- **Monitored period:** starts at T_e minus 395 days and ends at T_e plus 179 days. The e-detector replaces the plain e-process used in Part B for the event study only. Features use only past days (lag and 7-day roll use days up to t-1).
- **Quiet period:** the first 365 monitored days (T_e minus 395 to T_e minus 31). Any head or whole-forecast alarm here is a false alarm.
- **Buffer:** T_e minus 30 to T_e minus 1 (30 days). An alarm here is reported as a buffer alarm and counted as neither a detection nor a false alarm, because reporting can move slightly before a public communication.
- **Detection window:** T_e to T_e plus 179 days (180 days). The first alarm here is a detection and the delay is the number of days from T_e. A first alarm that occurred earlier (quiet or buffer) is classified by that earlier alarm. No alarm anywhere is a miss, right-censored at 180 days.
- **Alarm rules (e-detector):** each head (occurrence, magnitude) and the whole forecast is monitored with an e-detector in the sense of Saha and Ramdas: D_t is the sum over start days s of the betting mixture restarted at s (bets, features and tolerances as in `R/eprocess.R`: 0.02 occurrence, 0.03 magnitude and whole). A plain e-process was checked in a simulation on 2026-11-08 (occurrence shift of +1.5 on the logit after 395 quiet days) and lost power after a long quiet period (detected 2 of 20 simulated shifts, against 20 of 20 for the e-detector). An alarm is D_t at or above 2n/0.05 for a head and n/0.05 for the whole forecast, where n is the number of monitored days (575); by the e-detector average-run-length bound this keeps the chance of a false alarm over the monitored period near or below 0.05 for each head. "Any head" is the earlier crossing of the two heads.
- **Comparators** (CUSUM k=0.5, EWMA lambda=0.1, one-sided Poisson-type MaxSPRT): same frozen forecast, same windows. Each threshold is the 95th percentile of the maximum of its statistic over the whole monitored period, from B=300 series simulated from the frozen forecast with the observed covariates. This gives all methods the same nominal 5% false-alarm rate per series.
- **Randomization:** randomized PIT seeds fixed (2025 positive part, 2024 whole forecast); comparator simulation seed derived from the run identifier by `run_seed()`.
- **Placebo windows:** for each event that ran, the same calendar windows are applied to each panel series in `data/panel.csv` (and every other series that passes the screen) that is not the event's own drug and reaction. Any alarm there is a false alarm. These series were not individually checked for events in the window, so the placebo false-alarm rate is a conservative (upper) estimate.
- **Expected head:** fixed before any run and recorded in `events/events.csv`. A communication that adds a new reaction is expected to raise the occurrence head first; one that changes dose or use is expected to raise the magnitude head first. Head attribution counts a detection as correct if the first head to cross equals the expected head; same-day crossings are not counted correct.
- **Series screen:** applied to each event's own training window: weekday zero share between 0.15 and 0.85, more than 30 positive weekdays, mean positive weekday count at least 1.3, training at least 1500 days, and the follow-up inside the data (to 2024-12-31). Events failing it are excluded and listed with the reason.
- **Fewer than 10 events passing:** the study is reported as a small case series with every candidate listed and its outcome. There is no change of design and no switch to weekly counts.
- **Events:** the list is `events/events.csv`, frozen with this protocol. One event per drug-reaction pair, first communication only. Overlapping events sharing drugs (for example the two fluoroquinolone rows) are reported with a note and not treated as independent in any summary rate beyond the descriptive table.
- **What counts as verified:** date, URL and preferred term checked by the author against the FDA page and openFDA before the freeze. Only verified rows are run. Rows are never added or dropped after the freeze; a mistake is logged in the deviations log below.

## Deviations log

| Date | What changed | Why |
|---|---|---|
| (none yet) | | |
