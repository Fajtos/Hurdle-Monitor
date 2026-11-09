# Event study (Track A)

Head-wise e-detector monitor and comparators on FAERS daily report counts around dated FDA safety communications. The design is fixed in `data/protocol.md` (Parts D and D2). The event list is `events/events.csv`.

## Layout

| Path | Content |
|---|---|
| `R/event_study.R` | windows, series screen, forecast, e-detector, comparators, classification, openFDA download |
| `application/07_event_selftest.R` | check of the code on simulated series, no network |
| `application/08_event_fetch.R` | download of the daily series for the events marked `verified = yes` |
| `application/09_event_screen.R` | series screen per event, on that event's own training window |
| `application/10_event_study.R` | event runs and placebo runs, tables and traces |
| `events/events.csv` | events: drug names, reaction term, date, expected head, source, verification flag |
| `data/faers_events/` | downloaded series, `checksums.txt`, `retrieved.txt` |
| `out/events/` | results |

## Requirements

R with the packages in `renv.lock` (`renv::restore()`), and an openFDA API key stored as `OPENFDA_API_KEY` in `.Renviron`. The key file is not committed.

## Run order

All commands are run from the project root. In an R console use `source("application/NN_name.R")`, in a terminal use `Rscript application/NN_name.R`.

1. `07_event_selftest.R`
2. `08_event_fetch.R`
3. `09_event_screen.R`
4. `10_event_study.R`

## Outputs in `out/events/`

- `tables/event_screen.csv`: screen result and reason per event.
- `tables/event_runs.csv`: one row per run with the outcome class and date for every method.
- `tables/event_summary.csv` and `.tex`: false alarm, buffer alarm, detected, missed and median delay by method, for events and placebo windows.
- `tables/head_attribution.csv`: detections in which the first head to cross matched the expected head.
- `raw/event_traces.rds`: monitor and comparator paths with thresholds.
- `session_info.txt` and `constants.csv`: package versions and constants.

## Conventions

- Dates in `events/events.csv` use the form `YYYY-MM-DD`.
- A reaction term ending in `*` is used as an unquoted wildcard in the openFDA query.
- Drugs separated by `|` are summed per day.
- Events with fewer than 10 passing the screen are reported as a case series with every candidate listed.
