# FAERS panel data

Retrieved: 2026-11-05
Source: openFDA drug-event endpoint, https://api.fda.gov/drug/event.json

## Files

- r2.csv
- r3.csv
- r4.csv
- r6.csv
- r8.csv

## Columns

- `date`: calendar date, daily from 2012-01-01 to 2024-12-31 (4749 rows).
- `count`: integer count of adverse-event reports for the series on that date.

## Panel definition

Each series is defined by the generic drug name and, where present, a
MedDRA reaction term. See `data/panel.csv` for the exact definitions.

## Retrieval protocol

The series were retrieved by `application/01_retrieve_faers.R` using the
openFDA drug-event endpoint with `count=receivedate`, one request per
two-year window, and days without any reports filled as zero. Retrieval
used an authenticated api.data.gov key stored outside the repository.

## Verification

SHA-256 checksums of the CSVs are in `checksums.txt`. If the retrieval is
re-run and the checksums differ, the upstream data have changed or the
panel definition has been modified.
