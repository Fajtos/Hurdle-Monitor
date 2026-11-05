# Pre-specified simulation protocol

- Simulation scenarios are defined in 05_simulation.R and 06_stress_test.R.
- Replications per scenario: 400 for the main simulation, 400 for the stress test.
- Training window: 2000 days.
- Monitoring window: 3000 days.
- Shift time: 1000 days after the start of monitoring.
- E-process margins: strict (0, 0) and tolerant (0.02, 0.03).
- Bet grid: ±{0.05, 0.1, 0.2, 0.4}.
- Feature grid: bias, slope, serial dependence.
- Alpha: 0.05.
- Per-head threshold: 2/alpha. Overall threshold: 1/alpha.
- RNG seeds: id × 1000 + 7 (main simulation), id × 1000 + 23 (stress test).
- PIT randomisation uses seeds id + 1 and id + 2.

The protocol was fixed before any replication was run. The test
configuration in TEST_MODE is for pipeline verification only and is
not used in the paper.