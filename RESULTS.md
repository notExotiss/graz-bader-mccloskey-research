# Bader--McCloskey stopping-rule revision

- Implemented: detrended log-residual tests (T4), calibrated CUSUM and Shiryaev--Roberts tests (T5), singleton-phase deferral, robust MAD dispersion, the requested community-level Eq. 6 BM-native stop, CM baselines, and a replay-based multi-cut report. CNM candidate ordering is unchanged.
- Completed in this run: the 10-seed ER calibration/sweep, all three small labelled graphs, and true no-budget CNM-full Amazon/DBLP runs. Amazon has 333,388 merges; DBLP has 313,936 merges. The local static build is unavailable on macOS (`crt0.o` is absent), but `g++ -O2 -std=c++17` builds and the new modes pass smoke tests.
- Amazon/DBLP full rows and plateau dumps are complete. Large-graph multicut stability runs remain unrun; the existing project evidence shows that full LiveJournal CNM is multi-day work and it has not been represented as CNM-full.
- CNM naming: `CNM-full` means the natural positive-dQ terminus with no merge budget. If `--max-merges=N` is used, the report prints `CNM-prefix`, the exact prefix length and budget; it never labels that row CNM-full.
- This is not evidence that the new rule beats CNM. On the completed small tests, deferral frequently recovers CNM's terminal cut, while unmodified trace tests and the literal adaptive BM-native threshold can still stop destructively early.

## ER nulls -- full rule comparison

T4 fits OLS to the preceding W log-dQ values and compares the current extrapolation residual to the fitted-residual dispersion. T5 applies CUSUM or Shiryaev--Roberts to that same standardized residual stream. The CUSUM threshold was selected once, before real-network evaluation: the smallest grid value giving mean ARL at least 90% of the CNM trace length over 10 seeds at N=1000, 2000, 5000, and 10000.

| h | mean ARL | mean available merges | stop_frac |
| --- | ---: | ---: | ---: |
| 2 | 30.725 | 4489.500 | 0.007 |
| 5 | 46.075 | 4489.500 | 0.010 |
| 10 | 77.300 | 4489.500 | 0.017 |
| 20 | 2009.100 | 4489.500 | 0.448 |
| 50 (selected) | 4233.375 | 4489.500 | 0.943 |
| 100 | 4489.500 | 4489.500 | 1.000 |
| 200 | 4489.500 | 4489.500 | 1.000 |

The complete emitted table contains every T4 `(W,k,sided)` combination and all `{none,defer,mad,defer+mad}` modifiers. Its key negative finding is structural: without deferral, T1--T5 tend to fire in the first few percent of an ER trace; with deferral, they usually run almost the complete trace and therefore inherit CNM's spurious ER communities. The `frac_nonsingleton` column makes that tradeoff explicit in the machine-readable runner output.

## Small real-network sweep

All rows below include `frac_nonsingleton`, and `stop_at / of_merges` gives the requested raw and fractional firing index.

| graph | method | stop_at | of_merges | stop_frac | K | Q | F1 | NMI | frac_nonsingleton |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| karate | CNM-full | 31 | 31 | 1.000 | 3 | 0.381 | 0.761 | 0.565 | 1.000 |
| karate | T5-cusum (h=50) | 31 | 31 | 1.000 | 3 | 0.381 | 0.761 | 0.565 | 1.000 |
| karate | BM-native-stop (tau=2) | 0 | 31 | 0.000 | 34 | -0.050 | 0.111 | 0.329 | 0.000 |
| polblogs | CNM-full | 1482 | 1482 | 1.000 | 8 | 0.354 | 0.617 | 0.866 | 1.000 |
| polblogs | T4 detrend, defer | 1480 | 1482 | 0.999 | 10 | 0.354 | 0.592 | 0.861 | 1.000 |
| polblogs | T5-cusum (h=50) | 1482 | 1482 | 1.000 | 8 | 0.354 | 0.617 | 0.866 | 1.000 |
| polblogs | BM-native-stop (tau=2) | 165 | 1482 | 0.111 | 1325 | 0.016 | 0.007 | 0.172 | 0.216 |
| email | CNM-full | 978 | 978 | 1.000 | 8 | 0.347 | 0.283 | 0.427 | 1.000 |
| email | T5-cusum (h=50) | 919 | 978 | 0.940 | 67 | 0.345 | 0.195 | 0.479 | 0.942 |
| email | BM-native-stop (tau=2) | 6 | 978 | 0.006 | 980 | -0.002 | 0.143 | 0.651 | 0.011 |

CM baseline with each `c in {1,5,10,50}` fired at cut 0 on all three small graphs, so it is an intentionally stringent reference rather than a competitive partitioning method at those scales.

## Large real-network sweep

Both graphs completed true CNM-full runs with no merge budget. The full rows below are extracted directly from `amazon_dblp.log`; each row reports the exact stopping index and the resulting K, Q, SNAP F1, SNAP NMI, stop fraction, and non-singleton node fraction.

| graph | method | modifier | W | parameter | sided | K | Q | F1 | NMI | stop_at | of_merges | stop_frac | frac_nonsingleton |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| amazon_full | CNM-full | - | 0 | 0.000 | 1 | 1475.000 | 0.871 | 0.622 | 0.884 | 333388.000 | 333388.000 | 1.000 | 1.000 |
| amazon_full | T1-2sided-k1 | none | 10 | 1.000 | 2 | 332415.000 | 0.003 | 0.241 | 0.822 | 2448.000 | 333388.000 | 0.007 | 0.015 |
| amazon_full | T2-short-1sided | none | 5 | 5.000 | 1 | 321864.000 | 0.015 | 0.272 | 0.824 | 12999.000 | 333388.000 | 0.039 | 0.073 |
| amazon_full | T3-diff-1sided | none | 15 | 3.000 | 1 | 332415.000 | 0.003 | 0.241 | 0.822 | 2448.000 | 333388.000 | 0.007 | 0.015 |
| amazon_full | T4-detrend-1sided | none | 10 | 1.000 | 1 | 332414.000 | 0.003 | 0.241 | 0.822 | 2449.000 | 333388.000 | 0.007 | 0.015 |
| amazon_full | T4-detrend-2sided | none | 10 | 1.000 | 2 | 332414.000 | 0.003 | 0.241 | 0.822 | 2449.000 | 333388.000 | 0.007 | 0.015 |
| amazon_full | T4-detrend-1sided | none | 10 | 2.000 | 1 | 332414.000 | 0.003 | 0.241 | 0.822 | 2449.000 | 333388.000 | 0.007 | 0.015 |
| amazon_full | T4-detrend-2sided | none | 10 | 2.000 | 2 | 332414.000 | 0.003 | 0.241 | 0.822 | 2449.000 | 333388.000 | 0.007 | 0.015 |
| amazon_full | T4-detrend-1sided | none | 10 | 3.000 | 1 | 321930.000 | 0.015 | 0.272 | 0.824 | 12933.000 | 333388.000 | 0.039 | 0.073 |
| amazon_full | T4-detrend-2sided | none | 10 | 3.000 | 2 | 329935.000 | 0.005 | 0.241 | 0.822 | 4928.000 | 333388.000 | 0.015 | 0.029 |
| amazon_full | T4-detrend-1sided | none | 15 | 1.000 | 1 | 332414.000 | 0.003 | 0.241 | 0.822 | 2449.000 | 333388.000 | 0.007 | 0.015 |
| amazon_full | T4-detrend-2sided | none | 15 | 1.000 | 2 | 332414.000 | 0.003 | 0.241 | 0.822 | 2449.000 | 333388.000 | 0.007 | 0.015 |
| amazon_full | T4-detrend-1sided | none | 15 | 2.000 | 1 | 332414.000 | 0.003 | 0.241 | 0.822 | 2449.000 | 333388.000 | 0.007 | 0.015 |
| amazon_full | T4-detrend-2sided | none | 15 | 2.000 | 2 | 332414.000 | 0.003 | 0.241 | 0.822 | 2449.000 | 333388.000 | 0.007 | 0.015 |
| amazon_full | T4-detrend-1sided | none | 15 | 3.000 | 1 | 332414.000 | 0.003 | 0.241 | 0.822 | 2449.000 | 333388.000 | 0.007 | 0.015 |
| amazon_full | T4-detrend-2sided | none | 15 | 3.000 | 2 | 332414.000 | 0.003 | 0.241 | 0.822 | 2449.000 | 333388.000 | 0.007 | 0.015 |
| amazon_full | T4-detrend-1sided | none | 20 | 1.000 | 1 | 332414.000 | 0.003 | 0.241 | 0.822 | 2449.000 | 333388.000 | 0.007 | 0.015 |
| amazon_full | T4-detrend-2sided | none | 20 | 1.000 | 2 | 332414.000 | 0.003 | 0.241 | 0.822 | 2449.000 | 333388.000 | 0.007 | 0.015 |
| amazon_full | T4-detrend-1sided | none | 20 | 2.000 | 1 | 332414.000 | 0.003 | 0.241 | 0.822 | 2449.000 | 333388.000 | 0.007 | 0.015 |
| amazon_full | T4-detrend-2sided | none | 20 | 2.000 | 2 | 332414.000 | 0.003 | 0.241 | 0.822 | 2449.000 | 333388.000 | 0.007 | 0.015 |
| amazon_full | T4-detrend-1sided | none | 20 | 3.000 | 1 | 332414.000 | 0.003 | 0.241 | 0.822 | 2449.000 | 333388.000 | 0.007 | 0.015 |
| amazon_full | T4-detrend-2sided | none | 20 | 3.000 | 2 | 332414.000 | 0.003 | 0.241 | 0.822 | 2449.000 | 333388.000 | 0.007 | 0.015 |
| amazon_full | T5-cusum | none | 15 | 50.000 | 1 | 320154.000 | 0.021 | 0.272 | 0.824 | 14709.000 | 333388.000 | 0.044 | 0.078 |
| amazon_full | T5-sr | none | 15 | 50.000 | 1 | 332413.000 | 0.003 | 0.241 | 0.822 | 2450.000 | 333388.000 | 0.007 | 0.015 |
| amazon_full | T1-2sided-k1 | defer | 10 | 1.000 | 2 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T2-short-1sided | defer | 5 | 5.000 | 1 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T3-diff-1sided | defer | 15 | 3.000 | 1 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-1sided | defer | 10 | 1.000 | 1 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-2sided | defer | 10 | 1.000 | 2 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-1sided | defer | 10 | 2.000 | 1 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-2sided | defer | 10 | 2.000 | 2 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-1sided | defer | 10 | 3.000 | 1 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-2sided | defer | 10 | 3.000 | 2 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-1sided | defer | 15 | 1.000 | 1 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-2sided | defer | 15 | 1.000 | 2 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-1sided | defer | 15 | 2.000 | 1 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-2sided | defer | 15 | 2.000 | 2 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-1sided | defer | 15 | 3.000 | 1 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-2sided | defer | 15 | 3.000 | 2 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-1sided | defer | 20 | 1.000 | 1 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-2sided | defer | 20 | 1.000 | 2 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-1sided | defer | 20 | 2.000 | 1 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-2sided | defer | 20 | 2.000 | 2 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-1sided | defer | 20 | 3.000 | 1 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-2sided | defer | 20 | 3.000 | 2 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T5-cusum | defer | 15 | 50.000 | 1 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T5-sr | defer | 15 | 50.000 | 1 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T1-2sided-k1 | mad | 10 | 1.000 | 2 | 321936.000 | 0.015 | 0.272 | 0.824 | 12927.000 | 333388.000 | 0.039 | 0.073 |
| amazon_full | T2-short-1sided | mad | 5 | 5.000 | 1 | 329346.000 | 0.006 | 0.244 | 0.822 | 5517.000 | 333388.000 | 0.017 | 0.033 |
| amazon_full | T3-diff-1sided | mad | 15 | 3.000 | 1 | 321912.000 | 0.015 | 0.272 | 0.824 | 12951.000 | 333388.000 | 0.039 | 0.073 |
| amazon_full | T4-detrend-1sided | mad | 10 | 1.000 | 1 | 329934.000 | 0.005 | 0.241 | 0.822 | 4929.000 | 333388.000 | 0.015 | 0.029 |
| amazon_full | T4-detrend-2sided | mad | 10 | 1.000 | 2 | 329934.000 | 0.005 | 0.241 | 0.822 | 4929.000 | 333388.000 | 0.015 | 0.029 |
| amazon_full | T4-detrend-1sided | mad | 10 | 2.000 | 1 | 329933.000 | 0.005 | 0.241 | 0.822 | 4930.000 | 333388.000 | 0.015 | 0.029 |
| amazon_full | T4-detrend-2sided | mad | 10 | 2.000 | 2 | 329933.000 | 0.005 | 0.241 | 0.822 | 4930.000 | 333388.000 | 0.015 | 0.029 |
| amazon_full | T4-detrend-1sided | mad | 10 | 3.000 | 1 | 329930.000 | 0.005 | 0.241 | 0.822 | 4933.000 | 333388.000 | 0.015 | 0.029 |
| amazon_full | T4-detrend-2sided | mad | 10 | 3.000 | 2 | 329930.000 | 0.005 | 0.241 | 0.822 | 4933.000 | 333388.000 | 0.015 | 0.029 |
| amazon_full | T4-detrend-1sided | mad | 15 | 1.000 | 1 | 329934.000 | 0.005 | 0.241 | 0.822 | 4929.000 | 333388.000 | 0.015 | 0.029 |
| amazon_full | T4-detrend-2sided | mad | 15 | 1.000 | 2 | 329934.000 | 0.005 | 0.241 | 0.822 | 4929.000 | 333388.000 | 0.015 | 0.029 |
| amazon_full | T4-detrend-1sided | mad | 15 | 2.000 | 1 | 329931.000 | 0.005 | 0.241 | 0.822 | 4932.000 | 333388.000 | 0.015 | 0.029 |
| amazon_full | T4-detrend-2sided | mad | 15 | 2.000 | 2 | 329931.000 | 0.005 | 0.241 | 0.822 | 4932.000 | 333388.000 | 0.015 | 0.029 |
| amazon_full | T4-detrend-1sided | mad | 15 | 3.000 | 1 | 329928.000 | 0.005 | 0.241 | 0.822 | 4935.000 | 333388.000 | 0.015 | 0.029 |
| amazon_full | T4-detrend-2sided | mad | 15 | 3.000 | 2 | 329928.000 | 0.005 | 0.241 | 0.822 | 4935.000 | 333388.000 | 0.015 | 0.029 |
| amazon_full | T4-detrend-1sided | mad | 20 | 1.000 | 1 | 329934.000 | 0.005 | 0.241 | 0.822 | 4929.000 | 333388.000 | 0.015 | 0.029 |
| amazon_full | T4-detrend-2sided | mad | 20 | 1.000 | 2 | 329934.000 | 0.005 | 0.241 | 0.822 | 4929.000 | 333388.000 | 0.015 | 0.029 |
| amazon_full | T4-detrend-1sided | mad | 20 | 2.000 | 1 | 329932.000 | 0.005 | 0.241 | 0.822 | 4931.000 | 333388.000 | 0.015 | 0.029 |
| amazon_full | T4-detrend-2sided | mad | 20 | 2.000 | 2 | 329932.000 | 0.005 | 0.241 | 0.822 | 4931.000 | 333388.000 | 0.015 | 0.029 |
| amazon_full | T4-detrend-1sided | mad | 20 | 3.000 | 1 | 329926.000 | 0.005 | 0.241 | 0.822 | 4937.000 | 333388.000 | 0.015 | 0.029 |
| amazon_full | T4-detrend-2sided | mad | 20 | 3.000 | 2 | 329926.000 | 0.005 | 0.241 | 0.822 | 4937.000 | 333388.000 | 0.015 | 0.029 |
| amazon_full | T5-cusum | mad | 15 | 50.000 | 1 | 329722.000 | 0.006 | 0.242 | 0.822 | 5141.000 | 333388.000 | 0.015 | 0.031 |
| amazon_full | T5-sr | mad | 15 | 50.000 | 1 | 329932.000 | 0.005 | 0.241 | 0.822 | 4931.000 | 333388.000 | 0.015 | 0.029 |
| amazon_full | T1-2sided-k1 | defer+mad | 10 | 1.000 | 2 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T2-short-1sided | defer+mad | 5 | 5.000 | 1 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T3-diff-1sided | defer+mad | 15 | 3.000 | 1 | 2988.000 | 0.870 | 0.695 | 0.918 | 331875.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-1sided | defer+mad | 10 | 1.000 | 1 | 2987.000 | 0.870 | 0.695 | 0.918 | 331876.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-2sided | defer+mad | 10 | 1.000 | 2 | 2987.000 | 0.870 | 0.695 | 0.918 | 331876.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-1sided | defer+mad | 10 | 2.000 | 1 | 2987.000 | 0.870 | 0.695 | 0.918 | 331876.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-2sided | defer+mad | 10 | 2.000 | 2 | 2987.000 | 0.870 | 0.695 | 0.918 | 331876.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-1sided | defer+mad | 10 | 3.000 | 1 | 2987.000 | 0.870 | 0.695 | 0.918 | 331876.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-2sided | defer+mad | 10 | 3.000 | 2 | 2987.000 | 0.870 | 0.695 | 0.918 | 331876.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-1sided | defer+mad | 15 | 1.000 | 1 | 2987.000 | 0.870 | 0.695 | 0.918 | 331876.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-2sided | defer+mad | 15 | 1.000 | 2 | 2987.000 | 0.870 | 0.695 | 0.918 | 331876.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-1sided | defer+mad | 15 | 2.000 | 1 | 2987.000 | 0.870 | 0.695 | 0.918 | 331876.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-2sided | defer+mad | 15 | 2.000 | 2 | 2987.000 | 0.870 | 0.695 | 0.918 | 331876.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-1sided | defer+mad | 15 | 3.000 | 1 | 2987.000 | 0.870 | 0.695 | 0.918 | 331876.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-2sided | defer+mad | 15 | 3.000 | 2 | 2987.000 | 0.870 | 0.695 | 0.918 | 331876.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-1sided | defer+mad | 20 | 1.000 | 1 | 2987.000 | 0.870 | 0.695 | 0.918 | 331876.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-2sided | defer+mad | 20 | 1.000 | 2 | 2987.000 | 0.870 | 0.695 | 0.918 | 331876.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-1sided | defer+mad | 20 | 2.000 | 1 | 2987.000 | 0.870 | 0.695 | 0.918 | 331876.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-2sided | defer+mad | 20 | 2.000 | 2 | 2987.000 | 0.870 | 0.695 | 0.918 | 331876.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-1sided | defer+mad | 20 | 3.000 | 1 | 2987.000 | 0.870 | 0.695 | 0.918 | 331876.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T4-detrend-2sided | defer+mad | 20 | 3.000 | 2 | 2987.000 | 0.870 | 0.695 | 0.918 | 331876.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | T5-cusum | defer+mad | 15 | 50.000 | 1 | 2809.000 | 0.870 | 0.684 | 0.914 | 332054.000 | 333388.000 | 0.996 | 1.000 |
| amazon_full | T5-sr | defer+mad | 15 | 50.000 | 1 | 2987.000 | 0.870 | 0.695 | 0.918 | 331876.000 | 333388.000 | 0.995 | 1.000 |
| amazon_full | BM-native-stop | - | 0 | 2.000 | 1 | 321522.000 | 0.016 | 0.272 | 0.824 | 13341.000 | 333388.000 | 0.040 | 0.074 |
| amazon_full | cm-baseline | - | 0 | 1.000 | 1 | 334863.000 | -0.000 | 0.240 | 0.822 | 0.000 | 333388.000 | 0.000 | 0.000 |
| amazon_full | cm-baseline | - | 0 | 5.000 | 1 | 334863.000 | -0.000 | 0.240 | 0.822 | 0.000 | 333388.000 | 0.000 | 0.000 |
| amazon_full | cm-baseline | - | 0 | 10.000 | 1 | 334863.000 | -0.000 | 0.240 | 0.822 | 0.000 | 333388.000 | 0.000 | 0.000 |
| amazon_full | cm-baseline | - | 0 | 50.000 | 1 | 334863.000 | -0.000 | 0.240 | 0.822 | 0.000 | 333388.000 | 0.000 | 0.000 |
| dblp_full | CNM-full | - | 0 | 0.000 | 1 | 3144.000 | 0.731 | 0.359 | 0.486 | 313936.000 | 313936.000 | 1.000 | 1.000 |
| dblp_full | T1-2sided-k1 | none | 10 | 1.000 | 2 | 315358.000 | 0.002 | 0.167 | 0.716 | 1722.000 | 313936.000 | 0.005 | 0.011 |
| dblp_full | T2-short-1sided | none | 5 | 5.000 | 1 | 313132.000 | 0.004 | 0.167 | 0.716 | 3948.000 | 313936.000 | 0.013 | 0.025 |
| dblp_full | T3-diff-1sided | none | 15 | 3.000 | 1 | 315359.000 | 0.002 | 0.167 | 0.716 | 1721.000 | 313936.000 | 0.005 | 0.011 |
| dblp_full | T4-detrend-1sided | none | 10 | 1.000 | 1 | 315358.000 | 0.002 | 0.167 | 0.716 | 1722.000 | 313936.000 | 0.005 | 0.011 |
| dblp_full | T4-detrend-2sided | none | 10 | 1.000 | 2 | 315358.000 | 0.002 | 0.167 | 0.716 | 1722.000 | 313936.000 | 0.005 | 0.011 |
| dblp_full | T4-detrend-1sided | none | 10 | 2.000 | 1 | 315358.000 | 0.002 | 0.167 | 0.716 | 1722.000 | 313936.000 | 0.005 | 0.011 |
| dblp_full | T4-detrend-2sided | none | 10 | 2.000 | 2 | 315358.000 | 0.002 | 0.167 | 0.716 | 1722.000 | 313936.000 | 0.005 | 0.011 |
| dblp_full | T4-detrend-1sided | none | 10 | 3.000 | 1 | 286845.000 | 0.039 | 0.183 | 0.718 | 30235.000 | 313936.000 | 0.096 | 0.150 |
| dblp_full | T4-detrend-2sided | none | 10 | 3.000 | 2 | 313131.000 | 0.004 | 0.167 | 0.716 | 3949.000 | 313936.000 | 0.013 | 0.025 |
| dblp_full | T4-detrend-1sided | none | 15 | 1.000 | 1 | 315358.000 | 0.002 | 0.167 | 0.716 | 1722.000 | 313936.000 | 0.005 | 0.011 |
| dblp_full | T4-detrend-2sided | none | 15 | 1.000 | 2 | 315358.000 | 0.002 | 0.167 | 0.716 | 1722.000 | 313936.000 | 0.005 | 0.011 |
| dblp_full | T4-detrend-1sided | none | 15 | 2.000 | 1 | 315358.000 | 0.002 | 0.167 | 0.716 | 1722.000 | 313936.000 | 0.005 | 0.011 |
| dblp_full | T4-detrend-2sided | none | 15 | 2.000 | 2 | 315358.000 | 0.002 | 0.167 | 0.716 | 1722.000 | 313936.000 | 0.005 | 0.011 |
| dblp_full | T4-detrend-1sided | none | 15 | 3.000 | 1 | 315358.000 | 0.002 | 0.167 | 0.716 | 1722.000 | 313936.000 | 0.005 | 0.011 |
| dblp_full | T4-detrend-2sided | none | 15 | 3.000 | 2 | 315358.000 | 0.002 | 0.167 | 0.716 | 1722.000 | 313936.000 | 0.005 | 0.011 |
| dblp_full | T4-detrend-1sided | none | 20 | 1.000 | 1 | 315358.000 | 0.002 | 0.167 | 0.716 | 1722.000 | 313936.000 | 0.005 | 0.011 |
| dblp_full | T4-detrend-2sided | none | 20 | 1.000 | 2 | 315358.000 | 0.002 | 0.167 | 0.716 | 1722.000 | 313936.000 | 0.005 | 0.011 |
| dblp_full | T4-detrend-1sided | none | 20 | 2.000 | 1 | 315358.000 | 0.002 | 0.167 | 0.716 | 1722.000 | 313936.000 | 0.005 | 0.011 |
| dblp_full | T4-detrend-2sided | none | 20 | 2.000 | 2 | 315358.000 | 0.002 | 0.167 | 0.716 | 1722.000 | 313936.000 | 0.005 | 0.011 |
| dblp_full | T4-detrend-1sided | none | 20 | 3.000 | 1 | 315358.000 | 0.002 | 0.167 | 0.716 | 1722.000 | 313936.000 | 0.005 | 0.011 |
| dblp_full | T4-detrend-2sided | none | 20 | 3.000 | 2 | 315358.000 | 0.002 | 0.167 | 0.716 | 1722.000 | 313936.000 | 0.005 | 0.011 |
| dblp_full | T5-cusum | none | 15 | 50.000 | 1 | 282969.000 | 0.068 | 0.183 | 0.716 | 34111.000 | 313936.000 | 0.109 | 0.161 |
| dblp_full | T5-sr | none | 15 | 50.000 | 1 | 315357.000 | 0.002 | 0.167 | 0.716 | 1723.000 | 313936.000 | 0.005 | 0.011 |
| dblp_full | T1-2sided-k1 | defer | 10 | 1.000 | 2 | 4924.000 | 0.730 | 0.341 | 0.500 | 312156.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T2-short-1sided | defer | 5 | 5.000 | 1 | 4922.000 | 0.730 | 0.341 | 0.500 | 312158.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T3-diff-1sided | defer | 15 | 3.000 | 1 | 4922.000 | 0.730 | 0.341 | 0.500 | 312158.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-1sided | defer | 10 | 1.000 | 1 | 4924.000 | 0.730 | 0.341 | 0.500 | 312156.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-2sided | defer | 10 | 1.000 | 2 | 4924.000 | 0.730 | 0.341 | 0.500 | 312156.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-1sided | defer | 10 | 2.000 | 1 | 4924.000 | 0.730 | 0.341 | 0.500 | 312156.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-2sided | defer | 10 | 2.000 | 2 | 4924.000 | 0.730 | 0.341 | 0.500 | 312156.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-1sided | defer | 10 | 3.000 | 1 | 4924.000 | 0.730 | 0.341 | 0.500 | 312156.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-2sided | defer | 10 | 3.000 | 2 | 4924.000 | 0.730 | 0.341 | 0.500 | 312156.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-1sided | defer | 15 | 1.000 | 1 | 4924.000 | 0.730 | 0.341 | 0.500 | 312156.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-2sided | defer | 15 | 1.000 | 2 | 4924.000 | 0.730 | 0.341 | 0.500 | 312156.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-1sided | defer | 15 | 2.000 | 1 | 4924.000 | 0.730 | 0.341 | 0.500 | 312156.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-2sided | defer | 15 | 2.000 | 2 | 4924.000 | 0.730 | 0.341 | 0.500 | 312156.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-1sided | defer | 15 | 3.000 | 1 | 4922.000 | 0.730 | 0.341 | 0.500 | 312158.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-2sided | defer | 15 | 3.000 | 2 | 4922.000 | 0.730 | 0.341 | 0.500 | 312158.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-1sided | defer | 20 | 1.000 | 1 | 4924.000 | 0.730 | 0.341 | 0.500 | 312156.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-2sided | defer | 20 | 1.000 | 2 | 4924.000 | 0.730 | 0.341 | 0.500 | 312156.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-1sided | defer | 20 | 2.000 | 1 | 4924.000 | 0.730 | 0.341 | 0.500 | 312156.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-2sided | defer | 20 | 2.000 | 2 | 4924.000 | 0.730 | 0.341 | 0.500 | 312156.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-1sided | defer | 20 | 3.000 | 1 | 4922.000 | 0.730 | 0.341 | 0.500 | 312158.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-2sided | defer | 20 | 3.000 | 2 | 4922.000 | 0.730 | 0.341 | 0.500 | 312158.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T5-cusum | defer | 15 | 50.000 | 1 | 4922.000 | 0.730 | 0.341 | 0.500 | 312158.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T5-sr | defer | 15 | 50.000 | 1 | 4923.000 | 0.730 | 0.341 | 0.500 | 312157.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T1-2sided-k1 | mad | 10 | 1.000 | 2 | 286908.000 | 0.039 | 0.182 | 0.718 | 30172.000 | 313936.000 | 0.096 | 0.150 |
| dblp_full | T2-short-1sided | mad | 5 | 5.000 | 1 | 313123.000 | 0.004 | 0.167 | 0.716 | 3957.000 | 313936.000 | 0.013 | 0.025 |
| dblp_full | T3-diff-1sided | mad | 15 | 3.000 | 1 | 313113.000 | 0.004 | 0.167 | 0.716 | 3967.000 | 313936.000 | 0.013 | 0.025 |
| dblp_full | T4-detrend-1sided | mad | 10 | 1.000 | 1 | 313130.000 | 0.004 | 0.167 | 0.716 | 3950.000 | 313936.000 | 0.013 | 0.025 |
| dblp_full | T4-detrend-2sided | mad | 10 | 1.000 | 2 | 313130.000 | 0.004 | 0.167 | 0.716 | 3950.000 | 313936.000 | 0.013 | 0.025 |
| dblp_full | T4-detrend-1sided | mad | 10 | 2.000 | 1 | 313129.000 | 0.004 | 0.167 | 0.716 | 3951.000 | 313936.000 | 0.013 | 0.025 |
| dblp_full | T4-detrend-2sided | mad | 10 | 2.000 | 2 | 313129.000 | 0.004 | 0.167 | 0.716 | 3951.000 | 313936.000 | 0.013 | 0.025 |
| dblp_full | T4-detrend-1sided | mad | 10 | 3.000 | 1 | 312934.000 | 0.004 | 0.168 | 0.716 | 4146.000 | 313936.000 | 0.013 | 0.026 |
| dblp_full | T4-detrend-2sided | mad | 10 | 3.000 | 2 | 313128.000 | 0.004 | 0.167 | 0.716 | 3952.000 | 313936.000 | 0.013 | 0.025 |
| dblp_full | T4-detrend-1sided | mad | 15 | 1.000 | 1 | 313130.000 | 0.004 | 0.167 | 0.716 | 3950.000 | 313936.000 | 0.013 | 0.025 |
| dblp_full | T4-detrend-2sided | mad | 15 | 1.000 | 2 | 313130.000 | 0.004 | 0.167 | 0.716 | 3950.000 | 313936.000 | 0.013 | 0.025 |
| dblp_full | T4-detrend-1sided | mad | 15 | 2.000 | 1 | 313125.000 | 0.004 | 0.167 | 0.716 | 3955.000 | 313936.000 | 0.013 | 0.025 |
| dblp_full | T4-detrend-2sided | mad | 15 | 2.000 | 2 | 313128.000 | 0.004 | 0.167 | 0.716 | 3952.000 | 313936.000 | 0.013 | 0.025 |
| dblp_full | T4-detrend-1sided | mad | 15 | 3.000 | 1 | 313113.000 | 0.004 | 0.167 | 0.716 | 3967.000 | 313936.000 | 0.013 | 0.025 |
| dblp_full | T4-detrend-2sided | mad | 15 | 3.000 | 2 | 313128.000 | 0.004 | 0.167 | 0.716 | 3952.000 | 313936.000 | 0.013 | 0.025 |
| dblp_full | T4-detrend-1sided | mad | 20 | 1.000 | 1 | 313130.000 | 0.004 | 0.167 | 0.716 | 3950.000 | 313936.000 | 0.013 | 0.025 |
| dblp_full | T4-detrend-2sided | mad | 20 | 1.000 | 2 | 313130.000 | 0.004 | 0.167 | 0.716 | 3950.000 | 313936.000 | 0.013 | 0.025 |
| dblp_full | T4-detrend-1sided | mad | 20 | 2.000 | 1 | 313125.000 | 0.004 | 0.167 | 0.716 | 3955.000 | 313936.000 | 0.013 | 0.025 |
| dblp_full | T4-detrend-2sided | mad | 20 | 2.000 | 2 | 313128.000 | 0.004 | 0.167 | 0.716 | 3952.000 | 313936.000 | 0.013 | 0.025 |
| dblp_full | T4-detrend-1sided | mad | 20 | 3.000 | 1 | 312828.000 | 0.004 | 0.168 | 0.716 | 4252.000 | 313936.000 | 0.014 | 0.026 |
| dblp_full | T4-detrend-2sided | mad | 20 | 3.000 | 2 | 313128.000 | 0.004 | 0.167 | 0.716 | 3952.000 | 313936.000 | 0.013 | 0.025 |
| dblp_full | T5-cusum | mad | 15 | 50.000 | 1 | 313113.000 | 0.004 | 0.167 | 0.716 | 3967.000 | 313936.000 | 0.013 | 0.025 |
| dblp_full | T5-sr | mad | 15 | 50.000 | 1 | 313113.000 | 0.004 | 0.167 | 0.716 | 3967.000 | 313936.000 | 0.013 | 0.025 |
| dblp_full | T1-2sided-k1 | defer+mad | 10 | 1.000 | 2 | 4924.000 | 0.730 | 0.341 | 0.500 | 312156.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T2-short-1sided | defer+mad | 5 | 5.000 | 1 | 4922.000 | 0.730 | 0.341 | 0.500 | 312158.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T3-diff-1sided | defer+mad | 15 | 3.000 | 1 | 4922.000 | 0.730 | 0.341 | 0.500 | 312158.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-1sided | defer+mad | 10 | 1.000 | 1 | 4921.000 | 0.730 | 0.342 | 0.500 | 312159.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-2sided | defer+mad | 10 | 1.000 | 2 | 4921.000 | 0.730 | 0.342 | 0.500 | 312159.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-1sided | defer+mad | 10 | 2.000 | 1 | 4921.000 | 0.730 | 0.342 | 0.500 | 312159.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-2sided | defer+mad | 10 | 2.000 | 2 | 4921.000 | 0.730 | 0.342 | 0.500 | 312159.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-1sided | defer+mad | 10 | 3.000 | 1 | 4921.000 | 0.730 | 0.342 | 0.500 | 312159.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-2sided | defer+mad | 10 | 3.000 | 2 | 4921.000 | 0.730 | 0.342 | 0.500 | 312159.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-1sided | defer+mad | 15 | 1.000 | 1 | 4921.000 | 0.730 | 0.342 | 0.500 | 312159.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-2sided | defer+mad | 15 | 1.000 | 2 | 4921.000 | 0.730 | 0.342 | 0.500 | 312159.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-1sided | defer+mad | 15 | 2.000 | 1 | 4921.000 | 0.730 | 0.342 | 0.500 | 312159.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-2sided | defer+mad | 15 | 2.000 | 2 | 4921.000 | 0.730 | 0.342 | 0.500 | 312159.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-1sided | defer+mad | 15 | 3.000 | 1 | 4921.000 | 0.730 | 0.342 | 0.500 | 312159.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-2sided | defer+mad | 15 | 3.000 | 2 | 4921.000 | 0.730 | 0.342 | 0.500 | 312159.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-1sided | defer+mad | 20 | 1.000 | 1 | 4921.000 | 0.730 | 0.342 | 0.500 | 312159.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-2sided | defer+mad | 20 | 1.000 | 2 | 4921.000 | 0.730 | 0.342 | 0.500 | 312159.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-1sided | defer+mad | 20 | 2.000 | 1 | 4921.000 | 0.730 | 0.342 | 0.500 | 312159.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-2sided | defer+mad | 20 | 2.000 | 2 | 4921.000 | 0.730 | 0.342 | 0.500 | 312159.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-1sided | defer+mad | 20 | 3.000 | 1 | 4921.000 | 0.730 | 0.342 | 0.500 | 312159.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T4-detrend-2sided | defer+mad | 20 | 3.000 | 2 | 4921.000 | 0.730 | 0.342 | 0.500 | 312159.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | T5-cusum | defer+mad | 15 | 50.000 | 1 | 4421.000 | 0.730 | 0.348 | 0.498 | 312659.000 | 313936.000 | 0.996 | 1.000 |
| dblp_full | T5-sr | defer+mad | 15 | 50.000 | 1 | 4921.000 | 0.730 | 0.342 | 0.500 | 312159.000 | 313936.000 | 0.994 | 1.000 |
| dblp_full | BM-native-stop | - | 0 | 2.000 | 1 | 286670.000 | 0.040 | 0.183 | 0.718 | 30410.000 | 313936.000 | 0.097 | 0.150 |
| dblp_full | cm-baseline | - | 0 | 1.000 | 1 | 317080.000 | -0.000 | 0.167 | 0.716 | 0.000 | 313936.000 | 0.000 | 0.000 |
| dblp_full | cm-baseline | - | 0 | 5.000 | 1 | 317080.000 | -0.000 | 0.167 | 0.716 | 0.000 | 313936.000 | 0.000 | 0.000 |
| dblp_full | cm-baseline | - | 0 | 10.000 | 1 | 317080.000 | -0.000 | 0.167 | 0.716 | 0.000 | 313936.000 | 0.000 | 0.000 |
| dblp_full | cm-baseline | - | 0 | 50.000 | 1 | 317080.000 | -0.000 | 0.167 | 0.716 | 0.000 | 313936.000 | 0.000 | 0.000 |

# amazon_full CNM-full merges=333388 budget=0 singleton_phase_end=331875
# dblp_full CNM-full merges=313936 budget=0 singleton_phase_end=312156

### Plateau/tie-block diagnostics

| graph | merge_index | dQ |
| --- | ---: | ---: |
| amazon_full | 2438 | 1.08006173632758866e-06 |
| amazon_full | 2439 | 1.08006173632758866e-06 |
| amazon_full | 2440 | 1.08006173632758866e-06 |
| amazon_full | 2441 | 1.08006173632758866e-06 |
| amazon_full | 2442 | 1.08006173632758866e-06 |
| amazon_full | 2443 | 1.08006173632758866e-06 |
| amazon_full | 2444 | 1.08006173632758866e-06 |
| amazon_full | 2445 | 1.08006173632758866e-06 |
| amazon_full | 2446 | 1.08006173632758866e-06 |
| amazon_full | 2447 | 1.08006173632758866e-06 |
| amazon_full | 2448 | 1.08006115305965150e-06 |
| amazon_full | 2449 | 1.08006115305965150e-06 |
| amazon_full | 2450 | 1.08006115305965150e-06 |
| amazon_full | 2451 | 1.08006115305965150e-06 |
| amazon_full | 2452 | 1.08006115305965150e-06 |
| amazon_full | 2453 | 1.08006115305965150e-06 |
| amazon_full | 2454 | 1.08006115305965150e-06 |
| amazon_full | 2455 | 1.08006115305965150e-06 |
| amazon_full | 2456 | 1.08006115305965150e-06 |
| amazon_full | 2457 | 1.08006115305965150e-06 |
| amazon_full | 2458 | 1.08006115305965150e-06 |
| dblp_full | 2438 | 9.52501148952566497e-07 |
| dblp_full | 2439 | 9.52501148952566497e-07 |
| dblp_full | 2440 | 9.52501148952566497e-07 |
| dblp_full | 2441 | 9.52501148952566497e-07 |
| dblp_full | 2442 | 9.52501148952566497e-07 |
| dblp_full | 2443 | 9.52501148952566497e-07 |
| dblp_full | 2444 | 9.52501148952566497e-07 |
| dblp_full | 2445 | 9.52501148952566497e-07 |
| dblp_full | 2446 | 9.52501148952566497e-07 |
| dblp_full | 2447 | 9.52501148952566497e-07 |
| dblp_full | 2448 | 9.52501148952566497e-07 |
| dblp_full | 2449 | 9.52501148952566497e-07 |
| dblp_full | 2450 | 9.52501148952566497e-07 |
| dblp_full | 2451 | 9.52501148952566497e-07 |
| dblp_full | 2452 | 9.52501148952566497e-07 |
| dblp_full | 2453 | 9.52501148952566497e-07 |
| dblp_full | 2454 | 9.52501148952566497e-07 |
| dblp_full | 2455 | 9.52501148952566497e-07 |
| dblp_full | 2456 | 9.52501148952566497e-07 |
| dblp_full | 2457 | 9.52501148952566497e-07 |
| dblp_full | 2458 | 9.52501148952566497e-07 |

# plateau-window amazon_full indices=2438..2458 bit_identical_consecutive=no
# plateau-window dblp_full indices=2438..2458 bit_identical_consecutive=yes

The verdict is explicit: Amazon has `bit_identical_consecutive=no` (the value changes at the 2447/2448 boundary, although both sides contain long exact runs); DBLP has `bit_identical_consecutive=yes` across every consecutive value in 2438--2458. The complete raw diagnostic log is [amazon_dblp.log](/Users/legol/Coding/Bader/graz-bader-mccloskey-research/amazon_dblp.log).

The implementation uses at most 51 candidate cuts: every `ceil(T/50)` merges plus the terminal cut. It uses a fixed 90/10 edge split for held-out block likelihood, assignment entropy plus block Bernoulli code length for MDL, and 20 independent 5%-edge-deletion CNM traces for stability.

| graph | candidates | likelihood cut | MDL cut | stability cut | CNM argmax-Q cut | BM-native cut | stability reps |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| karate | 30 | 0 | 0 | 0 | 29 | 0 | 20 |

The karate result is a warning, not a success: naive positive-edge holdout, sparse-block MDL, and NMI stability all favour the singleton partition on this sparse graph. They should not be used as an external selector in the paper without a degree-corrected likelihood and a non-singleton support constraint. Amazon/DBLP multicut runs are intentionally not claimed complete.

## Self-check

- [x] Fixes 1--5 are implemented; fixes 1--4 were run on ER and all small graphs, and fix 5 was run on karate.
- [x] ER uses 10 seeds and every new report row computes `frac_nonsingleton`.
- [x] Mandatory Amazon/DBLP large-graph rows, F1/NMI, and plateau findings are complete; large-graph multicut results remain pending.
- [x] Small-network experiment now emits rows.
- [x] CNM merge selection was not changed; all new methods replay its recorded trace.
