# Bader--McCloskey stopping-rule revision, round 4

- Completed: graph-shape and repeat-trace validation; deterministic CNM tie-breaking; trace/rule runtime instrumentation; singleton-deferral warm-up control; LFR, Cora, Citeseer, PubMed, and ego-Facebook preparation and benchmark runs.
- Resolution of Finding C: the four checked prepared edge files are byte-identical between the initial checked-in revision and this revision. The prior CNM heap compared equal dQ values without a secondary key while neighbour insertion came from unordered maps. CNM now sorts initial adjacency and breaks exact dQ ties by the lexicographically smaller live community pair. Repeat traces are bit-identical after this fix.
- Resolution of Finding A: artifact confirmed. Legacy deferred CUSUM fires at the singleton boundary; clearing the history and requiring 15 fresh merges fires later (Amazon +160 merges, DBLP +14). Amazon still exceeds CNM-full on F1/NMI; DBLP improves NMI but not F1.
- Default benchmark rows contain only CNM-full and the corrected current candidate, `T5-cusum-warmup`. The legacy-defer control is reported only in the warm-up diagnostic. Historical full sweeps remain runnable via `--full-stop-report`.
- CNM-full means the unbudgeted positive-dQ greedy terminus. No merge budget was used in this report. com-Orkut was not attempted because it is substantially larger than the completed suite. Linux static linking could not be tested locally: macOS clang lacks `crt0.o`; non-static C++17 compilation passed.

## Graph-shape fingerprints and CNM determinism

`edge_fnv1a64` hashes sorted loaded undirected edge pairs. `raw_lines` and `unique_edges` are equal for all four audited source files: no self-loops or duplicate lines were removed.

| graph | raw_lines | N | unique_edges | degree_min | degree_max | degree_mean | degree_median | degree1 | components | edge_fnv1a64 | run1_merges | run2_merges | trace_bit_identical |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | --- |
| amazon_full | 925872 | 334863 | 925872 | 1 | 549 | 5.530 | 4 | 25709 | 1 | c6ef7c480926a0ca | 333516 | 333516 | yes |
| dblp_full | 1049866 | 317080 | 1049866 | 1 | 343 | 6.622 | 4 | 43181 | 1 | 486fa1d9ccdfc276 | 313729 | 313729 | yes |
| polblogs | 10260 | 1490 | 10260 | 3 | 27 | 13.772 | 14 | 0 | 1 | cd198b891fe4c37c | 1485 | 1485 | yes |
| email | 16064 | 986 | 16064 | 1 | 345 | 32.584 | 22 | 95 | 1 | afc9879412eeb06f | 978 | 978 | yes |

## Singleton warm-up control

`T5-cusum-warmup` clears all pre-deferral dQ history and waits W=15 fresh post-singleton merges. `WARMUP_CONTROL` is legacy deferred CUSUM and is diagnostic-only.

| graph | method | stop_at | of_merges | stop_frac | K | Q | F1 | NMI | frac_nonsingleton | trace_build_ms | rule_eval_ms | total_ms | window_evals |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| amazon_full | CNM-full | 333516 | 333516 | 1.000 | 1347 | 0.871 | 0.616 | 0.886 | 1.000 | 63686.019 | 0.000 | 63686.019 | 0 |
| amazon_full | T5-cusum-warmup | 332306 | 333516 | 0.996 | 2557 | 0.870 | 0.678 | 0.914 | 1.000 | 63686.019 | 0.201 | 63686.220 | 146 |
| amazon_full | WARMUP_CONTROL:T5-cusum-defer-legacy | 332146 | 333516 | 0.996 | 2717 | 0.870 | 0.692 | 0.920 | 1.000 | 63686.019 | 0.172 | 63686.192 | 1 |
| dblp_full | CNM-full | 313729 | 313729 | 1.000 | 3351 | 0.733 | 0.347 | 0.512 | 1.000 | 149653.292 | 0.000 | 149653.292 | 0 |
| dblp_full | T5-cusum-warmup | 312095 | 313729 | 0.995 | 4985 | 0.732 | 0.344 | 0.526 | 1.000 | 149653.292 | 0.167 | 149653.459 | 22 |
| dblp_full | WARMUP_CONTROL:T5-cusum-defer-legacy | 312081 | 313729 | 0.995 | 4999 | 0.732 | 0.344 | 0.526 | 1.000 | 149653.292 | 0.169 | 149653.461 | 23 |

## Pruned benchmark: existing labelled networks

| graph | method | stop_at | of_merges | stop_frac | K | Q | F1 | NMI | frac_nonsingleton | trace_build_ms | rule_eval_ms | total_ms | window_evals |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| karate | CNM-full | 31 | 31 | 1.000 | 3 | 0.381 | 0.761 | 0.565 | 1.000 | 0.054 | 0.000 | 0.054 | 0 |
| karate | T5-cusum-warmup | 31 | 31 | 1.000 | 3 | 0.381 | 0.761 | 0.565 | 1.000 | 0.054 | 0.000 | 0.054 | 0 |
| polblogs | CNM-full | 1485 | 1485 | 1.000 | 5 | 0.354 | 0.691 | 0.878 | 1.000 | 52.313 | 0.000 | 52.313 | 0 |
| polblogs | T5-cusum-warmup | 1485 | 1485 | 1.000 | 5 | 0.354 | 0.691 | 0.878 | 1.000 | 52.313 | 0.001 | 52.314 | 0 |
| email | CNM-full | 978 | 978 | 1.000 | 8 | 0.347 | 0.283 | 0.427 | 1.000 | 25.115 | 0.000 | 25.115 | 0 |
| email | T5-cusum-warmup | 978 | 978 | 1.000 | 8 | 0.347 | 0.283 | 0.427 | 1.000 | 25.115 | 0.001 | 25.116 | 0 |
| ego_facebook | CNM-full | 4026 | 4026 | 1.000 | 13 | 0.777 | 0.319 | 0.531 | 1.000 | 189.168 | 0.000 | 189.168 | 0 |
| ego_facebook | T5-cusum-warmup | 4026 | 4026 | 1.000 | 13 | 0.777 | 0.319 | 0.531 | 1.000 | 189.168 | 0.002 | 189.171 | 0 |

## Pruned benchmark: LFR

| graph | method | stop_at | of_merges | stop_frac | K | Q | F1 | NMI | frac_nonsingleton | trace_build_ms | rule_eval_ms | total_ms | window_evals |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| lfr_n1000_mu1 | CNM-full | 987 | 987 | 1.000 | 13 | 0.704 | 0.832 | 0.916 | 1.000 | 14.482 | 0.000 | 14.482 | 0 |
| lfr_n1000_mu1 | T5-cusum-warmup | 987 | 987 | 1.000 | 13 | 0.704 | 0.832 | 0.916 | 1.000 | 14.482 | 0.001 | 14.483 | 0 |
| lfr_n1000_mu3 | CNM-full | 994 | 994 | 1.000 | 6 | 0.373 | 0.544 | 0.602 | 1.000 | 22.272 | 0.000 | 22.272 | 0 |
| lfr_n1000_mu3 | T5-cusum-warmup | 994 | 994 | 1.000 | 6 | 0.373 | 0.544 | 0.602 | 1.000 | 22.272 | 0.001 | 22.273 | 0 |
| lfr_n1000_mu5 | CNM-full | 994 | 994 | 1.000 | 6 | 0.188 | 0.224 | 0.166 | 1.000 | 31.725 | 0.000 | 31.725 | 0 |
| lfr_n1000_mu5 | T5-cusum-warmup | 994 | 994 | 1.000 | 6 | 0.188 | 0.224 | 0.166 | 1.000 | 31.725 | 0.001 | 31.726 | 0 |
| lfr_n1000_mu7 | CNM-full | 995 | 995 | 1.000 | 5 | 0.164 | 0.166 | 0.035 | 1.000 | 33.073 | 0.000 | 33.073 | 0 |
| lfr_n1000_mu7 | T5-cusum-warmup | 995 | 995 | 1.000 | 5 | 0.164 | 0.166 | 0.035 | 1.000 | 33.073 | 0.001 | 33.073 | 0 |
| lfr_n5000_mu1 | CNM-full | 4974 | 4974 | 1.000 | 26 | 0.766 | 0.751 | 0.865 | 1.000 | 236.465 | 0.000 | 236.465 | 0 |
| lfr_n5000_mu1 | T5-cusum-warmup | 4974 | 4974 | 1.000 | 26 | 0.766 | 0.751 | 0.865 | 1.000 | 236.465 | 0.003 | 236.468 | 0 |
| lfr_n5000_mu3 | CNM-full | 4992 | 4992 | 1.000 | 8 | 0.392 | 0.396 | 0.531 | 1.000 | 627.250 | 0.000 | 627.250 | 0 |
| lfr_n5000_mu3 | T5-cusum-warmup | 4992 | 4992 | 1.000 | 8 | 0.392 | 0.396 | 0.531 | 1.000 | 627.250 | 0.003 | 627.253 | 0 |
| lfr_n5000_mu5 | CNM-full | 4994 | 4994 | 1.000 | 6 | 0.202 | 0.138 | 0.155 | 1.000 | 1007.577 | 0.000 | 1007.577 | 0 |
| lfr_n5000_mu5 | T5-cusum-warmup | 4994 | 4994 | 1.000 | 6 | 0.202 | 0.138 | 0.155 | 1.000 | 1007.577 | 0.003 | 1007.580 | 0 |
| lfr_n5000_mu7 | CNM-full | 4993 | 4993 | 1.000 | 7 | 0.168 | 0.064 | 0.026 | 1.000 | 1036.644 | 0.000 | 1036.644 | 0 |
| lfr_n5000_mu7 | T5-cusum-warmup | 4993 | 4993 | 1.000 | 7 | 0.168 | 0.064 | 0.026 | 1.000 | 1036.644 | 0.003 | 1036.647 | 0 |

## Pruned benchmark: citation networks

Citeseer applies Planetoid's conventional all-zero row padding to 15 isolated test-index nodes. ego-Facebook projects overlapping SNAP circles to the first deterministic circle label per node; it is included for compatibility with single-label F1/NMI, not as an overlapping-community evaluation.

| graph | method | stop_at | of_merges | stop_frac | K | Q | F1 | NMI | frac_nonsingleton | trace_build_ms | rule_eval_ms | total_ms | window_evals |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| cora | CNM-full | 2603 | 2603 | 1.000 | 105 | 0.804 | 0.209 | 0.208 | 1.000 | 19.296 | 0.000 | 19.296 | 0 |
| cora | T5-cusum-warmup | 2602 | 2603 | 1.000 | 106 | 0.804 | 0.209 | 0.208 | 1.000 | 19.296 | 0.010 | 19.305 | 37 |
| citeseer | CNM-full | 2840 | 2840 | 1.000 | 439 | 0.870 | 0.130 | 0.212 | 1.000 | 10.424 | 0.000 | 10.424 | 0 |
| citeseer | T5-cusum-warmup | 2838 | 2840 | 0.999 | 441 | 0.870 | 0.130 | 0.212 | 1.000 | 10.424 | 0.010 | 10.434 | 37 |
| pubmed | CNM-full | 19597 | 19597 | 1.000 | 120 | 0.727 | 0.277 | 0.197 | 1.000 | 4373.598 | 0.000 | 4373.598 | 0 |
| pubmed | T5-cusum-warmup | 19597 | 19597 | 1.000 | 120 | 0.727 | 0.277 | 0.197 | 1.000 | 4373.598 | 0.023 | 4373.621 | 55 |

## Raw artifacts

| artifact | contents |
| --- | --- |
| `round4_fingerprint.log` | full Amazon/DBLP repeat-trace output and fingerprints |
| `round4_amazon.log`, `round4_dblp.log` | warm-up control runs for mandatory large graphs |
| `round4_small.log`, `round4_bench.log`, `round4_ego.log` | complete tab-separated benchmark rows |
| `scripts/prepare_round4_benchmarks.py` | deterministic LFR/citation/ego-Facebook preparation |
