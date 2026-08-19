# Rolling-window dQ stopping rule — implementation status and findings

Working notes for the CNM-style stopping criterion (keep the last W dQ scores,
stop when the current dQ is more than k standard deviations from their mean).
Written so the work survives a context reset.

## What is implemented (cpp/graz.cpp)

- `struct RollStop` — the rule itself. Window W, multiplier k, one/two-sided.
  The tripping value is NOT added to the window (it is the anomaly under test).
- `cnm_greedy_trace(g, max_pushes)` — CNM greedy agglomeration recording every
  accepted dQ and merge pair. Vector-backed max-heap with lazy deletion by
  version stamps, plus **stale-entry compaction** when the heap exceeds the
  ceiling. Compaction drops exactly the entries the pop loop would have skipped,
  so it is result-preserving; it is what makes the big graphs fit in RAM.
- `cnm_labels_at(trace, cut)` — replays the first `cut` merges via union-find.
- `rollstop_cut(trace, W, k, two_sided)` — where the rule cuts a recorded trace.
- `rollq(...)` / `Rule::ROLLQ` — the same rule ported onto the Leiden move loop,
  to check results are not an artifact of the agglomerative schedule.
- `dump_trace_diagnostics(...)` — prints the dQ trace head, quartiles, first-20
  window mean/sd/cv, and the (mu-dq)/sd ratio at the firing point.
- Experiments 9 (large graphs), 10 (small labelled), 11 (ER nulls).
- Flags: `--only-roll`, `--roll-small`, `--roll-er`, `--max-pushes=N`.

Design point that made the sweep affordable: the rule is a pure function of the
dQ sequence, so ONE agglomeration answers every (W, k, sided) combination. The
tuning sweep costs one CNM run, not one run per setting.

## Verification

CNM-full on karate reproduces the published Clauset-Newman-Moore result exactly:
**K=3, Q=0.381** (paper: Q=0.3807, 3 communities). Merge counts karate 31,
polblogs 1482, email 977 are stable across every refactor.

## Finding 1 — the rule PASSES the Corollary 1 test that Eq. 6 fails

On Erdos-Renyi graphs (no communities; true K = 0), K counting size >= 3:

| N | Louvain | GRAZ t=2 | CNM-full | rolling rule |
|---|---|---|---|---|
| 1000 | 14 | 58 | 8 | **0.0 – 1.7** |
| 2000 | 14 | 58 | 9.7 | **0.0** |
| 5000 | 14 | 79 | 12 | **0.0** |

The rolling rule reports essentially zero spurious communities where GRAZ's Eq. 6
reports 58–79 and plain Louvain reports 14. It also moves in the CORRECT
direction: raising k makes it more permissive (0.33 -> 6.67 at N=1000), whereas
raising tau in Eq. 6 made fragmentation *worse*. This is the test that mattered.

## Finding 2 — but it fires far too early on real graphs, and why

dQ scales as ~1/m, so the spread available to the test collapses as graphs grow:

| graph | m | typical dQ | fires at | of merges |
|---|---|---|---|---|
| karate | 78 | ~1e-2 | 30 | 31 |
| email | 25,571 | ~1e-5 | 42 | 977 |
| polblogs | 16,714 | ~1e-5 | 17 | 1482 |

On karate (dQ well separated) the rule stops at 30/31 — the natural terminus, and
NMI 0.537–0.565 beats both Louvain (0.498) and GRAZ t=2 (0.369, which
over-segments to K=16). On the 1000+ node graphs the early merges have
near-identical dQ, so the window sd collapses toward zero, any small dip reads as
"many sd away", and the run is cut at ~1–3% of its merges (Q ~ 0.000, K ~ 1480).

This is a property of the dQ SEQUENCE, not of the implementation. The mechanism is
measured, not assumed: `dump_trace_diagnostics` prints the first-20 window
mean/sd/cv and the (mu-dq)/sd ratio at the firing point.

Implication: the rule needs a scale-free formulation (test in log-dQ space,
require a minimum relative drop, or use a robust/relative dispersion floor)
before it can serve as a stopping criterion on large graphs. The two-sided
variant is additionally wrong in principle — an upward outlier is an unusually
GOOD merge and halting on it is backwards — and it does cut earlier than the
one-sided variant on karate (K=13 vs K=4 at W=10).

## Finding 3 — the scale-free fix works, but there is a real trade-off

Two scale-free variants were added (both are replays of the same trace, so nearly
free):

- `roll-log` — the same mean/sd test on log(dQ).
- `roll-relXX` — require the drop to exceed k sd's AND be at least XX% of the
  window mean, so a merely flat window cannot fire on numerical noise.

`roll-log` does NOT fix it (polblogs still cut at 20/1482). Log space rescales the
values but the *relative* spread is what collapses, so the ratio still blows up.
That hypothesis was wrong and the measurement said so.

`roll-rel50` fixes the SMALL real graphs only. On polblogs (1,490 nodes) it cuts at
1479/1482, giving K=11, Q=0.349, NMI=0.818 versus CNM-full's K=8, Q=0.349,
NMI=0.826 — essentially full quality.

On Amazon (334,863 nodes) it does NOT hold up: it cuts at 13,367 of 333,311 merges
(4.0%), leaving K=321,496 and Q=0.017 against CNM-full's K=1,552 and Q=0.868.
Identical for every (W, k) in the sweep, so it is not a tuning failure. The
relative floor buys roughly a 5x later cut than the raw rule (0.7% -> 4.0%) and
that is all; it does not survive two more orders of magnitude of graph size.

But it loses the property that made the rule interesting. On ER nulls
(true K = 0, K counting size >= 3):

| N | Louvain | GRAZ t=2 | CNM-full | roll-1sided | roll-rel25 | roll-rel50 |
|---|---|---|---|---|---|---|
| 1000 | 14 | 58 | 8.7 | **1.0** | 29.7 | 8.7 |
| 2000 | 14.3 | 58.3 | 9.7 | **0.0** | 89.3 | 11.7 |
| 5000 | 14 | 79.3 | 10.3 | **0.0** | 444.7 | 12.0 |
| 10000 | 24.3 | 93.3 | 13.3 | **0.3** | 1572.3 | 15.0 |

So:
- the literal rule (`roll-1sided`) is an excellent noise filter and a bad
  community detector;
- `roll-rel50` is a good community detector and no longer filters noise (it
  matches CNM-full, i.e. the stop stops mattering);
- `roll-rel25` is worse than either — the relative floor and the sd test fight
  each other, and it fires in the middle of the run.

No single (W, k, min_rel) setting does both jobs on both graph classes. That is the
honest headline: the criterion is sound and cheap, and it decisively beats Eq. 6 at
rejecting structureless graphs, but as a general-purpose stopping rule it needs a
dispersion estimate that is stable across the whole dQ range, not merely rescaled.

## Reproducibility note

Tie-breaking among exactly-equal dQ merges depends on heap insertion order, so
`std::make_heap` vs incremental push shifts polblogs NMI between 0.826 and 0.834
(K=8, Q=0.349, and 1482 merges identical in both). Both are valid CNM runs.
Storing dq as float was tried and rejected for the same reason: it reordered
near-ties (email 977 -> 978 merges).

## State of runs

- Small graphs (karate, polblogs, email) + ER nulls: DONE.
  `results/roll_small_er.txt`, `results/roll_small_v2.txt`
- Large graphs: ALL FOUR DONE. amazon_full, dblp_full, youtube ran to the greedy
  terminus; livejournal ran under a 200,000-merge budget (see below).
  `roll_amazon.txt`, `roll_dblp.txt`, `roll_youtube.txt`, `roll_lj_budget.txt`

### LiveJournal: why a merge budget, and why the results are still exact

A full CNM on livejournal (N=3,997,962, M=34,681,189) is not feasible. Measured:
the first 105,000 merges take **8.9 seconds**, but by merge 155,000 one hub
community reaches degree 494,254, and CNM re-pushes the survivor's whole
neighbour list every merge. Merges 155k-355k took 7 hours. The run reached
403,429 merges in 9.1 h at a still-decaying 5.5 merges/s; the ~3.93M merges a
full run needs project to **~7.5 days**, and that is an underestimate.

Fix: `--max-merges=N`, a BUDGET rather than a guard. CNM is deterministic, so the
recorded prefix is identical to the prefix of an uncapped run; the budget only
decides when to stop appending. With a 200,000 budget the run finished in 85 min,
and the stopping rule fires at **26,396 of 200,000**, well inside the prefix, so
every roll row is exact. Only the CNM-full reference is lost, and it is reported
as `CNM-prefix` with an explicit note rather than passed off as argmax Q.

Regression after adding the budget: karate 31 merges / K=3 / Q=0.381, polblogs
1482 merges / K=8 / Q=0.349 / NMI=0.826 — unchanged.

LiveJournal confirms the pattern at the largest scale tested: all 24 literal-rule
settings cut within one merge of each other (26,396-26,397), and the Leiden port
fires after 18 moves.

### Why the large-graph runs kept dying (two separate causes)

**Cause 1 — heap growth from hub re-pushing.** Measured, and it corrected a wrong
first guess. I assumed the neighbour maps were the memory hog; instrumenting them
disproved it — `nbr_entries` *shrinks* over the run (1.84M -> 1.62M at merge 45k).
The actual driver is `max deg`, the largest community's neighbour count, which
explodes 549 -> 17,773 by merge 45k. CNM re-pushes the survivor's entire
neighbour list on every merge, so one hub refills a 20M-entry heap every ~1,400
merges. Stale-entry compaction fixes it: ~95% of entries are dead
(20,001,093 -> 890,938 live), and compaction is result-preserving because it drops
exactly the entries the pop loop would have skipped. With it, memory stays flat
and the run progresses steadily.

**Cause 2 — the runs were being killed, not crashing.** `nohup ... &` inside the
tool call meant the process died when its parent shell exited, always mid-run.
The same binary reached 45,000 merges when run in the foreground. Use the
harness's own background mechanism (`run_in_background`), not `nohup &`.

Neither cause was a defect in the stopping rule or in CNM correctness.

**Cause 3 (LiveJournal only) — the seed heap alone exceeds the ceiling.** With
34.7M edges, every edge seeds a positive-dQ candidate, so the heap holds
34,681,188 entries after ONE merge. The guard fired correctly and refused to
report a partial run, but the 20M ceiling was simply too small for this graph.
These entries are all genuinely live, so compaction cannot help; the fix is a
larger ceiling (`--max-pushes=120000000`, about 2.9 GB at 24 B/entry).
- Notion page to receive results: "Improved Statistical Functionality"
  https://app.notion.com/p/Improved-Statistical-Functionality-3af5f411773980c7a09dd11fefea842b
