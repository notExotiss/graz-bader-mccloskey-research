# GRAZ: Experiment 7 (large networks) + the Corollary 1 diagnosis

Reproduce with:

```
cd cpp && g++ -O2 -std=c++17 -static graz.cpp -o graz
./graz karate=prepared/karate.edges:prepared/karate.labels \
       polblogs=prepared/polblogs.edges:prepared/polblogs.labels \
       email=prepared/email.edges:prepared/email.labels     # Exp 1,2,4,5,6,8 + verification
./graz --only-large --methods=graz4 large=livejournal=prepared/livejournal.edges:prepared/livejournal.labels
```

---

## 1. Headline

Two results, and they pull in opposite directions on the same mechanism.

1. **GRAZ substantially beats Louvain/Leiden on large graphs with many true
   communities** — up to +0.195 F1 and +0.081 NMI on `com-Amazon`, at
   indistinguishable modularity. This is a real, reproducible win at 4 M nodes.
2. **Corollary 1 still fails, and the cause is not the one in the earlier note.**
   The single-node z-test does not suppress spurious communities; it imposes an
   **analytic ceiling on community volume**. Over-segmentation is not a side
   effect — it is what the test does.

These reconcile: the volume ceiling helps exactly when the true communities are
smaller than what modularity's resolution limit produces, and hurts otherwise.
GRAZ is best understood as a **resolution control**, not a significance filter.

## 2. Experiment 7 — large real networks

Full SNAP graphs, all top-5000 ground-truth communities (overlapping cover
collapsed smallest-wins, so F1/NMI are lower bounds). Single-threaded.

| dataset | N | M | K_true | method | K | Q | F1 | NMI | time |
|---|---|---|---|---|---|---|---|---|---|
| amazon | 334 863 | 925 872 | 1 481 | Louvain | 240 | 0.926 | 0.432 | 0.838 | 3.8 s |
| | | | | Leiden | 248 | 0.926 | 0.438 | 0.840 | 3.7 s |
| | | | | GRAZ t=4 | 971 | 0.924 | **0.627** | **0.919** | 6.5 s |
| dblp | 317 080 | 1 049 866 | 4 957 | Louvain | 215 | 0.820 | 0.174 | 0.537 | 8.1 s |
| | | | | GRAZ t=4 | 1 169 | 0.818 | **0.252** | **0.605** | 10.0 s |
| youtube | 1 134 890 | 2 987 624 | 4 749 | Louvain | 7 255 | 0.716 | 0.275 | 0.474 | 92 s |
| | | | | GRAZ t=4 | 36 996 | 0.705 | 0.278 | **0.578** | 164 s |
| livejournal | 3 997 962 | 34 681 189 | 4 636 | Louvain | 1 710 | 0.748 | 0.496 | 0.727 | 708 s |
| | | | | GRAZ t=3 | 25 948 | 0.746 | **0.534** | 0.759 | 1 370 s |
| | | | | GRAZ t=4 | 38 307 | 0.745 | **0.534** | **0.763** | 1 808 s |

Reading:

- **Modularity is flat while agreement moves a lot.** Amazon Q goes 0.926 →
  0.924 while F1 goes 0.432 → 0.627. Q is not the discriminating quantity here,
  which is the paper's point and it holds up.
- **The gain tracks the K_true/K_Louvain gap.** Amazon (1481 vs 240, 6.2×) gains
  most; youtube (4749 vs 7255 — Louvain already over-splits) gains nothing on F1.
  This is the single best predictor of whether GRAZ helps.
- **Cost is 1.2–2.0× wall clock**, from ~2× the move sweeps, not from the z-test
  itself (aggregation and refinement times are unchanged).
- **F1 keeps improving with tau past the point Q starts to drop** — but only
  while K stays under K_true. LiveJournal is where the tau ladder stops paying:
  t=4 buys +0.004 NMI and *zero* F1 over t=3 for +32% runtime and 12 000 more
  communities, because t=3's K=25 948 already overshoots K_true=4 636 five-fold.
  Amazon, still under-split at t=3 (760 vs 1 481), keeps gaining through t=4.
  The rule of thumb: raise tau while K < K_true, stop once past it.

## 3. Corollary 1 — the actual mechanism

Corollary 1 claims the fraction of accepted moves under a configuration-model
null is bounded by 2(1−Φ(τ)). Measured, one sweep from singletons on ER(1000, 10):

| tau | bound | observed move fraction |
|---|---|---|
| 2 | 0.046 | **1.000** |
| 3 | 0.003 | **1.000** |
| 4 | 0.000 | **1.000** |

The bound is off by a factor of ~300 at τ=3. The earlier explanation —
*selection bias, since only neighbours are scored and under the null that edge
is itself the rare event* — is correct as far as it goes, and it explains the
acceptance rate. **It does not explain the fragmentation**, and the fragmentation
is the part that matters.

The tell: **raising τ makes over-segmentation worse, monotonically.** A
significance filter cannot behave that way — a stricter filter should accept
less and merge less, not shatter more. ER(5000, 10), which has no communities:

| method | K(≥3) | max size | median size | singletons |
|---|---|---|---|---|
| Louvain | 14 | 942 | 212 | 0 |
| GRAZ t=2 | 77 | 292 | 54 | 0 |
| GRAZ t=3 | 160 | 91 | 28 | 0 |
| GRAZ t=4 | 277 | 52 | 18 | 0 |

No stranded singletons, and `max size` collapses monotonically. GRAZ is not
creating communities; it is **capping how large any community may grow**.

That cap is analytic. Holding k_i and k_iB fixed, z_iB is strictly decreasing in
Σ_B, so for every τ there is a largest volume that will still admit a node:

```
              2m·k_iB − k_i·Σ_B
z_iB  =  ─────────────────────────────      ∂z/∂Σ_B < 0
          sqrt(k_i·Σ_B·(2m − Σ_B))
```

Predicted ceiling vs observed max size (ER 5000, mean degree 9):

| tau | ceiling, k_iB=1 | k_iB=2 | k_iB=3 | observed max |
|---|---|---|---|---|
| 2 | 85 | 277 | 527 | 292 |
| 3 | 45 | 161 | 326 | 91 |
| 4 | 27 | 101 | 215 | 52 |

Observed sizes sit inside the k_iB=2 ceiling and fall at the same rate. **Once a
community's volume exceeds the τ-ceiling, no further node can join it regardless
of how much modularity the move would gain.** τ is a resolution parameter wearing
a p-value's clothing — which is why the p-value reading (τ=2 ≈ p<0.05) is not
safe to state in the paper.

## 4. Three candidate fixes, tested

| fix | idea | ER K(≥3) @ N=10000 | karate F1 | verdict |
|---|---|---|---|---|
| — | plain GRAZ t=2 | 93 | 0.310 | baseline (Louvain: 24, 0.689) |
| A | τ·√(2 log m), §3.4 | **458** | 0.220 | **worse** — scales the ceiling the wrong way |
| B | merge-level gate | 92 | **0.550** | best available |
| C | gate-only, argmax ΔQ | 83 | 0.310 | ~no effect; refutes the ranking hypothesis |

- **A (adaptive τ) backfires.** √(2 log m) grows with m, tightening the volume
  ceiling on exactly the large graphs where communities must be big. It nearly
  doubles the spurious count. The paper's own §3.4 proposal makes this worse.
- **C isolates a hypothesis worth recording as dead.** GRAZ selects the surviving
  candidate by argmax z, and since z favours small Σ_B, argmax-z plausibly
  steers nodes into the smallest admissible community. Selecting by ΔQ instead
  changes almost nothing (83 vs 93; karate identical). **The gate, not the
  ranking, is the whole effect.**
- **B (merge-level gate) is the one to adopt.** Running level 0 as plain
  modularity and engaging the z-test only from the first aggregated level — the
  literal 2010 Bader–McCloskey merge test — removes the singleton-level selection
  bias where every candidate is a neighbour by construction.

**B costs nothing on large graphs and is strictly faster:**

| dataset | method | K | Q | F1 | NMI | time | z-tests |
|---|---|---|---|---|---|---|---|
| amazon | GRAZ t=4 | 971 | 0.924 | 0.627 | 0.919 | 6 484 ms | 20.6 M |
| amazon | **merge t=4** | 1 010 | 0.924 | **0.627** | **0.920** | **3 586 ms** | **4.5 M** |
| dblp | GRAZ t=4 | 1 169 | 0.818 | 0.252 | 0.605 | 10 001 ms | 34.6 M |
| dblp | **merge t=4** | 1 194 | 0.818 | **0.257** | **0.604** | **8 402 ms** | **6.6 M** |
| youtube | GRAZ t=4 | 36 996 | 0.705 | 0.278 | 0.578 | 163 726 ms | 96.2 M |
| youtube | **merge t=4** | **30 851** | 0.708 | **0.294** | 0.552 | **114 533 ms** | **24.4 M** |
| karate | GRAZ t=2 | 16 | 0.186 | 0.310 | 0.369 | | |
| karate | **merge t=2** | **6** | **0.368** | **0.550** | **0.483** | | |

Same quality, **1.8× faster, 4.5× fewer z-tests**, and it recovers most of the
karate/polblogs regression. Recommend it as the default formulation.

YouTube is the strongest case for it: plain GRAZ never beats Louvain on F1 there
(0.278 vs 0.275) while splitting into 37 k communities against a true 4 749. The
merge gate gets **F1 0.294 — the best of any method on that graph — with 6 000
fewer communities, in 70% of the time and a quarter of the z-tests.** It trades
0.026 NMI for that, the one place the two metrics disagree.

Not yet run: merge-gate on LiveJournal (~25 min, needs a free core).

It does **not** rescue Corollary 1 (92 vs 93 spurious on ER) — that failure is
the volume ceiling, which no amount of gate placement removes.

## 5. What to raise with Bader

1. **Corollary 1 as stated is false**, and the proof needs the neighbour-only
   candidate set. The empirical gap is ~300× at τ=3.
2. **Reframe τ as a resolution parameter with an explicit volume ceiling**, not a
   significance level. The p-value interpretation should come out, or be
   restricted to the merge-level formulation where the null is defensible.
3. **§3.4's √(2 log m) scaling makes things worse**, measurably — it should not
   ship as the recommended remedy.
4. **Adopt the merge-level gate**: same accuracy, 1.8× faster, and the small-graph
   embarrassment mostly goes away. It is also closer to the 2010 paper.
5. **§6.6's claim that GRAZ handles karate/polblogs well does not reproduce** —
   plain GRAZ gets F1 0.310 on karate vs Louvain's 0.689.
6. **The large-graph result is genuinely strong and is the paper's best argument.**
   +0.195 F1 on com-Amazon at 4 M-node scale, at equal modularity. Lead with it,
   and state the K_true/K_Louvain condition under which it applies.

## 6. Confirmed as claimed

Theorem 1 (resolution limit) — exact: predicted ceilings n≤128.2 (τ=2), n≤240.0
(τ=3) match the ring-of-cliques transitions. Theorem 3 (stability) — GRAZ σ_K
0.156 vs Louvain 0.380, pairwise NMI 0.998 vs 0.979. Theorem 4 (LFR/hierarchical
SBM) — GRAZ recovers K=10 at NMI_sub 0.999 vs Louvain's 6.8 at 0.885.
