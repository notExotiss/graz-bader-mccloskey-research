# GRAZ — Bader–McCloskey Correction for Louvain / Leiden

An open-source implementation of **GRAZ (Greedy Residual-Adjusted Z-test)**: the
Bader–McCloskey statistical significance correction applied to modularity-based
community detection. GRAZ is a drop-in replacement for the move-acceptance rule
of Louvain and Leiden that gates every candidate node move by the standardized
Pearson residual of its move under the configuration-model null.

> Based on David A. Bader, *"GRAZ: A Statistically Calibrated Resolution-Limit-Free
> Replacement for Louvain and Leiden"* (NJIT), which extends the 2010 sketch of
> D. A. Bader and J. McCloskey, *"Modularity and Graph Algorithms"* (SIAM AN10).

## The one-line idea

Louvain/Leiden accept a node move if **modularity goes up (ΔQ > 0)** — which is
secretly a statistical test with the significance bar set to **zero**. GRAZ adds a
real z-test, so a move must *also* be statistically significant:

```
accept move i -> B   iff   z_iB > tau   AND   dQ > 0
```

where (Eq. 6 of the paper)

```
            2m * k_iB  -  k_i * Sigma_B
z_iB  =  ---------------------------------------
          sqrt( k_i * Sigma_B * (2m - Sigma_B) )
```

`tau` is the single knob. The paper reads it as a p-value (`tau = 2` ≈ p < 0.05,
`tau = 3` ≈ p < 0.001), but **that interpretation does not survive measurement** —
see [results/FINDINGS.md](results/FINDINGS.md). Because `z_iB` is strictly
decreasing in `Sigma_B`, `tau` acts as a ceiling on community *volume*, so it
behaves as a resolution parameter rather than a significance level. Corollary 1's
non-imposition bound fails by ~300× at `tau = 3`.

## Install

```bash
pip install networkx python-igraph scikit-learn matplotlib numpy scipy
```

## Use

```python
from graz import graz_communities, louvain_communities, leiden_communities

# edges: list of (u, v) integer pairs; nodes 0..N-1
labels = graz_communities(n_nodes, edges, tau=2.0, seed=0)   # dict node -> community
```

## Reproduce the Section 6 experiments

```bash
python experiments/run_all.py          # full suite (Exp 1-6)
python experiments/run_all.py 1 2 5    # selected experiments
```

Outputs land in `results/` (CSV + `summary.json`) and `results/figures/` (PNG).

## Experiments

| # | Name | What it shows |
|---|------|---------------|
| 1 | Ring of cliques | Resolution-limit mitigation (Theorem 1) — exact prediction |
| 2 | Erdős–Rényi nulls | Spurious-community suppression (Corollary 1) |
| 3 | LFR benchmarks | Locally adaptive resolution (Theorem 4) |
| 4 | Hierarchical SBM | Multi-scale resolution |
| 5 | Stability | Run-to-run stability under node reordering (Theorem 3) |
| 6 | Real networks | Karate, PolBlogs, and 3 SNAP datasets with ground truth |
| 7 | Large real networks | Full com-Amazon/DBLP/YouTube/LiveJournal (to 4M nodes), quality + timing (C++ only) |
| 8 | Corollary 1 fixes | Adaptive tau, merge-level gate, gate-only ranking (C++ only) |

Experiments 7 and 8 live in `cpp/graz.cpp`; see
[results/FINDINGS.md](results/FINDINGS.md) for what they show.

## Layout

```
graz/
  algorithm.py   core multilevel optimizer + Louvain/Leiden/CPM/GRAZ move rules
  datasets.py    synthetic generators + SNAP downloaders
  metrics.py     NMI, F1, community-count, size distributions
experiments/
  run_all.py     the six Section 6 experiments
  common.py      shared plotting / CSV helpers
data/            cached SNAP downloads
results/         CSVs, summary.json, figures/
```
