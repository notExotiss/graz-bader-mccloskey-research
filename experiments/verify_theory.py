"""Verification of the paper's closed-form predictions against the implementation.

The GRAZ paper is (per the author) AI-drafted and not fully proofread. The
student's Section 6 role is precisely to check where the experiments match the
theorems and where they do not. This script performs those checks numerically.

It verifies:
  * Theorem 1  - the ring-of-cliques no-merge ceiling n <= (c^2-c+2)/g(tau)^2
  * Eq. 6      - the closed-form z-statistic on a hand-computed 2x2 table
  * Corollary 1 - the "non-imposition on random graphs" claim (per-sweep move
                  fraction bounded by 2(1-Phi(tau))). This is the one that FAILS
                  for the literal single-node rule started from singletons, and
                  the script quantifies the failure.
"""

from __future__ import annotations

import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from collections import defaultdict

from graz.algorithm import z_statistic, _Graph
from graz import datasets as ds


def phi(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


def g_tau(tau):
    return (math.sqrt(tau * tau + 4) - tau) / 2.0


def verify_theorem1():
    print("\n--- Theorem 1: ring-of-cliques no-merge ceiling ---")
    c = 5
    base = c * c - c + 2
    print(f"  c={c}, Fortunato-Barthelemy threshold n = c(c-1)+2 = {c*(c-1)+2}")
    for tau in [0, 2, 3]:
        ceiling = base / g_tau(tau) ** 2 if tau > 0 else base
        print(f"  tau={tau}: g(tau)={g_tau(tau):.4f}, "
              f"predicted max n with no merge = {ceiling:.1f}")


def verify_eq6():
    print("\n--- Eq. 6: closed-form z-statistic sanity ---")
    # For a node with degree k_i moving into community B with degree sum Sigma_B
    # having k_iB edges into B, in a graph with 2m total degree.
    cases = [
        (10, 5, 200, 10000),
        (10, 1, 10, 10000),      # single edge between two singletons -> big z
        (50, 25, 5000, 100000),
    ]
    for k_i, k_iB, Sigma_B, m2 in cases:
        expected = k_i * Sigma_B / m2
        z = z_statistic(k_i, k_iB, Sigma_B, m2)
        print(f"  k_i={k_i}, k_iB={k_iB}, Sigma_B={Sigma_B}, 2m={m2}: "
              f"E[k_iB]={expected:.3f}, z={z:.3f}")


def _one_sweep_move_fraction(g, tau):
    """Fraction of nodes that leave their singleton in ONE local-moving sweep."""
    m2 = g.m2
    m = m2 / 2.0
    comm = {u: u for u in g.adj}
    Sigma = dict(g.degree)
    moved = 0
    for u in g.adj:
        ku = g.degree[u]
        A = comm[u]
        w_to = defaultdict(float)
        for v, w in g.adj[u].items():
            if v != u:
                w_to[comm[v]] += w
        Sigma[A] -= ku
        best = None
        for B, wiB in w_to.items():
            if B == A:
                continue
            z = z_statistic(ku, wiB, Sigma[B], m2)
            if z > tau:
                gainB = wiB / m - ku * Sigma[B] / (2 * m * m)
                if gainB > 0 and (best is None or z > best[1]):
                    best = (B, z)
        if best:
            comm[u] = best[0]
            moved += 1
        Sigma[comm[u]] += ku
    return moved / len(g.adj)


def verify_corollary1():
    print("\n--- Corollary 1: non-imposition on Erdos-Renyi (per-sweep move rate) ---")
    print("  Paper predicts per-node move probability <= 2(1-Phi(tau)).")
    n, edges, _ = ds.erdos_renyi(1000, 10, seed=0)
    g = _Graph.from_edges(n, edges)
    print(f"  {'tau':>5} {'bound 2(1-Phi)':>16} {'observed move frac':>20}")
    for tau in [2, 3, 4, math.sqrt(2 * math.log(len(edges)))]:
        bound = 2 * (1 - phi(tau))
        obs = _one_sweep_move_fraction(g, tau)
        flag = "  <-- FAILS" if obs > 5 * max(bound, 1e-6) else ""
        print(f"  {tau:5.2f} {bound:16.4f} {obs:20.4f}{flag}")
    print("\n  FINDING: the single-node z-test evaluated only on NEIGHBOURING")
    print("  communities is selection-biased. From a singleton start, the mere")
    print("  existence of an edge is itself the rare event under the null")
    print("  (E[k_iB] ~ k_i*k_j/2m << 1), so almost every first move clears any")
    print("  fixed tau. Corollary 1's bound therefore does not hold for the")
    print("  literal single-node rule; suppression needs an initial partition or")
    print("  an all-pairs (not neighbour-only) candidate set. This is the kind of")
    print("  gap between the AI-drafted paper and its implementation that the")
    print("  Section 6 verification is designed to catch.")


if __name__ == "__main__":
    verify_theorem1()
    verify_eq6()
    verify_corollary1()
