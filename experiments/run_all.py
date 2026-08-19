"""Section 6 experiments for GRAZ (Bader-McCloskey correction on Louvain/Leiden).

Run:  python experiments/run_all.py            # full suite
      python experiments/run_all.py 1 2 5      # selected experiments

Each experiment writes a CSV to results/ and one or more figures to
results/figures/. A machine-readable results/summary.json aggregates the
headline numbers used to compile the report.
"""

from __future__ import annotations

import json
import math
import os
import statistics
import sys

import numpy as np

from common import COLORS, FIG_DIR, RESULTS_DIR, save_csv, save_fig, Timer
import matplotlib.pyplot as plt

from graz.algorithm import (
    louvain_communities,
    leiden_communities,
    cpm_leiden_communities,
    graz_communities,
    modularity,
)
from graz import datasets as ds
from graz import metrics as mx

SUMMARY = {}


def _adj_from_edges(n, edges):
    adj = {i: {} for i in range(n)}
    for u, v in edges:
        if u == v:
            continue
        adj[u][v] = 1.0
        adj[v][u] = 1.0
    return adj


def g_tau(tau):
    return (math.sqrt(tau * tau + 4) - tau) / 2.0


# ==========================================================================
#  Experiment 1: Ring of cliques (exact resolution-limit demonstration)
# ==========================================================================
def experiment_1():
    print("\n=== Experiment 1: Ring of cliques (resolution limit) ===")
    c = 5
    ns = [20, 25, 30, 50, 100, 200, 400]
    n_seeds = 5
    rows = []
    series = {m: [] for m in ["Louvain", "Leiden", "CPM-Leiden", "GRAZ (t=2)", "GRAZ (t=3)"]}
    nmi_series = {m: [] for m in series}

    for n in ns:
        N, edges, gt = ds.ring_of_cliques(n, c)
        counts = {m: [] for m in series}
        nmis = {m: [] for m in series}
        for seed in range(n_seeds):
            methods = {
                "Louvain": louvain_communities(N, edges, seed=seed),
                "Leiden": leiden_communities(N, edges, seed=seed),
                "CPM-Leiden": cpm_leiden_communities(N, edges, gamma=0.1, seed=seed),
                "GRAZ (t=2)": graz_communities(N, edges, tau=2.0, seed=seed),
                "GRAZ (t=3)": graz_communities(N, edges, tau=3.0, seed=seed),
            }
            for m, lab in methods.items():
                counts[m].append(mx.num_communities(lab))
                nmis[m].append(mx.nmi(lab, gt))
        # predicted GRAZ thresholds from Theorem 1: no-merge iff n <= 22/g(tau)^2
        thr2 = (c * c - c + 2) / g_tau(2.0) ** 2
        thr3 = (c * c - c + 2) / g_tau(3.0) ** 2
        row = [n, n]
        for m in series:
            mean_k = statistics.mean(counts[m])
            std_k = statistics.pstdev(counts[m])
            mean_nmi = statistics.mean(nmis[m])
            series[m].append((mean_k, std_k))
            nmi_series[m].append(mean_nmi)
            row.extend([round(mean_k, 1), round(mean_nmi, 3)])
        row.extend([round(thr2, 1), round(thr3, 1)])
        rows.append(row)
        print(f"  n={n:3d}  " + "  ".join(
            f"{m}={statistics.mean(counts[m]):.0f}" for m in series))

    header = ["n", "ground_truth_K"]
    for m in series:
        header += [f"{m}_K", f"{m}_NMI"]
    header += ["GRAZ_t2_predicted_max_n", "GRAZ_t3_predicted_max_n"]
    save_csv("exp1_ring_of_cliques.csv", header, rows)

    # figure: recovered K vs n
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(ns, ns, "--", color=COLORS["Ground truth"], label="Ground truth (K=n)")
    for m in series:
        means = [v[0] for v in series[m]]
        stds = [v[1] for v in series[m]]
        ax.errorbar(ns, means, yerr=stds, marker="o", capsize=3,
                    color=COLORS[m], label=m)
    ax.set_xlabel("Number of cliques n (size c=5)")
    ax.set_ylabel("Recovered community count K")
    ax.set_title("Experiment 1: Ring of cliques resolution limit")
    ax.legend()
    ax.grid(alpha=0.3)
    save_fig(fig, "exp1_ring_of_cliques.png")

    SUMMARY["exp1"] = {
        "c": c, "n_values": ns, "n_seeds": n_seeds,
        "recovered_K": {m: [v[0] for v in series[m]] for m in series},
        "nmi": {m: nmi_series[m] for m in series},
        "graz_t2_predicted_ceiling": round((c * c - c + 2) / g_tau(2.0) ** 2, 1),
        "graz_t3_predicted_ceiling": round((c * c - c + 2) / g_tau(3.0) ** 2, 1),
        "fb_threshold": c * (c - 1) + 2,
    }


# ==========================================================================
#  Experiment 2: Erdos-Renyi nulls (spurious-community suppression)
# ==========================================================================
def experiment_2():
    print("\n=== Experiment 2: Erdos-Renyi nulls (spurious communities) ===")
    Ns = [1000, 2000, 5000, 10000]
    avg_degree = 10
    n_seeds = 3
    rows = []
    series = {m: [] for m in ["Louvain", "Leiden", "GRAZ (t=2)", "GRAZ (t=3)"]}
    guimera = []

    for N in Ns:
        counts = {m: [] for m in series}
        m_edges_list = []
        for seed in range(n_seeds):
            n, edges, _ = ds.erdos_renyi(N, avg_degree, seed=seed)
            m_edges = len(edges)
            m_edges_list.append(m_edges)
            methods = {
                "Louvain": louvain_communities(n, edges, seed=seed),
                "Leiden": leiden_communities(n, edges, seed=seed),
                "GRAZ (t=2)": graz_communities(n, edges, tau=2.0, seed=seed),
                "GRAZ (t=3)": graz_communities(n, edges, tau=3.0, seed=seed),
            }
            for mth, lab in methods.items():
                # count only non-trivial communities (size >= 3) as "found structure"
                sizes = mx.community_sizes(lab)
                counts[mth].append(sum(1 for s in sizes if s >= 3))
        mbar = statistics.mean(m_edges_list)
        gpred = math.sqrt(mbar / 2.0)
        guimera.append(gpred)
        row = [N, round(mbar)]
        for mth in series:
            mk = statistics.mean(counts[mth])
            series[mth].append(mk)
            row.append(round(mk, 1))
        row.append(round(gpred, 1))
        rows.append(row)
        print(f"  N={N:5d} m={mbar:.0f}  " +
              "  ".join(f"{mth}={statistics.mean(counts[mth]):.0f}" for mth in series) +
              f"  Guimera~{gpred:.0f}")

    header = ["N", "m_edges"] + [f"{m}_K" for m in series] + ["Guimera_sqrt_m_over_2"]
    save_csv("exp2_erdos_renyi.csv", header, rows)

    # log-log figure
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(Ns, guimera, "--", color=COLORS["Ground truth"],
            label="Guimera prediction sqrt(m/2)")
    for m in series:
        ax.plot(Ns, [max(v, 0.5) for v in series[m]], marker="o",
                color=COLORS[m], label=m)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of nodes N")
    ax.set_ylabel("Recovered communities (size >= 3)")
    ax.set_title("Experiment 2: Spurious communities on random graphs")
    ax.legend()
    ax.grid(alpha=0.3, which="both")
    save_fig(fig, "exp2_erdos_renyi.png")

    SUMMARY["exp2"] = {
        "N_values": Ns, "avg_degree": avg_degree,
        "recovered_K": {m: series[m] for m in series},
        "guimera_prediction": guimera,
    }


# ==========================================================================
#  Experiment 3: Heterogeneous LFR (locally adaptive resolution)
# ==========================================================================
def experiment_3():
    print("\n=== Experiment 3: LFR benchmarks (locally adaptive resolution) ===")
    mus = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    N = 1000
    rows = []
    series = {m: [] for m in ["Louvain", "Leiden", "CPM-Leiden", "GRAZ (t=2)"]}

    for mu in mus:
        try:
            n, edges, gt = ds.lfr_benchmark(N=N, mu=mu, seed=1)
        except Exception as e:
            print(f"  mu={mu}: LFR failed ({e}); skipping")
            continue
        methods = {
            "Louvain": louvain_communities(n, edges, seed=1),
            "Leiden": leiden_communities(n, edges, seed=1),
            "CPM-Leiden": cpm_leiden_communities(n, edges, gamma=0.1, seed=1),
            "GRAZ (t=2)": graz_communities(n, edges, tau=2.0, seed=1),
        }
        row = [mu, len(set(gt.values()))]
        for m, lab in methods.items():
            val = mx.nmi(lab, gt)
            series[m].append(val)
            row.append(round(val, 3))
        rows.append(row)
        print(f"  mu={mu:.2f}  " + "  ".join(
            f"{m}={mx.nmi(lab,gt):.2f}" for m, lab in methods.items()))

    header = ["mu", "ground_truth_K"] + [f"{m}_NMI" for m in series]
    save_csv("exp3_lfr.csv", header, rows)

    fig, ax = plt.subplots(figsize=(7, 5))
    valid_mus = mus[:len(next(iter(series.values())))]
    for m in series:
        ax.plot(valid_mus, series[m], marker="o", color=COLORS[m], label=m)
    ax.set_xlabel("Mixing parameter mu")
    ax.set_ylabel("NMI vs ground truth")
    ax.set_title("Experiment 3: LFR benchmark (N=1000)")
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_ylim(0, 1.05)
    save_fig(fig, "exp3_lfr.png")

    SUMMARY["exp3"] = {"mu_values": valid_mus, "nmi": {m: series[m] for m in series}}


# ==========================================================================
#  Experiment 4: Hierarchical SBM (multi-scale resolution)
# ==========================================================================
def experiment_4():
    print("\n=== Experiment 4: Hierarchical SBM (multi-scale) ===")
    n_seeds = 5
    rows = []
    stats = {m: {"K": [], "nmi_sub": [], "nmi_super": []}
             for m in ["Louvain", "Leiden", "CPM-Leiden", "GRAZ (t=2)", "GRAZ (t=3)"]}

    for seed in range(n_seeds):
        N, edges, gt_sub, gt_super = ds.hierarchical_sbm(seed=seed)
        methods = {
            "Louvain": louvain_communities(N, edges, seed=seed),
            "Leiden": leiden_communities(N, edges, seed=seed),
            "CPM-Leiden": cpm_leiden_communities(N, edges, gamma=0.1, seed=seed),
            "GRAZ (t=2)": graz_communities(N, edges, tau=2.0, seed=seed),
            "GRAZ (t=3)": graz_communities(N, edges, tau=3.0, seed=seed),
        }
        for m, lab in methods.items():
            stats[m]["K"].append(mx.num_communities(lab))
            stats[m]["nmi_sub"].append(mx.nmi(lab, gt_sub))
            stats[m]["nmi_super"].append(mx.nmi(lab, gt_super))

    for m in stats:
        row = [m,
               round(statistics.mean(stats[m]["K"]), 1),
               round(statistics.pstdev(stats[m]["K"]), 2),
               round(statistics.mean(stats[m]["nmi_sub"]), 3),
               round(statistics.mean(stats[m]["nmi_super"]), 3)]
        rows.append(row)
        print(f"  {m:12s} K={row[1]:5}  NMI_sub={row[3]:.2f}  NMI_super={row[4]:.2f}")

    header = ["method", "mean_K", "std_K", "mean_NMI_subblocks", "mean_NMI_superblocks"]
    save_csv("exp4_hierarchical_sbm.csv", header, rows)

    fig, ax = plt.subplots(figsize=(7, 5))
    labels = list(stats.keys())
    xs = np.arange(len(labels))
    sub = [statistics.mean(stats[m]["nmi_sub"]) for m in labels]
    sup = [statistics.mean(stats[m]["nmi_super"]) for m in labels]
    ax.bar(xs - 0.2, sub, 0.4, label="NMI vs 10 sub-blocks", color="#1f77b4")
    ax.bar(xs + 0.2, sup, 0.4, label="NMI vs 2 super-blocks", color="#ff7f0e")
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel("NMI")
    ax.set_title("Experiment 4: Hierarchical SBM (2 x 5 x 50 nodes)")
    ax.legend()
    ax.grid(alpha=0.3, axis="y")
    save_fig(fig, "exp4_hierarchical_sbm.png")

    SUMMARY["exp4"] = {m: {k: round(statistics.mean(v), 3) for k, v in stats[m].items()}
                       for m in stats}


# ==========================================================================
#  Experiment 5: Stability under node-visit reordering
# ==========================================================================
def experiment_5():
    print("\n=== Experiment 5: Stability under node reordering ===")
    n_runs = 40
    # one LFR instance (per paper's minimum-viable recommendation)
    n, edges, gt = ds.lfr_benchmark(N=1000, mu=0.3, seed=3)
    method_fns = {
        "Louvain": lambda s: louvain_communities(n, edges, seed=s),
        "Leiden": lambda s: leiden_communities(n, edges, seed=s),
        "GRAZ (t=2)": lambda s: graz_communities(n, edges, tau=2.0, seed=s),
        "GRAZ (t=3)": lambda s: graz_communities(n, edges, tau=3.0, seed=s),
    }
    rows = []
    stab = {}
    for m, fn in method_fns.items():
        runs = [fn(s) for s in range(n_runs)]
        Ks = [mx.num_communities(r) for r in runs]
        # mean pairwise NMI across run pairs (sample to keep it cheap)
        import itertools
        pair_nmis = []
        for a, b in itertools.combinations(range(n_runs), 2):
            pred_a = [runs[a][u] for u in gt]
            pred_b = [runs[b][u] for u in gt]
            from sklearn.metrics import normalized_mutual_info_score
            pair_nmis.append(normalized_mutual_info_score(pred_a, pred_b))
        mean_pair = statistics.mean(pair_nmis)
        std_K = statistics.pstdev(Ks)
        stab[m] = {"mean_pairwise_nmi": mean_pair, "std_K": std_K,
                   "mean_K": statistics.mean(Ks)}
        rows.append([m, round(statistics.mean(Ks), 1), round(std_K, 3),
                     round(mean_pair, 4)])
        print(f"  {m:12s} meanK={statistics.mean(Ks):.1f}  stdK={std_K:.3f}  "
              f"meanPairwiseNMI={mean_pair:.4f}")

    save_csv("exp5_stability.csv",
             ["method", "mean_K", "std_K", "mean_pairwise_NMI"], rows)

    # stability ratio headline
    ratio_lei = stab["Leiden"]["std_K"] / max(stab["GRAZ (t=2)"]["std_K"], 1e-9)
    print(f"  >> std(K) reduction Leiden / GRAZ(t=2) = {ratio_lei:.2f}x")

    fig, ax = plt.subplots(figsize=(7, 5))
    labels = list(stab.keys())
    xs = np.arange(len(labels))
    ax.bar(xs, [stab[m]["std_K"] for m in labels],
           color=[COLORS[m] for m in labels])
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel("std of recovered community count over runs")
    ax.set_title(f"Experiment 5: Run-to-run instability ({n_runs} shuffles)")
    ax.grid(alpha=0.3, axis="y")
    save_fig(fig, "exp5_stability.png")

    SUMMARY["exp5"] = {"n_runs": n_runs, "stability": stab,
                       "std_K_ratio_leiden_over_graz2": round(ratio_lei, 2)}


# ==========================================================================
#  Experiment 6: Real networks with metadata ground truth
# ==========================================================================
def experiment_6():
    print("\n=== Experiment 6: Real networks ===")
    rows = []
    detail = {}

    datasets = []
    # small canonical
    datasets.append(("Karate", ds.karate_club(), 2))
    datasets.append(("PolBlogs", ds.polblogs(), 2))
    # SNAP (with ground truth)
    try:
        datasets.append(("SNAP email-Eu-core", ds.load_snap_email_eu_core(), None))
    except Exception as e:
        print(f"  email-Eu-core skipped: {e}")
    try:
        datasets.append(("SNAP com-Amazon", ds.load_snap_com_amazon(max_nodes=12000, top_communities=1500), None))
    except Exception as e:
        print(f"  com-Amazon skipped: {e}")
    try:
        datasets.append(("SNAP com-DBLP", ds.load_snap_com_dblp(max_nodes=12000, top_communities=1500), None))
    except Exception as e:
        print(f"  com-DBLP skipped: {e}")

    for name, data, truth_K in datasets:
        n, edges, gt = data
        adj = _adj_from_edges(n, edges)
        print(f"\n  -- {name} (N={n}, m={len(edges)}, "
              f"truth_K={len(set(gt.values())) if gt else '?'}) --")
        methods = {
            "Louvain": louvain_communities(n, edges, seed=1),
            "Leiden": leiden_communities(n, edges, seed=1),
            "GRAZ (t=2)": graz_communities(n, edges, tau=2.0, seed=1),
            "GRAZ (t=3)": graz_communities(n, edges, tau=3.0, seed=1),
        }
        detail[name] = {}
        for m, lab in methods.items():
            K = mx.num_communities(lab)
            Q = modularity(adj, lab)
            f1 = mx.community_f1(lab, gt) if gt else float("nan")
            nmi_v = mx.nmi(lab, gt) if gt else float("nan")
            rows.append([name, m, K, round(Q, 4),
                         round(f1, 4) if not math.isnan(f1) else "",
                         round(nmi_v, 4) if not math.isnan(nmi_v) else ""])
            detail[name][m] = {"K": K, "Q": round(Q, 4),
                               "F1": None if math.isnan(f1) else round(f1, 4),
                               "NMI": None if math.isnan(nmi_v) else round(nmi_v, 4)}
            print(f"     {m:12s} K={K:5d}  Q={Q:.3f}  "
                  f"F1={'n/a' if math.isnan(f1) else f'{f1:.3f}'}  "
                  f"NMI={'n/a' if math.isnan(nmi_v) else f'{nmi_v:.3f}'}")

    save_csv("exp6_real_networks.csv",
             ["dataset", "method", "K", "modularity_Q", "F1", "NMI"], rows)

    # figure: community-size distribution for the largest SNAP dataset (tail comparison)
    big = None
    for name, data, _ in datasets:
        if "DBLP" in name:
            big = (name, data)
    if big:
        name, (n, edges, gt) = big
        fig, ax = plt.subplots(figsize=(7, 5))
        for m, fn in [("Leiden", lambda: leiden_communities(n, edges, seed=1)),
                      ("GRAZ (t=2)", lambda: graz_communities(n, edges, tau=2.0, seed=1))]:
            sizes = sorted(mx.community_sizes(fn()), reverse=True)
            ax.plot(range(1, len(sizes) + 1), sizes, marker=".",
                    color=COLORS[m], label=m)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Community rank")
        ax.set_ylabel("Community size")
        ax.set_title(f"Experiment 6: Community-size distribution ({name})")
        ax.legend()
        ax.grid(alpha=0.3, which="both")
        save_fig(fig, "exp6_size_distribution.png")

    SUMMARY["exp6"] = detail


# ==========================================================================
#  Driver
# ==========================================================================
EXPERIMENTS = {
    1: experiment_1,
    2: experiment_2,
    3: experiment_3,
    4: experiment_4,
    5: experiment_5,
    6: experiment_6,
}


def main():
    which = [int(x) for x in sys.argv[1:]] or list(EXPERIMENTS)
    for i in which:
        with Timer() as t:
            EXPERIMENTS[i]()
        print(f"  [experiment {i} took {t.dt:.1f}s]")
    with open(os.path.join(RESULTS_DIR, "summary.json"), "w") as f:
        json.dump(SUMMARY, f, indent=2)
    print(f"\nAll done. Summary -> {os.path.join(RESULTS_DIR, 'summary.json')}")


if __name__ == "__main__":
    main()
