"""Shared helpers for the Section 6 experiments."""

from __future__ import annotations

import csv
import os
import sys
import time

# make the graz package importable when running scripts directly
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib

matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt  # noqa: E402

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
FIG_DIR = os.path.join(RESULTS_DIR, "figures")
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

# consistent colours per method across all figures
COLORS = {
    "Louvain": "#d62728",
    "Leiden": "#ff7f0e",
    "CPM-Leiden": "#2ca02c",
    "GRAZ (t=2)": "#1f77b4",
    "GRAZ (t=3)": "#9467bd",
    "Ground truth": "#333333",
}


def save_csv(name, header, rows):
    path = os.path.join(RESULTS_DIR, name)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    print(f"  wrote {path}")
    return path


def save_fig(fig, name):
    path = os.path.join(FIG_DIR, name)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {path}")
    return path


class Timer:
    def __enter__(self):
        self.t = time.perf_counter()
        return self

    def __exit__(self, *a):
        self.dt = time.perf_counter() - self.t
