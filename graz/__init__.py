"""GRAZ: Greedy Residual-Adjusted Z-test for community detection.

An open-source implementation of the Bader-McCloskey statistical correction
applied to modularity-based community detection (Louvain / Leiden).

Reference:
    David A. Bader. "GRAZ: A Statistically Calibrated Resolution-Limit-Free
    Replacement for Louvain and Leiden." NJIT.
    Builds on: D. A. Bader and J. McCloskey, "Modularity and Graph Algorithms,"
    SIAM AN10 Minisymposium on Analyzing Massive Real-World Graphs, 12 July 2010.
"""

from .algorithm import (
    graz_communities,
    louvain_communities,
    leiden_communities,
    cpm_leiden_communities,
    z_statistic,
    modularity,
)

__all__ = [
    "graz_communities",
    "louvain_communities",
    "leiden_communities",
    "cpm_leiden_communities",
    "z_statistic",
    "modularity",
]

__version__ = "1.0.0"
