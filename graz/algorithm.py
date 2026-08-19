"""Core community-detection algorithms.

This module implements a single multilevel (Louvain-style) optimizer whose
move-acceptance rule is pluggable, giving four methods that run on identical
machinery so comparisons are apples-to-apples:

    * louvain_communities      - accept move iff dQ > 0                (Blondel 2008)
    * leiden_communities       - Louvain + connectivity refinement      (Traag 2019)
    * cpm_leiden_communities   - Constant Potts Model objective         (Traag 2011)
    * graz_communities         - accept move iff z > tau AND dQ > 0     (Bader, GRAZ)

The GRAZ acceptance rule is the *only* structural change from Louvain. Every
candidate move of node i into community B is gated by the standardized Pearson
residual of the move under the configuration-model null (Eq. 6 of the paper):

    z_iB = (2m * k_iB  -  k_i * Sigma_B)
           -----------------------------------------
             sqrt( k_i * Sigma_B * (2m - Sigma_B) )

where k_i is the (weighted) degree of i, k_iB is the weight of edges from i into
B, Sigma_B is the degree sum of B excluding i, and 2m is the sum of all degrees.
This is the single-node-move generalization of the 2010 Bader-McCloskey
standardized-residual test S_ij = R_ij / std(R_ij).

Graphs are handled as weighted, undirected adjacency dicts so that the
aggregation phase (super-nodes with self-loops) works uniformly.
"""

from __future__ import annotations

import math
import random
from collections import defaultdict


# --------------------------------------------------------------------------- #
#  Internal weighted-graph representation
# --------------------------------------------------------------------------- #
class _Graph:
    """Weighted undirected graph with self-loops (used for aggregated levels).

    adj[u] maps neighbour -> weight. A self-loop adj[u][u] = w contributes 2w to
    the degree of u (an internal edge of an aggregated super-node).
    """

    __slots__ = ("adj", "degree", "m2", "size")

    def __init__(self, adj, size=None):
        self.adj = adj
        # size[u] = number of ORIGINAL nodes that super-node u represents (for CPM)
        self.size = size if size is not None else {u: 1 for u in adj}
        self.degree = {}
        for u, nbrs in adj.items():
            d = 0.0
            for v, w in nbrs.items():
                d += 2.0 * w if v == u else w
            self.degree[u] = d
        self.m2 = sum(self.degree.values())  # 2m = sum of degrees

    @classmethod
    def from_edges(cls, n_nodes, edges, weights=None):
        adj = {i: {} for i in range(n_nodes)}
        if weights is None:
            weights = [1.0] * len(edges)
        for (u, v), w in zip(edges, weights):
            if u == v:
                adj[u][u] = adj[u].get(u, 0.0) + w
            else:
                adj[u][v] = adj[u].get(v, 0.0) + w
                adj[v][u] = adj[v].get(u, 0.0) + w
        return cls(adj)


# --------------------------------------------------------------------------- #
#  Standalone z-statistic (Eq. 6) - also used by the verification experiments
# --------------------------------------------------------------------------- #
def z_statistic(k_i, k_iB, Sigma_B, m2):
    """Standardized Pearson residual of the single-node move (Eq. 6).

    Parameters
    ----------
    k_i     : (weighted) degree of the moving node i
    k_iB    : weight of edges from i into community B
    Sigma_B : degree sum of community B, EXCLUDING node i
    m2      : 2m, the sum of all degrees in the graph

    Returns z_iB, asymptotically N(0,1) under the configuration-model null.
    """
    denom = k_i * Sigma_B * (m2 - Sigma_B)
    if denom <= 0.0:
        return 0.0
    return (m2 * k_iB - k_i * Sigma_B) / math.sqrt(denom)


# --------------------------------------------------------------------------- #
#  Modularity of a partition (Eq. 1 / community-additive Eq. 2)
# --------------------------------------------------------------------------- #
def modularity(graph_adj, labels, resolution=1.0):
    """Newman-Girvan modularity of a labelling on a raw (unweighted) adjacency.

    graph_adj : dict u -> {v: w}
    labels    : dict u -> community id  (or list indexed by node)
    """
    g = graph_adj if isinstance(graph_adj, _Graph) else _Graph(graph_adj)
    m2 = g.m2
    if m2 == 0:
        return 0.0
    L_in = defaultdict(float)   # 2 * internal edges of community
    D = defaultdict(float)      # degree sum of community
    for u, nbrs in g.adj.items():
        cu = labels[u]
        D[cu] += g.degree[u]
        for v, w in nbrs.items():
            if labels[v] == cu:
                L_in[cu] += 2.0 * w if v == u else w
    q = 0.0
    for c in D:
        q += L_in[c] / m2 - resolution * (D[c] / m2) ** 2
    return q


# --------------------------------------------------------------------------- #
#  One level of local moving (the pluggable heart of every method)
# --------------------------------------------------------------------------- #
def _local_moving(g, rng, rule, tau, gamma):
    """Run local-moving sweeps to convergence on graph g.

    rule in {"modularity", "graz", "cpm"}.
    Returns comm: dict node -> community id.
    """
    m2 = g.m2
    m = m2 / 2.0
    comm = {u: u for u in g.adj}          # start: every node its own community
    Sigma = dict(g.degree)                # community degree sums (volume)
    n_nodes = {u: g.size[u] for u in g.adj}  # community node counts (for CPM)

    nodes = list(g.adj)
    improved_any = True
    while improved_any:
        improved_any = False
        rng.shuffle(nodes)
        for u in nodes:
            ku = g.degree[u]
            su = g.size[u]
            A = comm[u]

            # weight from u into each neighbouring community (self-loop excluded)
            w_to = defaultdict(float)
            for v, w in g.adj[u].items():
                if v == u:
                    continue
                w_to[comm[v]] += w

            # remove u from A
            Sigma[A] -= ku
            n_nodes[A] -= su

            w_to_A = w_to.get(A, 0.0)

            if rule == "cpm":
                # CPM gain of (re)inserting u (size su) into community C:
                #   w_uC - gamma * su * n_C   (n_C excludes u after removal)
                best_c, best_gain = A, w_to_A - gamma * su * n_nodes[A]
                for c, wc in w_to.items():
                    gain = wc - gamma * su * n_nodes[c]
                    if gain > best_gain + 1e-15:
                        best_c, best_gain = c, gain
                target = best_c

            elif rule == "graz":
                # Baseline: gain of staying in A (reinserting into A)
                gain_A = w_to_A / m - ku * Sigma[A] / (2.0 * m * m) if m > 0 else 0.0
                best_c, best_z = A, None
                for B, w_iB in w_to.items():
                    if B == A:
                        continue
                    SB = Sigma[B]
                    z = z_statistic(ku, w_iB, SB, m2)
                    if z <= tau:
                        continue
                    gain_B = w_iB / m - ku * SB / (2.0 * m * m) if m > 0 else 0.0
                    dQ = gain_B - gain_A                 # full move: leave A, join B
                    if dQ > 0.0 and (best_z is None or z > best_z):
                        best_c, best_z = B, z
                target = best_c

            else:  # "modularity"  (Louvain / Leiden local moving)
                best_c = A
                best_gain = w_to_A / m - ku * Sigma[A] / (2.0 * m * m) if m > 0 else 0.0
                for c, wc in w_to.items():
                    gain = wc / m - ku * Sigma[c] / (2.0 * m * m) if m > 0 else 0.0
                    if gain > best_gain + 1e-15:
                        best_c, best_gain = c, gain
                target = best_c

            # (re)insert u into target
            comm[u] = target
            Sigma[target] += ku
            n_nodes[target] += su
            if target != A:
                improved_any = True
    return comm


# --------------------------------------------------------------------------- #
#  Aggregation: collapse communities into super-nodes
# --------------------------------------------------------------------------- #
def _aggregate(g, comm):
    """Build the aggregated graph; return (new_graph, mapping old_comm -> new_id)."""
    ids = sorted(set(comm.values()))
    remap = {c: i for i, c in enumerate(ids)}
    new_adj = {i: {} for i in range(len(ids))}
    new_size = {i: 0 for i in range(len(ids))}
    for u in g.adj:
        new_size[remap[comm[u]]] += g.size[u]
    for u, nbrs in g.adj.items():
        cu = remap[comm[u]]
        for v, w in nbrs.items():
            cv = remap[comm[v]]
            if v == u:  # existing self-loop: internal weight preserved
                new_adj[cu][cu] = new_adj[cu].get(cu, 0.0) + w
            elif cu == cv:  # edge internal to a community -> half becomes self-loop
                new_adj[cu][cu] = new_adj[cu].get(cu, 0.0) + w / 2.0
            else:
                new_adj[cu][cv] = new_adj[cu].get(cv, 0.0) + w
    return _Graph(new_adj, size=new_size), remap


# --------------------------------------------------------------------------- #
#  Connectivity refinement (Leiden's headline guarantee, simplified)
# --------------------------------------------------------------------------- #
def _split_disconnected(base_adj, labels):
    """Split any internally-disconnected community into its connected components.

    This is the defect Leiden fixed in Louvain: a community whose induced
    subgraph is disconnected is broken into connected pieces.
    """
    from collections import deque

    members = defaultdict(list)
    for u, c in labels.items():
        members[c].append(u)

    new_labels = {}
    next_id = 0
    for c, nodes in members.items():
        nodeset = set(nodes)
        seen = set()
        for start in nodes:
            if start in seen:
                continue
            # BFS within the community's induced subgraph
            comp = []
            dq = deque([start])
            seen.add(start)
            while dq:
                x = dq.popleft()
                comp.append(x)
                for y in base_adj.get(x, {}):
                    if y in nodeset and y not in seen:
                        seen.add(y)
                        dq.append(y)
            for x in comp:
                new_labels[x] = next_id
            next_id += 1
    return new_labels


# --------------------------------------------------------------------------- #
#  Multilevel driver
# --------------------------------------------------------------------------- #
def _run_multilevel(base_graph, rule, tau, gamma, seed, refine):
    rng = random.Random(seed)
    g = base_graph
    # node_to_cur[original node] -> node id at the current (aggregated) level
    node_to_cur = {u: u for u in base_graph.adj}
    level = 0
    while True:
        comm = _local_moving(g, rng, rule, tau, gamma)
        n_comm = len(set(comm.values()))
        if n_comm == len(g.adj):
            break  # no coarsening happened -> converged
        g, remap = _aggregate(g, comm)
        # push each original node to its new super-node id: cur -> comm[cur] -> remap
        node_to_cur = {u: remap[comm[cur]] for u, cur in node_to_cur.items()}
        level += 1
        if level > 100:
            break

    if refine:
        node_to_cur = _split_disconnected(base_graph.adj, node_to_cur)

    # normalise labels to 0..K-1
    ids = sorted(set(node_to_cur.values()))
    remap = {c: i for i, c in enumerate(ids)}
    return {u: remap[c] for u, c in node_to_cur.items()}


# --------------------------------------------------------------------------- #
#  Public API
# --------------------------------------------------------------------------- #
def _to_graph(n_nodes, edges, weights=None):
    return _Graph.from_edges(n_nodes, edges, weights)


def louvain_communities(n_nodes, edges, weights=None, seed=0):
    """Vanilla Louvain (Blondel et al. 2008): accept a move iff dQ > 0."""
    g = _to_graph(n_nodes, edges, weights)
    return _run_multilevel(g, "modularity", tau=0.0, gamma=1.0, seed=seed, refine=False)


def leiden_communities(n_nodes, edges, weights=None, seed=0):
    """Leiden (Traag et al. 2019): Louvain + connected-community refinement."""
    g = _to_graph(n_nodes, edges, weights)
    return _run_multilevel(g, "modularity", tau=0.0, gamma=1.0, seed=seed, refine=True)


def cpm_leiden_communities(n_nodes, edges, weights=None, gamma=0.05, seed=0):
    """CPM-Leiden (Traag et al. 2011): Constant Potts Model objective.

    gamma is the density-threshold resolution parameter (tuned by grid search).
    """
    g = _to_graph(n_nodes, edges, weights)
    return _run_multilevel(g, "cpm", tau=0.0, gamma=gamma, seed=seed, refine=True)


def graz_communities(n_nodes, edges, weights=None, tau=2.0, seed=0, refine=True):
    """GRAZ (Bader): accept a move iff z_iB > tau AND dQ > 0.

    tau is the single hyperparameter, interpretable as a significance level:
        tau = 2  ->  nominal p < 0.05   (default)
        tau = 3  ->  nominal p < 0.0013 (strict)
    """
    g = _to_graph(n_nodes, edges, weights)
    return _run_multilevel(g, "graz", tau=tau, gamma=1.0, seed=seed, refine=refine)
