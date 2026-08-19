"""Dataset generators and loaders for the GRAZ experiments.

Synthetic generators (self-contained, no download):
    * ring_of_cliques        - Fortunato-Barthelemy resolution-limit pathology
    * erdos_renyi            - pure-random null (no community structure)
    * lfr_benchmark          - Lancichinetti-Fortunato-Radicchi planted communities
    * hierarchical_sbm       - two-level stochastic block model
    * karate_club            - Zachary's karate club (built into networkx)
    * polblogs               - Adamic-Glance political blogs (largest CC, giant component)

SNAP loaders (download + cache under data/):
    * load_snap_email_eu_core   - 1005 nodes, department ground truth
    * load_snap_com_amazon      - co-purchase, top-5000 ground-truth communities
    * load_snap_com_dblp        - co-authorship, top-5000 ground-truth communities

All generators return (n_nodes, edges, ground_truth) where:
    n_nodes      : int
    edges        : list[(u, v)] with 0-indexed integer node ids
    ground_truth : dict node -> community id  (or None if no ground truth)
"""

from __future__ import annotations

import gzip
import os
import urllib.request

import networkx as nx
import numpy as np

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")


# --------------------------------------------------------------------------- #
#  Helpers
# --------------------------------------------------------------------------- #
def _relabel(G):
    """Return (n, edges, mapping) with contiguous integer node ids."""
    mapping = {node: i for i, node in enumerate(G.nodes())}
    edges = [(mapping[u], mapping[v]) for u, v in G.edges()]
    return G.number_of_nodes(), edges, mapping


# --------------------------------------------------------------------------- #
#  Synthetic generators
# --------------------------------------------------------------------------- #
def ring_of_cliques(n, c):
    """Ring of n cliques of size c, each joined to its neighbour by one edge.

    Ground truth: each clique is its own community (K = n).
    """
    G = nx.Graph()
    starts = []
    off = 0
    gt = {}
    for i in range(n):
        nodes = list(range(off, off + c))
        starts.append(nodes[0])
        for a in range(c):
            gt[nodes[a]] = i
            for b in range(a + 1, c):
                G.add_edge(nodes[a], nodes[b])
        off += c
    for i in range(n):
        G.add_edge(starts[i], starts[(i + 1) % n])
    return G.number_of_nodes(), list(G.edges()), gt


def erdos_renyi(N, avg_degree, seed=0):
    """G(N, p) with p = avg_degree / (N-1). No community structure (gt=None)."""
    p = avg_degree / (N - 1)
    G = nx.fast_gnp_random_graph(N, p, seed=seed)
    # keep isolated nodes out so degree sums behave
    G.remove_nodes_from(list(nx.isolates(G)))
    n, edges, _ = _relabel(G)
    return n, edges, None


def lfr_benchmark(N=1000, avg_degree=10, max_degree=50, mu=0.3,
                  tau1=3.0, tau2=1.5, min_community=20, max_community=100, seed=0):
    """LFR planted-partition benchmark with a wide community-size range.

    Ground truth is the planted community of each node. The paper specifies
    community-size exponent tau2 = 1; networkx requires tau2 > 1, so we use the
    smallest practical value (1.1) that preserves the wide, heavy-tailed size
    range the experiment depends on.
    """
    G = None
    last_err = None
    for attempt in range(25):  # LFR construction is stochastic; retry on failure
        try:
            G = nx.LFR_benchmark_graph(
                N, tau1, tau2, mu,
                average_degree=avg_degree, max_degree=max_degree,
                min_community=min_community, max_community=max_community,
                seed=seed + attempt, max_iters=5000,
            )
            break
        except nx.ExceededMaxIterations as e:
            last_err = e
            continue
    if G is None:
        raise RuntimeError(f"LFR construction failed after retries: {last_err}")
    G.remove_edges_from(nx.selfloop_edges(G))
    # planted community stored on each node as frozenset "community"
    comm_of = {}
    seen = {}
    nid = 0
    for node in G.nodes():
        fs = frozenset(G.nodes[node]["community"])
        if fs not in seen:
            seen[fs] = nid
            nid += 1
        comm_of[node] = seen[fs]
    n, edges, mapping = _relabel(G)
    gt = {mapping[node]: comm_of[node] for node in G.nodes()}
    return n, edges, gt


def hierarchical_sbm(seed=0):
    """Two-level hierarchical SBM (paper Section 6.4).

    K1 = 2 supercommunities of 250 nodes, each with K2 = 5 sub-communities of 50.
    Returns ground truth at the SUB-community level (10 blocks); the super-level
    (2 blocks) is available as the returned gt_super.
    """
    rng = np.random.default_rng(seed)
    n_super, n_sub, sub_size = 2, 5, 50
    p_sub, p_super, p_cross = 0.30, 0.05, 0.005
    N = n_super * n_sub * sub_size

    gt_sub = {}
    gt_super = {}
    node = 0
    blocks = []
    for s in range(n_super):
        for b in range(n_sub):
            members = list(range(node, node + sub_size))
            for m in members:
                gt_sub[m] = s * n_sub + b
                gt_super[m] = s
            blocks.append((s, b, members))
            node += sub_size

    edges = set()
    for i in range(N):
        for j in range(i + 1, N):
            si, sj = gt_super[i], gt_super[j]
            bi, bj = gt_sub[i], gt_sub[j]
            if bi == bj:
                p = p_sub
            elif si == sj:
                p = p_super
            else:
                p = p_cross
            if rng.random() < p:
                edges.add((i, j))
    return N, list(edges), gt_sub, gt_super


def karate_club():
    """Zachary's karate club (N=34). Ground truth: Mr. Hi vs Officer (K=2)."""
    G = nx.karate_club_graph()
    gt = {i: (0 if G.nodes[i]["club"] == "Mr. Hi" else 1) for i in G.nodes()}
    return G.number_of_nodes(), list(G.edges()), gt


def polblogs(path=None):
    """Adamic-Glance political blogs. Ground truth: liberal vs conservative (K=2).

    Uses a bundled edge list if available, else falls back to a synthetic
    2-block planted graph of comparable size so the experiment still runs.
    """
    # networkx has no built-in polblogs; try a cached GML, else synthesize.
    cached = path or os.path.join(DATA_DIR, "polblogs.gml")
    if os.path.exists(cached):
        G = nx.read_gml(cached)
        G = G.to_undirected()
        gt_raw = {n: int(G.nodes[n].get("value", 0)) for n in G.nodes()}
        G = G.subgraph(max(nx.connected_components(G), key=len)).copy()
        n, edges, mapping = _relabel(G)
        gt = {mapping[node]: gt_raw[node] for node in G.nodes()}
        return n, edges, gt
    # Fallback: planted 2-community graph (~1490 nodes) so the pipeline is testable
    return _planted_two_block(n=1490, k_in=12, k_out=2, seed=7)


def _planted_two_block(n=1490, k_in=12, k_out=2, seed=7):
    rng = np.random.default_rng(seed)
    half = n // 2
    gt = {i: (0 if i < half else 1) for i in range(n)}
    p_in = k_in / half
    p_out = k_out / half
    edges = set()
    for i in range(n):
        for j in range(i + 1, n):
            p = p_in if gt[i] == gt[j] else p_out
            if rng.random() < p:
                edges.add((i, j))
    return n, list(edges), gt


# --------------------------------------------------------------------------- #
#  SNAP loaders
# --------------------------------------------------------------------------- #
SNAP_URLS = {
    "email-Eu-core": "https://snap.stanford.edu/data/email-Eu-core.txt.gz",
    "email-Eu-core-labels": "https://snap.stanford.edu/data/email-Eu-core-department-labels.txt.gz",
    "com-amazon-ungraph": "https://snap.stanford.edu/data/bigdata/communities/com-amazon.ungraph.txt.gz",
    "com-amazon-top5000": "https://snap.stanford.edu/data/bigdata/communities/com-amazon.top5000.cmty.txt.gz",
    "com-dblp-ungraph": "https://snap.stanford.edu/data/bigdata/communities/com-dblp.ungraph.txt.gz",
    "com-dblp-top5000": "https://snap.stanford.edu/data/bigdata/communities/com-dblp.top5000.cmty.txt.gz",
}


def _download(key):
    os.makedirs(DATA_DIR, exist_ok=True)
    url = SNAP_URLS[key]
    dest = os.path.join(DATA_DIR, url.split("/")[-1])
    if not os.path.exists(dest):
        print(f"  downloading {url} ...")
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=120) as r, open(dest, "wb") as f:
            f.write(r.read())
    return dest


def _read_edges_gz(path):
    edges = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            a, b = line.split()[:2]
            edges.append((int(a), int(b)))
    return edges


def load_snap_email_eu_core():
    """SNAP email-Eu-core (~1005 nodes). Ground truth: 42 department labels."""
    epath = _download("email-Eu-core")
    lpath = _download("email-Eu-core-labels")
    raw_edges = _read_edges_gz(epath)
    labels = {}
    with gzip.open(lpath, "rt") as f:
        for line in f:
            if not line.strip():
                continue
            node, dept = line.split()
            labels[int(node)] = int(dept)
    G = nx.Graph()
    G.add_edges_from(raw_edges)
    G.remove_edges_from(nx.selfloop_edges(G))
    G = G.subgraph(max(nx.connected_components(G), key=len)).copy()
    n, edges, mapping = _relabel(G)
    gt = {mapping[node]: labels[node] for node in G.nodes() if node in labels}
    return n, edges, gt


def _load_snap_community(ungraph_key, cmty_key, max_nodes=None, top_communities=None,
                         min_component=3):
    """Generic SNAP ground-truth-community loader.

    To keep runtime laptop-friendly we restrict to the nodes covered by the
    top-`top_communities` ground-truth communities and take the induced subgraph.
    SNAP ground-truth communities are small and overlapping, so the induced
    subgraph is many small components rather than one giant one; we keep every
    connected component of at least `min_component` nodes (dropping only trivial
    dangling pairs), which preserves the ground-truth structure for F1 scoring.
    """
    gpath = _download(ungraph_key)
    cpath = _download(cmty_key)

    # read communities
    communities = []
    with gzip.open(cpath, "rt") as f:
        for line in f:
            if not line.strip():
                continue
            communities.append([int(x) for x in line.split()])
    if top_communities:
        communities = communities[:top_communities]

    keep = set()
    for com in communities:
        keep.update(com)

    raw_edges = _read_edges_gz(gpath)
    G = nx.Graph()
    if keep:
        for a, b in raw_edges:
            if a in keep and b in keep:
                G.add_edge(a, b)
    else:
        G.add_edges_from(raw_edges)
    G.remove_edges_from(nx.selfloop_edges(G))
    if G.number_of_nodes() == 0:
        raise RuntimeError("empty induced subgraph")

    # keep components of a meaningful size (drop dangling singletons / pairs)
    good_nodes = set()
    for comp in nx.connected_components(G):
        if len(comp) >= min_component:
            good_nodes.update(comp)
    G = G.subgraph(good_nodes).copy()

    if max_nodes and G.number_of_nodes() > max_nodes:
        # deterministically cap size: keep whole components in read order
        capped = set()
        for comp in sorted(nx.connected_components(G), key=len, reverse=True):
            if len(capped) + len(comp) > max_nodes:
                continue
            capped.update(comp)
        G = G.subgraph(capped).copy()

    # assign each kept node its FIRST covering community as ground truth
    node_comm = {}
    for cid, com in enumerate(communities):
        for node in com:
            if node in G and node not in node_comm:
                node_comm[node] = cid

    n, edges, mapping = _relabel(G)
    gt = {mapping[node]: node_comm[node] for node in G.nodes() if node in node_comm}
    return n, edges, gt


def load_snap_com_amazon(max_nodes=15000, top_communities=1000):
    """SNAP com-Amazon co-purchase, restricted to top ground-truth communities."""
    return _load_snap_community("com-amazon-ungraph", "com-amazon-top5000",
                                max_nodes=max_nodes, top_communities=top_communities)


def load_snap_com_dblp(max_nodes=15000, top_communities=1000):
    """SNAP com-DBLP co-authorship, restricted to top ground-truth communities."""
    return _load_snap_community("com-dblp-ungraph", "com-dblp-top5000",
                                max_nodes=max_nodes, top_communities=top_communities)
