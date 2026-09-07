#!/usr/bin/env python3
"""Prepare reproducible round-4 LFR and Planetoid citation benchmarks.

Citation archives are the public DGL mirrors of the Planetoid-formatted Cora,
Citeseer, and PubMed datasets.  Output is the project's node-id edge list plus
one integer ground-truth label per node format.
"""
from __future__ import annotations

import argparse
import gzip
import pickle
import tarfile
import urllib.request
import zipfile
from pathlib import Path

import networkx as nx
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PREPARED = ROOT / "cpp" / "prepared"
CACHE = ROOT / ".round4_downloads"


def write_graph(name: str, n: int, edges, labels) -> None:
    PREPARED.mkdir(parents=True, exist_ok=True)
    clean = sorted({(min(int(u), int(v)), max(int(u), int(v)))
                    for u, v in edges if u != v})
    with (PREPARED / f"{name}.edges").open("w") as out:
        for u, v in clean:
            out.write(f"{u} {v}\n")
    with (PREPARED / f"{name}.labels").open("w") as out:
        for u in range(n):
            out.write(f"{u} {int(labels[u])}\n")
    print(f"{name}: N={n} M={len(clean)} K={len(set(map(int, labels)))}")


def prepare_lfr() -> None:
    for n in (1000, 5000):
        for mu in (0.1, 0.3, 0.5, 0.7):
            # LFR may reject a particular random degree/community draw.  The
            # retry is deterministic and leaves the successful seed in stdout.
            for attempt in range(20):
                seed = 44000 + n + int(mu * 1000) + attempt
                try:
                    g = nx.LFR_benchmark_graph(
                        n, tau1=3.0, tau2=1.5, mu=mu,
                        average_degree=20, max_degree=80,
                        min_community=20, max_community=min(500, n // 4),
                        seed=seed,
                    )
                    break
                except (nx.ExceededMaxIterations, nx.NetworkXError):
                    if attempt == 19:
                        raise
            communities = {}
            for node in sorted(g):
                c = tuple(sorted(g.nodes[node]["community"]))
                communities.setdefault(c, len(communities))
            labels = [communities[tuple(sorted(g.nodes[u]["community"]))] for u in range(n)]
            write_graph(f"lfr_n{n}_mu{int(mu * 10)}", n, g.edges(), labels)


def read_pickle(zf: zipfile.ZipFile, member: str):
    return pickle.loads(zf.read(member), encoding="latin1")


def prepare_planetoid(name: str, archive: str, prefix: str) -> None:
    path = CACHE / archive
    if not path.exists():
        CACHE.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(f"https://data.dgl.ai/dataset/{archive}", path)
    with zipfile.ZipFile(path) as zf:
        names = zf.namelist()
        graph = read_pickle(zf, next(x for x in names if x.endswith(".graph")))
        ally = read_pickle(zf, next(x for x in names if x.endswith(".ally")))
        ty = read_pickle(zf, next(x for x in names if x.endswith(".ty")))
        test_file = next(x for x in names if x.endswith(".test.index"))
        test_reorder = [int(x) for x in zf.read(test_file).decode().split()]
    # The familiar Planetoid index correction: allx/ally is contiguous except
    # for the test range, which must be restored to its original node ids.
    test_range = sorted(test_reorder)
    n = max(max(graph), max(test_range)) + 1
    labels = np.full(n, -1, dtype=int)
    for i, row in enumerate(ally):
        if i < n:
            labels[i] = int(np.argmax(row))
    for pos, row in zip(test_range, ty):
        labels[pos] = int(np.argmax(row))
    # Citeseer has 15 isolated test-index holes.  The original Planetoid loader
    # pads their missing `ty` rows with zeros; argmax of that conventional
    # all-zero class vector is 0.  Reproduce that documented preprocessing so
    # node ids and labels stay aligned rather than silently dropping nodes.
    missing = int(np.sum(labels < 0))
    if missing:
        if name != "citeseer":
            raise RuntimeError(f"{name}: {missing} nodes lack a class label")
        labels[labels < 0] = 0
        print(f"{name}: Planetoid zero-row padding applied to {missing} isolated nodes")
    edges = ((u, v) for u, nbrs in graph.items() for v in nbrs)
    write_graph(name, n, edges, labels)


def prepare_ego_facebook() -> None:
    """Project SNAP's overlapping circles to one deterministic label/node.

    The project evaluator takes one label per node.  For a node in more than
    one ego circle, choose the numerically first circle in archive order; a
    node in no released circle receives its own singleton label.  The original
    overlapping circle files remain cached, so this lossy evaluation projection
    is explicit and reproducible.
    """
    edge_path = CACHE / "facebook_combined.txt.gz"
    tar_path = CACHE / "facebook.tar.gz"
    if not edge_path.exists():
        urllib.request.urlretrieve("https://snap.stanford.edu/data/facebook_combined.txt.gz", edge_path)
    if not tar_path.exists():
        urllib.request.urlretrieve("https://snap.stanford.edu/data/facebook.tar.gz", tar_path)
    edges = []
    nodes = set()
    with gzip.open(edge_path, "rt") as src:
        for line in src:
            u, v = map(int, line.split())
            edges.append((u, v)); nodes.update((u, v))
    circle_members = {}
    with tarfile.open(tar_path, "r:gz") as archive:
        for member in sorted((x for x in archive.getmembers() if x.name.endswith(".circles")),
                             key=lambda x: x.name):
            for line in archive.extractfile(member).read().decode().splitlines():
                parts = line.split()
                for node in map(int, parts[1:]):
                    circle_members.setdefault(node, []).append(member.name + ":" + parts[0])
    ids = sorted(nodes)
    remap = {node: i for i, node in enumerate(ids)}
    class_id = {}
    labels = []
    next_id = 0
    for node in ids:
        key = min(circle_members[node]) if node in circle_members else f"singleton:{node}"
        if key not in class_id:
            class_id[key] = next_id; next_id += 1
        labels.append(class_id[key])
    write_graph("ego_facebook", len(ids), ((remap[u], remap[v]) for u, v in edges), labels)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--lfr-only", action="store_true")
    args = ap.parse_args()
    prepare_lfr()
    if not args.lfr_only:
        prepare_planetoid("cora", "cora_v2.zip", "cora_v2")
        prepare_planetoid("citeseer", "citeseer.zip", "citeseer")
        prepare_planetoid("pubmed", "pubmed.zip", "pubmed")
        prepare_ego_facebook()


if __name__ == "__main__":
    main()
