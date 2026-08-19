"""Prepare FULL-SIZE SNAP graphs with ground-truth communities for the C++ benchmark.

Unlike graz/datasets.py::_load_snap_community -- which caps the graph at 15k
nodes / 1000 communities and takes an induced subgraph for the laptop-scale
Python experiments -- this writes the ENTIRE graph plus all 5000 ground-truth
communities, so the reference partition has communities in the thousands.

Writes  cpp/prepared/<name>.edges   "u v"          0-based contiguous ids
        cpp/prepared/<name>.labels  "node label"   only nodes covered by a community

SNAP com-* ground truth is OVERLAPPING. For a single-label (partition)
evaluation each covered node is assigned the SMALLEST community containing it,
which is the most specific/cohesive choice. The multi-membership fraction is
reported so the resulting NMI/F1 can be read as a lower bound: no partition can
reproduce an overlapping cover exactly.
"""
import gzip
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "cpp", "prepared")

# name -> (edge gz, community gz)
SETS = {
    "amazon_full": ("com-amazon.ungraph.txt.gz", "com-amazon.top5000.cmty.txt.gz"),
    "dblp_full": ("com-dblp.ungraph.txt.gz", "com-dblp.top5000.cmty.txt.gz"),
    "youtube": ("com-youtube.ungraph.txt.gz", "com-youtube.top5000.cmty.txt.gz"),
    "livejournal": ("com-lj.ungraph.txt.gz", "com-lj.top5000.cmty.txt.gz"),
}


def prep(name, edge_gz, cmty_gz):
    epath = os.path.join(DATA, edge_gz)
    cpath = os.path.join(DATA, cmty_gz)
    if not os.path.exists(epath) or not os.path.exists(cpath):
        print(f"[{name}] SKIP - missing source files")
        return

    # ---- communities first ------------------------------------------------
    comms = []
    with gzip.open(cpath, "rt") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            comms.append([int(x) for x in line.split()])

    member_count = defaultdict(int)
    for c in comms:
        for v in c:
            member_count[v] += 1

    # largest-first so smaller communities overwrite -> smallest wins
    label_of = {}
    for ci in sorted(range(len(comms)), key=lambda i: -len(comms[i])):
        for v in comms[ci]:
            label_of[v] = ci

    covered = len(label_of)
    overlapped = sum(1 for v in label_of if member_count[v] > 1)

    # ---- edges ------------------------------------------------------------
    remap = {}

    def rid(x):
        r = remap.get(x)
        if r is None:
            r = len(remap)
            remap[x] = r
        return r

    os.makedirs(OUT, exist_ok=True)
    n_edges = 0
    with gzip.open(epath, "rt") as f, open(
        os.path.join(OUT, name + ".edges"), "w", buffering=1 << 20
    ) as out:
        for line in f:
            if not line or line[0] == "#":
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            ia = rid(int(parts[0]))
            ib = rid(int(parts[1]))
            out.write(f"{ia} {ib}\n")
            n_edges += 1

    n_nodes = len(remap)

    # ---- labels for covered nodes present in the graph --------------------
    written = 0
    used_labels = set()
    with open(os.path.join(OUT, name + ".labels"), "w", buffering=1 << 20) as out:
        for raw, lab in label_of.items():
            nid = remap.get(raw)
            if nid is not None:
                out.write(f"{nid} {lab}\n")
                used_labels.add(lab)
                written += 1

    sizes = sorted((len(c) for c in comms), reverse=True)
    print(
        f"[{name}] nodes={n_nodes:,} edges={n_edges:,} "
        f"gt_communities={len(comms):,} (distinct after smallest-wins: {len(used_labels):,}) "
        f"labeled_nodes={written:,} of {covered:,} covered; "
        f"multi-membership={overlapped:,} ({100.0*overlapped/max(covered,1):.1f}%); "
        f"comm_size max={sizes[0]} median={sizes[len(sizes)//2]} min={sizes[-1]}"
    )
    sys.stdout.flush()


if __name__ == "__main__":
    want = sys.argv[1:] or list(SETS)
    for name in want:
        if name not in SETS:
            print(f"unknown set {name}")
            continue
        prep(name, *SETS[name])
