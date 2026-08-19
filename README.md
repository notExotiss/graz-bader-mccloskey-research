# GRAZ: Bader–McCloskey correction for community detection

This repository contains the complete research implementation, experiments, raw outputs, figures, datasets used by the experiments, and the relevant development history for the GRAZ project.

## Contents

- `cpp/` — standalone C++17 implementation and prepared edge-list inputs.
- `graz/` — Python implementation of Louvain, Leiden, CPM, and GRAZ.
- `experiments/` — Python experiment runners and theory checks.
- `results/` — CSV tables, raw console output, figures, and rolling-stop studies.
- `data/` — compressed SNAP source data used to prepare the test inputs.
- `research/` — source notes and extracted research context.
- `history/` — readable transcripts and the consolidated historical research record.
- `MANIFEST.md` — file-level provenance and packaging decisions.

## Build the C++ program

From `cpp/`:

```powershell
g++ -O2 -std=c++17 -static -static-libgcc -static-libstdc++ graz.cpp -o graz.exe
```

Run the built-in experiments:

```powershell
./graz.exe
```

Run rolling-window experiments on prepared graphs:

```powershell
./graz.exe --only-roll --roll-small
./graz.exe --only-roll --roll-er
./graz.exe --only-roll large=amazon_full=prepared/amazon_full.edges:prepared/amazon_full.labels --max-merges=333311
```

The large LiveJournal run was intentionally capped at 200,000 CNM merges because of its measured runtime. The stopping-rule cut occurs before that budget, so the reported rolling-window prefix is reproducible without completing the full agglomeration.

## Python environment

```powershell
pip install networkx python-igraph scikit-learn matplotlib numpy scipy
python experiments/run_all.py
```

## Main findings

The raw results are authoritative; the short summary is in `results/FINDINGS.md`. The work reproduces the ring-of-cliques resolution-limit result, tests GRAZ against Erdős–Rényi nulls, evaluates LFR and hierarchical SBM benchmarks, measures stability, and compares rolling ΔQ stopping rules across labelled and large SNAP networks.

## Provenance

The research began from the Bader–McCloskey material and the GRAZ paper context recorded in `research/`. The readable Claude Code transcripts in `history/` preserve the implementation and verification process. They are historical records, not executable instructions.

## License

No license was specified in the source project. Add a license before treating this repository as an externally redistributed open-source package.
