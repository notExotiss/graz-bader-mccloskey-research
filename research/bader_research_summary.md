# Research with Prof. David Bader (NJIT) — Consolidated Summary

> Single-file summary of everything on record about Aarit Malhotra's network-science
> research with Prof. David A. Bader, centered on the **GRAZ** community-detection paper.
> Compiled from `chat_log_full.md`, `GRAZ_cheatsheet.md`, and `graz_extracted.txt`.

---

## 1. The people & the setup

- **Student:** Aarit Malhotra — rising junior (Class of 2028), John P. Stevens High School, Edison, NJ. 4.0 GPA, CS / web-dev focus.
- **Mentor:** **Prof. David A. Bader**, Department of Data Science, **New Jersey Institute of Technology (NJIT)**. Email: `bader@njit.edu`.
- **Project:** Network science / community-detection algorithms, built around Bader's paper introducing **GRAZ**.
- **Aarit's role:** Run the **Section 6 experiments** that empirically confirm the paper's four theorems. The theory is done; the student reproduces the known failures of existing tools and then shows GRAZ removes them.
- **Rapport note (from call prep):** Bader opened with a math-competition story (AHSME); good to bond briefly on math competitions, then pivot to the work. He noted the paper is AI-drafted and not yet fully proofread, so demonstrating independent verification reads as mature.

---

## 2. The paper in one sentence

Louvain/Leiden accept a node move "if modularity goes up," which is secretly a statistical test with the significance bar set to **zero**. **GRAZ adds a real z-test** so a move must *also* be **statistically significant** (z > τ). That roughly 4-line change fixes three long-standing flaws of modularity.

**Full title:** *GRAZ: A Statistically Calibrated Resolution-Limit-Free Replacement for Louvain and Leiden* (David A. Bader, NJIT).

**GRAZ = Greedy Residual-Adjusted Z-test.**

---

## 3. Background field concepts (plain English)

- **Graph** = dots (nodes) + lines (edges). **Community detection** = automatically finding the dense clusters.
- **Modularity (Q)** = the standard score: "are there more edges inside a group than you'd expect by random chance?"
- **Configuration model** = the "random chance" baseline (the *null*). Key phrase.
- **Louvain (2008)** = fast, popular algorithm. **Leiden (2019)** = the improved default used in most pipelines today. Leiden fixed only Louvain's *disconnected-community* defect and left two deeper flaws of modularity untouched.

---

## 4. The three problems GRAZ targets

1. **Resolution limit** (Fortunato–Barthélemy, 2007): modularity is blind to small communities and merges them. Canonical test = **ring of cliques**; past ~22 cliques of size 5, standard tools report about half the true count.
2. **Spurious communities** (Guimerá, 2004): on a *purely random* graph, modularity still hallucinates ~√(m/2) fake communities (m = number of edges).
3. **Instability**: Leiden gives different answers when node processing order is shuffled.

---

## 5. The fix — how GRAZ works

- **New accept rule:** accept a move only if **modularity rises (ΔQ > 0) AND z-score > τ (tau)**.
- Every candidate move is gated by the **standardized Pearson residual** of its 2×2 contingency table under the configuration-model null (a closed-form z-statistic derived for the single-node move).
- **z-score** = how many standard deviations a move is from random noise. **τ = the significance dial**, the only knob:
  - τ = 2 → p < 0.05 (default)
  - τ = 3 → p < 0.001 (strict)
- **Selling points:** drop-in replacement, ~4 lines of code, nearly free (about one square root per move).
- **Lineage:** builds on a **2010 Bader–McCloskey** sketch; GRAZ finishes that idea.

---

## 6. The four theorems (claims)

1. **Resolution-limit mitigation** — on a ring of n cliques of size c, GRAZ refuses the spurious merge whenever `n ≤ (c² − c + 2)/g(τ)²`, where `g(τ) = (√(τ²+4) − τ)/2`. Pushes the threshold back ~6× at τ=2, ~11× at τ=3, and arbitrarily far as τ grows.
2. **Uniformly bounded false-positive rate** — per-move false-merge probability is bounded by **1 − Φ(τ)**, independent of graph size or community-size heterogeneity. (Standard modularity: ~50%, with no bound.)
3. **Convergence** — under the conjunction rule (z > τ ∧ ΔQ > 0), modified Louvain/Leiden still converges in finite time, with Q as a Lyapunov function; lands on a cleaner, narrower answer set → less instability.
4. **Regularized / locally adaptive modularity** — GRAZ is the greedy optimizer of a heteroscedastic penalty on modularity where the local regularization scale is set by the null variance at each move. This **per-move** adaptivity is the real conceptual edge over the **global** resolution parameters of Reichardt–Bornholdt (γ) and CPM.

---

## 7. Aarit's job — Section 6 experiments (the work plan)

The pattern: **first reproduce the failures with standard tools, then show GRAZ removes them.**

- **Week 1 — Ring of cliques (Exp 1 / Table 1):** run Louvain + Leiden, n = 20 → 400, count communities, show the collapse past n ≈ 22. Deliverable: matplotlib plot with error bars.
- **Week 2 — Erdős–Rényi random graphs (Exp 2 / Table 2):** show standard modularity hallucinates ~√(m/2) communities. Deliverable: log-log plot.
- **Weeks 3–4 — Implement GRAZ (~4 lines) together:** show it fixes both failures; then **stability (Exp 5):** 100 shuffled runs, demonstrate GRAZ varies less than Leiden.

**Environment:** Python 3.10+, `pip install networkx python-igraph scikit-learn matplotlib numpy`.

---

## 8. Talking points & smart questions (call prep)

**Lines that prove he read it:**
- "The core idea is that Louvain's accept rule is a statistical test with the threshold at zero, and GRAZ calibrates it with a z-test, right?"
- "I like that τ has an actual p-value meaning, instead of γ, which you grid-search."
- "The cleanest result to me was the random-graph one: modularity *invents* communities, and GRAZ provably doesn't."

**Smart questions to ask:**
1. "Table 1 predicts exact integer counts. Should I treat any deviation as a code bug, or is there real run-to-run randomness even on the ring of cliques?"
2. "Is the 4-line GRAZ edit going into networkx's Louvain or igraph's Leiden? Which codebase?"
3. "Should I verify my z-score formula against the Eq. 5 variance before trusting the experiments?"
4. "Is there a risk τ is too conservative and GRAZ misses *real* small communities? How do we tell a good refusal from a bad one?"

**If asked "what would you check first?":** "I'd verify the experiments match the theorems' exact predictions. If Table 1's GRAZ numbers don't line up with Theorem 1, that's exactly what the experiments are meant to catch."

---

## 9. Glossary

- **null / configuration model** — random baseline for comparison
- **ER graph** — Erdős–Rényi, pure-random graph with no built-in structure
- **NMI / F1** — metrics for "how close to the true answer"
- **LFR** — synthetic benchmark graphs with known ground-truth communities
- **SBM** — stochastic block model (a rival approach)
- **CPM** — Constant Potts Model; a rival fix that drops modularity's meaning
- **γ (gamma)** — the global resolution knob used by rival methods
- **τ (tau)** — GRAZ's significance knob (the p-value threshold)
- **Q** — modularity score
- **Φ** — standard normal CDF (so 1 − Φ(τ) is the upper-tail probability)

---

## 10. Source files (in this folder)

- `graz_extracted.txt` — full UTF-8 text extraction of the GRAZ paper (`graz.pdf`).
- `GRAZ_cheatsheet.md` — one-page call cheat sheet (condensed).
- `chat_log_full.md` — full chronological session log (context for how this work arose).

**Extraction gotcha on record:** the PDF read tooling initially failed (pdftoppm error; no pypdf; cp1252 `UnicodeEncodeError` on Greek τ). Fixed via `python -m pip install pypdf` and writing the extraction to a UTF-8 file (`graz_extracted.txt`) to avoid Windows cp1252 charmap errors on non-ASCII characters.
