# GRAZ Call Cheat-Sheet — David Bader

## ONE-SENTENCE SUMMARY
Louvain/Leiden accept a move "if modularity goes up," which is secretly a statistical test with the
significance bar set to **zero** — GRAZ adds a real z-test so moves must also be **statistically significant**.
That ~4-line change fixes three old flaws of modularity.

## THE FIELD (30 sec version)
- **Graph** = dots (nodes) + lines (edges). **Community detection** = auto-find the dense clusters.
- **Modularity (Q)** = the score everyone uses: "more edges inside a group than random chance?"
- The "random chance" baseline = **configuration model** (the *null*). Memorize this phrase.
- **Louvain** (2008) = fast popular algorithm. **Leiden** (2019) = improved default everyone uses now.

## THE 3 PROBLEMS GRAZ TARGETS
1. **Resolution limit** (Fortunato–Barthélemy 2007): modularity is blind to small communities and
   merges them. Test case = **ring of cliques**: past ~22 cliques (size 5), it reports half as many.
2. **Spurious communities** (Guimerá 2004): on a *purely random* graph, modularity still hallucinates
   ~√(m/2) fake communities.
3. **Instability**: Leiden gives different answers when you shuffle node order.

## THE FIX
- **GRAZ = Greedy Residual-Adjusted Z-test.**
- New rule: accept a move only if **modularity goes up AND z-score > τ (tau)**.
- **z-score** = "how many std-devs from random noise." **τ = significance dial** (the ONLY knob).
  - τ = 2 → p < 0.05 (default).  τ = 3 → p < 0.001 (strict).
- Selling points: **drop-in, ~4 lines of code, ~free (one square root per move).**
- Roots in a **2010 Bader–McCloskey** sketch (the NSA collaborator) — GRAZ finishes that idea.

## THE 4 THEOREMS (claims only)
1. Pushes resolution limit back ~6× at τ=2, ~11× at τ=3, arbitrarily far as τ grows.
2. False-positive rate capped at 1−Φ(τ), any graph size. (Standard modularity: ~50%, no bound.)
3. Still converges; lands on a cleaner/narrower answer set → less instability.
4. **Locally adaptive resolution**: penalty is *per-move*, not one global knob (γ/CPM are global).
   This is the real conceptual edge on mixed big+small communities.

## MY JOB (Section 6 — know this cold)
The theory is done; **I run the experiments that confirm it.**
- **Week 1** — Ring of cliques (Exp 1 / Table 1): run Louvain+Leiden, n=20→400, count communities,
  show collapse past n=22. Deliver matplotlib plot w/ error bars by Friday.
- **Week 2** — Erdős–Rényi random graphs (Exp 2 / Table 2): show they hallucinate ~√(m/2). Log-log plot.
- **Week 3–4** — implement GRAZ together (the ~4 lines), show it fixes both; then stability (Exp 5:
  100 shuffled runs, GRAZ varies less).
- Logic: first reproduce the *failures* with standard tools, then show GRAZ removes them.
- Setup: Python 3.10+, pip install `networkx python-igraph scikit-learn matplotlib numpy`.

## LINES TO SAY (prove I read it)
- "The core idea is Louvain's accept rule is a stat test with the threshold at zero, and GRAZ calibrates
  it with a z-test — right?"
- "I like that τ has an actual p-value meaning instead of γ, which you grid-search."
- "The cleanest result to me was the random-graph one — modularity *invents* communities, GRAZ provably doesn't."
- Rapport: bond briefly on math competitions (he led with the AHSME story), then pivot to the work.

## SMART QUESTIONS TO ASK
1. "Table 1 predicts exact integer counts — treat any deviation as a code bug, or is there real
   run-to-run randomness even on the ring of cliques?"
2. "Is the 4-line GRAZ edit going into networkx's Louvain or igraph's Leiden — which codebase?"
3. "Should I verify my z-score formula against the Eq. 5 variance before trusting the experiments?"
4. "Risk that τ is too conservative and GRAZ misses *real* small communities — how do we tell a good
   refusal from a bad one?"

## IF HE ASKS "WHAT WOULD YOU CHECK FIRST?"
"I'd verify the experiments match the theorems' exact predictions — and if Table 1's GRAZ numbers don't
line up with Theorem 1, that's exactly what the experiments are meant to catch." (He noted the paper is
Claude-drafted and unproofread, so this shows maturity.)

## GLOSSARY (don't get caught)
null / configuration model = random baseline · ER graph = pure-random, no structure ·
NMI / F1 = "how close to the true answer" · LFR = synthetic test graphs w/ known communities ·
SBM = stochastic block model (rival approach) · CPM = rival fix that drops modularity's meaning ·
γ (gamma) = the global knob rivals use · τ (tau) = GRAZ's significance knob.
