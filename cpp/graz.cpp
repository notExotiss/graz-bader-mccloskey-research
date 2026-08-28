// GRAZ: Greedy Residual-Adjusted Z-test for community detection.
//
// C++ implementation of the Bader-McCloskey statistical correction applied to
// modularity-based community detection (Louvain / Leiden).
//
// One multilevel optimizer with four pluggable move rules:
//   Louvain     : accept move iff dQ > 0                 (Blondel 2008)
//   Leiden      : Louvain + connectivity refinement       (Traag 2019)
//   CPM-Leiden  : Constant Potts Model objective          (Traag 2011)
//   GRAZ        : accept move iff z > tau AND dQ > 0       (Bader, GRAZ)
//
// The GRAZ acceptance rule is the only structural change from Louvain. Every
// candidate move of node i into community B is gated by the standardized
// Pearson residual of the move under the configuration-model null (Eq. 6):
//
//     z_iB = (2m*k_iB - k_i*Sigma_B) / sqrt(k_i*Sigma_B*(2m - Sigma_B))
//
// MEASURED CAVEAT (see ../results/FINDINGS.md): z_iB is strictly decreasing in
// Sigma_B, so tau imposes a hard ceiling on community VOLUME -- it behaves as a
// resolution parameter, not a significance level, and Corollary 1's
// non-imposition bound fails by ~300x at tau=3. GRAZ_MERGE (the literal 2010
// merge-level test) is the better default: same accuracy on large graphs, ~1.8x
// faster, ~4.5x fewer z-tests, and far less over-segmentation on small graphs.
//
// Build:  g++ -O2 -std=c++17 -static graz.cpp -o graz
//         (-static bundles libstdc++/libgcc so the exe runs without the
//          MinGW bin directory on PATH)
// Run:    ./graz                 (synthetic experiments 1,2,4,5 + verification)
//         ./graz name=edges.txt[:labels.txt] ...   (adds experiment 6)
//         ./graz large=name=edges[:labels] ...     (adds experiment 7: large
//                                                   graphs, quality + timing)
//         ./graz --only-large large=...            (skip the synthetic suite)
//         ./graz --methods=graz4 ...               (exp 7: only these methods;
//                                                   keys: louvain leiden graz2
//                                                   graz3 graz4)
//         ./graz --only-roll --roll-er --roll-small name=... large=...
//                                                  (exps 9/10/11: the rolling-
//                                                   window dQ stopping rule)
//         ./graz --full-stop-report name=... large=...
//                                                  (exp 13: T1--T5, BM-native,
//                                                   CM baseline, and ER-calibrated h)
//         ./graz --multicut-report name=... large=...
//                                                  (exp 14: held-out/MDL/stability cuts)
//
// A SECOND statistical criterion is implemented alongside Eq. 6: a rolling-window
// test on the sequence of accepted dQ values (CNM-style). Keep the last W dQ
// scores, and stop when the current dQ lies more than k standard deviations from
// their mean. Unlike Eq. 6 this makes no null-model assumption and never changes
// which move is chosen -- it only decides when to stop -- so it cannot degenerate
// into a resolution parameter. See struct RollStop and experiments 9/10/11.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <deque>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using std::size_t;
using std::vector;

// --------------------------------------------------------------------------- //
//  Timing instrumentation for the core algorithm phases.
//
//  A multilevel run decomposes into exactly three kinds of work:
//    local_move  : the move loop (this is where the GRAZ z-gate is evaluated)
//    aggregate   : collapsing communities into super-nodes between levels
//    refine      : the Leiden connectivity split at the end
//  Counters are process-global; call g_t.reset() before each measured run and
//  snapshot g_t immediately after, before computing metrics.
// --------------------------------------------------------------------------- //
struct Timings {
    double local_move_ms = 0.0;
    double aggregate_ms = 0.0;
    double refine_ms = 0.0;
    double total_ms = 0.0;
    long long move_sweeps = 0;   // local_moving passes over all nodes
    long long node_visits = 0;   // individual node move evaluations
    long long z_tests = 0;       // z_statistic calls in the GRAZ gate
    int levels = 0;              // multilevel coarsening levels

    void reset() { *this = Timings(); }
};
static Timings g_t;

// Method keys selected by --methods=a,b,c for experiment 7. Empty = run all.
static std::set<std::string> g_methods;

// --------------------------------------------------------------------------- //
//  Rolling-window stopping rule on the sequence of accepted dQ values.
//
//  This is a DIFFERENT kind of statistical test from GRAZ's Eq. 6. Eq. 6 scores
//  one candidate move against an analytic configuration-model null, and because
//  z_iB is strictly decreasing in Sigma_B it degenerates into a ceiling on
//  community VOLUME (see ../results/FINDINGS.md). The rule here instead scores
//  the current dQ against the EMPIRICAL distribution of recent accepted dQ:
//
//      keep the last W scores, take their mean mu and population sd s;
//      stop when  mu - dQ > k*s   (one-sided)   or   |dQ - mu| > k*s  (two-sided)
//
//  No null model, no closed-form assumption, and no dependence on community
//  volume -- so it cannot silently become a resolution parameter. W ~ 10-20 and
//  k ~ 1.75-2.5 are the ranges worth tuning.
//
//  One-sided vs two-sided is not cosmetic. A greedy dQ trace trends downward, so
//  the informative event is a dQ falling FAR BELOW recent history: merges have
//  stopped paying and the real structure is exhausted. An upward outlier is an
//  unusually GOOD merge and halting on it would be backwards -- but "more than
//  2 s.d. away from the mean" reads as two-sided, so both are measured.
// --------------------------------------------------------------------------- //
struct RollStop {
    int W = 15;              // window length
    double kmult = 2.0;      // how many sd's count as "far"
    bool two_sided = false;  // literal reading of "away from the mean"

    std::deque<double> win;
    long long seen = 0;      // dQ values offered to the test
    long long fired_at = -1; // 1-based index of the value that tripped it
    double trip_dq = 0.0, trip_mu = 0.0, trip_sd = 0.0;

    // True when the rule says stop. The tripping value is deliberately NOT added
    // to the window: it is the anomaly under test, not part of its own baseline.
    bool test(double dq) {
        ++seen;
        if ((int)win.size() >= W) {
            double mu = 0.0;
            for (double x : win) mu += x;
            mu /= (double)win.size();
            double ss = 0.0;
            for (double x : win) ss += (x - mu) * (x - mu);
            double sd = std::sqrt(ss / (double)win.size());
            if (sd > 0.0) {
                double dev = two_sided ? std::fabs(dq - mu) : (mu - dq);
                if (dev > kmult * sd) {
                    fired_at = seen; trip_dq = dq; trip_mu = mu; trip_sd = sd;
                    return true;
                }
            }
        }
        win.push_back(dq);
        if ((int)win.size() > W) win.pop_front();
        return false;
    }
};

// Active rolling stop for the multilevel port (Rule::ROLLQ). Null = inactive.
static RollStop *g_roll = nullptr;
static bool g_roll_fired = false;

using Clock = std::chrono::steady_clock;
inline double ms_since(const Clock::time_point &t0) {
    return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

// --------------------------------------------------------------------------- //
//  Weighted undirected graph with self-loops (used for aggregated levels).
//  size[u] = number of ORIGINAL nodes the super-node u represents (for CPM).
// --------------------------------------------------------------------------- //
struct Graph {
    // adjacency as parallel arrays per node: nbr[u] -> {(v, w), ...}
    vector<vector<std::pair<int, double>>> adj;
    vector<double> degree;  // weighted degree (self-loop counts twice)
    vector<long long> size; // node multiplicity for CPM
    double m2 = 0.0;        // 2m = sum of degrees

    int n() const { return (int)adj.size(); }

    void finalize() {
        int N = n();
        degree.assign(N, 0.0);
        m2 = 0.0;
        for (int u = 0; u < N; ++u) {
            double d = 0.0;
            for (auto &e : adj[u]) d += (e.first == u) ? 2.0 * e.second : e.second;
            degree[u] = d;
            m2 += d;
        }
        if (size.empty()) size.assign(N, 1);
    }

    static Graph from_edges(int n_nodes, const vector<std::pair<int, int>> &edges) {
        Graph g;
        g.adj.assign(n_nodes, {});
        // accumulate weights (dedup parallel edges) via per-node maps
        vector<std::unordered_map<int, double>> tmp(n_nodes);
        for (auto &e : edges) {
            int u = e.first, v = e.second;
            if (u == v) { tmp[u][u] += 1.0; }
            else { tmp[u][v] += 1.0; tmp[v][u] += 1.0; }
        }
        for (int u = 0; u < n_nodes; ++u)
            for (auto &kv : tmp[u]) g.adj[u].push_back({kv.first, kv.second});
        g.size.assign(n_nodes, 1);
        g.finalize();
        return g;
    }
};

// Shared graph-plus-label container used by the synthetic generators, prepared
// SNAP loaders, and the statistically independent edge split below.
struct Dataset { int n; vector<std::pair<int, int>> edges; vector<int> gt; };

// --------------------------------------------------------------------------- //
//  Standardized Pearson residual of the single-node move (Eq. 6).
// --------------------------------------------------------------------------- //
inline double z_statistic(double k_i, double k_iB, double Sigma_B, double m2) {
    double denom = k_i * Sigma_B * (m2 - Sigma_B);
    if (denom <= 0.0) return 0.0;
    return (m2 * k_iB - k_i * Sigma_B) / std::sqrt(denom);
}

// --------------------------------------------------------------------------- //
//  Newman-Girvan modularity of a labelling on a raw adjacency graph.
// --------------------------------------------------------------------------- //
double modularity(const Graph &g, const vector<int> &comm) {
    double m2 = g.m2;
    if (m2 == 0) return 0.0;
    std::unordered_map<int, double> Lin, D;
    for (int u = 0; u < g.n(); ++u) {
        int cu = comm[u];
        D[cu] += g.degree[u];
        for (auto &e : g.adj[u]) {
            if (comm[e.first] == cu)
                Lin[cu] += (e.first == u) ? 2.0 * e.second : e.second;
        }
    }
    double q = 0.0;
    for (auto &kv : D) {
        double din = Lin.count(kv.first) ? Lin[kv.first] : 0.0;
        q += din / m2 - (kv.second / m2) * (kv.second / m2);
    }
    return q;
}

//  GRAZ_MERGE is the literal 2010 Bader-McCloskey formulation: the z-gate is a
//  test on MERGING two communities, so level 0 (where every node is a singleton
//  and every candidate target is a neighbour by construction) runs plain
//  modularity, and the gate engages only from the first aggregated level up.
//  This removes the selection bias that makes Corollary 1 fail.
//  GRAZ_GATE reads Eq. 6 as a pure GATE rather than an objective: among the
//  candidates that clear tau, it takes the one maximising dQ, not the one
//  maximising z. This matters because z is inflated when Sigma_B is SMALL, so
//  argmax-z steers every node into the smallest admissible community -- which
//  is why raising tau makes fragmentation worse instead of better.
//  ROLLQ is the rolling-window stopping rule ported into the Louvain-family move
//  loop: moves are chosen by plain modularity, and the rule only decides WHEN to
//  stop, never which move to take. That separation is the point -- GRAZ's Eq. 6
//  changes the move choice and so distorts what is found, while a stopping rule
//  can only truncate.
enum class Rule { Modularity, GRAZ, CPM, GRAZ_MERGE, GRAZ_GATE, ROLLQ };

// --------------------------------------------------------------------------- //
//  One level of local moving to convergence on graph g.
// --------------------------------------------------------------------------- //
//  adaptive_tau applies the paper's own Section 3.4 scaling: the threshold
//  becomes tau * sqrt(2 log m) rather than a fixed constant, so it grows with
//  the number of candidate tests instead of staying pinned at tau = 2 or 3.
vector<int> local_moving(const Graph &g, std::mt19937 &rng, Rule rule,
                         double tau, double gamma, bool adaptive_tau = false) {
    auto t_start = Clock::now();
    int N = g.n();
    double m2 = g.m2, m = m2 / 2.0;
    if (adaptive_tau) tau *= std::sqrt(2.0 * std::log(std::max(m, 3.0)));
    vector<int> comm(N);
    std::iota(comm.begin(), comm.end(), 0);
    vector<double> Sigma = g.degree;
    vector<long long> nnodes = g.size;

    vector<int> order(N);
    std::iota(order.begin(), order.end(), 0);

    // reusable scratch for weights into neighbouring communities
    std::unordered_map<int, double> w_to;

    bool improved_any = true;
    while (improved_any) {
        improved_any = false;
        ++g_t.move_sweeps;
        std::shuffle(order.begin(), order.end(), rng);
        for (int u : order) {
            ++g_t.node_visits;
            double ku = g.degree[u];
            long long su = g.size[u];
            int A = comm[u];

            w_to.clear();
            for (auto &e : g.adj[u]) {
                if (e.first == u) continue;
                w_to[comm[e.first]] += e.second;
            }

            // remove u from A
            Sigma[A] -= ku;
            nnodes[A] -= su;
            double w_to_A = w_to.count(A) ? w_to[A] : 0.0;

            int target = A;
            if (rule == Rule::CPM) {
                double best_gain = w_to_A - gamma * (double)su * (double)nnodes[A];
                for (auto &kv : w_to) {
                    double gain = kv.second - gamma * (double)su * (double)nnodes[kv.first];
                    if (gain > best_gain + 1e-15) { best_gain = gain; target = kv.first; }
                }
            } else if (rule == Rule::GRAZ || rule == Rule::GRAZ_MERGE ||
                       rule == Rule::GRAZ_GATE) {
                bool by_gain = (rule == Rule::GRAZ_GATE);
                double gain_A = (m > 0) ? w_to_A / m - ku * Sigma[A] / (2.0 * m * m) : 0.0;
                double best_z = -1e300, best_gain = -1e300; bool found = false;
                for (auto &kv : w_to) {
                    int B = kv.first; double w_iB = kv.second;
                    if (B == A) continue;
                    double SB = Sigma[B];
                    ++g_t.z_tests;
                    double z = z_statistic(ku, w_iB, SB, m2);
                    if (z <= tau) continue;
                    double gain_B = (m > 0) ? w_iB / m - ku * SB / (2.0 * m * m) : 0.0;
                    double dQ = gain_B - gain_A;
                    if (dQ <= 0.0) continue;
                    // GRAZ/GRAZ_MERGE rank by z (the paper's argmax-z);
                    // GRAZ_GATE ranks by modularity gain among gate survivors.
                    if (by_gain) {
                        if (!found || dQ > best_gain) { best_gain = dQ; target = B; found = true; }
                    } else {
                        if (!found || z > best_z) { best_z = z; target = B; found = true; }
                    }
                }
            } else { // Modularity (and ROLLQ, which differs only in when it stops)
                double stay_gain = (m > 0) ? w_to_A / m - ku * Sigma[A] / (2.0 * m * m) : 0.0;
                double best_gain = stay_gain;
                for (auto &kv : w_to) {
                    double gain = (m > 0) ? kv.second / m - ku * Sigma[kv.first] / (2.0 * m * m) : 0.0;
                    if (gain > best_gain + 1e-15) { best_gain = gain; target = kv.first; }
                }
                // The rolling rule sees the dQ of each ACCEPTED move, in the order
                // accepted, and can halt the sweep mid-pass. Moves already made
                // stand; this is a stopping rule, not a rollback.
                if (rule == Rule::ROLLQ && g_roll && target != A) {
                    if (g_roll->test(best_gain - stay_gain)) {
                        comm[u] = target;
                        Sigma[target] += ku;
                        nnodes[target] += su;
                        g_roll_fired = true;
                        g_t.local_move_ms += ms_since(t_start);
                        return comm;
                    }
                }
            }

            comm[u] = target;
            Sigma[target] += ku;
            nnodes[target] += su;
            if (target != A) improved_any = true;
        }
    }
    g_t.local_move_ms += ms_since(t_start);
    return comm;
}

// --------------------------------------------------------------------------- //
//  Aggregation: collapse communities into super-nodes.
// --------------------------------------------------------------------------- //
Graph aggregate(const Graph &g, const vector<int> &comm, vector<int> &remap_out) {
    auto t_start = Clock::now();
    // remap community ids to 0..K-1
    std::map<int, int> remap;
    for (int u = 0; u < g.n(); ++u)
        if (!remap.count(comm[u])) { int id = (int)remap.size(); remap[comm[u]] = id; }
    int K = (int)remap.size();

    remap_out.assign(g.n(), 0);
    Graph ng;
    ng.adj.assign(K, {});
    ng.size.assign(K, 0);
    vector<std::unordered_map<int, double>> tmp(K);
    for (int u = 0; u < g.n(); ++u) ng.size[remap[comm[u]]] += g.size[u];

    for (int u = 0; u < g.n(); ++u) {
        int cu = remap[comm[u]];
        for (auto &e : g.adj[u]) {
            int cv = remap[comm[e.first]];
            if (e.first == u) tmp[cu][cu] += e.second;          // existing self-loop
            else if (cu == cv) tmp[cu][cu] += e.second / 2.0;    // internal edge -> half self-loop
            else tmp[cu][cv] += e.second;
        }
    }
    for (int c = 0; c < K; ++c)
        for (auto &kv : tmp[c]) ng.adj[c].push_back({kv.first, kv.second});
    ng.finalize();

    for (int u = 0; u < g.n(); ++u) remap_out[u] = remap[comm[u]];
    g_t.aggregate_ms += ms_since(t_start);
    return ng;
}

// --------------------------------------------------------------------------- //
//  Split internally-disconnected communities (Leiden's guarantee, simplified).
// --------------------------------------------------------------------------- //
vector<int> split_disconnected(const Graph &base, const vector<int> &labels) {
    auto t_start = Clock::now();
    int N = base.n();
    vector<int> out(N, -1);
    int next_id = 0;
    vector<char> seen(N, 0);
    for (int s = 0; s < N; ++s) {
        if (seen[s]) continue;
        int c = labels[s];
        std::queue<int> q; q.push(s); seen[s] = 1;
        while (!q.empty()) {
            int x = q.front(); q.pop();
            out[x] = next_id;
            for (auto &e : base.adj[x]) {
                int y = e.first;
                if (!seen[y] && labels[y] == c) { seen[y] = 1; q.push(y); }
            }
        }
        ++next_id;
    }
    g_t.refine_ms += ms_since(t_start);
    return out;
}

// --------------------------------------------------------------------------- //
//  Multilevel driver.
// --------------------------------------------------------------------------- //
vector<int> run_multilevel(const Graph &base, Rule rule, double tau, double gamma,
                           unsigned seed, bool refine, bool adaptive_tau = false) {
    auto t_start = Clock::now();
    std::mt19937 rng(seed);
    Graph g = base;
    vector<int> node_to_cur(base.n());
    std::iota(node_to_cur.begin(), node_to_cur.end(), 0);

    for (int level = 0; level < 100; ++level) {
        // GRAZ_MERGE defers the z-gate past the singleton level, where every
        // candidate target is a neighbour and the test is selection-biased.
        Rule lvl_rule = (rule == Rule::GRAZ_MERGE && level == 0) ? Rule::Modularity : rule;
        vector<int> comm = local_moving(g, rng, lvl_rule, tau, gamma, adaptive_tau);
        int ncomm = (int)std::set<int>(comm.begin(), comm.end()).size();
        if (ncomm == g.n()) break;  // no coarsening -> converged
        ++g_t.levels;
        vector<int> remap;
        Graph ng = aggregate(g, comm, remap);
        for (int u = 0; u < base.n(); ++u) node_to_cur[u] = remap[node_to_cur[u]];
        g = ng;
        // The rolling rule tripped inside the sweep above. Keep the partition
        // reached so far (including this level's aggregation) and stop climbing.
        if (rule == Rule::ROLLQ && g_roll_fired) break;
    }

    if (refine) node_to_cur = split_disconnected(base, node_to_cur);

    // normalise labels to 0..K-1
    std::map<int, int> remap;
    for (int u = 0; u < base.n(); ++u)
        if (!remap.count(node_to_cur[u])) { int id = (int)remap.size(); remap[node_to_cur[u]] = id; }
    vector<int> out(base.n());
    for (int u = 0; u < base.n(); ++u) out[u] = remap[node_to_cur[u]];
    g_t.total_ms += ms_since(t_start);
    return out;
}

// convenience wrappers -------------------------------------------------------
vector<int> louvain(const Graph &g, unsigned seed = 0) {
    return run_multilevel(g, Rule::Modularity, 0.0, 1.0, seed, false);
}
vector<int> leiden(const Graph &g, unsigned seed = 0) {
    return run_multilevel(g, Rule::Modularity, 0.0, 1.0, seed, true);
}
vector<int> cpm_leiden(const Graph &g, double gamma = 0.1, unsigned seed = 0) {
    return run_multilevel(g, Rule::CPM, 0.0, gamma, seed, true);
}
vector<int> graz(const Graph &g, double tau = 2.0, unsigned seed = 0) {
    return run_multilevel(g, Rule::GRAZ, tau, 1.0, seed, true);
}
// Fix A: threshold scales as tau*sqrt(2 log m) (paper Section 3.4).
vector<int> graz_adaptive(const Graph &g, double tau = 1.0, unsigned seed = 0) {
    return run_multilevel(g, Rule::GRAZ, tau, 1.0, seed, true, true);
}
// Fix B: z-gate applied at the merge/aggregate level only (Bader-McCloskey 2010).
vector<int> graz_merge(const Graph &g, double tau = 2.0, unsigned seed = 0) {
    return run_multilevel(g, Rule::GRAZ_MERGE, tau, 1.0, seed, true);
}
// Fix C: Eq. 6 as a pure gate -- select by dQ among survivors, not by argmax z.
vector<int> graz_gate(const Graph &g, double tau = 2.0, unsigned seed = 0) {
    return run_multilevel(g, Rule::GRAZ_GATE, tau, 1.0, seed, true);
}
// Rolling-window dQ stopping rule on the Leiden move loop. `out` receives the
// firing diagnostics so a caller can report whether the rule ever engaged.
vector<int> rollq(const Graph &g, int W, double kmult, bool two_sided,
                  RollStop *out = nullptr, unsigned seed = 0) {
    RollStop rs; rs.W = W; rs.kmult = kmult; rs.two_sided = two_sided;
    g_roll = &rs; g_roll_fired = false;
    vector<int> lab = run_multilevel(g, Rule::ROLLQ, 0.0, 1.0, seed, true);
    g_roll = nullptr; g_roll_fired = false;
    if (out) *out = rs;
    return lab;
}

// --------------------------------------------------------------------------- //
//  CNM greedy modularity agglomeration (Clauset-Newman-Moore 2004) with a full
//  dQ trace.
//
//  Start from singletons; repeatedly merge the connected community pair with the
//  largest
//      dQ_ij = 2 * ( w_ij/2m  -  (Sigma_i/2m)(Sigma_j/2m) )
//  and stop at the natural greedy terminus, the first step where no positive-dQ
//  merge remains. Q is monotonically increasing over that run, so the natural
//  terminus IS the argmax-Q cut -- CNM-full and "CNM at best Q" are the same
//  partition and only need to be computed once.
//
//  Why record the trace instead of testing inline: the stopping rule is a pure
//  function of the dQ sequence, and the merge list replays through a union-find,
//  so ONE agglomeration answers every (window, k, one/two-sided) combination.
//  The tuning sweep the user asked for therefore costs one CNM run, not one run
//  per setting -- which is what makes the sweep affordable at this graph size.
//
//  Lazy-deletion heap: a stale entry is recognised by its version stamps rather
//  than being erased, and each merge folds the smaller neighbour map into the
//  larger so total work stays near O(m log n).
// --------------------------------------------------------------------------- //
struct CnmTrace {
    int n = 0;
    double q_start = 0.0;               // Q of the all-singletons partition
    vector<double> dq;                  // dQ of each accepted merge, in order
    vector<std::pair<int, int>> merges; // (survivor, absorbed) community ids
    double build_ms = 0.0;
    bool complete = false;              // false if a guard aborted the run
    bool truncated = false;             // true if stopped early by a merge budget
    // Diagnostics for stopping rules which operate on the *unchanged* CNM
    // merge order.  Entry t describes the candidate accepted at merge t,
    // immediately before that merge is applied.
    vector<double> merge_z;             // Eq. 6, evaluated at community scale
    vector<long long> candidate_pairs;  // adjacent live community pairs
    size_t singleton_phase_end = 0;     // first cut after no mergeable singleton
};

// max_merges > 0 stops the agglomeration after that many merges. This is a
// BUDGET, not a guard: the recorded prefix is byte-identical to the prefix of an
// uncapped run, because CNM is deterministic and the cap only decides when to
// stop appending. Any stopping rule that fires strictly before the cap therefore
// gets exactly the answer it would have got from a full run. What is lost is the
// CNM-full reference row, which is why `truncated` is reported rather than
// quietly folded into the results.
CnmTrace cnm_greedy_trace(const Graph &g, long long max_pushes = 0,
                          long long max_merges = 0) {
    auto t0 = Clock::now();
    CnmTrace tr;
    tr.n = g.n();
    double m2 = g.m2;
    if (m2 <= 0 || g.n() == 0) { tr.complete = true; return tr; }

    int N = g.n();
    vector<std::unordered_map<int, double>> nbr(N);
    for (int u = 0; u < N; ++u)
        for (auto &e : g.adj[u])
            if (e.first != u) nbr[u][e.first] += e.second;

    vector<double> Sigma = g.degree;
    vector<int> csize(N, 1);
    vector<char> alive(N, 1);
    vector<unsigned> ver(N, 0);

    for (int u = 0; u < N; ++u) {
        double a = Sigma[u] / m2;
        tr.q_start -= a * a;
    }

    // 24 bytes rather than 32, via 32-bit version stamps (per-community merge
    // counts, far below 2^32). dq stays DOUBLE on purpose: storing it as float was
    // tried and it reordered near-tied candidates, shifting polblogs NMI by 0.008
    // and email by one merge. Memory is bounded by compacting stale entries
    // instead, which cannot change the result at all.
    struct Ent {
        double dq; int i, j; unsigned vi, vj;
        bool operator<(const Ent &o) const { return dq < o.dq; }  // max-heap
    };
    auto dq_of = [&](int i, int j, double w) {
        return 2.0 * (w / m2 - (Sigma[i] / m2) * (Sigma[j] / m2));
    };

    // Only positive-dQ candidates are ever pushed. The loop stops at the first
    // non-positive maximum, so an entry with dQ <= 0 can never be acted on: it
    // would either be skipped as stale or trigger the terminus. Pruning it at
    // push time is therefore equivalent, and it is what keeps the heap in RAM --
    // as communities grow, most candidate merges go negative, and retaining them
    // exhausted memory on the 334k-node graphs.
    // A vector-backed heap rather than std::priority_queue, so stale entries can
    // be GARBAGE-COLLECTED. This is what actually bounds memory: entries are never
    // erased on merge, only superseded, and on the 334k-node graphs the dead ones
    // grow without limit and exhaust RAM. Compaction drops exactly the entries the
    // pop loop would have skipped, so it changes no result -- only footprint.
    vector<Ent> pq;
    long long pushes = 0, compactions = 0;
    auto hpush = [&](const Ent &e) {
        pq.push_back(e);
        std::push_heap(pq.begin(), pq.end());
        ++pushes;
    };
    auto is_live = [&](const Ent &e) {
        return alive[e.i] && alive[e.j] && ver[e.i] == e.vi && ver[e.j] == e.vj;
    };
    auto compact = [&]() {
        size_t before = pq.size();
        vector<Ent> keep;
        keep.reserve(before / 2 + 1);
        for (const Ent &e : pq) if (is_live(e)) keep.push_back(e);
        pq.swap(keep);
        std::make_heap(pq.begin(), pq.end());
        ++compactions;
        std::cerr << "    [compact] " << before << " -> " << pq.size()
                  << " live entries after " << tr.merges.size() << " merges ("
                  << (long long)ms_since(t0) << " ms)\n";
    };

    for (int u = 0; u < N; ++u)
        for (auto &kv : nbr[u])
            if (u < kv.first) {
                double dq0 = dq_of(u, kv.first, kv.second);
                if (dq0 <= 0.0) continue;
                pq.push_back({dq0, u, kv.first, 0u, 0u});
                ++pushes;
            }
    std::make_heap(pq.begin(), pq.end());

    // `singleton_mergeable` counts live singleton communities which still have
    // a candidate neighbour.  It gives the CNM-axis analogue of GRAZ_MERGE's
    // level-zero deferral without changing one heap comparison or one merge.
    long long candidate_pairs = 0, singleton_mergeable = 0;
    for (int u = 0; u < N; ++u) {
        candidate_pairs += (long long)nbr[u].size();
        if (!nbr[u].empty()) ++singleton_mergeable;
    }
    candidate_pairs /= 2;
    if (singleton_mergeable == 0) tr.singleton_phase_end = 0;

    while (!pq.empty()) {
        Ent e = pq.front();
        std::pop_heap(pq.begin(), pq.end());
        pq.pop_back();
        if (!alive[e.i] || !alive[e.j]) continue;
        if (ver[e.i] != e.vi || ver[e.j] != e.vj) continue;  // stale
        if (e.dq <= 0.0) break;  // natural greedy terminus

        // Recompute the accepted merge's dQ in full double precision from live
        // state, so the recorded trace -- which the stopping rule consumes -- is
        // exact and the float in Ent only ever affects heap ORDER.
        auto itw = nbr[e.i].find(e.j);
        double w_ij = (itw == nbr[e.i].end()) ? 0.0 : itw->second;
        double dq_exact = dq_of(e.i, e.j, w_ij);
        if (dq_exact <= 0.0) break;  // float ordering promoted a non-merge

        // Fold the smaller neighbour map into the larger one.
        int a = e.i, b = e.j;
        if (nbr[a].size() < nbr[b].size()) std::swap(a, b);

        tr.dq.push_back(dq_exact);
        tr.merges.push_back({a, b});
        tr.merge_z.push_back(z_statistic(Sigma[e.i], w_ij, Sigma[e.j], m2));
        tr.candidate_pairs.push_back(candidate_pairs);

        if (csize[a] == 1 && !nbr[a].empty()) --singleton_mergeable;
        if (csize[b] == 1 && !nbr[b].empty()) --singleton_mergeable;
        const long long deg_a = (long long)nbr[a].size();
        const long long deg_b = (long long)nbr[b].size();

        nbr[a].erase(b);
        for (auto &kv : nbr[b]) {
            int x = kv.first;
            if (x == a) continue;
            nbr[a][x] += kv.second;
            nbr[x].erase(b);
            nbr[x][a] = nbr[a][x];
        }
        nbr[b].clear();
        Sigma[a] += Sigma[b];
        Sigma[b] = 0.0;
        csize[a] += csize[b]; csize[b] = 0;
        alive[b] = 0;
        ++ver[a]; ++ver[b];

        for (auto &kv : nbr[a]) {
            int x = kv.first;
            if (!alive[x]) continue;
            double dqx = dq_of(a, x, kv.second);
            if (dqx <= 0.0) continue;  // can never be acted on (see note above)
            hpush({dqx, std::min(a, x), std::max(a, x),
                   ver[std::min(a, x)], ver[std::max(a, x)]});
        }
        // Remove the two old incident neighbourhoods, restore their mutual
        // edge once, then add the union of the new survivor's neighbours.
        candidate_pairs += (long long)nbr[a].size() - deg_a - deg_b + 1;
        if (singleton_mergeable == 0 && tr.singleton_phase_end == 0)
            tr.singleton_phase_end = tr.merges.size();
        // Memory control. What costs RAM is the heap's LIVE size, so bound
        // pq.size(): first try compaction, which is free in result terms, and only
        // abort if the heap is still oversized with all-live entries.
        // Compaction is cheap and result-preserving (on these graphs ~95% of
        // entries are dead), so trigger on a LOW ceiling and abort only if the
        // heap is still oversized when every remaining entry is live. A hub
        // community re-pushes its whole neighbour list on each merge -- ~1000
        // entries by merge 19k on Amazon -- so the dead fraction regrows fast and
        // frequent compaction is much cheaper than a bigger ceiling.
        if (max_pushes > 0 && (long long)pq.size() > max_pushes) {
            compact();
            if ((long long)pq.size() > max_pushes) {
                std::cerr << "    [guard] " << pq.size() << " LIVE entries after "
                          << tr.merges.size() << " merges -- aborting\n";
                tr.build_ms = ms_since(t0);
                tr.complete = false;
                return tr;
            }
        }
        if (max_merges > 0 && (long long)tr.merges.size() >= max_merges) {
            std::cerr << "    [budget] stopping after " << tr.merges.size()
                      << " merges (" << (long long)ms_since(t0) << " ms)\n";
            tr.build_ms = ms_since(t0);
            tr.complete = true;
            tr.truncated = true;
            return tr;
        }
        if (tr.merges.size() % 5000 == 0) {
            size_t nb_entries = 0, nb_max = 0;
            for (int u = 0; u < N; ++u)
                if (alive[u]) {
                    nb_entries += nbr[u].size();
                    if (nbr[u].size() > nb_max) nb_max = nbr[u].size();
                }
            std::cerr << "    ..." << tr.merges.size() << " merges, heap "
                      << pq.size() << ", nbr_entries " << nb_entries
                      << " (max deg " << nb_max << "), "
                      << (long long)ms_since(t0) << " ms\n";
        }
    }
    tr.build_ms = ms_since(t0);
    tr.complete = true;
    return tr;
}

// Replay the first `cut` merges of a trace into a 0..K-1 labelling.
vector<int> cnm_labels_at(const CnmTrace &tr, size_t cut) {
    vector<int> p(tr.n);
    std::iota(p.begin(), p.end(), 0);
    std::function<int(int)> find = [&](int x) {
        while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; }
        return x;
    };
    if (cut > tr.merges.size()) cut = tr.merges.size();
    for (size_t t = 0; t < cut; ++t) {
        int a = find(tr.merges[t].first), b = find(tr.merges[t].second);
        if (a != b) p[b] = a;
    }
    std::map<int, int> remap;
    vector<int> out(tr.n);
    for (int u = 0; u < tr.n; ++u) {
        int r = find(u);
        auto it = remap.find(r);
        if (it == remap.end()) { int id = (int)remap.size(); remap[r] = id; out[u] = id; }
        else out[u] = it->second;
    }
    return out;
}

// --------------------------------------------------------------------------- //
//  Bader--McCloskey merge path with a configuration-model calibration.
//
//  The 2010 proposal is a test on the off-diagonal residual R_ij of the
//  collapsed contingency table, not on a node selected from its neighbours.
//  For two current communities with stub volumes s_i and s_j, total N=2m, and
//  observed cross-edge count l_ij, independence gives
//
//       E[l_ij] = s_i s_j / N
//       Var[l_ij] = s_i s_j (N-s_i) (N-s_j) / (N^2 (N-1)).
//
//  The latter is the fixed-margin 2x2-table variance.  `bm_pair_z` is therefore
//  the symmetric standardized residual S_ij requested by Bader--McCloskey.
//  It is deliberately NOT Eq. 6, whose one-node/neighbour conditioning caused
//  the earlier selection-bias failure.
//
//  A fixed z cutoff is still invalid here because the algorithm chooses the
//  maximum among many pairs.  `bm_path_bootstrap` calibrates that maximum by
//  replaying the ENTIRE greedy Bader path on degree-preserving configuration
//  graphs.  At stage t, p_t is the Monte-Carlo tail probability of the observed
//  maximum S_ij against null paths' stage-t maxima.  This accounts for both the
//  candidate scan and the data-dependent greedy history.
// --------------------------------------------------------------------------- //
inline double bm_pair_z(double si, double sj, double lij, double m2) {
    if (m2 <= 1.0 || si <= 0.0 || sj <= 0.0 || si >= m2 || sj >= m2) return 0.0;
    double expected = si * sj / m2;
    double var = si * sj * (m2 - si) * (m2 - sj) / (m2 * m2 * (m2 - 1.0));
    if (var <= 0.0) return 0.0;
    return (lij - expected) / std::sqrt(var);
}

struct BmTrace {
    int n = 0;
    vector<double> z;                   // selected max standardized residual
    vector<double> dq;                  // its (positive) modularity increment
    vector<std::pair<int, int>> merges; // survivor, absorbed
    double build_ms = 0.0;
};

// Literal Bader--McCloskey greedy rule: merge the positive pair with largest
// standardized residual.  This compact version is intentionally used only for
// bootstrap-sized graphs; the production-scale CNM trace above remains separate.
BmTrace bm_greedy_trace(const Graph &g, size_t max_merges = 0) {
    auto t0 = Clock::now();
    BmTrace tr; tr.n = g.n();
    const int N = g.n();
    const double m2 = g.m2;
    if (N == 0 || m2 <= 1.0) { tr.build_ms = ms_since(t0); return tr; }

    vector<std::unordered_map<int, double>> nbr(N);
    for (int u = 0; u < N; ++u)
        for (auto &e : g.adj[u]) if (e.first != u) nbr[u][e.first] += e.second;
    vector<double> sigma = g.degree;
    vector<char> alive(N, 1);
    vector<unsigned> ver(N, 0);

    struct Ent {
        double z; int i, j; unsigned vi, vj;
        bool operator<(const Ent &o) const { return z < o.z; }
    };
    auto score = [&](int i, int j, double w) { return bm_pair_z(sigma[i], sigma[j], w, m2); };
    vector<Ent> pq;
    for (int u = 0; u < N; ++u)
        for (auto &kv : nbr[u]) if (u < kv.first) {
            double z = score(u, kv.first, kv.second);
            if (z > 0.0) pq.push_back({z, u, kv.first, 0u, 0u});
        }
    std::make_heap(pq.begin(), pq.end());

    while (!pq.empty() && (max_merges == 0 || tr.merges.size() < max_merges)) {
        Ent e = pq.front(); std::pop_heap(pq.begin(), pq.end()); pq.pop_back();
        if (!alive[e.i] || !alive[e.j] || ver[e.i] != e.vi || ver[e.j] != e.vj) continue;
        auto it = nbr[e.i].find(e.j);
        double w = it == nbr[e.i].end() ? 0.0 : it->second;
        double z = score(e.i, e.j, w);
        if (z <= 0.0) break;

        int a = e.i, b = e.j;
        if (nbr[a].size() < nbr[b].size()) std::swap(a, b);
        tr.z.push_back(z);
        tr.dq.push_back(2.0 * (w / m2 - (sigma[a] / m2) * (sigma[b] / m2)));
        tr.merges.push_back({a, b});

        nbr[a].erase(b);
        for (auto &kv : nbr[b]) {
            int x = kv.first;
            if (x == a) continue;
            nbr[a][x] += kv.second;
            nbr[x].erase(b);
            nbr[x][a] = nbr[a][x];
        }
        nbr[b].clear();
        sigma[a] += sigma[b]; sigma[b] = 0.0;
        alive[b] = 0; ++ver[a]; ++ver[b];
        for (auto &kv : nbr[a]) {
            int x = kv.first;
            if (!alive[x]) continue;
            int i = std::min(a, x), j = std::max(a, x);
            double zx = score(i, j, kv.second);
            if (zx > 0.0) { pq.push_back({zx, i, j, ver[i], ver[j]}); std::push_heap(pq.begin(), pq.end()); }
        }
    }
    tr.build_ms = ms_since(t0);
    return tr;
}

vector<int> bm_labels_at(const BmTrace &tr, size_t cut) {
    vector<int> p(tr.n); std::iota(p.begin(), p.end(), 0);
    std::function<int(int)> find = [&](int x) {
        while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; }
        return x;
    };
    cut = std::min(cut, tr.merges.size());
    for (size_t t = 0; t < cut; ++t) {
        int a = find(tr.merges[t].first), b = find(tr.merges[t].second);
        if (a != b) p[b] = a;
    }
    std::map<int, int> remap;
    vector<int> out(tr.n);
    for (int u = 0; u < tr.n; ++u) {
        int r = find(u);
        if (!remap.count(r)) remap[r] = (int)remap.size();
        out[u] = remap[r];
    }
    return out;
}

// Uniform stub matching: exact degree sequence, loops and parallel edges kept
// as weighted edges.  The null is the configuration model used by modularity.
Graph configuration_null(const Graph &g, std::mt19937 &rng) {
    vector<int> stubs;
    stubs.reserve((size_t)std::llround(g.m2));
    for (int u = 0; u < g.n(); ++u) {
        long long d = std::llround(g.degree[u]);
        for (long long q = 0; q < d; ++q) stubs.push_back(u);
    }
    std::shuffle(stubs.begin(), stubs.end(), rng);
    vector<std::pair<int, int>> edges; edges.reserve(stubs.size() / 2);
    for (size_t k = 1; k < stubs.size(); k += 2) edges.push_back({stubs[k - 1], stubs[k]});
    return Graph::from_edges(g.n(), edges);
}

struct BmBootstrapResult {
    BmTrace observed;
    vector<double> p_stage;
    size_t cut = 0;                    // significant merges retained
    int bootstraps = 0;
    double alpha = 0.05;
    double total_ms = 0.0;
};

BmBootstrapResult bm_path_bootstrap(const Graph &g, int B, double alpha, unsigned seed) {
    auto t0 = Clock::now();
    BmBootstrapResult out; out.bootstraps = B; out.alpha = alpha;
    out.observed = bm_greedy_trace(g);
    vector<vector<double>> nullz;
    nullz.reserve(B);
    std::mt19937 rng(seed);
    for (int b = 0; b < B; ++b) {
        Graph ng = configuration_null(g, rng);
        BmTrace nt = bm_greedy_trace(ng, out.observed.z.size());
        nullz.push_back(std::move(nt.z));
    }
    out.cut = out.observed.z.size();
    for (size_t t = 0; t < out.observed.z.size(); ++t) {
        int tail = 0;
        for (const auto &z : nullz) {
            double nz = t < z.size() ? z[t] : 0.0;
            if (nz >= out.observed.z[t] - 1e-12) ++tail;
        }
        double p = (double)(tail + 1) / (double)(B + 1);
        out.p_stage.push_back(p);
        // This is a stagewise parametric-bootstrap test.  The first non-significant
        // maximum is not merged, so the reported partition contains only evidence
        // exceeding the chosen configuration-null level at every retained stage.
        if (p > alpha) { out.cut = t; break; }
    }
    out.total_ms = ms_since(t0);
    return out;
}

// --------------------------------------------------------------------------- //
//  BM-CV: sample-split, familywise-calibrated Bader residual selection.
//
//  The literal path test above has no power at singleton scale: a selected edge
//  is as surprising in a real graph as in a configuration graph.  BM-CV fixes
//  that without pretending the selected edge is an ordinary z-test.  One edge
//  split constructs a Bader dendrogram; the independent split evaluates every
//  cut by the total within-module Bader residual
//
//      R(P) = sum_C [ l_CC - S_C^2/(4m) ].
//
//  A configuration bootstrap of the VALIDATION graph evaluates the same fixed
//  discovery tree.  Its maximum standardized residual over every candidate cut
//  is the null reference, so the reported `p_global` includes the cut search.
//  This is a real statistical denomination: conditional on the discovery tree,
//  it is a Monte-Carlo familywise p-value under the degree-preserving null.
// --------------------------------------------------------------------------- //
struct EdgeSplit { Graph discovery, validation; };

EdgeSplit split_edges(const Dataset &d, double discovery_fraction, unsigned seed) {
    std::mt19937 rng(seed);
    std::bernoulli_distribution keep(discovery_fraction);
    vector<std::pair<int, int>> a, b; a.reserve(d.edges.size()); b.reserve(d.edges.size());
    for (auto e : d.edges) (keep(rng) ? a : b).push_back(e);
    // Keep both sides nonempty on tiny graphs; this is deterministic given seed.
    if (b.empty() && !a.empty()) { b.push_back(a.back()); a.pop_back(); }
    if (a.empty() && !b.empty()) { a.push_back(b.back()); b.pop_back(); }
    return {Graph::from_edges(d.n, a), Graph::from_edges(d.n, b)};
}

// Evaluate the discovery tree against another graph.  Each entry is the sum of
// Bader's diagonal residuals after the corresponding number of discovery merges.
vector<double> bm_within_residual_trace(const Graph &g,
                                        const vector<std::pair<int, int>> &merges) {
    vector<double> out;
    if (g.m2 <= 0.0) { out.assign(merges.size(), 0.0); return out; }
    int N = g.n();
    vector<std::unordered_map<int, double>> nbr(N);
    for (int u = 0; u < N; ++u)
        for (auto &e : g.adj[u]) if (e.first != u) nbr[u][e.first] += e.second;
    vector<double> sigma = g.degree;
    double lin = 0.0, sumsq = 0.0;
    for (double s : sigma) sumsq += s * s;
    out.reserve(merges.size());
    for (auto ab : merges) {
        int a = ab.first, b = ab.second;
        double wab = 0.0;
        auto it = nbr[a].find(b); if (it != nbr[a].end()) wab = it->second;
        lin += wab;
        sumsq += 2.0 * sigma[a] * sigma[b];
        nbr[a].erase(b);
        for (auto &kv : nbr[b]) {
            int x = kv.first;
            if (x == a) continue;
            nbr[a][x] += kv.second;
            nbr[x].erase(b);
            nbr[x][a] = nbr[a][x];
        }
        nbr[b].clear(); sigma[a] += sigma[b]; sigma[b] = 0.0;
        out.push_back(lin - sumsq / (2.0 * g.m2));
    }
    return out;
}

struct BmCvResult {
    BmTrace discovery_tree;
    size_t cut = 0;
    bool significant = false;
    double z_best = 0.0;
    double p_global = 1.0;
    int bootstraps = 0;
    double total_ms = 0.0;
};

BmCvResult bm_cv_calibrate(const Graph &validation, BmTrace tree, int B, double alpha,
                            unsigned seed) {
    auto t0 = Clock::now();
    BmCvResult out; out.bootstraps = B; out.discovery_tree = std::move(tree);
    const size_t T = out.discovery_tree.merges.size();
    if (T == 0 || validation.m2 <= 0.0) { out.total_ms = ms_since(t0); return out; }

    vector<double> obs = bm_within_residual_trace(validation, out.discovery_tree.merges);
    vector<vector<double>> nul; nul.reserve(B);
    std::mt19937 rng(seed ^ 0x9E3779B9u);
    for (int b = 0; b < B; ++b)
        nul.push_back(bm_within_residual_trace(configuration_null(validation, rng), out.discovery_tree.merges));

    vector<double> mu(T, 0.0), sd(T, 0.0), zobs(T, -1e300);
    for (size_t t = 0; t < T; ++t) {
        for (int b = 0; b < B; ++b) mu[t] += nul[b][t];
        mu[t] /= B;
        for (int b = 0; b < B; ++b) sd[t] += (nul[b][t] - mu[t]) * (nul[b][t] - mu[t]);
        sd[t] = std::sqrt(sd[t] / B);
        if (sd[t] > 1e-12) zobs[t] = (obs[t] - mu[t]) / sd[t];
    }
    out.cut = (size_t)(std::max_element(zobs.begin(), zobs.end()) - zobs.begin()) + 1;
    out.z_best = zobs[out.cut - 1];

    int tail = 0;
    for (int b = 0; b < B; ++b) {
        double maxz = -1e300;
        for (size_t t = 0; t < T; ++t)
            if (sd[t] > 1e-12) maxz = std::max(maxz, (nul[b][t] - mu[t]) / sd[t]);
        if (maxz >= out.z_best - 1e-12) ++tail;
    }
    out.p_global = (double)(tail + 1) / (double)(B + 1);
    out.significant = out.p_global <= alpha;
    if (!out.significant) out.cut = 0;
    out.total_ms = ms_since(t0);
    return out;
}

// Literal Bader ranking for discovery, followed by independent calibration.
BmCvResult bm_cv_select(const Dataset &d, int B, double alpha, unsigned seed) {
    EdgeSplit sp = split_edges(d, 0.60, seed);
    return bm_cv_calibrate(sp.validation, bm_greedy_trace(sp.discovery), B, alpha, seed);
}

// The practical correction: CNM supplies a deterministic candidate dendrogram;
// Bader residuals and the held-out configuration bootstrap decide WHICH cut is
// statistically supported.  This avoids using a residual z-score as an unsafe
// local merge ranking while preserving Bader's actual residual criterion.
BmCvResult cnm_bm_cv_select(const Dataset &d, int B, double alpha, unsigned seed) {
    EdgeSplit sp = split_edges(d, 0.60, seed);
    CnmTrace cnm = cnm_greedy_trace(sp.discovery);
    BmTrace tree; tree.n = d.n; tree.merges = std::move(cnm.merges);
    return bm_cv_calibrate(sp.validation, std::move(tree), B, alpha, seed);
}

// Where the rolling rule cuts a recorded dQ trace (merge count to keep).
size_t rollstop_cut(const CnmTrace &tr, int W, double kmult, bool two_sided,
                    RollStop *out = nullptr) {
    RollStop rs; rs.W = W; rs.kmult = kmult; rs.two_sided = two_sided;
    size_t cut = tr.dq.size();
    for (size_t t = 0; t < tr.dq.size(); ++t)
        if (rs.test(tr.dq[t])) { cut = t; break; }  // exclude the tripping merge
    if (out) *out = rs;
    return cut;
}

// --------------------------------------------------------------------------- //
//  Scale-free variants of the same idea.
//
//  The raw rule fails on large graphs for a measurable reason: dQ ~ 1/m, so as m
//  grows the accepted dQ values crowd together, the window sd collapses toward
//  zero, and ANY small dip clears k*sd. On Amazon that cuts the run at 0.7% of
//  its merges. The test is measuring absolute spread, which is not scale-free.
//
//  LOG: run the identical mean/sd test on log(dQ). A multiplicative drop is then
//  a fixed distance regardless of the absolute magnitude of dQ, so the criterion
//  no longer depends on m. This is the smallest change that makes the statistic
//  scale-invariant while keeping the user's formulation intact.
//
//  REL: a floor on relative dispersion -- require the drop to exceed k sd's AND
//  be a real fraction of the window mean, so a window that is merely flat cannot
//  fire on numerical noise.
// --------------------------------------------------------------------------- //
size_t rollstop_cut_log(const CnmTrace &tr, int W, double kmult, RollStop *out = nullptr) {
    RollStop rs; rs.W = W; rs.kmult = kmult; rs.two_sided = false;
    size_t cut = tr.dq.size();
    for (size_t t = 0; t < tr.dq.size(); ++t) {
        double v = tr.dq[t];
        if (v <= 0.0) { cut = t; break; }
        if (rs.test(std::log(v))) { cut = t; break; }
    }
    if (out) *out = rs;
    return cut;
}

// DIFF: test the STEP, not the level. The window still supplies the dispersion
// estimate (sd over the last W dQ values), but the quantity compared against
// k*sd is |dQ_prev - dQ_now| -- the change from the immediately preceding merge.
//
// Why this can behave differently from the level test: a greedy dQ trace drifts
// downward, so the level test measures a value against a mean it has already
// drifted away from, and on large graphs that drift dominates the collapsed sd.
// The step test is drift-insensitive -- consecutive merges share the same local
// scale -- so it fires on an abrupt break rather than on accumulated drift.
//
// two_sided here means |dQ_prev - dQ_now| (any abrupt change); one-sided means
// (dQ_prev - dQ_now), i.e. only an abrupt DROP counts.
size_t rollstop_cut_diff(const CnmTrace &tr, int W, double kmult, bool two_sided,
                         RollStop *out = nullptr) {
    std::deque<double> win;
    size_t cut = tr.dq.size();
    RollStop diag; diag.W = W; diag.kmult = kmult; diag.two_sided = two_sided;
    for (size_t t = 0; t < tr.dq.size(); ++t) {
        double dq = tr.dq[t];
        if ((int)win.size() >= W && t > 0) {
            double mu = 0.0;
            for (double x : win) mu += x;
            mu /= (double)win.size();
            double ss = 0.0;
            for (double x : win) ss += (x - mu) * (x - mu);
            double sd = std::sqrt(ss / (double)win.size());
            double prev = tr.dq[t - 1];
            double step = two_sided ? std::fabs(prev - dq) : (prev - dq);
            if (sd > 0.0 && step > kmult * sd) {
                diag.fired_at = (long long)t + 1; diag.trip_dq = dq;
                diag.trip_mu = mu; diag.trip_sd = sd;
                cut = t; break;
            }
        }
        win.push_back(dq);
        if ((int)win.size() > W) win.pop_front();
    }
    if (out) *out = diag;
    return cut;
}

// min_rel = required drop as a fraction of the window mean (e.g. 0.5 = halved).
size_t rollstop_cut_rel(const CnmTrace &tr, int W, double kmult, double min_rel,
                        RollStop *out = nullptr) {
    std::deque<double> win;
    size_t cut = tr.dq.size();
    RollStop diag; diag.W = W; diag.kmult = kmult;
    for (size_t t = 0; t < tr.dq.size(); ++t) {
        double dq = tr.dq[t];
        if ((int)win.size() >= W) {
            double mu = 0.0;
            for (double x : win) mu += x;
            mu /= (double)win.size();
            double ss = 0.0;
            for (double x : win) ss += (x - mu) * (x - mu);
            double sd = std::sqrt(ss / (double)win.size());
            if (sd > 0.0 && mu > 0.0 && (mu - dq) > kmult * sd &&
                (mu - dq) > min_rel * mu) {
                diag.fired_at = (long long)t + 1; diag.trip_dq = dq;
                diag.trip_mu = mu; diag.trip_sd = sd;
                cut = t; break;
            }
        }
        win.push_back(dq);
        if ((int)win.size() > W) win.pop_front();
    }
    if (out) *out = diag;
    return cut;
}

// --------------------------------------------------------------------------- //
// David-review fixes 1--3: a common, replay-only stopping-rule layer.
//
// The modifiers deliberately live outside `cnm_greedy_trace`: they never alter
// the priority queue, candidate scores, or accepted merge list.  This makes it
// possible to compare an unmodified and a modified rule on exactly the same
// dendrogram.
// --------------------------------------------------------------------------- //
enum class Dispersion { SD, MAD };
struct StopModifiers {
    bool defer_singletons = false;
    Dispersion dispersion = Dispersion::SD;
    double epsilon = 1e-6;
};

double median_of(vector<double> x) {
    if (x.empty()) return 0.0;
    size_t k = x.size() / 2;
    std::nth_element(x.begin(), x.begin() + k, x.end());
    double hi = x[k];
    if (x.size() & 1) return hi;
    std::nth_element(x.begin(), x.begin() + k - 1, x.end());
    return 0.5 * (hi + x[k - 1]);
}

double dispersion_of(const vector<double> &x, Dispersion which) {
    if (x.empty()) return 0.0;
    if (which == Dispersion::MAD) {
        double med = median_of(x);
        vector<double> dev; dev.reserve(x.size());
        for (double v : x) dev.push_back(std::fabs(v - med));
        return 1.4826 * median_of(std::move(dev));
    }
    double mu = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
    double ss = 0.0;
    for (double v : x) ss += (v - mu) * (v - mu);
    return std::sqrt(ss / x.size());
}

double mean_of(const vector<double> &x) {
    return x.empty() ? 0.0 : std::accumulate(x.begin(), x.end(), 0.0) / x.size();
}

bool may_test_at(const CnmTrace &tr, size_t t, const StopModifiers &mod) {
    return !mod.defer_singletons || t >= tr.singleton_phase_end;
}

// MAD is a robust replacement, not an invitation to divide by a near-zero
// plateau.  The relative floor is intentionally enabled for MAD only: SD is
// retained as the exact historical baseline, while `--dispersion=mad` is the
// requested plateau-protected modifier.
bool degenerate_window(const vector<double> &window, double scale,
                       const StopModifiers &mod) {
    if (scale <= 0.0) return true;
    if (mod.dispersion != Dispersion::MAD) return false;
    return scale / std::max(std::fabs(mean_of(window)), 1e-300) < mod.epsilon;
}

// `kind` is level (`raw` / `log`) or a first-difference test (`diff`).
size_t rollstop_cut_mod(const CnmTrace &tr, int W, double kmult, bool two_sided,
                        const std::string &kind, const StopModifiers &mod) {
    vector<double> win;
    size_t cut = tr.dq.size();
    for (size_t t = 0; t < tr.dq.size(); ++t) {
        double value = tr.dq[t];
        if (kind == "log") {
            if (value <= 0.0) { if (may_test_at(tr, t, mod)) { cut = t; break; } }
            value = std::log(std::max(value, 1e-300));
        }
        if ((int)win.size() == W && may_test_at(tr, t, mod)) {
            double scale = dispersion_of(win, mod.dispersion);
            if (!degenerate_window(win, scale, mod)) {
                double statistic = kind == "diff"
                    ? (two_sided ? std::fabs(tr.dq[t - 1] - tr.dq[t])
                                 : tr.dq[t - 1] - tr.dq[t])
                    : (two_sided ? std::fabs(value - mean_of(win))
                                 : mean_of(win) - value);
                if (statistic > kmult * scale) { cut = t; break; }
            }
        }
        win.push_back(value);
        if ((int)win.size() > W) win.erase(win.begin());
    }
    return cut;
}

struct DetrendedPoint {
    bool ready = false;
    double residual = 0.0;
    double scale = 0.0;
    vector<double> log_window;
    vector<double> residual_window;
};

// Fit y_i=a+b*i to the W observations immediately before t.  The residual
// scale is calculated from the fitted residuals, never from the declining raw
// dQ level.  This is the explicit Fix-1 distinction from rollstop_cut_log.
DetrendedPoint detrended_log_point(const CnmTrace &tr, size_t t, int W,
                                   Dispersion dispersion) {
    DetrendedPoint out;
    if (t < (size_t)W || t >= tr.dq.size() || tr.dq[t] <= 0.0) return out;
    out.log_window.reserve(W);
    for (size_t i = t - W; i < t; ++i) {
        if (tr.dq[i] <= 0.0) return out;
        out.log_window.push_back(std::log(tr.dq[i]));
    }
    const double xbar = 0.5 * (W - 1), ybar = mean_of(out.log_window);
    double sxx = 0.0, sxy = 0.0;
    for (int i = 0; i < W; ++i) {
        double dx = i - xbar;
        sxx += dx * dx; sxy += dx * (out.log_window[i] - ybar);
    }
    double b = sxx > 0.0 ? sxy / sxx : 0.0;
    double a = ybar - b * xbar;
    out.residual_window.reserve(W);
    for (int i = 0; i < W; ++i)
        out.residual_window.push_back(out.log_window[i] - (a + b * i));
    out.residual = std::log(tr.dq[t]) - (a + b * W);
    out.scale = dispersion_of(out.residual_window, dispersion);
    out.ready = true;
    return out;
}

size_t rollstop_cut_detrend(const CnmTrace &tr, int W, double kmult, bool two_sided,
                            const StopModifiers &mod) {
    for (size_t t = 0; t < tr.dq.size(); ++t) {
        if (!may_test_at(tr, t, mod)) continue;
        DetrendedPoint p = detrended_log_point(tr, t, W, mod.dispersion);
        if (!p.ready || degenerate_window(p.log_window, p.scale, mod)) continue;
        double v = two_sided ? std::fabs(p.residual) : -p.residual;
        if (v > kmult * p.scale) return t;
    }
    return tr.dq.size();
}

// Fix 2.  Work in residual-standard-deviation units; otherwise delta=1 would
// depend on the arbitrary log-dQ scale.  `h` is therefore the single,
// calibratable control parameter for both procedures.
size_t rollstop_cut_cusum(const CnmTrace &tr, int W, double h, double delta,
                          const StopModifiers &mod, bool shiryaev_roberts) {
    double stat = 0.0;
    for (size_t t = 0; t < tr.dq.size(); ++t) {
        if (!may_test_at(tr, t, mod)) continue;
        DetrendedPoint p = detrended_log_point(tr, t, W, mod.dispersion);
        if (!p.ready || degenerate_window(p.log_window, p.scale, mod)) continue;
        double x = p.residual / p.scale;
        if (shiryaev_roberts) {
            double lr = std::exp(std::min(700.0, -delta * x - 0.5 * delta * delta));
            stat = std::min(1e300, (1.0 + stat) * lr);
        } else {
            stat = std::max(0.0, stat + (-x) - delta / 2.0);
        }
        if (stat > h) return t;
    }
    return tr.dq.size();
}

// Fix 4.  These values were recorded from the exact candidate CNM selected;
// adaptive_tau is recomputed with the current number of adjacent pairs.
size_t bm_native_stop_cut(const CnmTrace &tr, double tau) {
    for (size_t t = 0; t < tr.merge_z.size(); ++t) {
        double m = std::max(3LL, tr.candidate_pairs[t]);
        double tau_eff = tau * std::sqrt(2.0 * std::log(m));
        if (tr.merge_z[t] <= tau_eff) return t;
    }
    return tr.merges.size();
}

size_t cm_baseline_cut(const CnmTrace &tr, double m2, double c) {
    if (m2 <= 0.0) return tr.merges.size();
    double threshold = 2.0 * c / m2; // c / m, where m=m2/2.
    for (size_t t = 0; t < tr.dq.size(); ++t)
        if (tr.dq[t] < threshold) return t;
    return tr.merges.size();
}

// --------------------------------------------------------------------------- //
//  Why the rule fires where it does: dump the head of the dQ trace and the window
//  statistics at the firing point. This is the diagnostic that shows the failure
//  mode is a property of the dQ SEQUENCE, not of the implementation.
// --------------------------------------------------------------------------- //
void dump_trace_diagnostics(const std::string &name, const CnmTrace &tr) {
    std::cout << "\n-- dQ trace diagnostics: " << name << " (" << tr.merges.size()
              << " merges) --\n";
    // Scientific notation on purpose: dQ scales as ~1/m, so on the larger graphs
    // every value is ~1e-5 and the default 3-decimal format prints them all as
    // "0.000", which hides the very scale gap that drives the failure mode.
    std::cout.unsetf(std::ios::fixed);
    std::cout << std::scientific;
    std::cout.precision(3);
    std::cout << "  first 15 dQ:";
    for (size_t i = 0; i < tr.dq.size() && i < 15; ++i) std::cout << " " << tr.dq[i];
    std::cout << "\n";
    if (tr.dq.size() > 15) {
        std::cout << "  quartile dQ: ";
        for (int q = 1; q <= 4; ++q) {
            size_t idx = (size_t)((double)q / 4.0 * (tr.dq.size() - 1));
            std::cout << " [" << idx << "]=" << tr.dq[idx];
        }
        std::cout << "\n";
    }
    // Coefficient of variation of the first window is the diagnostic number: it
    // says how much RELATIVE spread the test has to work with at the start.
    if (tr.dq.size() >= 20) {
        double mu = 0.0;
        for (int i = 0; i < 20; ++i) mu += tr.dq[i];
        mu /= 20.0;
        double ss = 0.0;
        for (int i = 0; i < 20; ++i) ss += (tr.dq[i] - mu) * (tr.dq[i] - mu);
        double sd = std::sqrt(ss / 20.0);
        std::cout << "  first-20 window: mu=" << mu << " sd=" << sd
                  << " cv=" << (mu != 0.0 ? sd / mu : 0.0) << "\n";
    }
    std::cout.unsetf(std::ios::scientific);
    std::cout << std::fixed;
    std::cout.precision(3);
    for (int W : {10, 15, 20}) {
        RollStop rs;
        size_t cut = rollstop_cut(tr, W, 2.0, false, &rs);
        std::cout << "  W=" << W << " k=2 1sided: cut=" << cut << "/" << tr.merges.size();
        if (rs.fired_at > 0)
            std::cout << "  window mu=" << rs.trip_mu << " sd=" << rs.trip_sd
                      << " dq=" << rs.trip_dq
                      << "  (mu-dq)/sd=" << (rs.trip_sd > 0 ? (rs.trip_mu - rs.trip_dq) / rs.trip_sd : 0.0);
        else
            std::cout << "  never fired";
        std::cout << "\n";
    }
    std::cout.flush();
}

// --------------------------------------------------------------------------- //
//  Metrics.
// --------------------------------------------------------------------------- //
int num_communities(const vector<int> &c) {
    return (int)std::set<int>(c.begin(), c.end()).size();
}

// counts of community sizes (descending)
vector<int> community_sizes(const vector<int> &c) {
    std::unordered_map<int, int> cnt;
    for (int x : c) cnt[x]++;
    vector<int> s;
    for (auto &kv : cnt) s.push_back(kv.second);
    std::sort(s.rbegin(), s.rend());
    return s;
}

double frac_nonsingleton(const vector<int> &c) {
    std::unordered_map<int, int> cnt;
    for (int x : c) ++cnt[x];
    long long kept = 0;
    for (int x : c) if (cnt[x] >= 2) ++kept;
    return c.empty() ? 0.0 : (double)kept / c.size();
}

// NMI over the subset of nodes present in ground truth (gt[u] < 0 = absent)
double nmi(const vector<int> &pred, const vector<int> &gt) {
    vector<int> a, b;
    for (size_t u = 0; u < gt.size(); ++u)
        if (gt[u] >= 0) { a.push_back(pred[u]); b.push_back(gt[u]); }
    if (a.empty()) return 0.0;
    double N = (double)a.size();
    std::map<int, double> ca, cb;
    std::map<std::pair<int, int>, double> cab;
    for (size_t i = 0; i < a.size(); ++i) { ca[a[i]]++; cb[b[i]]++; cab[{a[i], b[i]}]++; }
    double Ha = 0, Hb = 0, I = 0;
    for (auto &kv : ca) { double p = kv.second / N; Ha -= p * std::log(p); }
    for (auto &kv : cb) { double p = kv.second / N; Hb -= p * std::log(p); }
    for (auto &kv : cab) {
        double pij = kv.second / N;
        double pi = ca[kv.first.first] / N, pj = cb[kv.first.second] / N;
        I += pij * std::log(pij / (pi * pj));
    }
    if (Ha == 0 && Hb == 0) return 1.0;
    return I / ((Ha + Hb) / 2.0);
}

// best-match averaged F1 (symmetric), over nodes present in ground truth
double community_f1(const vector<int> &pred, const vector<int> &gt) {
    vector<int> a, b;
    for (size_t u = 0; u < gt.size(); ++u)
        if (gt[u] >= 0) { a.push_back(pred[u]); b.push_back(gt[u]); }
    if (a.empty()) return 0.0;
    // Build only the sparse contingency table.  The previous implementation
    // materialized node sets and compared every predicted community to every
    // truth community, which is quadratic for the singleton-heavy early CNM
    // cuts (and made DBLP metric generation look like a hang).
    std::unordered_map<int, long long> np, nt;
    std::unordered_map<unsigned long long, long long> nij;
    for (size_t i = 0; i < a.size(); ++i) {
        ++np[a[i]]; ++nt[b[i]];
        unsigned long long key = (static_cast<unsigned long long>(static_cast<unsigned>(a[i])) << 32) |
                                 static_cast<unsigned>(b[i]);
        ++nij[key];
    }
    double pred_total = 0.0, truth_total = 0.0;
    std::unordered_map<int, double> best_pred, best_truth;
    for (auto &kv : nij) {
        int p = (int)(kv.first >> 32), t = (int)(kv.first & 0xffffffffu);
        double inter = kv.second;
        double f = 2.0 * inter / (np[p] + nt[t]);
        best_pred[p] = std::max(best_pred[p], f);
        best_truth[t] = std::max(best_truth[t], f);
    }
    for (auto &kv : np) pred_total += best_pred[kv.first];
    for (auto &kv : nt) truth_total += best_truth[kv.first];
    return 0.5 * (pred_total / np.size() + truth_total / nt.size());
}

double phi(double x) { return 0.5 * (1 + std::erf(x / std::sqrt(2.0))); }
double g_tau(double tau) { return (std::sqrt(tau * tau + 4) - tau) / 2.0; }

// --------------------------------------------------------------------------- //
//  Dataset generators (synthetic; self-contained).
// --------------------------------------------------------------------------- //
Dataset ring_of_cliques(int n, int c) {
    Dataset d; d.n = n * c; d.gt.assign(n * c, 0);
    vector<int> starts;
    int off = 0;
    for (int i = 0; i < n; ++i) {
        int s0 = off;
        starts.push_back(s0);
        for (int a = 0; a < c; ++a) {
            d.gt[off + a] = i;
            for (int b = a + 1; b < c; ++b) d.edges.push_back({off + a, off + b});
        }
        off += c;
    }
    for (int i = 0; i < n; ++i) d.edges.push_back({starts[i], starts[(i + 1) % n]});
    return d;
}

Dataset erdos_renyi(int N, double avg_degree, unsigned seed) {
    Dataset d; d.n = N; d.gt.assign(N, -1);  // no ground truth
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> U(0, 1);
    double p = avg_degree / (N - 1);
    for (int i = 0; i < N; ++i)
        for (int j = i + 1; j < N; ++j)
            if (U(rng) < p) d.edges.push_back({i, j});
    return d;
}

// two-level hierarchical SBM: 2 super x 5 sub x 50 nodes (10 sub-blocks)
Dataset hierarchical_sbm(unsigned seed, vector<int> &gt_super) {
    Dataset d;
    int n_super = 2, n_sub = 5, sub_size = 50;
    int N = n_super * n_sub * sub_size;
    d.n = N; d.gt.assign(N, 0); gt_super.assign(N, 0);
    double p_sub = 0.30, p_super = 0.05, p_cross = 0.005;
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> U(0, 1);
    int node = 0;
    for (int s = 0; s < n_super; ++s)
        for (int b = 0; b < n_sub; ++b)
            for (int k = 0; k < sub_size; ++k) {
                d.gt[node] = s * n_sub + b; gt_super[node] = s; node++;
            }
    for (int i = 0; i < N; ++i)
        for (int j = i + 1; j < N; ++j) {
            double p = (d.gt[i] == d.gt[j]) ? p_sub
                       : (gt_super[i] == gt_super[j]) ? p_super : p_cross;
            if (U(rng) < p) d.edges.push_back({i, j});
        }
    return d;
}

// --------------------------------------------------------------------------- //
//  Edge-list + label loaders (for SNAP / karate / polblogs via TSV files).
//  Node ids are remapped to a contiguous 0..N-1 range; keeps the giant-or-all
//  induced structure exactly as produced by the Python data-prep step.
// --------------------------------------------------------------------------- //
Dataset load_edgelist(const std::string &edge_path, const std::string &label_path) {
    Dataset d;
    std::unordered_map<long long, int> remap;
    auto id = [&](long long x) {
        auto it = remap.find(x);
        if (it != remap.end()) return it->second;
        int nid = (int)remap.size(); remap[x] = nid; return nid;
    };
    std::ifstream ef(edge_path);
    std::string line;
    while (std::getline(ef, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        long long a, b; if (!(ss >> a >> b)) continue;
        int ia = id(a); int ib = id(b);  // evaluate in a defined order
        d.edges.push_back({ia, ib});
    }
    d.n = (int)remap.size();
    d.gt.assign(d.n, -1);
    if (!label_path.empty()) {
        std::ifstream lf(label_path);
        while (std::getline(lf, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            long long node; int lab; if (!(ss >> node >> lab)) continue;
            auto it = remap.find(node);
            if (it != remap.end()) d.gt[it->second] = lab;
        }
    }
    return d;
}

// --------------------------------------------------------------------------- //
//  Experiments.
// --------------------------------------------------------------------------- //
double mean(const vector<double> &v) {
    return v.empty() ? 0.0 : std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}
double pstdev(const vector<double> &v) {
    if (v.empty()) return 0.0;
    double mu = mean(v), s = 0;
    for (double x : v) s += (x - mu) * (x - mu);
    return std::sqrt(s / v.size());
}

void experiment_1() {
    std::cout << "\n=== Experiment 1: Ring of cliques (resolution limit) ===\n";
    int c = 5, n_seeds = 5;
    vector<int> ns = {20, 25, 30, 50, 100, 200, 400};
    std::cout << "n\ttruth\tLouvain\tLeiden\tCPM\tGRAZt2\tGRAZt3\n";
    for (int n : ns) {
        Dataset d = ring_of_cliques(n, c);
        Graph g = Graph::from_edges(d.n, d.edges);
        vector<double> kl, kle, kc, k2, k3;
        for (int s = 0; s < n_seeds; ++s) {
            kl.push_back(num_communities(louvain(g, s)));
            kle.push_back(num_communities(leiden(g, s)));
            kc.push_back(num_communities(cpm_leiden(g, 0.1, s)));
            k2.push_back(num_communities(graz(g, 2.0, s)));
            k3.push_back(num_communities(graz(g, 3.0, s)));
        }
        std::cout << n << "\t" << n << "\t" << (int)mean(kl) << "\t" << (int)mean(kle)
                  << "\t" << (int)mean(kc) << "\t" << (int)mean(k2) << "\t"
                  << (int)mean(k3) << "\n";
    }
    std::cout << "  (Theorem 1 ceilings: GRAZ t=2 -> n<=" << (c*c-c+2)/(g_tau(2)*g_tau(2))
              << ", t=3 -> n<=" << (c*c-c+2)/(g_tau(3)*g_tau(3)) << ")\n";
}

void experiment_2() {
    std::cout << "\n=== Experiment 2: Erdos-Renyi nulls (spurious communities) ===\n";
    vector<int> Ns = {1000, 2000, 5000, 10000};
    int n_seeds = 3;
    std::cout << "N\tm\tLouvain\tLeiden\tGRAZt2\tGRAZt3\tGuimera\n";
    for (int N : Ns) {
        vector<double> kl, kle, k2, k3, ms;
        for (int s = 0; s < n_seeds; ++s) {
            Dataset d = erdos_renyi(N, 10.0, s);
            Graph g = Graph::from_edges(d.n, d.edges);
            ms.push_back((double)d.edges.size());
            auto big = [](const vector<int> &c) {
                auto s = community_sizes(c); int k = 0;
                for (int x : s) if (x >= 3) k++; return (double)k; };
            kl.push_back(big(louvain(g, s)));
            kle.push_back(big(leiden(g, s)));
            k2.push_back(big(graz(g, 2.0, s)));
            k3.push_back(big(graz(g, 3.0, s)));
        }
        double mbar = mean(ms);
        std::cout << N << "\t" << (int)mbar << "\t" << (int)mean(kl) << "\t"
                  << (int)mean(kle) << "\t" << (int)mean(k2) << "\t" << (int)mean(k3)
                  << "\t" << (int)std::sqrt(mbar / 2.0) << "\n";
    }
}

void experiment_4() {
    std::cout << "\n=== Experiment 4: Hierarchical SBM (multi-scale) ===\n";
    int n_seeds = 5;
    std::cout << "method\tK\tNMI_sub\tNMI_super\n";
    struct Acc { vector<double> K, ns, np; };
    std::map<std::string, Acc> acc;
    for (int s = 0; s < n_seeds; ++s) {
        vector<int> gt_super;
        Dataset d = hierarchical_sbm(s, gt_super);
        Graph g = Graph::from_edges(d.n, d.edges);
        auto rec = [&](const std::string &name, const vector<int> &lab) {
            acc[name].K.push_back(num_communities(lab));
            acc[name].ns.push_back(nmi(lab, d.gt));
            acc[name].np.push_back(nmi(lab, gt_super));
        };
        rec("Louvain", louvain(g, s));
        rec("Leiden", leiden(g, s));
        rec("CPM-Leiden", cpm_leiden(g, 0.1, s));
        rec("GRAZ t=2", graz(g, 2.0, s));
        rec("GRAZ t=3", graz(g, 3.0, s));
    }
    for (const char *m : {"Louvain", "Leiden", "CPM-Leiden", "GRAZ t=2", "GRAZ t=3"}) {
        auto &a = acc[m];
        std::cout << m << "\t" << mean(a.K) << "\t" << mean(a.ns) << "\t"
                  << mean(a.np) << "\n";
    }
}

void experiment_5() {
    std::cout << "\n=== Experiment 5: Stability under node reordering ===\n";
    // hierarchical SBM instance (self-contained, needs no download)
    vector<int> gt_super;
    Dataset d = hierarchical_sbm(3, gt_super);
    Graph g = Graph::from_edges(d.n, d.edges);
    int runs = 40;
    std::cout << "method\tmeanK\tstdK\tmeanPairwiseNMI\n";
    auto eval = [&](const std::string &name, auto fn) {
        vector<vector<int>> R;
        vector<double> Ks;
        for (int s = 0; s < runs; ++s) { R.push_back(fn(s)); Ks.push_back(num_communities(R.back())); }
        vector<double> pn;
        for (int i = 0; i < runs; ++i)
            for (int j = i + 1; j < runs; ++j) {
                // pairwise NMI treating one run as "ground truth" of the other
                vector<int> gt = R[j];
                pn.push_back(nmi(R[i], gt));
            }
        std::cout << name << "\t" << mean(Ks) << "\t" << pstdev(Ks) << "\t" << mean(pn) << "\n";
    };
    eval("Louvain", [&](int s) { return louvain(g, s); });
    eval("Leiden", [&](int s) { return leiden(g, s); });
    eval("GRAZ t=2", [&](int s) { return graz(g, 2.0, s); });
    eval("GRAZ t=3", [&](int s) { return graz(g, 3.0, s); });
}

void experiment_6(int argc, char **argv) {
    std::cout << "\n=== Experiment 6: Real networks (from prepared edge lists) ===\n";
    // args: pairs of  name=edgefile[:labelfile]
    if (argc < 3) {
        std::cout << "  (skipped: pass datasets as name=edges.txt[:labels.txt] arguments)\n";
        return;
    }
    std::cout << "dataset\tmethod\tK\tQ\tF1\tNMI\n";
    for (int i = 1; i < argc; ++i) {
        std::string spec = argv[i];
        auto eq = spec.find('=');
        if (eq == std::string::npos) continue;  // not a dataset spec
        std::string name = spec.substr(0, eq);
        std::string rest = spec.substr(eq + 1);
        auto colon = rest.find(':');
        std::string ef = colon == std::string::npos ? rest : rest.substr(0, colon);
        std::string lf = colon == std::string::npos ? "" : rest.substr(colon + 1);
        Dataset d = load_edgelist(ef, lf);
        if (d.n == 0) continue;
        Graph g = Graph::from_edges(d.n, d.edges);
        auto row = [&](const std::string &m, const vector<int> &lab) {
            std::cout << name << "\t" << m << "\t" << num_communities(lab) << "\t"
                      << modularity(g, lab) << "\t" << community_f1(lab, d.gt) << "\t"
                      << nmi(lab, d.gt) << "\n";
        };
        row("Louvain", louvain(g, 1));
        row("Leiden", leiden(g, 1));
        row("GRAZ t=2", graz(g, 2.0, 1));
        row("GRAZ t=3", graz(g, 3.0, 1));
    }
}

// --------------------------------------------------------------------------- //
//  Experiment 7: large real networks with thousands of ground-truth communities.
//
//  Reports quality (K, Q, F1, NMI) AND the timing decomposition of every run,
//  so the cost of the GRAZ z-gate can be read directly against the identical
//  Louvain/Leiden move loop. Times are single-threaded wall clock.
//
//  Ground truth here is SNAP's top-5000 overlapping communities collapsed to a
//  partition (smallest-community-wins), and it covers only a subset of nodes;
//  F1/NMI are computed on covered nodes only and are a lower bound.
// --------------------------------------------------------------------------- //
void experiment_7(const vector<std::string> &specs) {
    std::cout << "\n=== Experiment 7: Large real networks (quality + timing) ===\n";
    if (specs.empty()) {
        std::cout << "  (skipped: pass large=edges[:labels] specs)\n";
        return;
    }
    std::cout << "dataset\tN\tM\tmethod\tK\tQ\tF1\tNMI\t"
                 "total_ms\tmove_ms\taggr_ms\trefine_ms\tlevels\tsweeps\tvisits\tz_tests\n";

    for (const std::string &spec : specs) {
        auto eq = spec.find('=');
        if (eq == std::string::npos) continue;
        std::string name = spec.substr(0, eq);
        std::string rest = spec.substr(eq + 1);
        auto colon = rest.find(':');
        std::string ef = colon == std::string::npos ? rest : rest.substr(0, colon);
        std::string lf = colon == std::string::npos ? "" : rest.substr(colon + 1);

        auto t_load = Clock::now();
        Dataset d = load_edgelist(ef, lf);
        if (d.n == 0) { std::cout << "  (" << name << ": empty/missing)\n"; continue; }
        Graph g = Graph::from_edges(d.n, d.edges);
        double load_ms = ms_since(t_load);
        long long M = (long long)d.edges.size();
        int gt_k = 0;
        { std::set<int> s; for (int x : d.gt) if (x >= 0) s.insert(x); gt_k = (int)s.size(); }
        std::cerr << "  [" << name << "] loaded N=" << d.n << " M=" << M
                  << " gt_communities=" << gt_k << " in " << (long long)load_ms << " ms\n";

        auto run = [&](const std::string &key, const std::string &mname, Rule rule,
                       double tau, double gamma, bool refine) {
            // --methods=... restricts which rows are produced, so a single
            // expensive configuration can be re-run without redoing the rest.
            if (!g_methods.empty() && !g_methods.count(key)) return;
            g_t.reset();
            vector<int> lab = run_multilevel(g, rule, tau, gamma, 1, refine);
            Timings t = g_t;  // snapshot before metric computation
            std::cout << name << "\t" << d.n << "\t" << M << "\t" << mname << "\t"
                      << num_communities(lab) << "\t" << modularity(g, lab) << "\t"
                      << community_f1(lab, d.gt) << "\t" << nmi(lab, d.gt) << "\t"
                      << (long long)t.total_ms << "\t" << (long long)t.local_move_ms << "\t"
                      << (long long)t.aggregate_ms << "\t" << (long long)t.refine_ms << "\t"
                      << t.levels << "\t" << t.move_sweeps << "\t" << t.node_visits << "\t"
                      << t.z_tests << "\n";
            std::cout.flush();
        };

        run("louvain", "Louvain", Rule::Modularity, 0.0, 1.0, false);
        run("leiden", "Leiden", Rule::Modularity, 0.0, 1.0, true);
        run("graz2", "GRAZ t=2", Rule::GRAZ, 2.0, 1.0, true);
        run("graz3", "GRAZ t=3", Rule::GRAZ, 3.0, 1.0, true);
        run("graz4", "GRAZ t=4", Rule::GRAZ, 4.0, 1.0, true);
        // The merge-level gate is the only Corollary 1 fix that helped on the
        // small networks; it has to be checked here too, because the large-graph
        // F1/NMI gains are the reason to keep the z-test at all.
        run("merge2", "GRAZ-merge t=2", Rule::GRAZ_MERGE, 2.0, 1.0, true);
        run("merge3", "GRAZ-merge t=3", Rule::GRAZ_MERGE, 3.0, 1.0, true);
        run("merge4", "GRAZ-merge t=4", Rule::GRAZ_MERGE, 4.0, 1.0, true);
    }
}

// --------------------------------------------------------------------------- //
//  Experiment 8: do the proposed Corollary 1 fixes actually work?
//
//  The literal single-node z-test (Eq. 6) started from singletons is selection-
//  biased: it only ever scores NEIGHBOURING communities, where an edge already
//  exists, and under the configuration-model null that edge is itself the rare
//  event. So nearly every first move clears any fixed tau, Corollary 1's
//  spurious-community bound fails, and GRAZ over-fragments random graphs.
//
//  Two fixes are evaluated against plain GRAZ on both failure sites:
//    A. adaptive tau = tau*sqrt(2 log m)   -- the paper's own Section 3.4
//    B. merge-level gate                   -- the literal 2010 formulation
//
//  Panel 1 is the Corollary 1 test proper: an ER graph has no communities, so
//  the correct answer is "few or none" and Louvain's count is the reference for
//  what modularity alone already imposes. Panel 2 checks the fixes did not just
//  trade the false positives for a loss of real structure.
// --------------------------------------------------------------------------- //
void experiment_8(int argc, char **argv) {
    std::cout << "\n=== Experiment 8: Corollary 1 fixes on ER nulls ===\n";
    std::cout << "N\tm\tLouvain\tGRAZt2\tGRAZt3\tAdaptM\tMergeT2\tMergeT3\tGateT2\tGateT3\n";
    int n_seeds = 3;
    for (int N : {1000, 2000, 5000, 10000}) {
        vector<double> kl, k2, k3, ka, km2, km3, kg2, kg3, ms;
        for (int s = 0; s < n_seeds; ++s) {
            Dataset d = erdos_renyi(N, 10.0, s);
            Graph g = Graph::from_edges(d.n, d.edges);
            ms.push_back((double)d.edges.size());
            auto big = [](const vector<int> &c) {
                auto sz = community_sizes(c); int k = 0;
                for (int x : sz) if (x >= 3) k++; return (double)k; };
            kl.push_back(big(louvain(g, s)));
            k2.push_back(big(graz(g, 2.0, s)));
            k3.push_back(big(graz(g, 3.0, s)));
            ka.push_back(big(graz_adaptive(g, 1.0, s)));
            km2.push_back(big(graz_merge(g, 2.0, s)));
            km3.push_back(big(graz_merge(g, 3.0, s)));
            kg2.push_back(big(graz_gate(g, 2.0, s)));
            kg3.push_back(big(graz_gate(g, 3.0, s)));
        }
        std::cout << N << "\t" << (int)mean(ms) << "\t" << (int)mean(kl) << "\t"
                  << (int)mean(k2) << "\t" << (int)mean(k3) << "\t" << (int)mean(ka)
                  << "\t" << (int)mean(km2) << "\t" << (int)mean(km3)
                  << "\t" << (int)mean(kg2) << "\t" << (int)mean(kg3) << "\n";
    }
    std::cout << "  (ER graphs have NO communities: lower is better; Louvain is the\n"
                 "   reference for what modularity alone already over-imposes.)\n";

    // Mechanism check. Raising tau makes fragmentation WORSE, which is backwards
    // if the gate merely suppresses spurious communities. The size distribution
    // separates the two candidate explanations: either the gate is inventing
    // extra mid-size communities, or it is BLOCKING the merges that would
    // consolidate the initial singletons, leaving stranded fragments behind.
    std::cout << "\n  Size distribution on ER N=5000 (why higher tau is worse):\n";
    std::cout << "  method\tK_all\tK_ge3\tmaxSize\tmedSize\tfrac_in_singletons\n";
    {
        Dataset d = erdos_renyi(5000, 10.0, 0);
        Graph g = Graph::from_edges(d.n, d.edges);
        auto describe = [&](const std::string &nm, const vector<int> &lab) {
            auto sz = community_sizes(lab);
            std::sort(sz.begin(), sz.end());
            long long singles = 0;
            for (int x : sz) if (x == 1) singles += 1;
            std::cout << "  " << nm << "\t" << sz.size() << "\t"
                      << std::count_if(sz.begin(), sz.end(), [](int x) { return x >= 3; })
                      << "\t" << sz.back() << "\t" << sz[sz.size() / 2] << "\t"
                      << (double)singles / (double)g.n() << "\n";
        };
        describe("Louvain  ", louvain(g, 0));
        describe("GRAZ t=2 ", graz(g, 2.0, 0));
        describe("GRAZ t=3 ", graz(g, 3.0, 0));
        describe("GRAZ t=4 ", graz(g, 4.0, 0));

        // The ceiling is analytic. Holding k_i and k_iB fixed, z_iB is strictly
        // decreasing in Sigma_B, so for each tau there is a LARGEST community
        // volume that will still admit a node carrying k_iB edges. Past it the
        // move is rejected however much modularity it would gain -- the growth
        // of a community is what the test blocks, not its spuriousness.
        std::cout << "\n  Analytic volume ceiling: max Sigma_B admitting a node with\n"
                     "  k_i=" << (int)(g.m2 / g.n()) << " (mean degree), 2m=" << (long long)g.m2
                  << ".  Volume/mean-degree ~ community size in nodes.\n";
        std::cout << "  tau\tk_iB=1\tk_iB=2\tk_iB=3\t(as node counts)\n";
        double ki = g.m2 / g.n();
        for (double tau : {2.0, 3.0, 4.0}) {
            std::cout << "  " << tau;
            for (double kiB : {1.0, 2.0, 3.0}) {
                double best = 0.0;
                for (double S = ki; S < g.m2; S += ki)  // scan in one-node steps
                    if (z_statistic(ki, kiB, S, g.m2) > tau) best = S;
                std::cout << "\t" << (int)(best / ki);
            }
            std::cout << "\n";
        }
    }

    std::cout << "\n=== Experiment 8b: the same fixes on real networks ===\n";
    std::cout << "dataset\tmethod\tK\tQ\tF1\tNMI\n";
    for (int i = 1; i < argc; ++i) {
        std::string spec = argv[i];
        auto eq = spec.find('=');
        if (eq == std::string::npos) continue;
        std::string name = spec.substr(0, eq);
        if (name == "large" || spec.rfind("--", 0) == 0) continue;
        std::string rest = spec.substr(eq + 1);
        auto colon = rest.find(':');
        std::string ef = colon == std::string::npos ? rest : rest.substr(0, colon);
        std::string lf = colon == std::string::npos ? "" : rest.substr(colon + 1);
        Dataset d = load_edgelist(ef, lf);
        if (d.n == 0) continue;
        Graph g = Graph::from_edges(d.n, d.edges);
        auto row = [&](const std::string &m, const vector<int> &lab) {
            std::cout << name << "\t" << m << "\t" << num_communities(lab) << "\t"
                      << modularity(g, lab) << "\t" << community_f1(lab, d.gt) << "\t"
                      << nmi(lab, d.gt) << "\n";
        };
        row("Louvain", louvain(g, 1));
        row("GRAZ t=2", graz(g, 2.0, 1));
        row("GRAZ-adapt", graz_adaptive(g, 1.0, 1));
        row("GRAZ-merge t=2", graz_merge(g, 2.0, 1));
        row("GRAZ-merge t=3", graz_merge(g, 3.0, 1));
        row("GRAZ-gate t=2", graz_gate(g, 2.0, 1));
        row("GRAZ-gate t=3", graz_gate(g, 3.0, 1));
    }
}

void verify_theory() {
    std::cout << "\n=== Verification: Corollary 1 (non-imposition on random graphs) ===\n";
    Dataset d = erdos_renyi(1000, 10.0, 0);
    Graph g = Graph::from_edges(d.n, d.edges);
    std::cout << "tau\tbound_2(1-Phi)\tobserved_move_frac\n";
    for (double tau : {2.0, 3.0, 4.0}) {
        // one sweep from singletons: fraction of nodes that move
        int N = g.n(); double m2 = g.m2, m = m2 / 2.0;
        vector<int> comm(N); std::iota(comm.begin(), comm.end(), 0);
        vector<double> Sigma = g.degree;
        int moved = 0;
        std::unordered_map<int, double> w_to;
        for (int u = 0; u < N; ++u) {
            double ku = g.degree[u]; int A = comm[u];
            w_to.clear();
            for (auto &e : g.adj[u]) if (e.first != u) w_to[comm[e.first]] += e.second;
            Sigma[A] -= ku;
            double best_z = -1e300; int target = A; bool found = false;
            for (auto &kv : w_to) {
                int B = kv.first; if (B == A) continue;
                double z = z_statistic(ku, kv.second, Sigma[B], m2);
                if (z > tau) {
                    double gainB = kv.second / m - ku * Sigma[B] / (2 * m * m);
                    if (gainB > 0 && (!found || z > best_z)) { best_z = z; target = B; found = true; }
                }
            }
            if (found) { comm[u] = target; moved++; }
            Sigma[comm[u]] += ku;
        }
        std::cout << tau << "\t" << 2 * (1 - phi(tau)) << "\t" << (double)moved / N << "\n";
    }
    std::cout << "  Single-node z-test only sees neighbours (where an edge exists);\n"
                 "  under the null that edge is itself the rare event, so almost every\n"
                 "  first move clears any fixed tau. Corollary 1's bound does not hold.\n";

    // Same one-sweep-from-singletons measurement, but with the adaptive
    // threshold: does tau*sqrt(2 log m) actually pull the acceptance rate down
    // toward the 2(1-Phi(tau)) bound Corollary 1 claims?
    std::cout << "\n  With adaptive tau = tau*sqrt(2 log m)  [m = " << (long long)(g.m2 / 2)
              << ", scale = " << std::sqrt(2.0 * std::log(std::max(g.m2 / 2.0, 3.0))) << "]:\n";
    std::cout << "  tau\ttau_eff\tbound_2(1-Phi)\tobserved_move_frac\n";
    for (double tau : {1.0, 2.0, 3.0}) {
        int N = g.n(); double m2 = g.m2, m = m2 / 2.0;
        double tau_eff = tau * std::sqrt(2.0 * std::log(std::max(m, 3.0)));
        vector<int> comm(N); std::iota(comm.begin(), comm.end(), 0);
        vector<double> Sigma = g.degree;
        int moved = 0;
        std::unordered_map<int, double> w_to;
        for (int u = 0; u < N; ++u) {
            double ku = g.degree[u]; int A = comm[u];
            w_to.clear();
            for (auto &e : g.adj[u]) if (e.first != u) w_to[comm[e.first]] += e.second;
            Sigma[A] -= ku;
            double best_z = -1e300; int target = A; bool found = false;
            for (auto &kv : w_to) {
                int B = kv.first; if (B == A) continue;
                double z = z_statistic(ku, kv.second, Sigma[B], m2);
                if (z > tau_eff) {
                    double gainB = kv.second / m - ku * Sigma[B] / (2 * m * m);
                    if (gainB > 0 && (!found || z > best_z)) { best_z = z; target = B; found = true; }
                }
            }
            if (found) { comm[u] = target; moved++; }
            Sigma[comm[u]] += ku;
        }
        std::cout << "  " << tau << "\t" << tau_eff << "\t" << 2 * (1 - phi(tau_eff))
                  << "\t" << (double)moved / N << "\n";
    }
}

// --------------------------------------------------------------------------- //
//  Experiment 9: the rolling-window dQ stopping rule on the large graphs.
//
//  Panel A: CNM run to its natural greedy terminus = the reference. Because Q
//           rises monotonically along a positive-dQ greedy run, that terminus is
//           also the argmax-Q partition, so it is the honest "no early stop" row.
//  Panel B: the (W, k, sided) sweep, every row a replay of the SAME trace.
//  Panel C: the rule ported onto the Leiden move loop, to check it is not an
//           artifact of the agglomerative schedule.
// --------------------------------------------------------------------------- //
void experiment_9(const vector<std::string> &specs, long long max_pushes,
                  long long max_merges = 0) {
    std::cout << "\n=== Experiment 9: rolling-window dQ stopping rule (large graphs) ===\n";
    if (specs.empty()) { std::cout << "  (skipped: pass large=edges[:labels])\n"; return; }

    for (const std::string &spec : specs) {
        auto eq = spec.find('=');
        if (eq == std::string::npos) continue;
        std::string name = spec.substr(0, eq);
        std::string rest = spec.substr(eq + 1);
        auto colon = rest.find(':');
        std::string ef = colon == std::string::npos ? rest : rest.substr(0, colon);
        std::string lf = colon == std::string::npos ? "" : rest.substr(colon + 1);

        Dataset d = load_edgelist(ef, lf);
        if (d.n == 0) { std::cout << "  (" << name << ": empty/missing)\n"; continue; }
        Graph g = Graph::from_edges(d.n, d.edges);
        long long M = (long long)d.edges.size();
        std::cerr << "  [" << name << "] N=" << d.n << " M=" << M << " -- CNM...\n";

        g_t.reset();
        CnmTrace tr = cnm_greedy_trace(g, max_pushes, max_merges);
        if (!tr.complete) {
            std::cout << "  (" << name << ": CNM aborted by memory guard after "
                      << tr.merges.size() << " merges -- not reported)\n";
            std::cout.flush();
            continue;
        }
        std::cerr << "    trace: " << tr.merges.size() << " merges in "
                  << (long long)tr.build_ms << " ms\n";

        std::cout << "\n-- " << name << " (N=" << d.n << ", M=" << M
                  << ", CNM merges=" << tr.merges.size() << ", trace_ms="
                  << (long long)tr.build_ms
                  << (tr.truncated ? ", TRUNCATED at merge budget" : "") << ")\n";
        if (tr.truncated)
            std::cout << "   NOTE: agglomeration stopped at a merge budget, so the "
                         "reference row is a PREFIX, not argmax Q. Stopping-rule rows "
                         "that cut before the budget are still exact.\n";
        std::cout << "method\tW\tk\tsided\tstop_at\tof_merges\tfrac\tK\tQ\tF1\tNMI\n";

        auto row = [&](const std::string &mname, int W, double k, const char *sided,
                       size_t cut) {
            vector<int> lab = cnm_labels_at(tr, cut);
            std::cout << mname << "\t" << (W ? std::to_string(W) : "-") << "\t"
                      << (k > 0 ? std::to_string(k).substr(0, 4) : "-") << "\t" << sided
                      << "\t" << cut << "\t" << tr.merges.size() << "\t"
                      << (tr.merges.empty() ? 0.0 : (double)cut / tr.merges.size()) << "\t"
                      << num_communities(lab) << "\t" << modularity(g, lab) << "\t"
                      << community_f1(lab, d.gt) << "\t" << nmi(lab, d.gt) << "\n";
            std::cout.flush();
        };

        // Panel A -- reference: the full greedy run (= argmax Q).
        row(tr.truncated ? "CNM-prefix" : "CNM-full", 0, 0, "-", tr.merges.size());

        // Panel B -- the tuning sweep the user described.
        for (int W : {10, 15, 20})
            for (double k : {1.75, 1.9, 2.0, 2.5}) {
                RollStop fired;
                size_t cut = rollstop_cut(tr, W, k, false, &fired);
                row("roll-1sided", W, k, "1", cut);
                cut = rollstop_cut(tr, W, k, true, &fired);
                row("roll-2sided", W, k, "2", cut);
            }

        // Panel B2 -- scale-free variants. The raw rule cuts these graphs at well
        // under 1% of their merges because dQ ~ 1/m makes the window sd collapse;
        // these two test the same "far from recent history" idea in a form that
        // does not depend on the absolute size of dQ.
        for (int W : {10, 15, 20})
            for (double k : {1.75, 2.0, 2.5})
                row("roll-log", W, k, "1", rollstop_cut_log(tr, W, k));
        for (int W : {15, 20})
            for (double k : {2.0, 2.5})
                for (double mr : {0.25, 0.50})
                    row("roll-rel" + std::to_string((int)(mr * 100)), W, k, "1",
                        rollstop_cut_rel(tr, W, k, mr));

        // Panel B3 -- the three requested tests.
        //   T1  two-sided at k=1 (a deliberately LOOSER threshold than 1.75-2.5)
        //   T2  shorter history + higher k (less averaging, stricter multiplier)
        //   T3  step test: |dQ_prev - dQ_now| > k*sd, sd still from the last W
        for (double k : {1.0})
            for (int W : {10, 15, 20})
                row("T1-2sided-k1", W, k, "2", rollstop_cut(tr, W, k, true));
        for (int W : {5, 7, 10})
            for (double k : {3.0, 4.0, 5.0}) {
                row("T2-short-1sided", W, k, "1", rollstop_cut(tr, W, k, false));
                row("T2-short-2sided", W, k, "2", rollstop_cut(tr, W, k, true));
            }
        for (int W : {10, 15, 20})
            for (double k : {1.0, 2.0, 2.5, 3.0}) {
                row("T3-diff-2sided", W, k, "2", rollstop_cut_diff(tr, W, k, true));
                row("T3-diff-1sided", W, k, "1", rollstop_cut_diff(tr, W, k, false));
            }

        // Panel C -- same rule on the Leiden move loop, W=15 k=2 only.
        for (bool ts : {false, true}) {
            RollStop fired;
            g_t.reset();
            vector<int> lab = rollq(g, 15, 2.0, ts, &fired, 1);
            Timings t = g_t;
            std::cout << (ts ? "leiden-roll-2sided" : "leiden-roll-1sided")
                      << "\t15\t2.00\t" << (ts ? "2" : "1") << "\t"
                      << fired.fired_at << "\t" << fired.seen << "\t-\t"
                      << num_communities(lab) << "\t" << modularity(g, lab) << "\t"
                      << community_f1(lab, d.gt) << "\t" << nmi(lab, d.gt)
                      << "\t(total_ms=" << (long long)t.total_ms << ")\n";
            std::cout.flush();
        }
    }
}

// --------------------------------------------------------------------------- //
//  Experiment 10: same rule on the small labelled networks, where GRAZ's Eq. 6
//  over-segmented badly (karate 4 -> 16 communities). A stopping rule should not
//  be able to do that, since it can only truncate a greedy run, never split.
// --------------------------------------------------------------------------- //
void experiment_10(int argc, char **argv) {
    std::cout << "\n=== Experiment 10: rolling-window rule on small networks ===\n";
    std::cout << "dataset\tmethod\tW\tk\tsided\tstop_at\tof_merges\tK\tQ\tF1\tNMI\n";
    for (int i = 1; i < argc; ++i) {
        std::string spec = argv[i];
        if (spec.rfind("large=", 0) == 0) continue;
        auto eq = spec.find('=');
        if (eq == std::string::npos) continue;
        std::string name = spec.substr(0, eq), rest = spec.substr(eq + 1);
        auto colon = rest.find(':');
        std::string ef = colon == std::string::npos ? rest : rest.substr(0, colon);
        std::string lf = colon == std::string::npos ? "" : rest.substr(colon + 1);
        Dataset d = load_edgelist(ef, lf);
        if (d.n == 0) continue;
        Graph g = Graph::from_edges(d.n, d.edges);
        CnmTrace tr = cnm_greedy_trace(g);
        dump_trace_diagnostics(name, tr);

        auto row = [&](const std::string &mn, int W, double k, const char *sd,
                       const vector<int> &lab, long long cut) {
            std::cout << name << "\t" << mn << "\t" << (W ? std::to_string(W) : "-")
                      << "\t" << (k > 0 ? std::to_string(k).substr(0, 4) : "-") << "\t" << sd
                      << "\t" << cut << "\t" << tr.merges.size() << "\t"
                      << num_communities(lab) << "\t" << modularity(g, lab) << "\t"
                      << community_f1(lab, d.gt) << "\t" << nmi(lab, d.gt) << "\n";
        };
        row("Louvain", 0, 0, "-", louvain(g, 1), -1);
        row("GRAZ t=2", 0, 0, "-", graz(g, 2.0, 1), -1);
        row("CNM-full", 0, 0, "-", cnm_labels_at(tr, tr.merges.size()),
            (long long)tr.merges.size());
        for (int W : {10, 15, 20})
            for (double k : {1.75, 2.0, 2.5}) {
                size_t c1 = rollstop_cut(tr, W, k, false);
                row("roll-1sided", W, k, "1", cnm_labels_at(tr, c1), (long long)c1);
                size_t c2 = rollstop_cut(tr, W, k, true);
                row("roll-2sided", W, k, "2", cnm_labels_at(tr, c2), (long long)c2);
                size_t c3 = rollstop_cut_log(tr, W, k);
                row("roll-log", W, k, "1", cnm_labels_at(tr, c3), (long long)c3);
            }
        for (int W : {15, 20})
            for (double k : {2.0, 2.5})
                for (double mr : {0.25, 0.50}) {
                    size_t c4 = rollstop_cut_rel(tr, W, k, mr);
                    row("roll-rel" + std::to_string((int)(mr * 100)), W, k, "1",
                        cnm_labels_at(tr, c4), (long long)c4);
                }
        // The three requested tests, same panel layout as experiment 9.
        for (int W : {10, 15, 20}) {
            size_t t1 = rollstop_cut(tr, W, 1.0, true);
            row("T1-2sided-k1", W, 1.0, "2", cnm_labels_at(tr, t1), (long long)t1);
        }
        for (int W : {5, 7, 10})
            for (double k : {3.0, 4.0, 5.0}) {
                size_t a = rollstop_cut(tr, W, k, false);
                row("T2-short-1sided", W, k, "1", cnm_labels_at(tr, a), (long long)a);
                size_t b = rollstop_cut(tr, W, k, true);
                row("T2-short-2sided", W, k, "2", cnm_labels_at(tr, b), (long long)b);
            }
        for (int W : {10, 15, 20})
            for (double k : {1.0, 2.0, 2.5, 3.0}) {
                size_t a = rollstop_cut_diff(tr, W, k, true);
                row("T3-diff-2sided", W, k, "2", cnm_labels_at(tr, a), (long long)a);
                size_t b = rollstop_cut_diff(tr, W, k, false);
                row("T3-diff-1sided", W, k, "1", cnm_labels_at(tr, b), (long long)b);
            }
        std::cout.flush();
    }
}

// --------------------------------------------------------------------------- //
//  Experiment 11: the Corollary 1 test, applied to the stopping rule.
//
//  This is the experiment that decides whether the rule is a real statistical
//  significance criterion or just another knob. An ER graph has NO communities,
//  so the correct K is 0 and every community reported is spurious. GRAZ's Eq. 6
//  fails here badly and gets WORSE as tau rises (58 -> 103 -> 277 at N=1000).
//  A stopping rule that keys on the empirical dQ distribution should instead cut
//  the greedy run early on a structureless graph and report FEWER communities
//  than the unstopped baseline -- and it must get monotonically stricter as k
//  falls, which is the direction-of-effect check Eq. 6 fails.
// --------------------------------------------------------------------------- //
void experiment_11() {
    std::cout << "\n=== Experiment 11: rolling-window rule on ER nulls (Corollary 1) ===\n";
    std::cout << "  (ER has NO communities: true K = 0, lower is better. K counts size>=3.)\n";
    std::cout << "N\tm\tmethod\tW\tk\tsided\tK_ge3\tK_all\tQ\tstop_frac\n";
    int n_seeds = 3;
    for (int N : {1000, 2000, 5000, 10000}) {
        struct Acc { vector<double> k3, kall, q, sf; };
        std::map<std::string, Acc> acc;
        vector<double> ms;
        std::map<std::string, std::pair<int, double>> cfg;  // label -> (W,k)
        std::map<std::string, const char *> side;

        for (int s = 0; s < n_seeds; ++s) {
            Dataset d = erdos_renyi(N, 10.0, s);
            Graph g = Graph::from_edges(d.n, d.edges);
            ms.push_back((double)d.edges.size());
            CnmTrace tr = cnm_greedy_trace(g);
            auto big = [](const vector<int> &c) {
                auto sz = community_sizes(c); int k = 0;
                for (int x : sz) if (x >= 3) k++; return (double)k; };
            auto rec = [&](const std::string &lbl, int W, double k, const char *sd,
                           const vector<int> &lab, double frac) {
                acc[lbl].k3.push_back(big(lab));
                acc[lbl].kall.push_back(num_communities(lab));
                acc[lbl].q.push_back(modularity(g, lab));
                acc[lbl].sf.push_back(frac);
                cfg[lbl] = {W, k}; side[lbl] = sd;
            };
            rec("Louvain", 0, 0, "-", louvain(g, s), 1.0);
            rec("GRAZ t=2", 0, 0, "-", graz(g, 2.0, s), 1.0);
            rec("CNM-full", 0, 0, "-", cnm_labels_at(tr, tr.merges.size()), 1.0);
            for (int W : {10, 15, 20})
                for (double k : {1.75, 2.0, 2.5}) {
                    char b1[64], b2[64];
                    snprintf(b1, sizeof b1, "roll-1sided W=%d k=%.2f", W, k);
                    snprintf(b2, sizeof b2, "roll-2sided W=%d k=%.2f", W, k);
                    double den = tr.merges.empty() ? 1.0 : (double)tr.merges.size();
                    size_t c1 = rollstop_cut(tr, W, k, false);
                    rec(b1, W, k, "1", cnm_labels_at(tr, c1), (double)c1 / den);
                    size_t c2 = rollstop_cut(tr, W, k, true);
                    rec(b2, W, k, "2", cnm_labels_at(tr, c2), (double)c2 / den);
                    // A scale-free variant must ALSO still reject ER nulls, or it
                    // has only traded one failure mode for the other.
                    char b3[64];
                    snprintf(b3, sizeof b3, "roll-log   W=%d k=%.2f", W, k);
                    size_t c3 = rollstop_cut_log(tr, W, k);
                    rec(b3, W, k, "1", cnm_labels_at(tr, c3), (double)c3 / den);
                    // roll-rel is the variant that recovered full quality on the
                    // real graphs, so it is the one that most needs checking here.
                    char b4[64];
                    snprintf(b4, sizeof b4, "roll-rel50 W=%d k=%.2f", W, k);
                    size_t c4 = rollstop_cut_rel(tr, W, k, 0.50);
                    rec(b4, W, k, "1", cnm_labels_at(tr, c4), (double)c4 / den);
                    char b5[64];
                    snprintf(b5, sizeof b5, "roll-rel25 W=%d k=%.2f", W, k);
                    size_t c5 = rollstop_cut_rel(tr, W, k, 0.25);
                    rec(b5, W, k, "1", cnm_labels_at(tr, c5), (double)c5 / den);
                }
            // The three requested tests. A variant that improves the real graphs
            // must still reject ER nulls, or it has only moved the failure.
            {
                double den = tr.merges.empty() ? 1.0 : (double)tr.merges.size();
                for (int W : {10, 15, 20}) {
                    char b[64];
                    snprintf(b, sizeof b, "T1-2sided-k1 W=%d k=1.00", W);
                    size_t c = rollstop_cut(tr, W, 1.0, true);
                    rec(b, W, 1.0, "2", cnm_labels_at(tr, c), (double)c / den);
                }
                for (int W : {5, 7, 10})
                    for (double k : {3.0, 4.0, 5.0}) {
                        char ba[72], bb[72];
                        snprintf(ba, sizeof ba, "T2-short-1sided W=%d k=%.2f", W, k);
                        snprintf(bb, sizeof bb, "T2-short-2sided W=%d k=%.2f", W, k);
                        size_t a = rollstop_cut(tr, W, k, false);
                        rec(ba, W, k, "1", cnm_labels_at(tr, a), (double)a / den);
                        size_t b2 = rollstop_cut(tr, W, k, true);
                        rec(bb, W, k, "2", cnm_labels_at(tr, b2), (double)b2 / den);
                    }
                for (int W : {10, 15, 20})
                    for (double k : {1.0, 2.0, 2.5, 3.0}) {
                        char ba[72], bb[72];
                        snprintf(ba, sizeof ba, "T3-diff-2sided W=%d k=%.2f", W, k);
                        snprintf(bb, sizeof bb, "T3-diff-1sided W=%d k=%.2f", W, k);
                        size_t a = rollstop_cut_diff(tr, W, k, true);
                        rec(ba, W, k, "2", cnm_labels_at(tr, a), (double)a / den);
                        size_t b2 = rollstop_cut_diff(tr, W, k, false);
                        rec(bb, W, k, "1", cnm_labels_at(tr, b2), (double)b2 / den);
                    }
            }
        }
        for (auto &kv : acc) {
            auto &a = kv.second;
            int W = cfg[kv.first].first; double k = cfg[kv.first].second;
            std::cout << N << "\t" << (int)mean(ms) << "\t" << kv.first << "\t"
                      << (W ? std::to_string(W) : "-") << "\t"
                      << (k > 0 ? std::to_string(k).substr(0, 4) : "-") << "\t"
                      << side[kv.first] << "\t" << mean(a.k3) << "\t" << mean(a.kall)
                      << "\t" << mean(a.q) << "\t" << mean(a.sf) << "\n";
        }
        std::cout.flush();
    }
}

// --------------------------------------------------------------------------- //
// Experiment 13: David-review full stopping-rule report.
//
// This is intentionally a tab-separated, rectangular report so it can be
// redirected into a raw .tsv and converted verbatim to the plain Markdown
// tables in RESULTS.md.  A run pays for one CNM trace per graph; every rule is
// then a replay of that immutable trace.
// --------------------------------------------------------------------------- //
struct StopRuleSpec {
    std::string method, modifier;
    int W = 0;
    double parameter = 0.0;       // k for T1--T4, h for T5, c for CM
    bool two_sided = false;
    int type = 0;                 // 0 raw, 1 diff, 2 detrend, 3 cusum, 4 SR, 5 BM, 6 CM
    StopModifiers mod;
};

vector<StopRuleSpec> full_rule_specs(double h) {
    vector<StopRuleSpec> out;
    const std::vector<std::pair<std::string, StopModifiers>> modifiers = {
        {"none", {}},
        {"defer", {true, Dispersion::SD, 1e-6}},
        {"mad", {false, Dispersion::MAD, 1e-6}},
        {"defer+mad", {true, Dispersion::MAD, 1e-6}}
    };
    // Carried-over best representatives.  The exhaustive legacy grids remain
    // available in experiments 9--11; this focused comparison avoids making a
    // 100-row ER table turn into an unreadable 1,000-row tuning table.
    for (const auto &mm : modifiers) {
        out.push_back({"T1-2sided-k1", mm.first, 10, 1.0, true, 0, mm.second});
        out.push_back({"T2-short-1sided", mm.first, 5, 5.0, false, 0, mm.second});
        out.push_back({"T3-diff-1sided", mm.first, 15, 3.0, false, 1, mm.second});
        for (int W : {10, 15, 20}) for (double k : {1.0, 2.0, 3.0}) {
            out.push_back({"T4-detrend-1sided", mm.first, W, k, false, 2, mm.second});
            out.push_back({"T4-detrend-2sided", mm.first, W, k, true, 2, mm.second});
        }
        out.push_back({"T5-cusum", mm.first, 15, h, false, 3, mm.second});
        out.push_back({"T5-sr", mm.first, 15, h, false, 4, mm.second});
    }
    out.push_back({"BM-native-stop", "-", 0, 2.0, false, 5, {}});
    for (double c : {1.0, 5.0, 10.0, 50.0})
        out.push_back({"cm-baseline", "-", 0, c, false, 6, {}});
    return out;
}

size_t apply_stop_rule(const CnmTrace &tr, double m2, const StopRuleSpec &s) {
    switch (s.type) {
    case 0: return rollstop_cut_mod(tr, s.W, s.parameter, s.two_sided, "raw", s.mod);
    case 1: return rollstop_cut_mod(tr, s.W, s.parameter, s.two_sided, "diff", s.mod);
    case 2: return rollstop_cut_detrend(tr, s.W, s.parameter, s.two_sided, s.mod);
    case 3: return rollstop_cut_cusum(tr, s.W, s.parameter, 1.0, s.mod, false);
    case 4: return rollstop_cut_cusum(tr, s.W, s.parameter, 1.0, s.mod, true);
    case 5: return bm_native_stop_cut(tr, s.parameter);
    case 6: return cm_baseline_cut(tr, m2, s.parameter);
    }
    return tr.merges.size();
}

struct RuleAggregate {
    double k_all = 0, k_ge3 = 0, q = 0, f1 = 0, nmi_v = 0, nonsingle = 0;
    double stop = 0, total = 0, edges = 0;
    int n = 0;
};

void add_rule_metrics(RuleAggregate &a, const Graph &g, const Dataset &d,
                      const CnmTrace &tr, size_t cut) {
    vector<int> lab = cnm_labels_at(tr, cut);
    int k3 = 0;
    for (int s : community_sizes(lab)) if (s >= 3) ++k3;
    a.k_all += num_communities(lab); a.k_ge3 += k3;
    a.q += modularity(g, lab); a.f1 += community_f1(lab, d.gt); a.nmi_v += nmi(lab, d.gt);
    a.nonsingle += frac_nonsingleton(lab); a.stop += cut;
    a.total += tr.merges.size(); a.edges += d.edges.size(); ++a.n;
}

double calibrate_cusum_h(const vector<CnmTrace> &traces) {
    // An empirical in-control calibration: select the smallest predeclared h
    // with mean ARL >= 90% of the actual available CNM trace length.  No graph
    // may choose a bespoke threshold after looking at its ground truth.
    for (double h : {2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0}) {
        double arl = 0.0, total = 0.0;
        for (const CnmTrace &tr : traces) {
            StopModifiers plain;
            arl += rollstop_cut_cusum(tr, 15, h, 1.0, plain, false);
            total += tr.merges.size();
        }
        if (total > 0.0 && arl / total >= 0.90) return h;
    }
    return 200.0;
}

void print_rule_aggregate_row(const std::string &dataset, const StopRuleSpec &s,
                              const RuleAggregate &a, bool labelled) {
    double d = std::max(1, a.n);
    std::cout << dataset << "\t" << s.method << "\t" << s.modifier << "\t" << s.W << "\t"
              << s.parameter << "\t" << (s.two_sided ? "2" : "1") << "\t"
              << a.stop / d << "\t" << a.total / d << "\t"
              << (a.total > 0 ? a.stop / a.total : 0.0) << "\t"
              << a.k_all / d << "\t" << a.k_ge3 / d << "\t" << a.q / d << "\t";
    if (labelled) std::cout << a.f1 / d << "\t" << a.nmi_v / d;
    else std::cout << "NA\tNA";
    std::cout << "\t" << a.nonsingle / d << "\n";
}

void experiment_13(int argc, char **argv, const vector<std::string> &large_specs,
                   long long max_pushes, long long max_merges, int n_seeds) {
    std::cout << "\n=== Experiment 13: full David-review stopping-rule comparison ===\n";
    std::cout << "# T4 fits OLS to the preceding log-dQ window and compares the current "
                 "trend-extrapolation residual with the fitted-residual scale. T5 applies "
                 "CUSUM/SR to the same standardized residual stream; h is calibrated once "
                 "on ER nulls, rather than swept per real graph.\n";
    std::cout << "# `mad` uses 1.4826*MAD and suppresses a near-flat window when its "
                 "relative scale is below 1e-6; `defer` begins testing only after no "
                 "mergeable singleton community remains. BM-native-stop is separate: it "
                 "tests Eq. 6 on CNM's selected community pair with adaptive tau.\n";
    std::cout << "dataset\tmethod\tmodifier\tW\tparameter\tsided\tstop_at\tof_merges\t"
                 "stop_frac\tK_all\tK_ge3\tQ\tF1\tNMI\tfrac_nonsingleton\n";

    // First pass creates exactly the requested 10-seed ER calibration corpus.
    struct ErCase { int N; Dataset d; CnmTrace tr; };
    vector<ErCase> er;
    vector<CnmTrace> null_traces;
    for (int N : {1000, 2000, 5000, 10000}) for (int seed = 0; seed < n_seeds; ++seed) {
        Dataset d = erdos_renyi(N, 10.0, seed);
        Graph g = Graph::from_edges(d.n, d.edges);
        CnmTrace tr = cnm_greedy_trace(g, max_pushes);
        if (!tr.complete) { std::cerr << "ER trace aborted; omitting N=" << N << " seed=" << seed << "\n"; continue; }
        null_traces.push_back(tr);
        er.push_back({N, std::move(d), std::move(tr)});
    }
    double h = calibrate_cusum_h(null_traces);
    std::cout << "# CUSUM calibration (delta=1, W=15): selected h=" << h
              << " by the >=90% mean-ARL rule.\n";
    for (double candidate_h : {2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0}) {
        double arl = 0.0, total = 0.0;
        for (const CnmTrace &tr : null_traces) {
            arl += rollstop_cut_cusum(tr, 15, candidate_h, 1.0, {}, false);
            total += tr.merges.size();
        }
        std::cout << "CALIBRATION\tT5-cusum\t-\t15\t" << candidate_h << "\t1\t"
                  << (null_traces.empty() ? 0.0 : arl / null_traces.size()) << "\t"
                  << (null_traces.empty() ? 0.0 : total / null_traces.size()) << "\t"
                  << (total > 0 ? arl / total : 0.0) << "\tNA\tNA\tNA\tNA\tNA\tNA\n";
    }
    vector<StopRuleSpec> specs = full_rule_specs(h);
    for (int N : {1000, 2000, 5000, 10000}) {
        std::map<std::string, RuleAggregate> acc;
        for (const ErCase &c : er) if (c.N == N) {
            Graph g = Graph::from_edges(c.d.n, c.d.edges);
            for (const StopRuleSpec &s : specs) {
                size_t cut = apply_stop_rule(c.tr, g.m2, s);
                add_rule_metrics(acc[s.method + "|" + s.modifier + "|" + std::to_string(s.W) + "|" +
                                     std::to_string(s.parameter) + "|" + std::to_string(s.two_sided)],
                                 g, c.d, c.tr, cut);
            }
        }
        for (const StopRuleSpec &s : specs) {
            std::string key = s.method + "|" + s.modifier + "|" + std::to_string(s.W) + "|" +
                              std::to_string(s.parameter) + "|" + std::to_string(s.two_sided);
            print_rule_aggregate_row("ER-" + std::to_string(N), s, acc[key], false);
        }
    }

    auto run_real = [&](const std::string &name, const Dataset &d, long long budget) {
        Graph g = Graph::from_edges(d.n, d.edges);
        CnmTrace tr = cnm_greedy_trace(g, max_pushes, budget);
        if (!tr.complete) { std::cout << "# " << name << " CNM trace aborted by memory guard.\n"; return; }
        std::cout << "# " << name << " CNM-" << (tr.truncated ? "prefix" : "full")
                  << " merges=" << tr.merges.size() << " budget=" << budget
                  << " singleton_phase_end=" << tr.singleton_phase_end << "\n";
        if ((name.find("amazon") != std::string::npos || name.find("dblp") != std::string::npos) &&
            tr.dq.size() >= 2458) {
            bool identical = true;
            for (size_t t = 2438; t < 2458; ++t) if (tr.dq[t] != tr.dq[t + 1]) identical = false;
            std::cout << "# plateau-window " << name << " indices=2438..2458 bit_identical_consecutive="
                      << (identical ? "yes" : "no") << "\n";
            std::streamsize old_precision = std::cout.precision();
            std::ios::fmtflags old_flags = std::cout.flags();
            std::cout << std::scientific << std::setprecision(17);
            for (size_t t = 2438; t <= 2458; ++t)
                std::cout << "TRACE\t" << name << "\t" << t << "\t" << tr.dq[t] << "\n";
            std::cout.flags(old_flags); std::cout.precision(old_precision);
        }
        StopRuleSpec ref{tr.truncated ? "CNM-prefix" : "CNM-full", "-", 0, 0, false, 6, {}};
        RuleAggregate ra; add_rule_metrics(ra, g, d, tr, tr.merges.size());
        print_rule_aggregate_row(name, ref, ra,
            std::any_of(d.gt.begin(), d.gt.end(), [](int x) { return x >= 0; }));
        for (const StopRuleSpec &s : specs) {
            RuleAggregate a; add_rule_metrics(a, g, d, tr, apply_stop_rule(tr, g.m2, s));
            print_rule_aggregate_row(name, s, a,
                std::any_of(d.gt.begin(), d.gt.end(), [](int x) { return x >= 0; }));
        }
    };
    // Small graph specs use the existing name=edges:labels convention.
    for (int i = 1; i < argc; ++i) {
        std::string spec = argv[i];
        if (spec.rfind("--", 0) == 0 || spec.rfind("large=", 0) == 0) continue;
        auto eq = spec.find('='); if (eq == std::string::npos) continue;
        auto colon = spec.find(':', eq + 1);
        Dataset d = load_edgelist(spec.substr(eq + 1, colon == std::string::npos ? std::string::npos : colon - eq - 1),
                                  colon == std::string::npos ? "" : spec.substr(colon + 1));
        if (d.n) run_real(spec.substr(0, eq), d, 0);
    }
    for (const std::string &spec : large_specs) {
        auto eq = spec.find('='); if (eq == std::string::npos) continue;
        auto colon = spec.find(':', eq + 1);
        Dataset d = load_edgelist(spec.substr(eq + 1, colon == std::string::npos ? std::string::npos : colon - eq - 1),
                                  colon == std::string::npos ? "" : spec.substr(colon + 1));
        if (d.n) run_real(spec.substr(0, eq), d, max_merges);
    }
}

// --------------------------------------------------------------------------- //
// Fix 5: post-hoc multi-cut selection.  Candidate cuts are every ceil(T/50)-th
// merge plus the terminal cut, so at most 51 partitions are evaluated.  This
// downsampling is stated in the report rather than hiding a search budget.
// --------------------------------------------------------------------------- //
using BlockKey = unsigned long long;
BlockKey block_key(int a, int b) {
    if (a > b) std::swap(a, b);
    return (static_cast<BlockKey>(static_cast<unsigned>(a)) << 32) |
           static_cast<unsigned>(b);
}

double entropy_bits(double p) {
    if (p <= 0.0 || p >= 1.0) return 0.0;
    return -p * std::log2(p) - (1.0 - p) * std::log2(1.0 - p);
}

struct BlockScores { double heldout_ll = -1e300, mdl = 1e300; };

BlockScores score_partition_blocks(const Dataset &train, const Dataset &test,
                                  const vector<int> &lab) {
    BlockScores out;
    std::unordered_map<int, long long> sizes;
    for (int x : lab) ++sizes[x];
    std::unordered_map<BlockKey, long long> etrain, etest;
    for (auto e : train.edges) ++etrain[block_key(lab[e.first], lab[e.second])];
    for (auto e : test.edges) ++etest[block_key(lab[e.first], lab[e.second])];
    std::unordered_set<BlockKey> keys;
    for (auto &x : etrain) keys.insert(x.first);
    for (auto &x : etest) keys.insert(x.first);
    // Include empty blocks implicitly only in MDL through their zero entropy.
    double graph_bits = 0.0;
    for (BlockKey key : keys) {
        int a = (int)(key >> 32), b = (int)(key & 0xffffffffu);
        double possible = a == b ? 0.5 * sizes[a] * (sizes[a] - 1.0)
                                 : (double)sizes[a] * sizes[b];
        if (possible <= 0.0) continue;
        // Jeffreys smoothing makes held-out likelihood finite for unseen blocks.
        double p = (etrain[key] + 0.5) / (possible + 1.0);
        double held_nonedges = std::max(0.0, possible - etrain[key] - etest[key]);
        out.heldout_ll += etest[key] * std::log(p) + held_nonedges * std::log(1.0 - p);
        double pall = (etrain[key] + etest[key] + 0.5) / (possible + 1.0);
        graph_bits += possible * entropy_bits(pall);
    }
    double assignment_bits = 0.0, n = lab.size();
    for (auto &kv : sizes) {
        double p = kv.second / n;
        assignment_bits -= kv.second * std::log2(p);
    }
    out.mdl = assignment_bits + graph_bits;
    return out;
}

void experiment_14_multicut(int argc, char **argv, const vector<std::string> &large_specs,
                            long long max_pushes, int stability_reps) {
    std::cout << "\n=== Experiment 14: multicut-report ===\n";
    std::cout << "# Candidate cuts: every ceil(T/50)-th CNM merge plus T (at most 51). "
                 "Held-out likelihood uses a fixed 90/10 edge split with Jeffreys-smoothed "
                 "block probabilities; MDL is assignment entropy plus block Bernoulli code length.\n";
    std::cout << "graph\tcandidates\tlikelihood_cut\tmdl_cut\tstability_cut\t"
                 "CNM_argmaxQ_cut\tBM_native_cut\tstability_reps\n";
    auto run = [&](const std::string &name, const Dataset &d) {
        EdgeSplit split = split_edges(d, 0.90, 0x51EC7u);
        CnmTrace tr = cnm_greedy_trace(split.discovery, max_pushes);
        if (!tr.complete) { std::cout << name << "\tABORTED\n"; return; }
        size_t stride = std::max<size_t>(1, (tr.merges.size() + 49) / 50);
        vector<size_t> cuts;
        for (size_t c = 0; c < tr.merges.size(); c += stride) cuts.push_back(c);
        if (cuts.empty() || cuts.back() != tr.merges.size()) cuts.push_back(tr.merges.size());
        vector<vector<int>> labs; labs.reserve(cuts.size());
        vector<BlockScores> scores(cuts.size());
        for (size_t cut : cuts) {
            labs.push_back(cnm_labels_at(tr, cut));
        }
        // score_partition_blocks consumes edge lists; reconstruct them from the
        // split graphs once, preserving weights as unweighted input edges.
        Dataset train{d.n, {}, {}}, test{d.n, {}, {}};
        for (int u = 0; u < split.discovery.n(); ++u)
            for (auto e : split.discovery.adj[u]) if (u < e.first)
                train.edges.push_back({u, e.first});
        for (int u = 0; u < split.validation.n(); ++u)
            for (auto e : split.validation.adj[u]) if (u < e.first)
                test.edges.push_back({u, e.first});
        for (size_t i = 0; i < cuts.size(); ++i) scores[i] = score_partition_blocks(train, test, labs[i]);
        size_t best_ll = 0, best_mdl = 0;
        for (size_t i = 1; i < cuts.size(); ++i) {
            if (scores[i].heldout_ll > scores[best_ll].heldout_ll) best_ll = i;
            if (scores[i].mdl < scores[best_mdl].mdl) best_mdl = i;
        }
        vector<double> stability(cuts.size(), 0.0);
        std::mt19937 rng(0x57AB1Eu);
        for (int rep = 0; rep < stability_reps; ++rep) {
            std::bernoulli_distribution keep(0.95);
            Dataset sub{d.n, {}, d.gt};
            for (auto e : d.edges) if (keep(rng)) sub.edges.push_back(e);
            CnmTrace st = cnm_greedy_trace(Graph::from_edges(sub.n, sub.edges), max_pushes);
            if (!st.complete) continue;
            for (size_t i = 0; i < cuts.size(); ++i)
                stability[i] += nmi(labs[i], cnm_labels_at(st, std::min(cuts[i], st.merges.size())));
        }
        size_t best_stability = 0;
        for (size_t i = 1; i < cuts.size(); ++i)
            if (stability[i] > stability[best_stability]) best_stability = i;
        std::cout << name << "\t" << cuts.size() << "\t" << cuts[best_ll] << "\t"
                  << cuts[best_mdl] << "\t" << cuts[best_stability] << "\t"
                  << tr.merges.size() << "\t" << bm_native_stop_cut(tr, 2.0) << "\t"
                  << stability_reps << "\n";
    };
    for (int i = 1; i < argc; ++i) {
        std::string spec = argv[i];
        if (spec.rfind("--", 0) == 0 || spec.rfind("large=", 0) == 0) continue;
        auto eq = spec.find('='); if (eq == std::string::npos) continue;
        auto colon = spec.find(':', eq + 1);
        Dataset d = load_edgelist(spec.substr(eq + 1, colon == std::string::npos ? std::string::npos : colon - eq - 1),
                                  colon == std::string::npos ? "" : spec.substr(colon + 1));
        if (d.n) run(spec.substr(0, eq), d);
    }
    for (const std::string &spec : large_specs) {
        auto eq = spec.find('='); if (eq == std::string::npos) continue;
        auto colon = spec.find(':', eq + 1);
        Dataset d = load_edgelist(spec.substr(eq + 1, colon == std::string::npos ? std::string::npos : colon - eq - 1),
                                  colon == std::string::npos ? "" : spec.substr(colon + 1));
        if (d.n) run(spec.substr(0, eq), d);
    }
}

// --------------------------------------------------------------------------- //
//  Experiment 12: BM-PATH, a calibrated Bader--McCloskey merge criterion.
//
//  This is deliberately opt-in (`--bm-path`).  A run uses B configuration-model
//  paths for every graph, so the small/medium benchmark is the development
//  target; it must earn a large-network run rather than silently charging one.
// --------------------------------------------------------------------------- //
void experiment_12(int argc, char **argv, int B, double alpha) {
    std::cout << "\n=== Experiment 12: BM-PATH configuration-bootstrap residual test ===\n";
    std::cout << "graph\tN\tM\tmethod\tB\talpha\tmerges\tK\tQ\tF1\tNMI\tlast_S\tp_next\tms\n";

    auto run = [&](const std::string &name, const Dataset &d) {
        Graph g = Graph::from_edges(d.n, d.edges);
        bool labelled = std::any_of(d.gt.begin(), d.gt.end(), [](int x) { return x >= 0; });
        auto print = [&](const std::string &method, const vector<int> &lab, size_t merges,
                         double last_s, double pnext, double ms) {
            std::cout << name << "\t" << d.n << "\t" << d.edges.size() << "\t" << method
                      << "\t" << B << "\t" << alpha << "\t" << merges << "\t"
                      << num_communities(lab) << "\t" << modularity(g, lab) << "\t";
            if (labelled) std::cout << community_f1(lab, d.gt) << "\t" << nmi(lab, d.gt);
            else std::cout << "NA\tNA";
            std::cout << "\t" << last_s << "\t";
            if (pnext < 0.0) std::cout << "NA"; else std::cout << pnext;
            std::cout << "\t" << ms << "\n";
        };

        auto c0 = Clock::now();
        CnmTrace cnm = cnm_greedy_trace(g);
        print("CNM-full", cnm_labels_at(cnm, cnm.merges.size()), cnm.merges.size(), 0.0, -1.0,
              ms_since(c0));

        auto b0 = Clock::now();
        BmTrace raw = bm_greedy_trace(g);
        double raw_s = raw.z.empty() ? 0.0 : raw.z.back();
        print("BM-residual-full", bm_labels_at(raw, raw.merges.size()), raw.merges.size(), raw_s,
              -1.0, ms_since(b0));

        BmBootstrapResult bm = bm_path_bootstrap(g, B, alpha, 0xBADC0DEu + (unsigned)d.n);
        double accepted_s = bm.cut == 0 ? 0.0 : bm.observed.z[bm.cut - 1];
        double pnext = bm.cut < bm.p_stage.size() ? bm.p_stage[bm.cut] : -1.0;
        print("BM-PATH", bm_labels_at(bm.observed, bm.cut), bm.cut, accepted_s, pnext, bm.total_ms);

        BmCvResult cv = bm_cv_select(d, B, alpha, 0xC0FFEEu + (unsigned)d.n);
        if (cv.significant) {
            print("BM-CV", bm_labels_at(cv.discovery_tree, cv.cut), cv.cut, cv.z_best,
                  cv.p_global, cv.total_ms);
        } else {
            // `K=0` means "no statistically supported partition", not a graph
            // labelling with zero nodes.  Avoid calling partition metrics for it.
            std::cout << name << "\t" << d.n << "\t" << d.edges.size() << "\tBM-CV-none\t"
                      << B << "\t" << alpha << "\t0\t0\tNA\tNA\tNA\t"
                      << cv.z_best << "\t" << cv.p_global << "\t" << cv.total_ms << "\n";
        }
    };

    // The resolution-limit instance and the degree-corrected hierarchical SBM
    // have known planted labels.  ER has no planted partition: a calibrated
    // method should retain zero statistically significant merges there.
    run("ring40-c5", ring_of_cliques(40, 5));
    vector<int> gt_super;
    run("hierarchical-sbm", hierarchical_sbm(17, gt_super));
    run("ER-250-d10", erdos_renyi(250, 10.0, 17));

    // Optional prepared labelled data: e.g.
    //   --bm-path karate=prepared/karate.edges:prepared/karate.labels
    for (int i = 1; i < argc; ++i) {
        std::string spec = argv[i];
        if (spec.rfind("--", 0) == 0 || spec.rfind("large=", 0) == 0) continue;
        auto eq = spec.find('=');
        if (eq == std::string::npos) continue;
        std::string name = spec.substr(0, eq), rest = spec.substr(eq + 1);
        auto colon = rest.find(':');
        Dataset d = load_edgelist(rest.substr(0, colon),
                                  colon == std::string::npos ? "" : rest.substr(colon + 1));
        if (d.n > 0) run(name, d);
    }
}

int main(int argc, char **argv) {
    std::cout.setf(std::ios::fixed);
    std::cout.precision(3);

    // Split args: "large=..." specs go to experiment 7, "--only-large" skips the
    // synthetic suite so the big benchmark can be re-run on its own.
    vector<std::string> large;
    bool only_large = false;
    // The rolling-stop experiments (9/10/11) are opt-in so the CNM agglomeration
    // is only paid for when asked. --only-roll runs 9 alone on the large specs.
    bool only_roll = false, want_roll_small = false, want_roll_er = false, want_bm_path = false;
    bool want_full_stop_report = false;
    bool want_multicut_report = false;
    int report_seeds = 10;  // David-review ER calibration requirement
    int stability_reps = 20;
    int bm_bootstraps = 99;
    double bm_alpha = 0.05;
    // CNM heap ceiling in entries. Ent is 24 B, so 20M entries is about 480 MB.
    // Crossing it triggers stale-entry compaction (result-preserving); only a heap
    // still oversized when fully live aborts the run.
    long long max_pushes = 20000000LL;
    // Optional merge budget. 0 = run to the greedy terminus. Used for graphs where
    // a single hub community makes the tail of the agglomeration cost days: the
    // recorded prefix is identical either way, so a budget far beyond where the
    // stopping rule fires costs nothing but the CNM-full reference row.
    long long max_merges = 0;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--only-large") { only_large = true; continue; }
        if (a == "--only-roll") { only_roll = true; continue; }
        if (a == "--roll-small") { want_roll_small = true; continue; }
        if (a == "--roll-er") { want_roll_er = true; continue; }
        if (a == "--bm-path") { want_bm_path = true; continue; }
        if (a == "--full-stop-report") { want_full_stop_report = true; continue; }
        if (a == "--multicut-report") { want_multicut_report = true; continue; }
        if (a.rfind("--report-seeds=", 0) == 0) {
            report_seeds = std::max(1, atoi(a.substr(15).c_str()));
            continue;
        }
        if (a.rfind("--stability-reps=", 0) == 0) {
            stability_reps = std::max(1, atoi(a.substr(17).c_str()));
            continue;
        }
        if (a.rfind("--bm-bootstrap=", 0) == 0) {
            bm_bootstraps = std::max(19, atoi(a.substr(15).c_str()));
            continue;
        }
        if (a.rfind("--bm-alpha=", 0) == 0) {
            bm_alpha = atof(a.substr(11).c_str());
            continue;
        }
        if (a.rfind("--max-pushes=", 0) == 0) {
            max_pushes = atoll(a.substr(13).c_str());
            continue;
        }
        if (a.rfind("--max-merges=", 0) == 0) {
            max_merges = atoll(a.substr(13).c_str());
            continue;
        }
        if (a.rfind("--methods=", 0) == 0) {
            std::stringstream ss(a.substr(10));
            std::string tok;
            while (std::getline(ss, tok, ',')) if (!tok.empty()) g_methods.insert(tok);
            continue;
        }
        if (a.rfind("large=", 0) == 0) large.push_back(a.substr(6));
    }

    if (only_roll) {                      // rolling-stop work only
        if (want_roll_er) experiment_11();
        if (want_roll_small) experiment_10(argc, argv);
        experiment_9(large, max_pushes, max_merges);
        return 0;
    }

    if (want_bm_path) {
        experiment_12(argc, argv, bm_bootstraps, bm_alpha);
        return 0;
    }

    if (want_full_stop_report) {
        experiment_13(argc, argv, large, max_pushes, max_merges, report_seeds);
        return 0;
    }

    if (want_multicut_report) {
        experiment_14_multicut(argc, argv, large, max_pushes, stability_reps);
        return 0;
    }

    if (!only_large) {
        experiment_1();
        experiment_2();
        experiment_4();
        experiment_5();
        experiment_6(argc, argv);
        experiment_8(argc, argv);
    }
    experiment_7(large);
    if (!only_large) verify_theory();
    if (want_roll_er) experiment_11();
    if (want_roll_small) experiment_10(argc, argv);
    experiment_9(large, max_pushes);
    return 0;
}
