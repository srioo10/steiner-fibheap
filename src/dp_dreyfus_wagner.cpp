// dp_dreyfus_wagner.cpp -- Exact paradigm 1: dynamic programming over terminal subsets
// (Dreyfus & Wagner, Networks 1, 1971) with the Dijkstra layer of
// Erickson, Monma & Veinott (Math. Oper. Res. 12, 1987).
//
// Fix a root terminal r and let Q = T \ {r}, q = |Q| = k - 1. For S subset of Q, v in V:
//     dp[S][v] = cost of a cheapest tree that spans S + {v}.
// Base:   dp[{t}][v] = d(t, v).
// Step:   (merge)  m[S][v]  = min over S1 with 0 < S1 < S of  dp[S1][v] + dp[S \ S1][v]
//         (grow)   dp[S][v] = min over u of  m[S][u] + d(u, v)
// Answer: dp[Q][r].
// The grow step is a shortest-path problem with initial labels m[S][.]; it is solved by
// ONE multi-source Dijkstra on the Fibonacci heap (n INSERT, <= 2m DECREASE-KEY, n
// EXTRACT-MIN), instead of the O(n^2) Floyd-style relaxation of the original paper.
//
// Two savings:
//   * only splits with lowbit(S) in S1 are enumerated, so each unordered split is
//     visited once: sum over S of 2^{|S|-1} = (3^q + 1)/2 - 2^q   merges per vertex;
//   * the root is not a subset bit, which halves the table (2^q instead of 2^k rows).
// Time   O(3^q n + 2^q (m + n log n)),  memory 2^q * n * (8 + 4) bytes.
// The int32 pred table stores the last edge of the grow step, so the optimal tree is
// read back without recomputation (report, Section 5).
//
// Usage: dp_dreyfus_wagner <instance.stp> [--heap=fib|bin] [--time-limit=s] [--mem-gb=g] [--max-k=K]
#include "../include/steiner_common.hpp"

template <class PQ>
Result solve(const Graph& g, const Args& a, PQ& pq) {
    Result r;
    r.algo = "DP_DreyfusWagner";
    const int k = g.k();
    if (k == 1) { r.tree = {}; return r; }
    const int root_idx = (a.root >= 0 && a.root < k) ? a.root : 0;
    const int root = g.terminals[size_t(root_idx)];
    std::vector<int> Q;
    for (int i = 0; i < k; ++i) if (i != root_idx) Q.push_back(g.terminals[size_t(i)]);
    const int q = int(Q.size());
    const size_t n = size_t(g.n);
    const size_t layers = size_t(1) << q;

    double mem_gb = a.extra.count("mem-gb") ? std::stod(a.extra.at("mem-gb")) : 8.0;
    int max_k = a.extra.count("max-k") ? std::stoi(a.extra.at("max-k")) : 30;
    double need_gb = double(layers) * double(n) * 12.0 / 1e9;
    if (k > max_k || need_gb > mem_gb) {
        r.status = "skipped_memory"; r.has_tree = false;
        r.notes = "table_gb=" + std::to_string(need_gb);
        return r;
    }

    std::vector<int64_t> dp(layers * n, INF);
    std::vector<int32_t> pr(layers * n, -1);
    std::vector<int64_t> lab(n);
    std::vector<int> pred(n);
    Timer clock;
    unsigned long long merges = 0;

    for (size_t S = 1; S < layers; ++S) {
        int64_t* cur = &dp[S * n];
        if ((S & (S - 1)) == 0) {
            // base case: a single terminal
            int i = __builtin_ctzll(S);
            cur[size_t(Q[size_t(i)])] = 0;
        } else {
            // merge step: lowbit(S) always in S1, so {S1, S\S1} is visited once
            size_t low = S & (~S + 1);
            for (size_t S1 = (S - 1) & S; S1 > 0; S1 = (S1 - 1) & S) {
                if (!(S1 & low)) continue;
                const int64_t* A = &dp[S1 * n];
                const int64_t* B = &dp[(S ^ S1) * n];
                for (size_t v = 0; v < n; ++v) {
                    int64_t s = A[v] + B[v];
                    cur[v] = s < cur[v] ? s : cur[v];
                }
                ++merges;
            }
        }
        // grow step: one Dijkstra with initial labels cur[]
        std::copy(cur, cur + n, lab.begin());
        dijkstra(g, lab, pred, pq);
        std::copy(lab.begin(), lab.end(), cur);
        std::copy(pred.begin(), pred.end(), &pr[S * n]);
        if ((S & 1023) == 0 && clock.s() > a.time_limit) {
            r.status = "timeout"; r.has_tree = false;
            r.notes = "layers_done=" + std::to_string(S) + "/" + std::to_string(layers - 1);
            return r;
        }
    }

    const size_t FULL = layers - 1;
    const int64_t best = dp[FULL * n + size_t(root)];
    if (best >= INF) { r.status = "disconnected"; r.has_tree = false; return r; }

    // Read the tree back: pred edge -> stay in S; no pred -> it came from a merge
    // (or is a base terminal), find a split whose two halves add up exactly.
    std::vector<int> edges;
    std::vector<char> used(size_t(g.m()), 0);
    std::vector<std::pair<size_t, int>> stack{{FULL, root}};
    while (!stack.empty()) {
        auto [S, v] = stack.back(); stack.pop_back();
        int e = pr[S * n + size_t(v)];
        if (e >= 0) {
            if (!used[size_t(e)]++) edges.push_back(e);
            const Edge& ed = g.E[size_t(e)];
            stack.push_back({S, ed.u == v ? ed.v : ed.u});
            continue;
        }
        if ((S & (S - 1)) == 0) continue;   // base: v is the terminal itself
        size_t low = S & (~S + 1);
        int64_t target = dp[S * n + size_t(v)];
        bool found = false;
        for (size_t S1 = (S - 1) & S; S1 > 0 && !found; S1 = (S1 - 1) & S) {
            if (!(S1 & low)) continue;
            if (dp[S1 * n + size_t(v)] + dp[(S ^ S1) * n + size_t(v)] == target) {
                stack.push_back({S1, v}); stack.push_back({S ^ S1, v});
                found = true;
            }
        }
        if (!found) { r.status = "reconstruction_failed"; break; }
    }
    // The union of the traced paths can share vertices; MST + prune never raises cost.
    r.tree = mst_prune(g, edges, pq);
    int64_t c = cost_of(g, r.tree);
    if (c != best) r.status = "cost_mismatch_dp=" + std::to_string(best);
    r.notes = "dp_value=" + std::to_string(best) + ";subsets=" + std::to_string(layers - 1) +
              ";merges=" + std::to_string(merges) + ";table_mb=" + std::to_string(need_gb * 1000.0);
    return r;
}

int main(int argc, char** argv) {
    Args a = parse_args(argc, argv, "usage: dp_dreyfus_wagner <instance.stp> [--heap=fib|bin] [--time-limit=s] [--mem-gb=g] [--max-k=K] [--root=i] [--tree=file] [--header]");
    Graph g = read_stp(a.path);
    Timer t;
    Result r;
    if (a.heap == "bin") { BinPQ pq; r = solve(g, a, pq); }
    else { FibPQ pq; r = solve(g, a, pq); r.hs = pq.stats(); }
    r.time_ms = t.ms();
    return report(g, a, r);
}
