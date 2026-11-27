// greedy3_zelikovsky.cpp -- Greedy strategy G3: best-gain triples
// (A. Zelikovsky, "An 11/6-approximation algorithm for the network Steiner problem",
//  Algorithmica 9, 1993).
//
// Greedy criterion: G1 and G2 only ever connect terminals PAIRWISE, so a Steiner
// vertex that serves three terminals is never seen as one unit. G3 looks at every
// triple z = {a, b, c} of terminals and its best centre v(z),
//     d(z) = min_v d(a,v) + d(b,v) + d(c,v),
// and greedily contracts the triple with the largest WIN
//     win(z) = mst(F) - mst(F[z]) - d(z),
// where F is the current distance network and F[z] is F with a and b, b and c joined
// by zero-cost edges. It stops when no triple has positive win.
//
// Fast win evaluation (Lemma 4.3 in the report): if s_ab, s_bc, s_ac are the largest
// edge weights on the MST paths between the three pairs, then
//     mst(F) - mst(F[z]) = max(s_ab, s_bc, s_ac) + min(s_ab, s_bc, s_ac),
// so every round costs O(k^2) for the "save" matrix plus O(k^3) for the scan.
// d(z) is computed lazily (O(n)) and only for triples that can still beat the best
// win found so far, using the lower bound d(z) >= (D_ab + D_bc + D_ac) / 2.
//
// Total O(k (m + n log n) + k^4 + k^3 n) worst case. Ratio 11/6 (report, Sec. 4.3).
//
// Usage: greedy3_zelikovsky <instance.stp> [--heap=fib|bin] [--tree=file] [--header]
#include "../include/steiner_common.hpp"
#include <unordered_map>

using Mat = std::vector<std::vector<int64_t>>;

// save[x][y] = heaviest edge on the tree path x..y (tree given by parent array).
static Mat bottleneck_matrix(int k, const std::vector<int>& par, const Mat& F) {
    std::vector<std::vector<std::pair<int, int64_t>>> adj(static_cast<size_t>(k));
    for (int v = 0; v < k; ++v)
        if (par[size_t(v)] >= 0) {
            int64_t w = F[size_t(v)][size_t(par[size_t(v)])];
            adj[size_t(v)].push_back({par[size_t(v)], w});
            adj[size_t(par[size_t(v)])].push_back({v, w});
        }
    Mat save(static_cast<size_t>(k), std::vector<int64_t>(size_t(k), 0));
    std::vector<int> stack;
    std::vector<char> seen(static_cast<size_t>(k), 0);
    for (int s = 0; s < k; ++s) {
        std::fill(seen.begin(), seen.end(), 0);
        stack.assign(1, s); seen[size_t(s)] = 1;
        while (!stack.empty()) {
            int u = stack.back(); stack.pop_back();
            for (auto [v, w] : adj[size_t(u)]) {
                if (seen[size_t(v)]) continue;
                seen[size_t(v)] = 1;
                save[size_t(s)][size_t(v)] = std::max(save[size_t(s)][size_t(u)], w);
                stack.push_back(v);
            }
        }
    }
    return save;
}

template <class PQ>
Result solve(const Graph& g, PQ& pq) {
    Result r;
    r.algo = "G3_Zelikovsky";
    const int k = g.k();
    if (k == 1) { r.tree = {}; return r; }

    std::vector<std::vector<int64_t>> dist(static_cast<size_t>(k));
    std::vector<std::vector<int>> pred(static_cast<size_t>(k));
    for (int i = 0; i < k; ++i) dijkstra_from(g, g.terminals[size_t(i)], dist[size_t(i)], pred[size_t(i)], pq);
    Mat D(static_cast<size_t>(k), std::vector<int64_t>(size_t(k)));
    for (int i = 0; i < k; ++i)
        for (int j = 0; j < k; ++j) D[size_t(i)][size_t(j)] = dist[size_t(i)][size_t(g.terminals[size_t(j)])];
    for (int j = 1; j < k; ++j)
        if (D[0][size_t(j)] >= INF) { r.status = "disconnected"; r.has_tree = false; return r; }

    // memo for d(z): key = (a*k + b)*k + c with a<b<c
    std::unordered_map<int64_t, std::pair<int64_t, int>> center;
    auto triple_cost = [&](int a, int b, int c) -> std::pair<int64_t, int> {
        int64_t key = (int64_t(a) * k + b) * k + c;
        auto it = center.find(key);
        if (it != center.end()) return it->second;
        int64_t best = INF; int arg = -1;
        const auto &da = dist[size_t(a)], &db = dist[size_t(b)], &dc = dist[size_t(c)];
        for (int v = 0; v < g.n; ++v) {
            int64_t s = da[size_t(v)] + db[size_t(v)] + dc[size_t(v)];
            if (s < best) { best = s; arg = v; }
        }
        return center[key] = {best, arg};
    };

    Mat F = D;   // working distance network (contracted triples get zero edges)
    struct Chosen { int a, b, c, v; int64_t cost, win; };
    std::vector<Chosen> chosen;
    long long evaluated = 0, scanned = 0;
    for (;;) {
        std::vector<int> par = prim_dense(k, F, pq);
        Mat save = bottleneck_matrix(k, par, F);
        int64_t best_win = 0;
        Chosen bc{-1, -1, -1, -1, 0, 0};
        for (int a = 0; a < k; ++a)
            for (int b = a + 1; b < k; ++b) {
                int64_t sab = save[size_t(a)][size_t(b)];
                for (int c = b + 1; c < k; ++c) {
                    ++scanned;
                    int64_t sbc = save[size_t(b)][size_t(c)], sac = save[size_t(a)][size_t(c)];
                    int64_t gain = std::max({sab, sbc, sac}) + std::min({sab, sbc, sac});
                    if (gain <= best_win) continue;
                    // lower bound on d(z): each centre-to-terminal path is used by two pairs
                    int64_t lb = (D[size_t(a)][size_t(b)] + D[size_t(b)][size_t(c)] + D[size_t(a)][size_t(c)] + 1) / 2;
                    if (gain - lb <= best_win) continue;
                    ++evaluated;
                    auto [dz, v] = triple_cost(a, b, c);
                    if (gain - dz > best_win) { best_win = gain - dz; bc = {a, b, c, v, dz, best_win}; }
                }
            }
        if (best_win <= 0) break;
        chosen.push_back(bc);
        // contract: zero-cost edges a-b and b-c
        F[size_t(bc.a)][size_t(bc.b)] = F[size_t(bc.b)][size_t(bc.a)] = 0;
        F[size_t(bc.b)][size_t(bc.c)] = F[size_t(bc.c)][size_t(bc.b)] = 0;
    }

    // Build the subgraph: stars of the chosen triples + remaining (non-zero) MST edges of F.
    std::vector<int> par = prim_dense(k, F, pq);
    std::vector<int> edges;
    std::vector<char> used(size_t(g.m()), 0);
    auto add_path = [&](int term_idx, int target) {
        std::vector<int> p;
        path_edges(g, pred[size_t(term_idx)], target, p);
        for (int e : p) if (!used[size_t(e)]++) edges.push_back(e);
    };
    for (const Chosen& z : chosen) { add_path(z.a, z.v); add_path(z.b, z.v); add_path(z.c, z.v); }
    for (int j = 1; j < k; ++j) {
        int i = par[size_t(j)];
        if (F[size_t(i)][size_t(j)] == 0 && D[size_t(i)][size_t(j)] > 0) continue;  // inside a contracted triple
        add_path(i, g.terminals[size_t(j)]);
    }
    r.tree = mst_prune(g, edges, pq);
    r.notes = "triples=" + std::to_string(chosen.size()) + ";triples_scanned=" + std::to_string(scanned) +
              ";centres_evaluated=" + std::to_string(evaluated);
    return r;
}

int main(int argc, char** argv) {
    Args a = parse_args(argc, argv, "usage: greedy3_zelikovsky <instance.stp> [--heap=fib|bin] [--tree=file] [--header]");
    Graph g = read_stp(a.path);
    Timer t;
    Result r;
    if (a.heap == "bin") { BinPQ pq; r = solve(g, pq); }
    else { FibPQ pq; r = solve(g, pq); r.hs = pq.stats(); }
    r.time_ms = t.ms();
    return report(g, a, r);
}
