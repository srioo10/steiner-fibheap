// greedy2_sph.cpp -- Greedy strategy G2: shortest-path heuristic
// (Takahashi & Matsuyama, Math. Japonica 24, 1980), a Prim-like growth.
//
// Greedy criterion: keep ONE growing tree and always attach the terminal that is
// closest to the CURRENT TREE (not to another terminal, as G1 does).
//
//   T := {r}
//   while some terminal is outside T:
//       multi-source Dijkstra from all vertices of T (label 0), stopped as soon as
//       the first outside terminal t is extracted                O(m + n log n)
//       add the shortest path from T to t
//   MST of the induced subgraph + prune non-terminal leaves
// One root: O(k (m + n log n)). Ratio 2 - 2/k (report, Section 4.2).
// Because the result depends on the root, the solver runs every terminal as root
// (--roots=all, the default) and keeps the best tree: O(k^2 (m + n log n)).
// The notes column records the root-0 cost and the worst root, which the report uses
// to show the root dependence.
//
// Usage: greedy2_sph <instance.stp> [--heap=fib|bin] [--roots=all|one] [--root=i]
#include "../include/steiner_common.hpp"

template <class PQ>
std::vector<int> sph_from_root(const Graph& g, int root, PQ& pq) {
    std::vector<char> in_tree(size_t(g.n), 0), done(size_t(g.n), 0);
    std::vector<int64_t> dist(size_t(g.n));
    std::vector<int> pred(size_t(g.n));
    std::vector<int> edges;
    in_tree[size_t(root)] = 1;
    int left = g.k() - 1;
    while (left > 0) {
        // multi-source Dijkstra from the current tree, early exit at the first terminal
        std::fill(dist.begin(), dist.end(), INF);
        std::fill(pred.begin(), pred.end(), -1);
        std::fill(done.begin(), done.end(), 0);
        pq.reset(g.n);
        for (int v = 0; v < g.n; ++v) if (in_tree[size_t(v)]) { dist[size_t(v)] = 0; pq.push(v, 0); }
        int hit = -1;
        while (!pq.empty()) {
            int u = pq.pop().second;
            done[size_t(u)] = 1;
            if (g.is_term[size_t(u)] && !in_tree[size_t(u)]) { hit = u; break; }
            for (const Arc* a = g.begin(u); a != g.end(u); ++a) {
                if (done[size_t(a->to)]) continue;
                int64_t nd = dist[size_t(u)] + a->w;
                if (nd < dist[size_t(a->to)]) {
                    dist[size_t(a->to)] = nd; pred[size_t(a->to)] = a->e;
                    pq.push_or_decrease(a->to, nd);
                }
            }
        }
        if (hit < 0) return {};   // graph disconnected
        // walk back from the terminal to the tree, adding the path
        int v = hit;
        while (!in_tree[size_t(v)]) {
            in_tree[size_t(v)] = 1;
            int e = pred[size_t(v)];
            edges.push_back(e);
            v = g.E[size_t(e)].u == v ? g.E[size_t(e)].v : g.E[size_t(e)].u;
            if (g.is_term[size_t(v)] && !in_tree[size_t(v)]) --left;  // cannot happen (hit was nearest), kept for safety
        }
        --left;
    }
    return mst_prune(g, edges, pq);
}

template <class PQ>
Result solve(const Graph& g, const Args& a, PQ& pq) {
    Result r;
    r.algo = "G2_SPH";
    bool all = !(a.extra.count("roots") && a.extra.at("roots") == "one");
    int64_t best = INF, worst = -1, root0 = -1;
    int best_root = -1;
    for (int i = 0; i < g.k(); ++i) {
        if (!all && i != a.root) continue;
        std::vector<int> t = sph_from_root(g, g.terminals[size_t(i)], pq);
        if (t.empty() && g.k() > 1) { r.status = "disconnected"; r.has_tree = false; return r; }
        int64_t c = cost_of(g, t);
        if (i == a.root) root0 = c;
        worst = std::max(worst, c);
        if (c < best) { best = c; best_root = i; r.tree = t; }
    }
    r.algo = all ? "G2_SPH" : "G2_SPH_root";
    r.notes = "root0=" + std::to_string(root0) + ";worst_root=" + std::to_string(worst) +
              ";best_root_index=" + std::to_string(best_root);
    return r;
}

int main(int argc, char** argv) {
    Args a = parse_args(argc, argv, "usage: greedy2_sph <instance.stp> [--heap=fib|bin] [--roots=all|one] [--root=i] [--tree=file] [--header]");
    Graph g = read_stp(a.path);
    Timer t;
    Result r;
    if (a.heap == "bin") { BinPQ pq; r = solve(g, a, pq); }
    else { FibPQ pq; r = solve(g, a, pq); r.hs = pq.stats(); }
    r.time_ms = t.ms();
    return report(g, a, r);
}
