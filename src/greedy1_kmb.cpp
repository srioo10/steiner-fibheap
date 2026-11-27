// greedy1_kmb.cpp -- Greedy strategy G1: distance-network MST
// (Kou, Markowsky & Berman, Acta Informatica 15, 1981).
//
// Greedy criterion: repeatedly join the two terminal components that are CLOSEST in
// the shortest-path metric (Prim on the terminal distance network).
//
//   1. For every terminal t_i run Dijkstra (Fibonacci heap)       k * O(m + n log n)
//   2. D[i][j] = d_G(t_i, t_j); MST of the complete graph (K_k, D)    O(k^2)
//   3. Replace each MST edge by its shortest path in G                 O(k n)
//   4. MST of the subgraph induced by those vertices (Prim, Fib heap)  O(m + n log n)
//   5. Delete non-terminal leaves repeatedly                           O(n)
// Total O(k (m + n log n)). Ratio 2 - 2/l, l = number of leaves of an optimal tree
// (proof and a family where it is tight: report, Section 4.1; instances/adversarial/).
//
// Usage: greedy1_kmb <instance.stp> [--heap=fib|bin] [--tree=out.txt] [--header]
#include "../include/steiner_common.hpp"

template <class PQ>
Result solve(const Graph& g, PQ& pq) {
    Result r;
    r.algo = "G1_KMB";
    const int k = g.k();
    if (k == 1) { r.tree = {}; return r; }

    // Step 1: shortest-path trees from every terminal.
    std::vector<std::vector<int64_t>> dist(static_cast<size_t>(k));
    std::vector<std::vector<int>> pred(static_cast<size_t>(k));
    for (int i = 0; i < k; ++i) dijkstra_from(g, g.terminals[size_t(i)], dist[size_t(i)], pred[size_t(i)], pq);

    // Step 2: distance network on the terminals and its MST.
    std::vector<std::vector<int64_t>> D(static_cast<size_t>(k), std::vector<int64_t>(size_t(k)));
    for (int i = 0; i < k; ++i)
        for (int j = 0; j < k; ++j) D[size_t(i)][size_t(j)] = dist[size_t(i)][size_t(g.terminals[size_t(j)])];
    for (int j = 1; j < k; ++j)
        if (D[0][size_t(j)] >= INF) { r.status = "disconnected"; r.has_tree = false; return r; }
    std::vector<int> par = prim_dense(k, D, pq);

    // Step 3: expand every distance-network edge into a shortest path of G.
    std::vector<int> edges;
    std::vector<char> used(size_t(g.m()), 0);
    for (int j = 1; j < k; ++j) {
        std::vector<int> p;
        path_edges(g, pred[size_t(par[size_t(j)])], g.terminals[size_t(j)], p);
        for (int e : p) if (!used[size_t(e)]++) edges.push_back(e);
    }
    int64_t mst_d = 0;
    for (int j = 1; j < k; ++j) mst_d += D[size_t(par[size_t(j)])][size_t(j)];

    // Steps 4-5: MST of the induced subgraph, then prune non-terminal leaves.
    r.tree = mst_prune(g, edges, pq);
    r.notes = "mst_distance_network=" + std::to_string(mst_d) +
              ";expanded=" + std::to_string(cost_of(g, edges));
    return r;
}

int main(int argc, char** argv) {
    Args a = parse_args(argc, argv, "usage: greedy1_kmb <instance.stp> [--heap=fib|bin] [--tree=file] [--header]");
    Graph g = read_stp(a.path);
    Timer t;
    Result r;
    if (a.heap == "bin") { BinPQ pq; r = solve(g, pq); }
    else { FibPQ pq; r = solve(g, pq); r.hs = pq.stats(); }
    r.time_ms = t.ms();
    return report(g, a, r);
}
