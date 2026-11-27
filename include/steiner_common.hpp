// steiner_common.hpp -- building blocks shared by the five solvers:
// Dijkstra and Prim templated on the priority queue (FibPQ or BinPQ),
// path extraction, the MST-and-prune cleanup step, timing, CLI and CSV output.
#pragma once
#include "binary_heap.hpp"
#include "fib_heap.hpp"
#include "graph.hpp"
#include <chrono>
#include <cstring>
#include <map>
#include <sys/stat.h>

// ---------------------------------------------------------------- heap statistics
inline FibStats heap_stats(const FibPQ& q) { return q.stats(); }
inline FibStats heap_stats(const BinPQ&) { return FibStats(); }

// ---------------------------------------------------------------- Dijkstra
// Multi-source Dijkstra. dist[] arrives initialised: dist[v] < INF for every source
// (its starting label), INF elsewhere. On return dist[v] is the shortest
// label-plus-path value and pred[v] the last edge on that path (-1 for sources).
// Every vertex is inserted at most once and every arc causes at most one
// DECREASE-KEY, so the run costs O(m + n log n) with FibPQ, O(m log n) with BinPQ.
template <class PQ>
void dijkstra(const Graph& g, std::vector<int64_t>& dist, std::vector<int>& pred, PQ& pq) {
    pq.reset(g.n);
    pred.assign(size_t(g.n), -1);
    std::vector<char> done(size_t(g.n), 0);
    for (int v = 0; v < g.n; ++v)
        if (dist[size_t(v)] < INF) pq.push(v, dist[size_t(v)]);
    while (!pq.empty()) {
        auto [d, u] = pq.pop();
        done[size_t(u)] = 1;
        for (const Arc* a = g.begin(u); a != g.end(u); ++a) {
            if (done[size_t(a->to)]) continue;
            int64_t nd = d + a->w;
            if (nd < dist[size_t(a->to)]) {
                dist[size_t(a->to)] = nd;
                pred[size_t(a->to)] = a->e;
                pq.push_or_decrease(a->to, nd);
            }
        }
    }
}

// Single-source convenience wrapper.
template <class PQ>
void dijkstra_from(const Graph& g, int s, std::vector<int64_t>& dist, std::vector<int>& pred, PQ& pq) {
    dist.assign(size_t(g.n), INF);
    dist[size_t(s)] = 0;
    dijkstra(g, dist, pred, pq);
}

// Edges of the shortest path ending at v, following pred[] back to a source.
inline void path_edges(const Graph& g, const std::vector<int>& pred, int v, std::vector<int>& out) {
    while (pred[size_t(v)] >= 0) {
        int e = pred[size_t(v)];
        out.push_back(e);
        v = g.E[size_t(e)].u == v ? g.E[size_t(e)].v : g.E[size_t(e)].u;
    }
}

// ---------------------------------------------------------------- Prim
// Minimum spanning forest of the subgraph induced by in_set (Prim, one tree per
// component). Returns edge ids. O(m + n log n) with FibPQ.
template <class PQ>
std::vector<int> prim_induced(const Graph& g, const std::vector<char>& in_set, PQ& pq) {
    std::vector<int> tree;
    std::vector<int64_t> key(size_t(g.n), INF);
    std::vector<int> via(size_t(g.n), -1);
    std::vector<char> done(size_t(g.n), 0);
    pq.reset(g.n);
    for (int s = 0; s < g.n; ++s) {
        if (!in_set[size_t(s)] || done[size_t(s)]) continue;
        key[size_t(s)] = 0;
        pq.push(s, 0);
        while (!pq.empty()) {
            int u = pq.pop().second;
            done[size_t(u)] = 1;
            if (via[size_t(u)] >= 0) tree.push_back(via[size_t(u)]);
            for (const Arc* a = g.begin(u); a != g.end(u); ++a) {
                int v = a->to;
                if (!in_set[size_t(v)] || done[size_t(v)] || a->w >= key[size_t(v)]) continue;
                key[size_t(v)] = a->w;
                via[size_t(v)] = a->e;
                pq.push_or_decrease(v, a->w);
            }
        }
    }
    return tree;
}

// Prim on a complete graph given as a k x k matrix (the terminal distance network).
// Returns parent[] (parent[root] = -1). k inserts, <= k^2 DECREASE-KEYs, k EXTRACT-MINs:
// O(k^2) with FibPQ (the DECREASE-KEYs are O(1) amortized), O(k^2 log k) with BinPQ.
template <class PQ>
std::vector<int> prim_dense(int k, const std::vector<std::vector<int64_t>>& D, PQ& pq) {
    std::vector<int> par(size_t(k), -1);
    std::vector<int64_t> key(size_t(k), INF);
    std::vector<char> done(size_t(k), 0);
    pq.reset(k);
    key[0] = 0;
    pq.push(0, 0);
    while (!pq.empty()) {
        int u = pq.pop().second;
        done[size_t(u)] = 1;
        for (int v = 0; v < k; ++v) {
            if (done[size_t(v)] || D[size_t(u)][size_t(v)] >= key[size_t(v)]) continue;
            key[size_t(v)] = D[size_t(u)][size_t(v)];
            par[size_t(v)] = u;
            pq.push_or_decrease(v, key[size_t(v)]);
        }
    }
    return par;
}

// ---------------------------------------------------------------- cleanup
// Repeatedly delete non-terminal leaves from an edge set (they only add cost).
inline std::vector<int> prune_leaves(const Graph& g, const std::vector<int>& edges) {
    std::vector<int> degree(size_t(g.n), 0);
    std::vector<std::vector<int>> inc(size_t(g.n));
    std::vector<char> alive(edges.size(), 1);
    for (size_t i = 0; i < edges.size(); ++i) {
        const Edge& e = g.E[size_t(edges[i])];
        ++degree[size_t(e.u)]; ++degree[size_t(e.v)];
        inc[size_t(e.u)].push_back(int(i)); inc[size_t(e.v)].push_back(int(i));
    }
    std::vector<int> stack;
    for (int v = 0; v < g.n; ++v)
        if (degree[size_t(v)] == 1 && !g.is_term[size_t(v)]) stack.push_back(v);
    while (!stack.empty()) {
        int v = stack.back(); stack.pop_back();
        if (degree[size_t(v)] != 1) continue;
        for (int i : inc[size_t(v)]) {
            if (!alive[size_t(i)]) continue;
            alive[size_t(i)] = 0;
            const Edge& e = g.E[size_t(edges[size_t(i)])];
            int w = e.u == v ? e.v : e.u;
            --degree[size_t(v)];
            if (--degree[size_t(w)] == 1 && !g.is_term[size_t(w)]) stack.push_back(w);
        }
    }
    std::vector<int> out;
    for (size_t i = 0; i < edges.size(); ++i) if (alive[i]) out.push_back(edges[i]);
    return out;
}

// The classical final step of KMB-type heuristics: take every vertex touched by a
// candidate edge set, compute the MST of the subgraph they induce, and prune
// non-terminal leaves. The cost never increases (Lemma in report, Sec. 4).
template <class PQ>
std::vector<int> mst_prune(const Graph& g, const std::vector<int>& edges, PQ& pq) {
    std::vector<char> in(size_t(g.n), 0);
    for (int id : edges) { in[size_t(g.E[size_t(id)].u)] = 1; in[size_t(g.E[size_t(id)].v)] = 1; }
    if (edges.empty()) for (int t : g.terminals) in[size_t(t)] = 1;
    return prune_leaves(g, prim_induced(g, in, pq));
}

inline int64_t cost_of(const Graph& g, const std::vector<int>& edges) {
    int64_t c = 0;
    for (int id : edges) c += g.E[size_t(id)].w;
    return c;
}

// ---------------------------------------------------------------- timing
struct Timer {
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    double ms() const { return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count(); }
    double s() const { return ms() / 1000.0; }
};

// ---------------------------------------------------------------- CLI and output
struct Args {
    std::string path;
    std::string heap = "fib";       // fib | bin
    double time_limit = 300.0;      // seconds
    int root = 0;                   // terminal index used as root by SPH / DP
    bool header = false;
    std::string tree_out;           // optional: write the tree's edges (1-based endpoints)
    std::map<std::string, std::string> extra;  // any other --key=value
};

inline Args parse_args(int argc, char** argv, const char* usage) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        std::string s = argv[i];
        auto val = [&](const char* k) { return s.substr(std::strlen(k)); };
        if (s == "--header") a.header = true;
        else if (s.rfind("--heap=", 0) == 0) a.heap = val("--heap=");
        else if (s.rfind("--time-limit=", 0) == 0) a.time_limit = std::stod(val("--time-limit="));
        else if (s.rfind("--root=", 0) == 0) a.root = std::stoi(val("--root="));
        else if (s.rfind("--tree=", 0) == 0) a.tree_out = val("--tree=");
        else if (s.rfind("--", 0) == 0 && s.find('=') != std::string::npos)
            a.extra[s.substr(2, s.find('=') - 2)] = s.substr(s.find('=') + 1);
        else if (s == "-h" || s == "--help") { std::printf("%s\n", usage); std::exit(0); }
        else a.path = s;
    }
    if (a.path.empty() && !a.header) { std::fprintf(stderr, "%s\n", usage); std::exit(2); }
    if (a.heap != "fib" && a.heap != "bin") { std::fprintf(stderr, "--heap must be fib or bin\n"); std::exit(2); }
    return a;
}

// optima.csv is looked up next to the instance, then one directory up.
inline int64_t lookup_optimum(const Graph& g, const std::string& path) {
    std::string dir = path.find('/') == std::string::npos ? "." : path.substr(0, path.find_last_of('/'));
    struct stat sb;
    for (std::string c : {dir + "/optima.csv", dir + "/../optima.csv"})
        if (stat(c.c_str(), &sb) == 0) { int64_t o = known_optimum(g.name, c); if (o >= 0) return o; }
    return -1;
}

inline const char* csv_header() {
    return "algo,instance,n,m,k,cost,opt,gap_pct,time_ms,valid,status,"
           "heap,heap_inserts,heap_extracts,heap_decreases,heap_links,heap_cuts,heap_actual,notes";
}

struct Result {
    std::string algo, status = "ok", notes;
    std::vector<int> tree;
    bool has_tree = true;
    double time_ms = 0;
    FibStats hs;
};

// Validates the tree independently, prints one CSV row and returns the exit code.
inline int report(const Graph& g, const Args& a, Result r) {
    if (a.header) std::printf("%s\n", csv_header());
    int64_t opt = lookup_optimum(g, a.path);
    Validation v{false, -1, "no tree"};
    if (r.has_tree) v = validate_tree(g, r.tree);
    int64_t cost = r.has_tree ? v.cost : -1;
    double gap = (opt > 0 && cost >= 0) ? 100.0 * double(cost - opt) / double(opt) : -1;
    if (r.has_tree && !v.ok) { r.status = "INVALID:" + v.msg; }
    std::printf("%s,%s,%d,%d,%d,%lld,%lld,%.4f,%.3f,%d,%s,%s,%llu,%llu,%llu,%llu,%llu,%llu,%s\n",
                r.algo.c_str(), g.name.c_str(), g.n, g.m(), g.k(), (long long)cost, (long long)opt, gap,
                r.time_ms, v.ok ? 1 : 0, r.status.c_str(), a.heap.c_str(),
                (unsigned long long)r.hs.inserts, (unsigned long long)r.hs.extracts,
                (unsigned long long)r.hs.decreases, (unsigned long long)r.hs.links,
                (unsigned long long)r.hs.cuts, (unsigned long long)r.hs.actual_cost, r.notes.c_str());
    if (!a.tree_out.empty() && v.ok) {
        std::FILE* f = std::fopen(a.tree_out.c_str(), "w");
        if (f) {
            for (int id : r.tree) std::fprintf(f, "%d %d %lld\n", g.E[size_t(id)].u + 1, g.E[size_t(id)].v + 1, (long long)g.E[size_t(id)].w);
            std::fclose(f);
        }
    }
    return (r.has_tree && !v.ok) ? 3 : 0;
}
