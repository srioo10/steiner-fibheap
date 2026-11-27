// graph.hpp -- undirected weighted graph in CSR form, SteinLib .stp reader,
// and an independent validator for Steiner trees.
#pragma once
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

constexpr int64_t INF = INT64_MAX / 4;

struct Edge { int u, v; int64_t w; };
struct Arc { int to; int64_t w; int e; };   // e = index into Graph::E

struct Graph {
    std::string name;
    int n = 0;
    std::vector<Edge> E;
    std::vector<int> off;          // CSR offsets, size n+1
    std::vector<Arc> adj;          // CSR arcs, size 2|E|
    std::vector<int> terminals;    // 0-based vertex ids
    std::vector<char> is_term;

    int m() const { return int(E.size()); }
    int k() const { return int(terminals.size()); }
    int deg(int v) const { return off[size_t(v) + 1] - off[size_t(v)]; }
    const Arc* begin(int v) const { return adj.data() + off[size_t(v)]; }
    const Arc* end(int v) const { return adj.data() + off[size_t(v) + 1]; }

    // Builds CSR adjacency from E. O(n + m).
    void build() {
        off.assign(size_t(n) + 1, 0);
        for (const Edge& e : E) { ++off[size_t(e.u) + 1]; ++off[size_t(e.v) + 1]; }
        for (int i = 0; i < n; ++i) off[size_t(i) + 1] += off[size_t(i)];
        adj.assign(2 * E.size(), Arc{});
        std::vector<int> pos(off.begin(), off.end() - 1);
        for (int i = 0; i < m(); ++i) {
            const Edge& e = E[size_t(i)];
            adj[size_t(pos[size_t(e.u)]++)] = {e.v, e.w, i};
            adj[size_t(pos[size_t(e.v)]++)] = {e.u, e.w, i};
        }
        is_term.assign(size_t(n), 0);
        for (int t : terminals) is_term[size_t(t)] = 1;
    }
};

// ---- SteinLib STP reader (format v1.0: SECTION Graph / E u v w / SECTION Terminals / T v).
inline Graph read_stp(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);
    Graph g;
    {   // instance name = file name without directory and extension, lower case
        std::string b = path.substr(path.find_last_of('/') + 1);
        b = b.substr(0, b.find('.'));
        for (char& c : b) c = char(std::tolower(static_cast<unsigned char>(c)));
        g.name = b;
    }
    std::string line, section;
    while (std::getline(in, line)) {
        std::istringstream ss(line);
        std::string tok;
        if (!(ss >> tok)) continue;
        std::string up = tok;
        for (char& c : up) c = char(std::toupper(static_cast<unsigned char>(c)));
        if (up == "SECTION") { ss >> section; for (char& c : section) c = char(std::toupper(static_cast<unsigned char>(c))); continue; }
        if (up == "END") { section.clear(); continue; }
        if (section == "GRAPH") {
            if (up == "NODES") ss >> g.n;
            else if (up == "E" || up == "A") {
                Edge e; ss >> e.u >> e.v >> e.w;
                if (!ss) throw std::runtime_error("bad edge line: " + line);
                if (e.u < 1 || e.v < 1 || e.u > g.n || e.v > g.n) throw std::runtime_error("edge endpoint out of range: " + line);
                if (e.w < 0) throw std::runtime_error("negative edge weight: " + line);
                --e.u; --e.v;
                if (e.u != e.v) g.E.push_back(e);      // self-loops never help a Steiner tree
            }
        } else if (section == "TERMINALS") {
            if (up == "T") {
                int t; ss >> t;
                if (!ss || t < 1 || t > g.n) throw std::runtime_error("bad terminal line: " + line);
                g.terminals.push_back(t - 1);
            }
        }
    }
    std::sort(g.terminals.begin(), g.terminals.end());
    g.terminals.erase(std::unique(g.terminals.begin(), g.terminals.end()), g.terminals.end());
    if (g.n <= 0 || g.terminals.empty()) throw std::runtime_error("empty instance " + path);
    g.build();
    return g;
}

// Writes a graph back in STP format (used for the generated adversarial instances).
inline void write_stp(const Graph& g, const std::string& path, const std::string& remark) {
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) throw std::runtime_error("cannot write " + path);
    std::fprintf(f, "33D32945 STP File, STP Format Version 1.0\n\nSECTION Comment\nName    \"%s\"\nRemark  \"%s\"\nEND\n\n",
                 g.name.c_str(), remark.c_str());
    std::fprintf(f, "SECTION Graph\nNodes %d\nEdges %d\n", g.n, g.m());
    for (const Edge& e : g.E) std::fprintf(f, "E %d %d %lld\n", e.u + 1, e.v + 1, (long long)e.w);
    std::fprintf(f, "END\n\nSECTION Terminals\nTerminals %d\n", g.k());
    for (int t : g.terminals) std::fprintf(f, "T %d\n", t + 1);
    std::fprintf(f, "END\n\nEOF\n");
    std::fclose(f);
}

// ---- Disjoint-set union (path halving + union by size), used by the validator and Kruskal steps.
struct DSU {
    std::vector<int> p, s;
    explicit DSU(int n) : p(size_t(n)), s(size_t(n), 1) { std::iota(p.begin(), p.end(), 0); }
    int find(int x) { while (p[size_t(x)] != x) { p[size_t(x)] = p[size_t(p[size_t(x)])]; x = p[size_t(x)]; } return x; }
    bool unite(int a, int b) {
        a = find(a); b = find(b);
        if (a == b) return false;
        if (s[size_t(a)] < s[size_t(b)]) std::swap(a, b);
        p[size_t(b)] = a; s[size_t(a)] += s[size_t(b)];
        return true;
    }
};

// ---- Validator: independent of every solver. A set of edge ids is a Steiner tree iff
// (1) ids are valid and distinct, (2) the edges contain no cycle, (3) the edges form a
// single connected component, and (4) every terminal lies in it.
struct Validation { bool ok; int64_t cost; std::string msg; };

inline Validation validate_tree(const Graph& g, const std::vector<int>& edges) {
    std::vector<char> seen(size_t(g.m()), 0), touched(size_t(g.n), 0);
    DSU d(g.n);
    int64_t cost = 0;
    int nv = 0;
    for (int id : edges) {
        if (id < 0 || id >= g.m()) return {false, 0, "edge id out of range"};
        if (seen[size_t(id)]++) return {false, 0, "duplicate edge"};
        const Edge& e = g.E[size_t(id)];
        if (!d.unite(e.u, e.v)) return {false, 0, "cycle"};
        cost += e.w;
        if (!touched[size_t(e.u)]++) ++nv;
        if (!touched[size_t(e.v)]++) ++nv;
    }
    if (g.k() == 1) return {edges.empty(), cost, edges.empty() ? "" : "single terminal needs no edges"};
    if (nv != int(edges.size()) + 1) return {false, cost, "not connected"};
    int r = d.find(g.terminals[0]);
    for (int t : g.terminals)
        if (!touched[size_t(t)] || d.find(t) != r) return {false, cost, "terminal not spanned"};
    return {true, cost, ""};
}

// Known optimum from instances/steinlib/optima.csv (name,n,m,k,opt); -1 if unknown.
inline int64_t known_optimum(const std::string& name, const std::string& csv) {
    std::ifstream in(csv);
    std::string line;
    while (std::getline(in, line)) {
        std::vector<std::string> f;
        std::stringstream ss(line);
        std::string x;
        while (std::getline(ss, x, ',')) f.push_back(x);
        if (f.size() >= 5 && f[0] == name) return std::stoll(f[4]);
    }
    return -1;
}
