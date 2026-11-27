// bnb_steiner.cpp -- Exact paradigm 2: best-first branch-and-bound.
//
// Pipeline (report, Section 6):
//   1. Exact reductions at the root (they never remove every optimal tree):
//        R1 non-terminal of degree <= 1 is deleted
//        R2 terminal of degree 1: its edge is in every solution; fix it, contract
//        R3 non-terminal of degree 2: its two edges are merged into one
//        R4 long edge: edge (u,v) with w(u,v) > d(u,v) is on no shortest path; delete
//   2. Lower bound at every node: Wong's dual ascent on the bidirected cut LP
//        min sum c_a x_a  s.t.  x(delta^-(W)) >= 1  for W with W cap T != {}, r not in W.
//      Any dual-feasible y gives LB = sum y_W  <=  OPT (weak duality), so the bound is
//      admissible. Terminals waiting for an ascent step sit in a Fibonacci heap keyed
//      by the size of their current cut (smallest cut first, Wong 1984).
//   3. Upper bound at every node: shortest-path heuristic on the node graph and on
//      the subgraph of arcs whose reduced cost is zero.
//   4. Reduced-cost elimination: with reduced costs c', LB and incumbent UB, a
//      non-terminal v (an arc u->v) can be in no tree cheaper than UB if
//        LB + d'(r, v) + d'(v, T)            >= UB      (vertex test)
//        LB + d'(r, u) + c'(u,v) + d'(v, T)  >= UB      (arc test)
//      d' are shortest paths under c' (two Dijkstras on the Fibonacci heap).
//   5. Branching on a Steiner vertex v: child 1 makes v a terminal, child 2 deletes v.
//      Open nodes live in a Fibonacci heap keyed by their lower bound (best-first);
//      a node is pruned when LB >= UB. Integer weights make LB integral.
// The search stops with a proof of optimality (heap empty), or at the time / node
// limit, in which case the best tree and the global lower bound are both reported.
//
// Usage: bnb_steiner <instance.stp> [--time-limit=s] [--node-limit=N] [--no-reduce=1]
#include "../include/steiner_common.hpp"
#include <map>

// ============================================================ 1. reductions
struct RedEdge { int u, v; int64_t w; std::vector<int> orig; bool alive; };

struct Reducer {
    int n;
    std::vector<RedEdge> E;
    std::vector<std::map<int, int>> nb;   // neighbour -> edge index (one edge per pair)
    std::vector<char> term, alive;
    int64_t fixed_cost = 0;
    std::vector<int> fixed_orig;
    int r1 = 0, r2 = 0, r3 = 0, r4 = 0;

    explicit Reducer(const Graph& g) : n(g.n), nb(size_t(g.n)), term(g.is_term), alive(size_t(g.n), 1) {
        for (int i = 0; i < g.m(); ++i) add_edge(g.E[size_t(i)].u, g.E[size_t(i)].v, g.E[size_t(i)].w, {i});
    }
    int nterm() const { int c = 0; for (int v = 0; v < n; ++v) c += alive[size_t(v)] && term[size_t(v)]; return c; }
    // keeps only the cheapest of parallel edges
    void add_edge(int u, int v, int64_t w, std::vector<int> orig) {
        auto it = nb[size_t(u)].find(v);
        if (it != nb[size_t(u)].end()) {
            RedEdge& e = E[size_t(it->second)];
            if (e.w <= w) return;
            e.w = w; e.orig = std::move(orig);
            return;
        }
        E.push_back({u, v, w, std::move(orig), true});
        nb[size_t(u)][v] = int(E.size()) - 1;
        nb[size_t(v)][u] = int(E.size()) - 1;
    }
    void del_edge(int id) {
        RedEdge& e = E[size_t(id)];
        e.alive = false;
        nb[size_t(e.u)].erase(e.v); nb[size_t(e.v)].erase(e.u);
    }
    bool degree_pass() {
        bool any = false, changed = true;
        while (changed) {
            changed = false;
            for (int v = 0; v < n; ++v) {
                if (!alive[size_t(v)]) continue;
                size_t d = nb[size_t(v)].size();
                if (!term[size_t(v)] && d <= 1) {                         // R1
                    if (d == 1) del_edge(nb[size_t(v)].begin()->second);
                    alive[size_t(v)] = 0; ++r1; changed = true;
                } else if (term[size_t(v)] && d == 1 && nterm() > 1) {     // R2
                    int id = nb[size_t(v)].begin()->second;
                    RedEdge e = E[size_t(id)];
                    int u = e.u == v ? e.v : e.u;
                    fixed_cost += e.w;
                    fixed_orig.insert(fixed_orig.end(), e.orig.begin(), e.orig.end());
                    del_edge(id);
                    alive[size_t(v)] = 0; term[size_t(u)] = 1; ++r2; changed = true;
                } else if (!term[size_t(v)] && d == 2) {                  // R3
                    auto it = nb[size_t(v)].begin();
                    int a = it->first, ea = it->second; ++it;
                    int b = it->first, eb = it->second;
                    RedEdge A = E[size_t(ea)], B = E[size_t(eb)];
                    del_edge(ea); del_edge(eb);
                    alive[size_t(v)] = 0;
                    std::vector<int> o = A.orig; o.insert(o.end(), B.orig.begin(), B.orig.end());
                    add_edge(a, b, A.w + B.w, std::move(o));
                    ++r3; changed = true;
                }
            }
            any |= changed;
        }
        return any;
    }
    Graph to_graph(std::vector<int>& red_id) const {
        Graph g;
        g.n = n;
        red_id.clear();
        for (int i = 0; i < int(E.size()); ++i)
            if (E[size_t(i)].alive) { g.E.push_back({E[size_t(i)].u, E[size_t(i)].v, E[size_t(i)].w}); red_id.push_back(i); }
        for (int v = 0; v < n; ++v) if (alive[size_t(v)] && term[size_t(v)]) g.terminals.push_back(v);
        g.build();
        return g;
    }
    bool long_edge_pass(FibPQ& pq) {
        std::vector<int> rid;
        Graph g = to_graph(rid);
        std::vector<int64_t> d; std::vector<int> p;
        std::vector<int> kill;
        for (int u = 0; u < n; ++u) {
            if (!alive[size_t(u)] || g.deg(u) == 0) continue;
            dijkstra_from(g, u, d, p, pq);
            for (const Arc* a = g.begin(u); a != g.end(u); ++a)
                if (a->to > u && d[size_t(a->to)] < a->w) kill.push_back(rid[size_t(a->e)]);   // R4
        }
        for (int id : kill) if (E[size_t(id)].alive) { del_edge(id); ++r4; }
        return !kill.empty();
    }
    void run(FibPQ& pq) {
        for (int round = 0; round < 10; ++round) {
            bool a = degree_pass();
            bool b = long_edge_pass(pq);
            if (!a && !b) break;
        }
    }
};

// ============================================================ 2. dual ascent
// arc ids: edge e = (u,v) gives arc 2e = u->v and arc 2e+1 = v->u
static inline int arc_out(const Graph& g, int x, const Arc& a) { return g.E[size_t(a.e)].u == x ? 2 * a.e : 2 * a.e + 1; }
static inline int arc_in(const Graph& g, int x, const Arc& a) { return g.E[size_t(a.e)].u == x ? 2 * a.e + 1 : 2 * a.e; }

struct DualAscent { bool feasible = true; int64_t lb = 0; std::vector<int64_t> rc; long long steps = 0; };

static DualAscent dual_ascent(const Graph& g, int root) {
    DualAscent da;
    da.rc.resize(2 * size_t(g.m()));
    for (int e = 0; e < g.m(); ++e) da.rc[2 * size_t(e)] = da.rc[2 * size_t(e) + 1] = g.E[size_t(e)].w;
    FibHeap queue;                                   // active terminals, key = |cut|
    for (int t : g.terminals) if (t != root) queue.insert(g.deg(t), t);
    std::vector<int> stamp(size_t(g.n), 0), stack, W;
    int cur = 0;
    while (!queue.empty()) {
        int t = queue.extract_min().second;
        // W = all vertices that reach t through arcs of reduced cost 0
        ++cur; W.clear(); stack.assign(1, t); stamp[size_t(t)] = cur;
        bool has_root = false;
        while (!stack.empty()) {
            int x = stack.back(); stack.pop_back();
            W.push_back(x);
            if (x == root) { has_root = true; break; }
            for (const Arc* a = g.begin(x); a != g.end(x); ++a)
                if (stamp[size_t(a->to)] != cur && da.rc[size_t(arc_in(g, x, *a))] == 0) {
                    stamp[size_t(a->to)] = cur; stack.push_back(a->to);
                }
        }
        if (has_root) continue;                      // t is connected to the root: inactive
        // delta = smallest reduced cost on the arcs entering W
        int64_t delta = INF; int cut = 0;
        for (int x : W)
            for (const Arc* a = g.begin(x); a != g.end(x); ++a)
                if (stamp[size_t(a->to)] != cur) { delta = std::min(delta, da.rc[size_t(arc_in(g, x, *a))]); ++cut; }
        if (delta >= INF) { da.feasible = false; return da; }   // no arc enters W: infeasible
        for (int x : W)
            for (const Arc* a = g.begin(x); a != g.end(x); ++a)
                if (stamp[size_t(a->to)] != cur) da.rc[size_t(arc_in(g, x, *a))] -= delta;
        da.lb += delta;
        ++da.steps;
        queue.insert(cut, t);
    }
    return da;
}

// Directed Dijkstra on reduced costs. forward: d(root, v); reverse: d(v, T).
static void rc_dijkstra(const Graph& g, const std::vector<int64_t>& rc, const std::vector<int>& src,
                        bool forward, std::vector<int64_t>& d, FibPQ& pq) {
    d.assign(size_t(g.n), INF);
    pq.reset(g.n);
    for (int s : src) { d[size_t(s)] = 0; pq.push_or_decrease(s, 0); }
    std::vector<char> done(size_t(g.n), 0);
    while (!pq.empty()) {
        auto [du, u] = pq.pop();
        done[size_t(u)] = 1;
        for (const Arc* a = g.begin(u); a != g.end(u); ++a) {
            if (done[size_t(a->to)]) continue;
            int64_t c = rc[size_t(forward ? arc_out(g, u, *a) : arc_in(g, u, *a))];
            if (du + c < d[size_t(a->to)]) { d[size_t(a->to)] = du + c; pq.push_or_decrease(a->to, du + c); }
        }
    }
}

// ============================================================ 3. upper bound
// Shortest-path heuristic from `root` restricted to edges with allowed[e] (all if empty).
static std::vector<int> sph(const Graph& g, int root, const std::vector<char>& allowed, FibPQ& pq) {
    std::vector<char> in_tree(size_t(g.n), 0), done(size_t(g.n));
    std::vector<int64_t> dist(size_t(g.n));
    std::vector<int> pred(size_t(g.n)), edges;
    in_tree[size_t(root)] = 1;
    int left = 0;
    for (int t : g.terminals) left += (t != root);
    while (left > 0) {
        std::fill(dist.begin(), dist.end(), INF); std::fill(pred.begin(), pred.end(), -1); std::fill(done.begin(), done.end(), 0);
        pq.reset(g.n);
        for (int v = 0; v < g.n; ++v) if (in_tree[size_t(v)]) { dist[size_t(v)] = 0; pq.push(v, 0); }
        int hit = -1;
        while (!pq.empty()) {
            int u = pq.pop().second;
            done[size_t(u)] = 1;
            if (g.is_term[size_t(u)] && !in_tree[size_t(u)]) { hit = u; break; }
            for (const Arc* a = g.begin(u); a != g.end(u); ++a) {
                if (done[size_t(a->to)] || (!allowed.empty() && !allowed[size_t(a->e)])) continue;
                if (dist[size_t(u)] + a->w < dist[size_t(a->to)]) {
                    dist[size_t(a->to)] = dist[size_t(u)] + a->w; pred[size_t(a->to)] = a->e;
                    pq.push_or_decrease(a->to, dist[size_t(a->to)]);
                }
            }
        }
        if (hit < 0) return {};
        for (int v = hit; !in_tree[size_t(v)];) {
            if (g.is_term[size_t(v)]) --left;
            in_tree[size_t(v)] = 1;
            int e = pred[size_t(v)];
            edges.push_back(e);
            v = g.E[size_t(e)].u == v ? g.E[size_t(e)].v : g.E[size_t(e)].u;
        }
    }
    return mst_prune(g, edges, pq);
}

// ============================================================ 4-5. search
struct BBNode {
    std::vector<int> incl, excl_v, excl_e;   // decisions relative to the reduced root graph
    int64_t lb;
    int depth;
};

int main(int argc, char** argv) {
    Args a = parse_args(argc, argv, "usage: bnb_steiner <instance.stp> [--time-limit=s] [--node-limit=N] [--no-reduce=1] [--tree=file] [--header]");
    if (a.heap != "fib") { std::fprintf(stderr, "bnb_steiner uses the Fibonacci heap only\n"); return 2; }
    Graph g = read_stp(a.path);
    long long node_limit = a.extra.count("node-limit") ? std::stoll(a.extra.at("node-limit")) : 2000000;
    bool reduce = !(a.extra.count("no-reduce") && a.extra.at("no-reduce") == "1");
    Timer clock;
    FibPQ pq;
    Result res;
    res.algo = reduce ? "BnB" : "BnB_noreduce";

    // ---- 1. reductions
    Reducer red(g);
    if (reduce) red.run(pq);
    std::vector<int> rid;
    Graph g0 = red.to_graph(rid);
    const int n = g0.n;

    auto to_original = [&](const std::vector<int>& g0_edges) {
        std::vector<int> out = red.fixed_orig;
        for (int e : g0_edges) for (int o : red.E[size_t(rid[size_t(e)])].orig) out.push_back(o);
        return out;
    };

    int64_t UB = INF;
    std::vector<int> best;                         // g0 edge ids
    long long nodes = 0, pruned = 0, elim_v_root = 0, elim_e_root = 0;
    int64_t root_lb = -1, global_lb = 0;
    std::string status = "optimal";

    if (g0.k() <= 1) {
        UB = 0; best.clear(); root_lb = 0;
    } else {
        FibHeap open;                               // open nodes, key = lower bound
        std::vector<BBNode> pool;
        pool.push_back({{}, {}, {}, 0, 0});
        open.insert(0, 0);
        std::vector<int64_t> dr, dt;
        while (!open.empty()) {
            if (open.minimum()->key >= UB) { pruned += long(open.size()); break; }
            if (clock.s() > a.time_limit || nodes >= node_limit) {
                status = clock.s() > a.time_limit ? "timeout" : "node_limit";
                global_lb = open.minimum()->key;
                break;
            }
            int id = open.extract_min().second;
            BBNode node = pool[size_t(id)];
            ++nodes;

            // node graph: g0 minus excluded vertices/edges, plus included terminals
            std::vector<char> xv(size_t(n), 0), xe(size_t(g0.m()), 0);
            for (int v : node.excl_v) xv[size_t(v)] = 1;
            for (int e : node.excl_e) xe[size_t(e)] = 1;
            Graph gn; gn.n = n;
            std::vector<int> map_e;                 // gn edge -> g0 edge
            for (int e = 0; e < g0.m(); ++e) {
                const Edge& ed = g0.E[size_t(e)];
                if (xe[size_t(e)] || xv[size_t(ed.u)] || xv[size_t(ed.v)]) continue;
                gn.E.push_back(ed); map_e.push_back(e);
            }
            gn.terminals = g0.terminals;
            for (int v : node.incl) gn.terminals.push_back(v);
            std::sort(gn.terminals.begin(), gn.terminals.end());
            gn.terminals.erase(std::unique(gn.terminals.begin(), gn.terminals.end()), gn.terminals.end());
            gn.build();
            const int root = gn.terminals[0];

            // ---- 2. lower bound
            DualAscent da = dual_ascent(gn, root);
            if (!da.feasible) { ++pruned; continue; }
            int64_t lb = std::max(da.lb, node.lb);
            if (nodes == 1) root_lb = da.lb;
            if (lb >= UB) { ++pruned; continue; }

            // ---- 3. upper bound: SPH on the node graph and on the zero-reduced-cost subgraph
            std::vector<char> zero(size_t(gn.m()), 0);
            for (int e = 0; e < gn.m(); ++e) zero[size_t(e)] = da.rc[2 * size_t(e)] == 0 || da.rc[2 * size_t(e) + 1] == 0;
            for (const auto& allowed : {std::vector<char>(), zero}) {
                std::vector<int> t = sph(gn, root, allowed, pq);
                if (t.empty() && gn.k() > 1) continue;
                int64_t c = cost_of(gn, t);
                if (c < UB) { UB = c; best.clear(); for (int e : t) best.push_back(map_e[size_t(e)]); }
            }
            if (lb >= UB) { ++pruned; continue; }

            // ---- 4. reduced-cost elimination (valid for the whole subtree: UB only falls)
            rc_dijkstra(gn, da.rc, {root}, true, dr, pq);
            rc_dijkstra(gn, da.rc, gn.terminals, false, dt, pq);
            BBNode child_base = node;
            child_base.lb = lb;
            child_base.depth = node.depth + 1;
            for (int v = 0; v < n; ++v) {
                if (gn.is_term[size_t(v)] || gn.deg(v) == 0) continue;
                if (da.lb + dr[size_t(v)] + dt[size_t(v)] >= UB) { child_base.excl_v.push_back(v); xv[size_t(v)] = 1; if (nodes == 1) ++elim_v_root; }
            }
            for (int e = 0; e < gn.m(); ++e) {
                const Edge& ed = gn.E[size_t(e)];
                if (xv[size_t(ed.u)] || xv[size_t(ed.v)]) continue;
                bool f = da.lb + dr[size_t(ed.u)] + da.rc[2 * size_t(e)] + dt[size_t(ed.v)] >= UB;
                bool b = da.lb + dr[size_t(ed.v)] + da.rc[2 * size_t(e) + 1] + dt[size_t(ed.u)] >= UB;
                if (f && b) { child_base.excl_e.push_back(map_e[size_t(e)]); if (nodes == 1) ++elim_e_root; }
            }

            // ---- 5. branching vertex: live non-terminal with the smallest d'(r,v) + d'(v,T)
            int bv = -1; int64_t bscore = INF;
            for (int v = 0; v < n; ++v) {
                if (gn.is_term[size_t(v)] || gn.deg(v) == 0 || xv[size_t(v)]) continue;
                int64_t s = dr[size_t(v)] + dt[size_t(v)];
                if (s < bscore || (s == bscore && gn.deg(v) > gn.deg(bv))) { bscore = s; bv = v; }
            }
            if (bv < 0) continue;      // only terminals left: SPH + MST above is optimal here
            BBNode in = child_base, out = child_base;
            in.incl.push_back(bv);
            out.excl_v.push_back(bv);
            pool.push_back(std::move(in));  open.insert(lb, int(pool.size()) - 1);
            pool.push_back(std::move(out)); open.insert(lb, int(pool.size()) - 1);
            pool[size_t(id)] = BBNode();   // free the parent's lists
        }
        if (status == "optimal") global_lb = UB;
        res.hs = pq.stats();
        res.hs.add(open.stats());
    }

    if (UB >= INF) { res.has_tree = false; res.status = "no_solution"; }
    else res.tree = to_original(best);
    res.status = res.has_tree ? status : res.status;
    res.time_ms = clock.ms();
    int64_t total_ub = UB >= INF ? -1 : UB + red.fixed_cost;
    int64_t total_lb = (root_lb < 0 ? 0 : global_lb) + red.fixed_cost;
    res.notes = "nodes=" + std::to_string(nodes) + ";pruned=" + std::to_string(pruned) +
                ";root_lb=" + std::to_string(root_lb < 0 ? -1 : root_lb + red.fixed_cost) +
                ";final_lb=" + std::to_string(total_lb) + ";ub=" + std::to_string(total_ub) +
                ";red_deg1=" + std::to_string(red.r1) + ";red_term1=" + std::to_string(red.r2) +
                ";red_deg2=" + std::to_string(red.r3) + ";red_long=" + std::to_string(red.r4) +
                ";n_after=" + std::to_string(std::count(red.alive.begin(), red.alive.end(), 1)) +
                ";m_after=" + std::to_string(g0.m()) + ";k_after=" + std::to_string(g0.k()) +
                ";elim_v_root=" + std::to_string(elim_v_root) + ";elim_e_root=" + std::to_string(elim_e_root);
    // safety net: if merged reduction paths overlap, MST + prune on the original graph
    if (res.has_tree && !validate_tree(g, res.tree).ok) res.tree = mst_prune(g, res.tree, pq);
    return report(g, a, res);
}
