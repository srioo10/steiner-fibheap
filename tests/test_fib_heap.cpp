// test_fib_heap.cpp -- randomized test of FibHeap against std::multiset, plus an
// operation-by-operation check of the amortized bounds proved in the report.
//
// For every operation i we measure the actual cost c_i (FibStats::actual_cost delta)
// and the potential change Phi(D_i) - Phi(D_{i-1}), and assert
//     c_i + Phi(D_i) - Phi(D_{i-1})  <=  bound(op)
// with bound(INSERT)=2, bound(UNION)=1, bound(DECREASE-KEY)=5,
// bound(EXTRACT-MIN)=3 D(n)+1. Since Phi(D_0)=0 and Phi>=0, summing gives
// sum c_i <= sum bound, the aggregate statement, which is also asserted.
//
// Build: make test     Writes results/amortized_trace.csv (one row per 1000 ops).
#include "../include/fib_heap.hpp"
#include "../include/binary_heap.hpp"
#include <cstdio>
#include <cstdlib>
#include <random>
#include <set>
#include <string>

static int fails = 0;
#define CHECK(c, msg) do { if (!(c)) { ++fails; std::fprintf(stderr, "FAIL %s:%d %s\n", __FILE__, __LINE__, msg); } } while (0)

struct Ref {  // reference: multiset of (key, item)
    std::multiset<std::pair<int64_t, int>> s;
};

static void random_ops(uint64_t seed, int nops, FILE* trace) {
    std::mt19937_64 rng(seed);
    FibHeap h;
    std::multiset<std::pair<int64_t, int>> ref;
    std::vector<FibHeap::Handle> hd;      // handle per item
    std::vector<int64_t> key;             // current key per item
    std::vector<int> live;                // items currently in the heap
    std::vector<int> where;               // index of item in live
    uint64_t sum_actual = 0, sum_bound = 0;
    int64_t max_deg_seen = 0;

    for (int i = 0; i < nops; ++i) {
        int r = int(rng() % 100);
        int64_t phi0 = h.potential();
        uint64_t c0 = h.stats().actual_cost;
        uint64_t bound = 0;
        size_t n_before = h.size();
        if (live.empty() || r < 40) {                       // INSERT
            int it = int(hd.size());
            int64_t k = int64_t(rng() % 1000000);
            hd.push_back(h.insert(k, it)); key.push_back(k);
            where.push_back(int(live.size())); live.push_back(it);
            ref.insert({k, it});
            bound = 2;
        } else if (r < 70) {                                // DECREASE-KEY
            int it = live[rng() % live.size()];
            int64_t k = key[size_t(it)] - int64_t(rng() % 5000);
            ref.erase(ref.find({key[size_t(it)], it}));
            key[size_t(it)] = k; ref.insert({k, it});
            h.decrease_key(hd[size_t(it)], k);
            bound = 5;
        } else if (r < 95) {                                // EXTRACT-MIN
            auto got = h.extract_min();
            auto exp = *ref.begin();
            CHECK(got.first == exp.first, "extract_min returned wrong key");
            ref.erase(ref.find({got.first, got.second}));
            int it = got.second, w = where[size_t(it)];
            live[size_t(w)] = live.back(); where[size_t(live[size_t(w)])] = w; live.pop_back();
            bound = 3 * uint64_t(FibHeap::max_degree_bound(n_before)) + 1;
        } else {                                            // DELETE (decrease + extract)
            int it = live[rng() % live.size()];
            ref.erase(ref.find({key[size_t(it)], it}));
            h.erase(hd[size_t(it)]);
            int w = where[size_t(it)];
            live[size_t(w)] = live.back(); where[size_t(live[size_t(w)])] = w; live.pop_back();
            bound = 5 + 3 * uint64_t(FibHeap::max_degree_bound(n_before)) + 1;
        }
        int64_t phi1 = h.potential();
        uint64_t ci = h.stats().actual_cost - c0;
        int64_t amort = int64_t(ci) + phi1 - phi0;
        CHECK(phi1 >= 0, "negative potential");
        CHECK(amort <= int64_t(bound), "amortized bound violated");
        sum_actual += ci; sum_bound += bound;
        CHECK(h.size() == ref.size(), "size mismatch");
        if (!ref.empty()) CHECK(h.minimum()->key == ref.begin()->first, "minimum mismatch");
        if (i % 4096 == 0) CHECK(h.check(), "structural invariant broken");
        if (trace && i % 1000 == 0)
            std::fprintf(trace, "%llu,%d,%llu,%llu,%lld,%zu\n", (unsigned long long)seed, i,
                         (unsigned long long)sum_actual, (unsigned long long)sum_bound,
                         (long long)phi1, h.size());
        (void)max_deg_seen;
    }
    CHECK(sum_actual <= sum_bound, "aggregate bound violated");
    // drain
    int64_t last = std::numeric_limits<int64_t>::min();
    while (!h.empty()) {
        auto got = h.extract_min();
        CHECK(got.first >= last, "drain not sorted");
        last = got.first;
    }
    CHECK(h.check(), "empty heap invalid");
    std::printf("seed %llu: %d ops, links=%llu cuts=%llu cascading=%llu, sum actual=%llu <= sum bound=%llu\n",
                (unsigned long long)seed, nops, (unsigned long long)h.stats().links,
                (unsigned long long)h.stats().cuts, (unsigned long long)h.stats().cascading_cuts,
                (unsigned long long)sum_actual, (unsigned long long)sum_bound);
}

static void test_union() {
    FibHeap a, b;
    std::multiset<int64_t> ref;
    std::mt19937 rng(7);
    std::vector<FibHeap::Handle> hb;
    for (int i = 0; i < 2000; ++i) { int64_t k = rng() % 10000; a.insert(k, i); ref.insert(k); }
    for (int i = 0; i < 2000; ++i) { int64_t k = rng() % 10000; hb.push_back(b.insert(k, 2000 + i)); ref.insert(k); }
    for (int i = 0; i < 500; ++i) ref.erase(ref.find(b.extract_min().first));  // gives b real trees
    int64_t phi = a.potential() + b.potential();
    a.merge(b);
    CHECK(a.potential() == phi, "union must not change total potential");
    CHECK(b.empty() && a.size() == ref.size() && a.check(), "union broke the heap");
    while (!a.empty()) { auto g = a.extract_min(); CHECK(g.first == *ref.begin(), "union order"); ref.erase(ref.begin()); }
    std::printf("union test: done\n");
}

static void test_binary_heap() {
    std::mt19937 rng(3);
    BinPQ q(5000);
    std::vector<int64_t> k(5000, -1);
    std::multiset<std::pair<int64_t, int>> ref;
    for (int i = 0; i < 5000; ++i) { k[size_t(i)] = rng() % 100000; q.push(i, k[size_t(i)]); ref.insert({k[size_t(i)], i}); }
    for (int i = 0; i < 3000; ++i) {
        int v = int(rng() % 5000);
        if (!q.contains(v)) continue;
        ref.erase({k[size_t(v)], v}); k[size_t(v)] -= rng() % 1000; ref.insert({k[size_t(v)], v});
        q.decrease(v, k[size_t(v)]);
    }
    while (!q.empty()) { auto g = q.pop(); CHECK(g.first == ref.begin()->first, "binary heap order"); ref.erase(ref.begin()); }
    std::printf("binary heap test: done\n");
}

int main(int argc, char** argv) {
    std::string out = argc > 1 ? argv[1] : "results/amortized_trace.csv";
    FILE* trace = std::fopen(out.c_str(), "w");
    if (trace) std::fprintf(trace, "seed,op,sum_actual,sum_amortized_bound,potential,size\n");
    for (uint64_t s = 1; s <= 5; ++s) random_ops(s, 200000, trace);
    if (trace) std::fclose(trace);
    test_union();
    test_binary_heap();
    if (fails) { std::printf("%d checks FAILED\n", fails); return 1; }
    std::printf("all heap tests passed\n");
    return 0;
}
