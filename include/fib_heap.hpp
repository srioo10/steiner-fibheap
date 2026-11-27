// fib_heap.hpp -- Fibonacci heap (Fredman & Tarjan, J. ACM 1987; CLRS 3e ch. 19).
//
// Every Steiner tree solver in this project runs on this heap:
//   * Dijkstra / Prim inside the three greedy heuristics (DECREASE-KEY heavy),
//   * one Dijkstra per terminal subset in the Dreyfus-Wagner DP,
//   * the open-node priority queue and reduced-cost Dijkstras in branch-and-bound.
//
// Besides the ADT, the heap keeps exact counters of the work it does and of the
// potential Phi(H) = t(H) + 2 m(H)  (t = number of roots, m = number of marked nodes).
// tests/test_fib_heap.cpp uses them to check, operation by operation, the amortized
// bounds proved in the report:
//     INSERT        actual 1               amortized <= 2
//     UNION         actual 1               amortized <= 1
//     DECREASE-KEY  actual 1 + c (c cuts)  amortized <= 5
//     EXTRACT-MIN   actual 2c' + t         amortized <= 3 D(n) + 1   (c' = children of min)
// where D(n) = floor(log_phi n) bounds the degree of any node (report, Lemma 3).
#pragma once
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

struct FibStats {
    uint64_t inserts = 0, extracts = 0, decreases = 0, unions = 0, deletes = 0;
    uint64_t links = 0;              // CONSOLIDATE links (each removes one root)
    uint64_t cuts = 0;               // all cuts (first cut + cascading cuts)
    uint64_t cascading_cuts = 0;     // cuts of an already marked parent
    uint64_t roots_scanned = 0;      // roots visited by CONSOLIDATE
    uint64_t children_promoted = 0;  // children of an extracted min moved to the root list
    uint64_t actual_cost = 0;        // sum of the unit costs defined above
    void add(const FibStats& o) {
        inserts += o.inserts; extracts += o.extracts; decreases += o.decreases;
        unions += o.unions; deletes += o.deletes; links += o.links; cuts += o.cuts;
        cascading_cuts += o.cascading_cuts; roots_scanned += o.roots_scanned;
        children_promoted += o.children_promoted; actual_cost += o.actual_cost;
    }
};

class FibHeap {
public:
    using Key = int64_t;
    struct Node {
        Key key;
        int item;                     // caller's payload (vertex id, B&B node id, ...)
        Node *parent, *child, *left, *right;
        int degree;
        bool mark;
        bool in_heap;
    };
    using Handle = Node*;

    FibHeap() = default;
    FibHeap(const FibHeap&) = delete;
    FibHeap& operator=(const FibHeap&) = delete;

    // ---- MAKE-HEAP is the constructor; clear() empties the heap and recycles memory.
    void clear() {
        min_ = nullptr; n_ = 0; roots_ = 0; marks_ = 0;
        chunk_ = 0; used_ = 0;
    }

    bool empty() const { return n_ == 0; }
    size_t size() const { return n_; }
    int64_t potential() const { return int64_t(roots_) + 2 * int64_t(marks_); }
    size_t roots() const { return roots_; }
    size_t marked() const { return marks_; }
    const FibStats& stats() const { return st_; }
    void reset_stats() { st_ = FibStats(); }

    // ---- INSERT: O(1) worst case. New singleton tree goes into the root list.
    Handle insert(Key key, int item) {
        Node* x = alloc();
        x->key = key; x->item = item;
        x->parent = x->child = nullptr;
        x->degree = 0; x->mark = false; x->in_heap = true;
        x->left = x->right = x;
        add_root(x);
        if (!min_ || key < min_->key) min_ = x;
        ++n_;
        ++st_.inserts; st_.actual_cost += 1;
        return x;
    }

    // ---- MINIMUM: O(1) worst case.
    Handle minimum() const { return min_; }

    // ---- UNION: O(1) worst case. Splices other's root list into ours; other becomes empty.
    void merge(FibHeap& o) {
        if (&o == this || !o.min_) { ++st_.unions; st_.actual_cost += 1; return; }
        if (!min_) {
            min_ = o.min_;
        } else {
            Node* a = min_->right; Node* b = o.min_->left;
            min_->right = o.min_; o.min_->left = min_;
            a->left = b; b->right = a;
            if (o.min_->key < min_->key) min_ = o.min_;
        }
        n_ += o.n_; roots_ += o.roots_; marks_ += o.marks_;
        // take ownership of the other heap's node memory
        for (size_t i = 0; i < o.chunks_.size(); ++i) foreign_.push_back(std::move(o.chunks_[i]));
        for (auto& c : o.foreign_) foreign_.push_back(std::move(c));
        o.chunks_.clear(); o.foreign_.clear();
        o.min_ = nullptr; o.n_ = o.roots_ = o.marks_ = 0; o.chunk_ = o.used_ = 0;
        ++st_.unions; st_.actual_cost += 1;
    }

    // ---- EXTRACT-MIN: O(n) worst case, O(log n) amortized.
    // Returns (key, item) of the removed minimum.
    std::pair<Key, int> extract_min() {
        Node* z = min_;
        assert(z && "extract_min on empty heap");
        // 1. promote every child of z to the root list
        uint64_t c = 0;
        if (Node* ch = z->child) {
            Node* x = ch;
            do {
                Node* next = x->right;
                x->parent = nullptr;
                if (x->mark) { x->mark = false; --marks_; }
                x->left = x->right = x;
                add_root(x);
                ++c;
                x = next;
            } while (x != ch);
            z->child = nullptr;
        }
        // 2. remove z from the root list
        remove_root(z);
        z->in_heap = false;
        --n_;
        uint64_t t_scanned = roots_;  // = t(H) + c - 1
        if (n_ == 0) {
            min_ = nullptr;
        } else {
            min_ = z->right == z ? nullptr : z->right;
            consolidate();
        }
        ++st_.extracts;
        st_.children_promoted += c;
        st_.roots_scanned += t_scanned;
        st_.actual_cost += 1 + c + t_scanned;  // 1 + c + (t + c - 1) = 2c + t
        return {z->key, z->item};
    }

    // ---- DECREASE-KEY: O(log n) worst case (cascade), O(1) amortized.
    void decrease_key(Handle x, Key k) {
        assert(x->in_heap && k <= x->key);
        x->key = k;
        uint64_t cuts = 0;
        Node* y = x->parent;
        if (y && x->key < y->key) {
            cut(x, y); ++cuts;
            // CASCADING-CUT, iterative
            Node* z = y;
            while (Node* p = z->parent) {
                if (!z->mark) { z->mark = true; ++marks_; break; }
                cut(z, p); ++cuts; ++st_.cascading_cuts;
                z = p;
            }
        }
        if (x->key < min_->key) min_ = x;
        ++st_.decreases; st_.cuts += cuts; st_.actual_cost += 1 + cuts;
    }

    // ---- DELETE: DECREASE-KEY to -infinity, then EXTRACT-MIN. O(log n) amortized.
    void erase(Handle x) {
        ++st_.deletes;
        decrease_key_to_neg_inf(x);
        extract_min();
    }

    // Upper bound on node degree, D(n) = floor(log_phi n) (report, Corollary 4).
    static int max_degree_bound(size_t n) {
        if (n <= 1) return 0;
        static const double inv_log_phi = 1.0 / std::log((1.0 + std::sqrt(5.0)) / 2.0);
        return int(std::floor(std::log(double(n)) * inv_log_phi + 1e-9));
    }

    // Full structural check (heap order, degrees, parent links, counters). Used by tests.
    bool check() const {
        if (n_ == 0) return min_ == nullptr && roots_ == 0 && marks_ == 0;
        size_t cnt = 0, rts = 0, mk = 0;
        Node* r = min_;
        do {
            if (r->parent || r->key < min_->key) return false;
            ++rts;
            if (!check_tree(r, cnt, mk)) return false;
            r = r->right;
        } while (r != min_);
        return cnt == n_ && rts == roots_ && mk == marks_;
    }

private:
    static constexpr size_t kChunk = 1024;
    std::vector<std::unique_ptr<Node[]>> chunks_, foreign_;
    size_t chunk_ = 0, used_ = 0;
    Node* min_ = nullptr;
    size_t n_ = 0, roots_ = 0, marks_ = 0;
    FibStats st_;

    Node* alloc() {
        if (chunk_ < chunks_.size() && used_ == kChunk) { ++chunk_; used_ = 0; }
        if (chunk_ == chunks_.size()) { chunks_.emplace_back(new Node[kChunk]); used_ = 0; }
        return &chunks_[chunk_][used_++];
    }
    void add_root(Node* x) {
        if (!min_) { x->left = x->right = x; min_ = x; }
        else { x->right = min_->right; x->left = min_; min_->right->left = x; min_->right = x; }
        ++roots_;
    }
    void remove_root(Node* x) {
        x->left->right = x->right; x->right->left = x->left;
        --roots_;
    }
    // HEAP-LINK: make y a child of x (both roots, x->key <= y->key).
    void link(Node* y, Node* x) {
        remove_root(y);
        y->parent = x;
        if (!x->child) { x->child = y; y->left = y->right = y; }
        else { y->right = x->child->right; y->left = x->child; x->child->right->left = y; x->child->right = y; }
        ++x->degree;
        if (y->mark) { y->mark = false; --marks_; }
        ++st_.links;
    }
    void consolidate() {
        // degree table; D(n) + 2 slots are always enough (Corollary 4)
        int D = max_degree_bound(n_) + 2;
        std::vector<Node*>& A = scratch_;
        A.assign(size_t(D) + 1, nullptr);
        // snapshot the root list (it is modified while linking)
        std::vector<Node*>& rl = roots_buf_;
        rl.clear();
        Node* w = min_;
        do { rl.push_back(w); w = w->right; } while (w != min_);
        for (Node* x : rl) {
            int d = x->degree;
            while (A[size_t(d)]) {
                Node* y = A[size_t(d)];
                if (y->key < x->key || (y->key == x->key && y->item < x->item)) std::swap(x, y);
                link(y, x);
                A[size_t(d)] = nullptr;
                ++d;
                assert(d <= D);
            }
            A[size_t(d)] = x;
        }
        min_ = nullptr;
        for (Node* x : A)
            if (x && (!min_ || x->key < min_->key)) min_ = x;
    }
    void cut(Node* x, Node* y) {
        // remove x from y's child list
        if (x->right == x) y->child = nullptr;
        else {
            if (y->child == x) y->child = x->right;
            x->left->right = x->right; x->right->left = x->left;
        }
        --y->degree;
        x->parent = nullptr;
        if (x->mark) { x->mark = false; --marks_; }
        x->left = x->right = x;
        // add to root list without touching min_ (caller updates it)
        x->right = min_->right; x->left = min_; min_->right->left = x; min_->right = x;
        ++roots_;
    }
    void decrease_key_to_neg_inf(Node* x) {
        x->key = std::numeric_limits<Key>::min();
        uint64_t cuts = 0;
        Node* y = x->parent;
        if (y) {
            cut(x, y); ++cuts;
            Node* z = y;
            while (Node* p = z->parent) {
                if (!z->mark) { z->mark = true; ++marks_; break; }
                cut(z, p); ++cuts; ++st_.cascading_cuts;
                z = p;
            }
        }
        min_ = x;
        ++st_.decreases; st_.cuts += cuts; st_.actual_cost += 1 + cuts;
    }
    bool check_tree(const Node* x, size_t& cnt, size_t& mk) const {
        ++cnt;
        if (x->mark) ++mk;
        if (!x->in_heap) return false;
        int deg = 0;
        if (const Node* c = x->child) {
            const Node* y = c;
            do {
                if (y->parent != x || y->key < x->key) return false;
                if (y->right->left != y) return false;
                ++deg;
                if (!check_tree(y, cnt, mk)) return false;
                y = y->right;
            } while (y != c);
        }
        return deg == x->degree && deg <= max_degree_bound(n_) + 1;
    }
    std::vector<Node*> scratch_, roots_buf_;
};

// Item-indexed priority queue on top of FibHeap: the interface Dijkstra/Prim need.
// push(v,k) / decrease(v,k) / pop() with v in [0, n).
class FibPQ {
public:
    explicit FibPQ(int n = 0) { reset(n); }
    void reset(int n) { h_.clear(); hd_.assign(size_t(n), nullptr); }
    bool empty() const { return h_.empty(); }
    bool contains(int v) const { return hd_[size_t(v)] && hd_[size_t(v)]->in_heap; }
    void push(int v, int64_t k) { hd_[size_t(v)] = h_.insert(k, v); }
    void decrease(int v, int64_t k) { h_.decrease_key(hd_[size_t(v)], k); }
    // insert-or-decrease, the common Dijkstra idiom
    void push_or_decrease(int v, int64_t k) {
        if (contains(v)) { if (k < hd_[size_t(v)]->key) decrease(v, k); }
        else push(v, k);
    }
    std::pair<int64_t, int> pop() { auto r = h_.extract_min(); hd_[size_t(r.second)] = nullptr; return r; }
    const FibStats& stats() const { return h_.stats(); }
    FibHeap& heap() { return h_; }
private:
    FibHeap h_;
    std::vector<FibHeap::Handle> hd_;
};
