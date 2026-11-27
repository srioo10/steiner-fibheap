// binary_heap.hpp -- indexed binary min-heap with DECREASE-KEY (O(log n) for every op).
// Used only as the baseline in the heap ablation (results/heap_ablation.csv);
// it has the same interface as FibPQ so Dijkstra/Prim/DP are templated on the queue.
#pragma once
#include <cstdint>
#include <utility>
#include <vector>

class BinPQ {
public:
    explicit BinPQ(int n = 0) { reset(n); }
    void reset(int n) { a_.clear(); pos_.assign(size_t(n), -1); }
    bool empty() const { return a_.empty(); }
    bool contains(int v) const { return pos_[size_t(v)] >= 0; }
    void push(int v, int64_t k) {
        pos_[size_t(v)] = int(a_.size());
        a_.push_back({k, v});
        up(a_.size() - 1);
    }
    void decrease(int v, int64_t k) {
        size_t i = size_t(pos_[size_t(v)]);
        a_[i].first = k;
        up(i);
    }
    void push_or_decrease(int v, int64_t k) {
        if (contains(v)) { if (k < a_[size_t(pos_[size_t(v)])].first) decrease(v, k); }
        else push(v, k);
    }
    std::pair<int64_t, int> pop() {
        auto top = a_[0];
        pos_[size_t(top.second)] = -1;
        a_[0] = a_.back();
        a_.pop_back();
        if (!a_.empty()) { pos_[size_t(a_[0].second)] = 0; down(0); }
        return top;
    }
private:
    std::vector<std::pair<int64_t, int>> a_;
    std::vector<int> pos_;
    void place(size_t i) { pos_[size_t(a_[i].second)] = int(i); }
    void up(size_t i) {
        while (i > 0) {
            size_t p = (i - 1) / 2;
            if (a_[p].first <= a_[i].first) break;
            std::swap(a_[p], a_[i]); place(p); place(i); i = p;
        }
    }
    void down(size_t i) {
        for (;;) {
            size_t l = 2 * i + 1, r = l + 1, s = i;
            if (l < a_.size() && a_[l].first < a_[s].first) s = l;
            if (r < a_.size() && a_[r].first < a_[s].first) s = r;
            if (s == i) break;
            std::swap(a_[s], a_[i]); place(s); place(i); i = s;
        }
    }
};
