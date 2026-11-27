# Steiner trees on Fibonacci heaps

ADSA course assignment, IIITDM Kancheepuram.
Team: Sri Sooraj S (CS23I1006), Nishaant Bhaskar (CS23I1002).

Problem: Steiner Tree Problem in graphs (optimisation version).
Data structure: Fibonacci heap, with a full amortized analysis (report, Sections 2-3).
Five algorithms, all running on the heap:

| File | Paradigm | Algorithm | Guarantee |
|---|---|---|---|
| `src/greedy1_kmb.cpp` | greedy | distance-network MST (Kou, Markowsky, Berman 1981) | 2 - 2/l |
| `src/greedy2_sph.cpp` | greedy | shortest-path heuristic, every terminal as root (Takahashi, Matsuyama 1980) | 2 - 2/l |
| `src/greedy3_zelikovsky.cpp` | greedy | best-gain triples (Zelikovsky 1993) | 11/6 |
| `src/dp_dreyfus_wagner.cpp` | dynamic programming | Dreyfus-Wagner with one heap-based Dijkstra per subset | exact, O(3^k n + 2^k (m + n log n)) |
| `src/bnb_steiner.cpp` | branch-and-bound | reductions, Wong dual ascent bound, reduced-cost elimination, best-first search | exact |

Shared code is in `include/`: `fib_heap.hpp` (the heap, with operation counters),
`binary_heap.hpp` (baseline for the heap comparison only), `graph.hpp` (SteinLib reader,
tree validator), `steiner_common.hpp` (Dijkstra and Prim templated on the heap, clean-up step,
CSV output). Only the C++17 standard library is used.

## Build and run

```
make                      # g++ -std=c++17 -O2, binaries in build/
make test                 # heap tests: 10^6 random operations against std::multiset,
                          # amortized bound asserted after every operation
./build/greedy3_zelikovsky instances/steinlib/B/b18.stp --header
./build/dp_dreyfus_wagner  instances/steinlib/B/b13.stp
./build/bnb_steiner        instances/steinlib/C/c08.stp --time-limit=300
```

Each solver prints one CSV row: algorithm, instance, n, m, k, cost, known optimum, gap %,
time (ms), validity, status, heap counters and algorithm-specific notes. The tree is checked by
an independent validator before the row is printed. Common options:

| Option | Meaning |
|---|---|
| `--header` | print the CSV header first |
| `--heap=fib\|bin` | Fibonacci (default) or binary heap (G1, G2, G3, DP) |
| `--tree=FILE` | write the tree's edges (`u v w`, 1-based) |
| `--time-limit=S` | DP and B&B time limit in seconds (default 300) |
| `--max-k=K`, `--mem-gb=G` | DP: skip instances with more terminals / larger tables |
| `--roots=one --root=I` | G2: single root instead of all terminals |
| `--node-limit=N`, `--no-reduce=1` | B&B: node limit, switch reductions off |

## Instances

* `instances/steinlib/` - SteinLib sets B (18), C (20), P6E (15), I080 (100), with the published
  optima in `optima.csv` (source: https://steinlib.zib.de).
* `instances/adversarial/` - the tight star family for G1/G2 and the worst small instances found by
  a hill-climbing search for each greedy strategy.
* `instances/synthetic/` - 10x10 grids with k = 4..20 and sparse random graphs with n = 100..3200.

`python3 scripts/gen_instances.py` regenerates the last two folders (seeded, deterministic);
their optima come from the DP.

## Reproducing the report

```
bash scripts/run_bench.sh          # results/main.csv, results/heap_ablation.csv (about an hour)
python3 scripts/make_tables.py     # report/generated/*: every table, plot and number in the report
make report                        # report/report.pdf
```

`bash scripts/run_bench.sh --small` runs only SteinLib B (about a minute).
