#!/usr/bin/env python3
"""Cross-check all five solvers on random small graphs (n <= 14, zero weights, ties, k = 1..n).
The DP and B&B must agree; every greedy tree must be valid and no cheaper than that optimum.
    python3 scripts/fuzz.py [count]
"""
import os, random, subprocess, sys, tempfile

ALGOS = ["greedy1_kmb", "greedy2_sph", "greedy3_zelikovsky", "dp_dreyfus_wagner", "bnb_steiner"]
count = int(sys.argv[1]) if len(sys.argv) > 1 else 300
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
tmp = tempfile.mkdtemp()
bad = 0
for s in range(count):
    r = random.Random(s)
    n = r.randint(2, 14)
    edges = {}
    for v in range(2, n + 1):  # random spanning tree keeps the graph connected
        edges[(r.randint(1, v - 1), v)] = r.choice([0, 1, 1, 2, 3, 5, 8])
    for _ in range(r.randint(0, 2 * n)):
        a, b = sorted(r.sample(range(1, n + 1), 2))
        edges[(a, b)] = r.choice([0, 1, 2, 3, 5, 8])
    T = r.sample(range(1, n + 1), r.randint(1, n))
    p = os.path.join(tmp, f"f{s}.stp")
    with open(p, "w") as f:
        f.write(f"SECTION Graph\nNodes {n}\nEdges {len(edges)}\n")
        f.write("".join(f"E {a} {b} {w}\n" for (a, b), w in edges.items()))
        f.write(f"END\nSECTION Terminals\nTerminals {len(T)}\n" + "".join(f"T {t}\n" for t in T) + "END\nEOF\n")
    res = {}
    for a in ALGOS:
        o = subprocess.run([os.path.join(root, "build", a), p], capture_output=True, text=True, timeout=60)
        row = o.stdout.strip().split(",")
        res[a] = (int(row[5]), row[9]) if o.returncode == 0 else (None, o.stderr[:80])
    opt = res["dp_dreyfus_wagner"][0]
    if opt is None or res["bnb_steiner"][0] != opt or any(c is None or c < opt or v != "1" for c, v in res.values()):
        bad += 1
        print("mismatch", s, res)
print(f"{count - bad}/{count} random instances agree")
sys.exit(1 if bad else 0)
