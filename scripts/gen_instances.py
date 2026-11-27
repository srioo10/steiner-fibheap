#!/usr/bin/env python3
"""Generate the synthetic and adversarial instances (SteinLib STP format).

  instances/adversarial/star_lXX.stp   tight family for G1/G2: ratio -> 2 - 2/l
  instances/adversarial/adv_g{1,2,3}.stp
        smallest instance found by a hill-climbing search that maximises
        cost(greedy)/OPT for each greedy strategy (OPT from the DP)
  instances/synthetic/grid_kXX_sY.stp  10x10 grids, k = 4..20, for DP time vs k
  instances/synthetic/rand_nXXXX.stp   sparse random graphs, n = 100..3200, k = 10

Every generated set gets an optima.csv whose values come from the exact DP
(build/dp_dreyfus_wagner), so the benchmark can report gaps for them too.
Deterministic: all randomness is seeded.

Usage: python3 scripts/gen_instances.py [--skip-search]
"""
import os
import random
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "build")


def write_stp(path, name, n, edges, terms, remark):
    with open(path, "w") as f:
        f.write("33D32945 STP File, STP Format Version 1.0\n\nSECTION Comment\n")
        f.write(f'Name    "{name}"\nCreator "gen_instances.py"\nRemark  "{remark}"\nEND\n\n')
        f.write(f"SECTION Graph\nNodes {n}\nEdges {len(edges)}\n")
        for u, v, w in edges:
            f.write(f"E {u + 1} {v + 1} {w}\n")
        f.write(f"END\n\nSECTION Terminals\nTerminals {len(terms)}\n")
        for t in sorted(terms):
            f.write(f"T {t + 1}\n")
        f.write("END\n\nEOF\n")


def run(binary, path, *extra):
    out = subprocess.run([os.path.join(BIN, binary), path, *extra], capture_output=True, text=True, check=True).stdout
    row = out.strip().split(",")
    return int(row[5])  # cost column


def connected(n, edges):
    adj = [[] for _ in range(n)]
    for u, v, _ in edges:
        adj[u].append(v)
        adj[v].append(u)
    seen, stack = {0}, [0]
    while stack:
        x = stack.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    return len(seen) == n


# ----------------------------------------------------------------- tight star family
def star_family(outdir):
    """Centre c (Steiner) joined to l terminals by weight 10; terminals pairwise
    joined by weight 19. OPT = 10 l (the star). Every terminal pair is at distance
    19 < 20, so G1 and G2 build a path of l-1 clique edges: 19 (l-1).
    Ratio 1.9 (1 - 1/l), approaching the bound 2 - 2/l as 19/10 -> 2."""
    rows = []
    for l in range(3, 13):
        n = l + 1
        edges = [(0, i, 10) for i in range(1, n)]
        edges += [(i, j, 19) for i in range(1, n) for j in range(i + 1, n)]
        name = f"star_l{l:02d}"
        write_stp(os.path.join(outdir, name + ".stp"), name, n, edges, range(1, n),
                  "tight family for KMB/SPH: OPT = 10 l, greedy = 19 (l-1)")
        rows.append((name, n, len(edges), l, 10 * l))
    return rows


# ----------------------------------------------------------------- adversarial search
def random_instance(rng, n, k):
    while True:
        edges = {}
        for v in range(1, n):  # random spanning tree first, so it is connected
            u = rng.randrange(v)
            edges[(u, v)] = rng.randint(1, 20)
        for _ in range(rng.randint(n, 2 * n)):
            u, v = rng.sample(range(n), 2)
            edges[(min(u, v), max(u, v))] = rng.randint(1, 20)
        e = [(u, v, w) for (u, v), w in edges.items()]
        if connected(n, e):
            return e, set(rng.sample(range(n), k))


def ratio(binary, n, edges, terms, tmp):
    write_stp(tmp, "search", n, edges, terms, "search")
    opt = run("dp_dreyfus_wagner", tmp)
    return run(binary, tmp) / opt


def mutate(rng, n, edges, terms):
    e = list(edges)
    t = set(terms)
    r = rng.random()
    if r < 0.6:
        i = rng.randrange(len(e))
        u, v, w = e[i]
        e[i] = (u, v, max(1, min(30, w + rng.choice([-3, -2, -1, 1, 2, 3]))))
    elif r < 0.75 and len(e) > n:
        e.pop(rng.randrange(len(e)))
    elif r < 0.9:
        u, v = sorted(rng.sample(range(n), 2))
        if all((a, b) != (u, v) for a, b, _ in e):
            e.append((u, v, rng.randint(1, 20)))
    else:
        a = rng.choice(sorted(t))
        b = rng.choice([v for v in range(n) if v not in t])
        t.remove(a)
        t.add(b)
    return e, t


def search(binary, seed, outdir, tag, iters=1200, restarts=6):
    """Hill climbing on (weights, edges, terminal set) of small graphs, maximising
    cost(greedy) / OPT. Ties are accepted so the walk can cross plateaus."""
    rng = random.Random(seed)
    fd, tmp = tempfile.mkstemp(suffix=".stp")
    os.close(fd)
    best = (0.0, None)
    for r in range(restarts):
        n, k = rng.choice([(7, 4), (8, 4), (8, 5), (9, 5)])
        e, t = random_instance(rng, n, k)
        cur = ratio(binary, n, e, t, tmp)
        for _ in range(iters):
            e2, t2 = mutate(rng, n, e, t)
            if not connected(n, e2):
                continue
            r2 = ratio(binary, n, e2, t2, tmp)
            if r2 >= cur:
                e, t, cur = e2, t2, r2
        # shrink: drop edges that do not lower the ratio (smaller drawing in the report)
        for edge in sorted(e, key=lambda x: -x[2]):
            e2 = [x for x in e if x != edge]
            if connected(n, e2) and ratio(binary, n, e2, t, tmp) >= cur:
                e = e2
        print(f"  {tag} restart {r}: ratio {cur:.4f} (n={n}, |E|={len(e)}, k={len(t)})", flush=True)
        if cur > best[0]:
            best = (cur, (n, e, t))
    os.remove(tmp)
    n, e, t = best[1]
    name = f"adv_{tag}"
    path = os.path.join(outdir, name + ".stp")
    write_stp(path, name, n, e, t, f"hill-climbing search maximising {tag}/OPT, ratio {best[0]:.4f}")
    return (name, n, len(e), len(t), run("dp_dreyfus_wagner", path))


# ----------------------------------------------------------------- scaling sets
def grid(rng, side, k):
    n = side * side
    edges = []
    for r in range(side):
        for c in range(side):
            v = r * side + c
            if c + 1 < side:
                edges.append((v, v + 1, rng.randint(1, 10)))
            if r + 1 < side:
                edges.append((v, v + side, rng.randint(1, 10)))
    return n, edges, set(rng.sample(range(n), k))


def sparse_random(rng, n, k, avg_deg=4):
    edges = {}
    for v in range(1, n):
        u = rng.randrange(max(0, v - 20), v)  # local spanning tree keeps distances long
        edges[(u, v)] = rng.randint(1, 100)
    while len(edges) < n * avg_deg // 2:
        u = rng.randrange(n)
        v = min(n - 1, u + rng.randint(1, 30))
        if u != v:
            edges[(u, v)] = rng.randint(1, 100)
    return n, [(u, v, w) for (u, v), w in edges.items()], set(rng.sample(range(n), k))


def write_optima(outdir, rows):
    with open(os.path.join(outdir, "optima.csv"), "w") as f:
        f.write("name,n,m,k,opt\n")
        for r in rows:
            f.write(",".join(map(str, r)) + "\n")


def main():
    skip_search = "--skip-search" in sys.argv
    adv = os.path.join(ROOT, "instances", "adversarial")
    syn = os.path.join(ROOT, "instances", "synthetic")
    os.makedirs(adv, exist_ok=True)
    os.makedirs(syn, exist_ok=True)

    rows = star_family(adv)
    if not skip_search:
        for tag, binary, seed in [("g1", "greedy1_kmb", 11), ("g2", "greedy2_sph", 22), ("g3", "greedy3_zelikovsky", 33)]:
            print(f"searching a bad instance for {binary} ...", flush=True)
            rows.append(search(binary, seed, adv, tag))
    else:  # keep previously found gadgets
        old = os.path.join(adv, "optima.csv")
        if os.path.exists(old):
            rows += [tuple(l.strip().split(",")) for l in open(old).readlines()[1:] if l.startswith("adv_")]
    write_optima(adv, rows)

    rng = random.Random(2026)
    rows = []
    for k in range(4, 21, 2):
        for s in range(3):
            n, e, t = grid(rng, 10, k)
            name = f"grid_k{k:02d}_s{s}"
            write_stp(os.path.join(syn, name + ".stp"), name, n, e, t, "10x10 grid, weights 1..10")
            rows.append((name, n, len(e), k, None))
    for n in [100, 200, 400, 800, 1600, 3200]:
        nn, e, t = sparse_random(rng, n, 10)
        name = f"rand_n{n:04d}"
        write_stp(os.path.join(syn, name + ".stp"), name, nn, e, t, "sparse random graph, average degree 4")
        rows.append((name, nn, len(e), 10, None))
    print("computing exact optima for the synthetic set with the DP ...", flush=True)
    rows = [(nm, n, m, k, run("dp_dreyfus_wagner", os.path.join(syn, nm + ".stp"))) for nm, n, m, k, _ in rows]
    write_optima(syn, rows)
    print("done")


if __name__ == "__main__":
    main()
