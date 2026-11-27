#!/usr/bin/env python3
"""Turn results/*.csv into the tables, plot data and numbers used by report/report.tex.

Every number printed in the report comes from here (report/generated/), so the
report cannot drift from the measurements. Run after scripts/run_bench.sh:

    python3 scripts/make_tables.py
"""
import csv
import math
import os
import platform
import statistics as st
import subprocess
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
GEN = os.path.join(ROOT, "report", "generated")
os.makedirs(GEN, exist_ok=True)

ALG = ["G1_KMB", "G2_SPH", "G3_Zelikovsky", "DP_DreyfusWagner", "BnB"]
SHORT = {"G1_KMB": "G1", "G2_SPH": "G2", "G3_Zelikovsky": "G3", "DP_DreyfusWagner": "DP", "BnB": "B\\&B"}
GREEDY = ALG[:3]
PHI = (1 + 5 ** 0.5) / 2

values = {}


def put(key, val):
    values[key] = val


def fmt(x, nd=2):
    if isinstance(x, int):
        return f"{x:,}".replace(",", "{,}")
    return f"{x:.{nd}f}"


def setname(inst):
    i = inst.lower()
    if i.startswith("b") and len(i) == 3:
        return "B"
    if i.startswith("c") and len(i) == 3:
        return "C"
    if i.startswith("p6"):
        return "P6E"
    if i.startswith("i080"):
        return "I080"
    if i.startswith("star") or i.startswith("adv"):
        return "Adversarial"
    if i.startswith("grid"):
        return "Grid"
    if i.startswith("rand"):
        return "Random"
    return "?"


SETS = ["B", "C", "P6E", "I080", "Grid", "Random", "Adversarial"]


def notes(row):
    out = {}
    for kv in row["notes"].split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def load(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        for c in ["n", "m", "k", "cost", "opt", "heap_inserts", "heap_extracts", "heap_decreases", "heap_links", "heap_cuts", "heap_actual"]:
            r[c] = int(r[c]) if r[c] not in ("", None) else -1
        for c in ["gap_pct", "time_ms"]:
            r[c] = float(r[c]) if r[c] not in ("", None) else -1.0
        r["valid"] = r["valid"] == "1"
        r["set"] = setname(r["instance"])
        r["nt"] = notes(r)
    return rows


main = load(os.path.join(RES, "main.csv"))
by = defaultdict(dict)  # by[instance][algo] = row
for r in main:
    by[r["instance"]][r["algo"]] = r
inst_set = {i: setname(i) for i in by}
info = {}
for i, d in by.items():
    any_row = next(iter(d.values()))
    info[i] = (any_row["n"], any_row["m"], any_row["k"], any_row["opt"])

# ---------------------------------------------------------------- sanity checks (hard failures)
bad = [r for r in main if r["status"].startswith("INVALID") or r["status"].startswith("crash") or r["status"].startswith("cost_mismatch") or r["status"] == "reconstruction_failed"]
assert not bad, f"invalid rows: {[(r['algo'], r['instance'], r['status']) for r in bad]}"
for r in main:
    if r["valid"] and r["opt"] > 0:
        assert r["cost"] >= r["opt"], f"{r['algo']} {r['instance']} below optimum"
    if r["algo"] in ("DP_DreyfusWagner", "BnB") and r["status"] == "optimal" or (r["algo"] == "DP_DreyfusWagner" and r["status"] == "ok"):
        assert r["cost"] == r["opt"], f"exact solver mismatch {r['algo']} {r['instance']} {r['cost']} vs {r['opt']}"

n_steinlib = sum(1 for i in by if inst_set[i] in ("B", "C", "P6E", "I080"))
put("inst.steinlib", n_steinlib)
put("inst.generated", len(by) - n_steinlib)
put("inst.total", len(by))
put("rows.total", len(main))


def solved(r):
    return r["status"] in ("optimal", "ok") and r["valid"]


# ---------------------------------------------------------------- environment
cpu = "unknown CPU"
try:
    for line in open("/proc/cpuinfo"):
        if line.startswith("model name"):
            cpu = line.split(":", 1)[1].strip()
            break
except OSError:
    pass
gxx = subprocess.run(["g++", "--version"], capture_output=True, text=True).stdout.splitlines()[0]
mem_gb = 0
try:
    for line in open("/proc/meminfo"):
        if line.startswith("MemTotal"):
            mem_gb = round(int(line.split()[1]) / 1024 / 1024)
except OSError:
    pass
put("env.cpu", cpu.replace("(R)", "").replace("(TM)", ""))
put("env.cores", os.cpu_count())
put("env.mem", mem_gb)
put("env.gxx", gxx.replace("(GCC)", "").replace("_", "\\_").strip())
put("env.os", platform.system() + " " + platform.release().split("-")[0])

# ---------------------------------------------------------------- headline numbers
bb = [r for r in main if r["algo"] == "BnB"]
bb_std = [r for r in bb if r["set"] in ("B", "C", "P6E", "I080")]
put("bb.total", len(bb_std))
put("bb.solved", sum(solved(r) for r in bb_std))
put("bb.solved.C", sum(solved(r) for r in bb_std if r["set"] == "C"))
put("bb.solved.all", sum(solved(r) for r in bb))
put("bb.total.all", len(bb))
put("bb.tl", 300)
dp = [r for r in main if r["algo"] == "DP_DreyfusWagner"]
dp_ok = [r for r in dp if solved(r)]
put("dp.solved", len(dp_ok))
put("dp.attempted", sum(1 for r in dp if r["k"] <= 20))
put("dp.maxtime", fmt(max(r["time_ms"] for r in dp_ok) / 1000, 1))
put("dp.mem.k17", fmt(12 * 2 ** 16 * 100 / 1e6, 1))
put("dp.mem.k20", fmt(12 * 2 ** 19 * 80 / 1e6, 0))
exact_checked = sum(1 for r in main if r["algo"] in ("BnB", "DP_DreyfusWagner") and solved(r) and r["set"] in ("B", "C", "P6E", "I080"))
put("exact.checked", exact_checked)

known = [i for i in by if info[i][3] > 0]
put("gap.count", len(known))
for a, key in zip(GREEDY, ["g1", "g2", "g3"]):
    gaps = [by[i][a]["gap_pct"] for i in known if a in by[i]]
    put(f"{key}.meangap", fmt(st.mean(gaps), 2))
    put(f"{key}.maxgap", fmt(max(gaps), 2))
    put(f"{key}.optcount", sum(1 for g in gaps if g == 0))
r18 = by["b18"]["G3_Zelikovsky"]["nt"]
put("g3.b18.scanned", fmt(int(r18["triples_scanned"])))
put("g3.b18.centres", fmt(int(r18["centres_evaluated"])))

# ---------------------------------------------------------------- amortized analysis
trace = list(csv.DictReader(open(os.path.join(RES, "amortized_trace.csv"))))
with open(os.path.join(GEN, "amort_seed1.csv"), "w") as f:
    f.write("op,sum_actual,sum_amortized_bound,potential\n")
    for t in trace:
        if t["seed"] == "1":
            f.write(f"{t['op']},{t['sum_actual']},{t['sum_amortized_bound']},{t['potential']}\n")
last = [t for t in trace if t["seed"] == "1"][-1]
put("amort.test.ratio", fmt(int(last["sum_actual"]) / int(last["sum_amortized_bound"]), 2))


def dbound(N):
    return math.floor(math.log(max(N, 2)) / math.log(PHI))


rows_h = []
tot_ops = tot_act = tot_bound = 0
for a in ALG:
    rs = [r for r in main if r["algo"] == a and r["heap_actual"] > 0]
    ins = sum(r["heap_inserts"] for r in rs)
    ext = sum(r["heap_extracts"] for r in rs)
    dec = sum(r["heap_decreases"] for r in rs)
    lnk = sum(r["heap_links"] for r in rs)
    cut = sum(r["heap_cuts"] for r in rs)
    act = sum(r["heap_actual"] for r in rs)
    # largest heap a run can hold: n vertices (Dijkstra/Prim) or the open-node queue in B&B
    bound = 0
    for r in rs:
        N = r["n"] + (int(r["nt"].get("nodes", "0")) * 2 if a == "BnB" else 0)
        bound += 2 * r["heap_inserts"] + (3 * dbound(N) + 1) * r["heap_extracts"] + 5 * r["heap_decreases"]
    rows_h.append((SHORT[a], len(rs), ins, ext, dec, lnk, cut, act, bound))
    tot_ops += ins + ext + dec
    tot_act += act
    tot_bound += bound
put("amort.ops", f"{tot_ops / 1e9:.2f}\\,billion" if tot_ops > 1e9 else f"{tot_ops / 1e6:.0f}\\,million")
put("amort.ratio", fmt(tot_act / tot_bound, 2))
with open(os.path.join(GEN, "heapcounts.tex"), "w") as f:
    f.write("\\begin{table}[H]\\centering\\small\n\\begin{tabular}{@{}lrrrrrrrrr@{}}\\toprule\n")
    f.write("Solver & runs & \\textsc{Insert} & \\textsc{Extr.-Min} & \\textsc{Decr.-Key} & links & cuts & actual & bound & ratio\\\\\\midrule\n")
    for s, nr, ins, ext, dec, lnk, cut, act, bound in rows_h:
        f.write(f"{s} & {nr} & {ins/1e6:.1f}M & {ext/1e6:.1f}M & {dec/1e6:.1f}M & {lnk/1e6:.1f}M & {cut/1e6:.2f}M & {act/1e6:.1f}M & {bound/1e6:.1f}M & {act/bound:.2f}\\\\\n")
    f.write("\\bottomrule\\end{tabular}\n")
    f.write("\\caption{Heap work summed over all runs of the main benchmark (M = million). \\emph{actual} is the measured cost in the units of "
            "\\cref{sec:amortized}; \\emph{bound} is \\(2\\,\\#\\textsc{Insert}+(3D(N)+1)\\,\\#\\textsc{Extract-Min}+5\\,\\#\\textsc{Decrease-Key}\\) "
            "from \\cref{thm:sequence}, with \\(N\\) the largest possible heap of the run. For B\\&B the counts cover the vertex queue and the "
            "open-node queue. On these sparse benchmark graphs there are fewer decrease-keys than extract-mins in every solver; "
            "\\cref{sec:ablation} shows what this means for running time.}\\label{tab:heapcounts}\n\\end{table}\n")

# ---------------------------------------------------------------- star family
with open(os.path.join(GEN, "star_table.tex"), "w") as f:
    f.write("\\small\\begin{tabular}{@{}rrrrrrr@{}}\\toprule\n\\(\\ell\\) & OPT & G1 & G2 & G3 & G1/OPT & \\(2-2/\\ell\\)\\\\\\midrule\n")
    for l in range(3, 13):
        i = f"star_l{l:02d}"
        o = info[i][3]
        c = [by[i][a]["cost"] for a in GREEDY]
        f.write(f"{l} & {o} & {c[0]} & {c[1]} & {c[2]} & {c[0]/o:.3f} & {2-2/l:.3f}\\\\\n")
    f.write("\\bottomrule\\end{tabular}\n")
    f.write("\\captionof{table}{Measured costs on the star family.}\\label{tab:star}\n")

# ---------------------------------------------------------------- adversarial gadgets (drawn from the .stp files)
def read_stp(path):
    n, E, T = 0, [], []
    for line in open(path):
        p = line.split()
        if not p:
            continue
        if p[0] == "Nodes":
            n = int(p[1])
        elif p[0] == "E":
            E.append((int(p[1]) - 1, int(p[2]) - 1, int(p[3])))
        elif p[0] == "T":
            T.append(int(p[1]) - 1)
    return n, E, T


def tree_of(binary, path):
    tmp = os.path.join(GEN, "_tree.txt")
    subprocess.run([os.path.join(ROOT, "build", binary), path, f"--tree={tmp}"], capture_output=True, check=True)
    edges = set()
    for line in open(tmp):
        u, v, _ = map(int, line.split())
        edges.add((min(u, v) - 1, max(u, v) - 1))
    os.remove(tmp)
    return edges


gadget_tex = ["\\begin{figure}[H]\\centering"]
for tag, binary, a in [("g1", "greedy1_kmb", "G1_KMB"), ("g2", "greedy2_sph", "G2_SPH"), ("g3", "greedy3_zelikovsky", "G3_Zelikovsky")]:
    path = os.path.join(ROOT, "instances", "adversarial", f"adv_{tag}.stp")
    n, E, T = read_stp(path)
    g_tree = tree_of(binary, path)
    o_tree = tree_of("dp_dreyfus_wagner", path)
    cost, opt = by[f"adv_{tag}"][a]["cost"], info[f"adv_{tag}"][3]
    put(f"adv.{tag}.cost", cost)
    put(f"adv.{tag}.opt", opt)
    put(f"adv.{tag}.ratio", fmt(cost / opt, 3))
    gadget_tex.append("\\begin{subfigure}{0.32\\linewidth}\\centering\\begin{tikzpicture}[scale=1.25]")
    for v in range(n):
        ang = 90 + 360 * v / n
        style = "fill=accent,text=white" if v in T else "fill=white"
        gadget_tex.append(f"\\node[circle,draw,{style},inner sep=1pt,minimum size=4.6mm,font=\\scriptsize] (v{v}) at ({ang:.1f}:1.35) {{{v+1}}};")
    for u, v, w in E:
        key = (min(u, v), max(u, v))
        styles = []
        if key in o_tree:
            styles.append("line width=2.6pt,cG3!45")
        else:
            styles.append("black!30")
        gadget_tex.append(f"\\draw[{styles[0]}] (v{u}) -- node[font=\\tiny,fill=white,inner sep=0.4pt]{{{w}}} (v{v});")
        if key in g_tree:
            gadget_tex.append(f"\\draw[cG1,thick,dashed] (v{u}) -- (v{v});")
    gadget_tex.append("\\end{tikzpicture}")
    gadget_tex.append(f"\\caption{{{SHORT[a]}: cost {cost}, OPT {opt}, ratio {cost/opt:.3f}.}}\\end{{subfigure}}\\hfill")
gadget_tex[-1] = gadget_tex[-1].rstrip("\\hfill")
gadget_tex.append("\\caption{Worst instances found by the search, one per strategy. Filled vertices are terminals. The wide light "
                  "band marks an optimal tree (from the DP); the dashed line is the tree the greedy strategy returns. "
                  "Numbers on edges are weights.}\\label{fig:gadgets}\\end{figure}")
with open(os.path.join(GEN, "gadgets.tex"), "w") as f:
    f.write("\n".join(gadget_tex) + "\n")

# ---------------------------------------------------------------- greedy gap table
with open(os.path.join(GEN, "greedy_gaps.tex"), "w") as f:
    f.write("\\begin{table}[H]\\centering\\small\n\\begin{tabular}{@{}lr" + "rrrr" * 3 + "@{}}\\toprule\n")
    f.write(" & & \\multicolumn{4}{c}{G1 KMB} & \\multicolumn{4}{c}{G2 SPH (best root)} & \\multicolumn{4}{c}{G3 Zelikovsky}\\\\\n")
    f.write("\\cmidrule(lr){3-6}\\cmidrule(lr){7-10}\\cmidrule(l){11-14}\n")
    f.write("Set & \\# & mean & max & opt & ms & mean & max & opt & ms & mean & max & opt & ms\\\\\\midrule\n")
    for s in SETS:
        ins = [i for i in known if inst_set[i] == s]
        if not ins:
            continue
        cells = []
        for a in GREEDY:
            g = [by[i][a]["gap_pct"] for i in ins]
            t = [by[i][a]["time_ms"] for i in ins]
            best_mean = min(st.mean([by[i][b]["gap_pct"] for i in ins]) for b in GREEDY)
            m = st.mean(g)
            mtxt = f"\\textbf{{{m:.2f}}}" if abs(m - best_mean) < 1e-9 else f"{m:.2f}"
            cells += [mtxt, f"{max(g):.2f}", str(sum(1 for x in g if x == 0)), f"{st.mean(t):.1f}"]
            put(f"gap.{s}.{SHORT[a]}", f"{m:.2f}")
        f.write(f"{s} & {len(ins)} & " + " & ".join(cells) + "\\\\\n")
    allg = []
    for a in GREEDY:
        g = [by[i][a]["gap_pct"] for i in known]
        t = [by[i][a]["time_ms"] for i in known]
        allg += [f"{st.mean(g):.2f}", f"{max(g):.2f}", str(sum(1 for x in g if x == 0)), f"{st.mean(t):.1f}"]
    f.write("\\midrule All & " + str(len(known)) + " & " + " & ".join(allg) + "\\\\\n")
    f.write("\\bottomrule\\end{tabular}\n\\caption{Gap to the optimum in percent (mean and maximum), number of instances solved "
            "optimally, and mean time in milliseconds. Bold: smallest mean gap in the row.}\\label{tab:gaps}\n\\end{table}\n")

std_known = [i for i in known if inst_set[i] in ("B", "C", "P6E", "I080")]
put("gap.max.steinlib", fmt(max(by[i][a]["gap_pct"] for i in std_known for a in GREEDY), 1))
for a in GREEDY:
    put(f"time.C.{SHORT[a]}", fmt(st.mean(by[i][a]["time_ms"] for i in known if inst_set[i] == "C") / 1000, 2))

# ---------------------------------------------------------------- SPH root dependence
with open(os.path.join(GEN, "sph_roots.tex"), "w") as f:
    f.write("\\begin{table}[H]\\centering\\small\n\\begin{tabular}{@{}lrrrrr@{}}\\toprule\n")
    f.write("Set & G1 & G2, first terminal & G2, worst root & G2, best root & worst \\(>\\) G1\\\\\\midrule\n")
    worse_total = 0
    for s in SETS:
        ins = [i for i in known if inst_set[i] == s]
        if not ins:
            continue
        g1 = [by[i]["G1_KMB"]["gap_pct"] for i in ins]
        o = {i: info[i][3] for i in ins}
        r0 = [100 * (int(by[i]["G2_SPH"]["nt"]["root0"]) - o[i]) / o[i] for i in ins]
        wr = [100 * (int(by[i]["G2_SPH"]["nt"]["worst_root"]) - o[i]) / o[i] for i in ins]
        bs = [by[i]["G2_SPH"]["gap_pct"] for i in ins]
        worse = sum(1 for i in ins if int(by[i]["G2_SPH"]["nt"]["worst_root"]) > by[i]["G1_KMB"]["cost"])
        worse_total += worse
        f.write(f"{s} & {st.mean(g1):.2f} & {st.mean(r0):.2f} & {st.mean(wr):.2f} & {st.mean(bs):.2f} & {worse}/{len(ins)}\\\\\n")
    f.write("\\bottomrule\\end{tabular}\n\\caption{Mean gap (\\%) of G2 depending on the root terminal, next to G1. The last column counts "
            "instances where the worst root gives a more expensive tree than G1.}\\label{tab:sphroots}\n\\end{table}\n")
    put("sph.worse_than_g1", worse_total)

# ---------------------------------------------------------------- DP results
with open(os.path.join(GEN, "dp_vs_k.csv"), "w") as f:
    f.write("k,grid_s,i080_s,grid_merges\n")
    for k in range(4, 21, 2):
        gt = [by[i]["DP_DreyfusWagner"]["time_ms"] / 1000 for i in by if inst_set[i] == "Grid" and info[i][2] == k and solved(by[i]["DP_DreyfusWagner"])]
        it = [by[i]["DP_DreyfusWagner"]["time_ms"] / 1000 for i in by if inst_set[i] == "I080" and info[i][2] == k and solved(by[i]["DP_DreyfusWagner"])]
        f.write(f"{k},{st.mean(gt) if gt else 'nan'},{st.mean(it) if it else 'nan'},{(3**(k-1)+1)//2 - 2**(k-1)}\n")
with open(os.path.join(GEN, "dp_vs_n.csv"), "w") as f:
    f.write("n,dp_s,bnb_s,g1_s,g3_s\n")
    for i in sorted([i for i in by if inst_set[i] == "Random"], key=lambda x: info[x][0]):
        d = by[i]
        f.write(f"{info[i][0]},{d['DP_DreyfusWagner']['time_ms']/1000},{d['BnB']['time_ms']/1000},{d['G1_KMB']['time_ms']/1000},{d['G3_Zelikovsky']['time_ms']/1000}\n")
# fitted constant: time / (3^(k-1) n) on grids, for the plot's reference line
grid_k20 = [by[i]["DP_DreyfusWagner"]["time_ms"] / 1000 for i in by if inst_set[i] == "Grid" and info[i][2] == 20]
c_fit = st.mean(grid_k20) / (3 ** 19 * 100) if grid_k20 else 1e-9
put("dp.cfit", f"{c_fit:.3e}")
put("dp.grid.k20", fmt(st.mean(grid_k20), 1) if grid_k20 else "--")
with open(os.path.join(GEN, "dp_results.tex"), "w") as f:
    f.write("\\begin{table}[H]\\centering\\small\n\\begin{tabular}{@{}lrrrrr@{}}\\toprule\n")
    f.write("Set & instances & with \\(k\\le20\\) & solved & equal to OPT & max time (s)\\\\\\midrule\n")
    for s in SETS:
        rs = [r for r in dp if r["set"] == s]
        if not rs:
            continue
        att = [r for r in rs if r["k"] <= 20]
        ok = [r for r in att if solved(r)]
        eq = sum(1 for r in ok if r["opt"] > 0 and r["cost"] == r["opt"])
        mt = f"{max(r['time_ms'] for r in ok)/1000:.2f}" if ok else "--"
        f.write(f"{s} & {len(rs)} & {len(att)} & {len(ok)} & {eq} & {mt}\\\\\n")
    f.write("\\bottomrule\\end{tabular}\n\\caption{DP results. Instances with \\(k>20\\) are not attempted (the table would need more than "
            "\\(12\\cdot2^{20}n\\) bytes and \\(3^{20}n\\) merges). Every solved instance matches the known optimum.}\\label{tab:dp}\n\\end{table}\n")
    f.write("""\\begin{figure}[H]\\centering
\\begin{tikzpicture}
\\begin{semilogyaxis}[width=0.48\\linewidth,height=5.4cm,xlabel={terminals \\(k\\)},ylabel={time (s)},
  legend pos=north west,legend style={font=\\scriptsize,draw=none},tick label style={font=\\footnotesize},label style={font=\\small},
  grid=major,grid style={black!8},xtick={4,6,...,20}]
\\addplot[cDP,thick,mark=*,mark size=1.6pt] table[x=k,y=grid_s,col sep=comma]{generated/dp_vs_k.csv};\\addlegendentry{grid, \\(n=100\\)}
\\addplot[cG3,thick,mark=square*,mark size=1.6pt] table[x=k,y=i080_s,col sep=comma]{generated/dp_vs_k.csv};\\addlegendentry{I080, \\(n=80\\)}
\\addplot[black!50,dashed,domain=4:20,samples=17] {""" + f"{c_fit:.6e}" + """*100*3^(x-1)};\\addlegendentry{\\(c\\cdot3^{k-1}n\\)}
\\end{semilogyaxis}
\\end{tikzpicture}\\hfill
\\begin{tikzpicture}
\\begin{loglogaxis}[width=0.48\\linewidth,height=5.4cm,xlabel={vertices \\(n\\) (\\(k=10\\))},ylabel={time (s)},
  legend pos=north west,legend style={font=\\scriptsize,draw=none},tick label style={font=\\footnotesize},label style={font=\\small},
  grid=major,grid style={black!8}]
\\addplot[cDP,thick,mark=*,mark size=1.6pt] table[x=n,y=dp_s,col sep=comma]{generated/dp_vs_n.csv};\\addlegendentry{DP}
\\addplot[cBB,thick,mark=triangle*,mark size=1.8pt] table[x=n,y=bnb_s,col sep=comma]{generated/dp_vs_n.csv};\\addlegendentry{B\\&B}
\\addplot[cG1,thick,mark=o,mark size=1.6pt] table[x=n,y=g1_s,col sep=comma]{generated/dp_vs_n.csv};\\addlegendentry{G1}
\\addplot[cG3,thick,mark=square,mark size=1.6pt] table[x=n,y=g3_s,col sep=comma]{generated/dp_vs_n.csv};\\addlegendentry{G3}
\\end{loglogaxis}
\\end{tikzpicture}
\\caption{Left: DP time against \\(k\\) (mean of the instances per \\(k\\)); the dashed line is \\(3^{k-1}n\\) scaled to the \\(k=20\\) grids,
and the measured curve follows it. Right: time against \\(n\\) on sparse random graphs with \\(k=10\\).}\\label{fig:dpscale}
\\end{figure}
""")

# ---------------------------------------------------------------- B&B results
with open(os.path.join(GEN, "bnb_results.tex"), "w") as f:
    f.write("\\begin{table}[H]\\centering\\small\\setlength{\\tabcolsep}{4pt}\n\\begin{tabular}{@{}lrrrrrrrrr@{}}\\toprule\n")
    f.write("Set & \\# & solved & \\(|V|\\) left & \\(|E|\\) left & root gap \\% & nodes: median & max & mean s & max s\\\\\\midrule\n")
    for s in SETS:
        rs = [r for r in bb if r["set"] == s]
        if not rs:
            continue
        ok = [r for r in rs if solved(r)]
        vleft = st.mean(int(r["nt"]["n_after"]) / r["n"] for r in rs) * 100
        eleft = st.mean(int(r["nt"]["m_after"]) / r["m"] for r in rs) * 100
        gaps = [100 * (r["opt"] - int(r["nt"]["root_lb"])) / r["opt"] for r in rs if r["opt"] > 0 and int(r["nt"]["root_lb"]) >= 0]
        nodes = [int(r["nt"]["nodes"]) for r in rs]
        t = [r["time_ms"] / 1000 for r in rs]
        f.write(f"{s} & {len(rs)} & {len(ok)} & {vleft:.0f}\\% & {eleft:.0f}\\% & {st.mean(gaps):.2f} & {int(st.median(nodes))} & {fmt(max(nodes))} & {st.mean(t):.2f} & {max(t):.2f}\\\\\n")
    f.write("\\bottomrule\\end{tabular}\n\\caption{Branch-and-bound with a 300\\,s limit. \\(|V|\\), \\(|E|\\) left: share of the graph that survives the "
            "reductions. Root gap: \\((\\OPT-\\mathrm{LB}_{\\text{root}})/\\OPT\\) for the dual ascent bound at the root.}\\label{tab:bnb}\n\\end{table}\n")
    uns = [r for r in bb if not solved(r)]
    put("bb.unsolved", len(uns))
    if uns:
        f.write("\\noindent Instances not solved within the limit, with the bounds proved when the search stopped:\n\\begin{center}\\small\\begin{tabular}{@{}lrrrrr@{}}\\toprule\n")
        f.write("Instance & \\(k\\) & lower bound & best tree & OPT & nodes\\\\\\midrule\n")
        for r in sorted(uns, key=lambda x: x["instance"]):
            f.write(f"{r['instance']} & {r['k']} & {r['nt']['final_lb']} & {r['cost']} & {r['opt'] if r['opt'] > 0 else '?'} & {fmt(int(r['nt']['nodes']))}\\\\\n")
        f.write("\\bottomrule\\end{tabular}\\end{center}\n")
    # DP vs B&B on instances both solved
    both = [i for i in by if "DP_DreyfusWagner" in by[i] and "BnB" in by[i] and solved(by[i]["DP_DreyfusWagner"]) and solved(by[i]["BnB"])]
    with open(os.path.join(GEN, "dp_vs_bb.csv"), "w") as g:
        g.write("k,dp,bb\n")
        for i in both:
            g.write(f"{info[i][2]},{max(by[i]['DP_DreyfusWagner']['time_ms'],0.01)/1000},{max(by[i]['BnB']['time_ms'],0.01)/1000}\n")
    bb_faster = sum(1 for i in both if by[i]["BnB"]["time_ms"] < by[i]["DP_DreyfusWagner"]["time_ms"])
    put("both.count", len(both))
    put("both.bbfaster", bb_faster)
    f.write("""\\begin{figure}[H]\\centering
\\begin{tikzpicture}
\\begin{loglogaxis}[width=0.42\\linewidth,height=5.6cm,xlabel={DP time (s)},ylabel={B\\&B time (s)},
  colorbar,colormap={kk}{color=(cG3) color=(cG1)},point meta=explicit,colorbar style={ylabel={\\(k\\)},font=\\footnotesize,width=2.5mm},
  tick label style={font=\\footnotesize},label style={font=\\small},grid=major,grid style={black!8}]
\\addplot[scatter,only marks,mark=*,mark size=1.5pt,fill opacity=0.7,draw opacity=0] table[x=dp,y=bb,meta=k,col sep=comma]{generated/dp_vs_bb.csv};
\\addplot[black!50,dashed,domain=1e-5:200] {x};
\\end{loglogaxis}
\\end{tikzpicture}
\\caption{The two exact methods on the """ + str(len(both)) + """ instances both solve. Below the diagonal B\\&B is faster; it is faster on """ + str(bb_faster) + """ of them.
Colour is the number of terminals: the DP time is set by \\(k\\), the B\\&B time by how much of the graph the bound settles.}\\label{fig:dpbb}
\\end{figure}
""")

# ---------------------------------------------------------------- heap ablation
abl = load(os.path.join(RES, "heap_ablation.csv"))
pair = defaultdict(dict)
for r in abl:
    pair[(r["algo"], r["instance"])][r["heap"]] = r
with open(os.path.join(GEN, "ablation.tex"), "w") as f:
    f.write("\\begin{table}[H]\\centering\\small\n\\begin{tabular}{@{}llrrrrr@{}}\\toprule\n")
    f.write("Solver & graphs & instances & Fibonacci (s) & binary (s) & Fib.\\,/\\,bin. & Fib. faster\\\\\\midrule\n")
    for a in ["G1_KMB", "G2_SPH", "G3_Zelikovsky", "DP_DreyfusWagner"]:
        for cls, test in [("sparse, \\(m/n<10\\)", lambda r: r["m"] < 10 * r["n"]), ("dense, \\(m/n\\ge10\\)", lambda r: r["m"] >= 10 * r["n"])]:
            ps = [p for (al, _), p in pair.items() if al == a and "fib" in p and "bin" in p and test(p["fib"])]
            if not ps:
                continue
            tf = sum(p["fib"]["time_ms"] for p in ps) / 1000
            tb = sum(p["bin"]["time_ms"] for p in ps) / 1000
            if a == "DP_DreyfusWagner":
                assert all(p["fib"]["cost"] == p["bin"]["cost"] for p in ps), "heap changed the DP optimum"
            faster = sum(1 for p in ps if p["fib"]["time_ms"] < p["bin"]["time_ms"])
            f.write(f"{SHORT[a]} & {cls} & {len(ps)} & {tf:.3f} & {tb:.3f} & {tf/tb:.2f} & {faster}/{len(ps)}\\\\\n")
            put(f"abl.{SHORT[a]}.{'dense' if 'dense' in cls else 'sparse'}.ratio", f"{tf/tb:.2f}")
    f.write("\\bottomrule\\end{tabular}\n\\caption{Total time with each heap (best of three serial runs per instance). Instances: SteinLib B, c01, c02, c06, c07, c11, c12, c16, c17, "
            "the 25 I080 instances i080-001 to i080-045, and the grids and random graphs (the DP only for \\(k\\le16\\)). The dense group is c16, c17 (\\(m/n=25\\)) and the "
            "five I080 graphs with \\(m/n=39.5\\). The DP returns the same optimum with both heaps; the greedy strategies can break ties between equal keys "
            "differently.}\\label{tab:ablation}\n\\end{table}\n")
    for inst in ["c17", "c12", "c16"]:
        p = pair.get(("DP_DreyfusWagner", inst), {})
        if "fib" in p and "bin" in p:
            put(f"abl.dp.{inst}.fib", fmt(p["fib"]["time_ms"], 1))
            put(f"abl.dp.{inst}.bin", fmt(p["bin"]["time_ms"], 1))
            put(f"abl.dp.{inst}.decr", fmt(p["fib"]["heap_decreases"]))
            put(f"abl.dp.{inst}.extr", fmt(p["fib"]["heap_extracts"]))
    dp_pairs = [p for (al, _), p in pair.items() if al == "DP_DreyfusWagner"]
    big = max(dp_pairs, key=lambda p: p["fib"]["time_ms"])
    put("abl.dp.big.inst", big["fib"]["instance"])
    put("abl.dp.big.fib", fmt(big["fib"]["time_ms"] / 1000, 2))
    put("abl.dp.big.bin", fmt(big["bin"]["time_ms"] / 1000, 2))
    put("abl.dp.big.decr", fmt(big["fib"]["heap_decreases"]))
    put("abl.dp.big.extr", fmt(big["fib"]["heap_extracts"]))

# ---------------------------------------------------------------- appendix tables
def cost_cell(r):
    if r is None:
        return "--"
    if not r["valid"]:
        return "--"
    c = str(r["cost"])
    if r["opt"] > 0 and r["cost"] == r["opt"]:
        c = f"\\textbf{{{c}}}"
    if r["algo"] in ("BnB",) and r["status"] != "optimal":
        c += "\\textsuperscript{t}"
    return c


def time_cell(r):
    if r is None or not r["valid"]:
        return "--"
    t = r["time_ms"] / 1000
    return f"{t:.3f}" if t < 10 else f"{t:.1f}"


with open(os.path.join(GEN, "appendix_tables.tex"), "w") as f:
    for s in ["B", "C", "P6E"]:
        ins = sorted([i for i in by if inst_set[i] == s])
        f.write(f"\\subsection*{{Set {s}}}\n\\begingroup\\scriptsize\\setlength{{\\tabcolsep}}{{3.2pt}}\n")
        f.write("\\begin{longtable}{@{}lrrrr" + "rr" * 5 + "@{}}\\toprule\n")
        f.write("Inst. & \\(n\\) & \\(m\\) & \\(k\\) & OPT & \\multicolumn{2}{c}{G1} & \\multicolumn{2}{c}{G2} & \\multicolumn{2}{c}{G3} & \\multicolumn{2}{c}{DP} & \\multicolumn{2}{c}{B\\&B}\\\\\n")
        f.write(" & & & & & cost & s & cost & s & cost & s & cost & s & cost & s\\\\\\midrule\\endhead\n")
        for i in ins:
            n, m, k, o = info[i]
            cells = []
            for a in ALG:
                r = by[i].get(a)
                cells += [cost_cell(r), time_cell(r)]
            f.write(f"{i.replace('_', '\\_')} & {n} & {m} & {k} & {o} & " + " & ".join(cells) + "\\\\\n")
        f.write("\\bottomrule\\end{longtable}\\endgroup\n")
    # I080 grouped by k
    f.write("\\subsection*{Set I080 (100 instances, grouped)}\n\\begin{table}[H]\\centering\\small\\begin{tabular}{@{}rr" + "rr" * 5 + "@{}}\\toprule\n")
    f.write("\\(k\\) & \\# & \\multicolumn{2}{c}{G1} & \\multicolumn{2}{c}{G2} & \\multicolumn{2}{c}{G3} & \\multicolumn{2}{c}{DP} & \\multicolumn{2}{c}{B\\&B}\\\\\n")
    f.write(" & & gap \\% & opt & gap \\% & opt & gap \\% & opt & solved & max s & solved & max s\\\\\\midrule\n")
    for k in sorted({info[i][2] for i in by if inst_set[i] == "I080"}):
        ins = [i for i in by if inst_set[i] == "I080" and info[i][2] == k]
        cells = []
        for a in GREEDY:
            g = [by[i][a]["gap_pct"] for i in ins]
            cells += [f"{st.mean(g):.2f}", str(sum(1 for x in g if x == 0))]
        for a in ["DP_DreyfusWagner", "BnB"]:
            ok = [by[i][a] for i in ins if solved(by[i][a])]
            cells += [str(len(ok)), f"{max(r['time_ms'] for r in ok)/1000:.2f}" if ok else "--"]
        f.write(f"{k} & {len(ins)} & " + " & ".join(cells) + "\\\\\n")
    f.write("\\bottomrule\\end{tabular}\\end{table}\n")

    for s in ["Adversarial", "Random", "Grid"]:
        ins = sorted([i for i in by if inst_set[i] == s])
        f.write(f"\\subsection*{{Set {s}}}\n\\begingroup\\scriptsize\\setlength{{\\tabcolsep}}{{3.2pt}}\n")
        f.write("\\begin{longtable}{@{}lrrrr" + "rr" * 5 + "@{}}\\toprule\n")
        f.write("Inst. & \\(n\\) & \\(m\\) & \\(k\\) & OPT & \\multicolumn{2}{c}{G1} & \\multicolumn{2}{c}{G2} & \\multicolumn{2}{c}{G3} & \\multicolumn{2}{c}{DP} & \\multicolumn{2}{c}{B\\&B}\\\\\n")
        f.write(" & & & & & cost & s & cost & s & cost & s & cost & s & cost & s\\\\\\midrule\\endhead\n")
        for i in ins:
            n, m, k, o = info[i]
            cells = []
            for a in ALG:
                r = by[i].get(a)
                cells += [cost_cell(r), time_cell(r)]
            f.write(f"{i.replace('_', '\\_')} & {n} & {m} & {k} & {o} & " + " & ".join(cells) + "\\\\\n")
        f.write("\\bottomrule\\end{longtable}\\endgroup\n")

# ---------------------------------------------------------------- overview table (one row per algorithm)
with open(os.path.join(GEN, "overview.tex"), "w") as f:
    f.write("\\begin{table}[H]\\centering\\small\n\\begin{tabular}{@{}llrrrr@{}}\\toprule\n")
    f.write("Solver & Paradigm & runs & optimal & mean gap (\\%) & total time (s)\\\\\\midrule\n")
    para = {"G1_KMB": "greedy", "G2_SPH": "greedy", "G3_Zelikovsky": "greedy", "DP_DreyfusWagner": "dynamic programming", "BnB": "branch-and-bound"}
    for a in ALG:
        rs = [by[i][a] for i in known if a in by[i]]
        valid = [r for r in rs if r["valid"]]
        opt = sum(1 for r in valid if r["cost"] == r["opt"])
        gaps = [r["gap_pct"] for r in valid]
        tt = sum(r["time_ms"] for r in rs) / 1000
        f.write(f"{SHORT[a]} & {para[a]} & {len(valid)} & {opt} & {st.mean(gaps):.3f} & {tt:.1f}\\\\\n")
    f.write("\\bottomrule\\end{tabular}\n\\caption{All instances with a known optimum (" + str(len(known)) + "). \\emph{runs}: runs that returned a "
            "tree (the DP skips \\(k>20\\)). For the exact methods a non-zero mean gap comes only from instances stopped at the time limit.}\\label{tab:overview}\n\\end{table}\n")

# ---------------------------------------------------------------- values.tex
with open(os.path.join(GEN, "values.tex"), "w") as f:
    f.write("% generated by scripts/make_tables.py, do not edit\n")
    for k, v in sorted(values.items()):
        f.write(f"\\expandafter\\def\\csname v@{k}\\endcsname{{{v}}}\n")
print(f"wrote {len(values)} values and tables to {GEN}")
