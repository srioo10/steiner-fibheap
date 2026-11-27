#!/usr/bin/env bash
# Runs the five solvers on every instance and writes
#   results/main.csv            one row per (algorithm, instance), Fibonacci heap
#   results/heap_ablation.csv   Fibonacci vs binary heap, same algorithms, run one at a time
# Options: --small   SteinLib B only (about a minute)
#          --ablation  heap comparison only (needs an idle machine)
#          --jobs N  parallel jobs for the main table (default 6; timings in the
#                    ablation are always measured serially)
# Every tree is checked by the validator in include/graph.hpp before its row is written.
set -euo pipefail
cd "$(dirname "$0")/.."
SMALL=0; JOBS=6; TL=300; ABL_ONLY=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --small) SMALL=1 ;;
    --ablation) ABL_ONLY=1 ;;
    --jobs) JOBS="$2"; shift ;;
    --time-limit) TL="$2"; shift ;;
  esac
  shift
done
make -s all
mkdir -p results

if [[ $SMALL == 1 ]]; then
  SETS="instances/steinlib/B"
else
  SETS="instances/steinlib/B instances/steinlib/C instances/steinlib/P6E instances/steinlib/I080 instances/adversarial instances/synthetic"
fi

HEADER=$(./build/greedy1_kmb instances/steinlib/B/b01.stp --header | head -1)
if [[ $ABL_ONLY == 0 ]]; then
jobs_file=$(mktemp)
for d in $SETS; do
  for f in "$d"/*.stp; do
    echo "greedy1_kmb $f" >> "$jobs_file"
    echo "greedy2_sph $f" >> "$jobs_file"
    echo "greedy3_zelikovsky $f" >> "$jobs_file"
    echo "dp_dreyfus_wagner $f --max-k=20 --time-limit=$TL" >> "$jobs_file"
    echo "bnb_steiner $f --time-limit=$TL" >> "$jobs_file"
  done
done

run_one() {
  # a crash or a hard kill still leaves a row, so nothing disappears from the table
  local algo=$1 file=$2; shift 2
  local out
  if out=$(timeout $(( ${TL%.*} + 120 )) "./build/$algo" "$file" "$@" 2>/dev/null); then
    echo "$out"
  else
    echo "$algo,$(basename "${file%.stp}"),,,,-1,-1,-1,-1,0,crash_or_killed,fib,0,0,0,0,0,0,"
  fi
}
export -f run_one
export TL

echo "running $(wc -l < "$jobs_file") jobs on $JOBS cores ..."
{
  echo "$HEADER"
  xargs -P "$JOBS" -L 1 bash -c 'run_one "$@"' _ < "$jobs_file"
} > results/main.csv.tmp
mv results/main.csv.tmp results/main.csv
rm -f "$jobs_file"
echo "wrote results/main.csv ($(($(wc -l < results/main.csv) - 1)) rows)"
fi

# ---- heap ablation: serial, best of 3 runs, instances where the DP finishes quickly
echo "heap ablation ..."
{
  echo "$HEADER"
  for f in instances/steinlib/B/*.stp instances/steinlib/C/c0[1267].stp instances/steinlib/C/c1[1267].stp instances/steinlib/I080/i080-0[0-4]*.stp instances/synthetic/*.stp; do
    k=$(grep -m1 -i '^Terminals' "$f" | awk '{print $2}' | tr -d '\r')
    for heap in fib bin; do
      for algo in greedy1_kmb greedy2_sph greedy3_zelikovsky dp_dreyfus_wagner; do
        [[ $algo == dp_dreyfus_wagner && $k -gt 16 ]] && continue
        best=""
        for rep in 1 2 3; do
          row=$("./build/$algo" "$f" --heap=$heap)
          t=$(echo "$row" | cut -d, -f9)
          if [[ -z $best ]] || awk "BEGIN{exit !($t < $(echo "$best" | cut -d, -f9))}"; then best=$row; fi
        done
        echo "$best"
      done
    done
  done
} > results/heap_ablation.csv
echo "wrote results/heap_ablation.csv"
