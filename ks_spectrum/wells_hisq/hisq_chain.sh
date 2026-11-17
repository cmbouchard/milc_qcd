#!/bin/bash
# hisq_chain.sh OUTDIR
# HISQ point-to-all propagators + Goldstone-pion 2pt (8 momenta) on Landau-fixed NERSC configs,
# one ks_spectrum_hisq job per config. OUTDIR holds hisq.conf and sites.tsv:
#   hisq.conf: ENS (ensemble dir) HL HT ERR CONFIGS ("1..1000" style: FIRST LAST)
#   sites.tsv: "config src t x y z" (B-spline convention, t first), 2 rows per config.
# Per config: in/cNNN.in, out/cNNN.out, corr/cNNN.txt, done/cNNN; progress.tsv row; RUNLOG.md.
# Restartable (done/ skipped; mkdir lock in lock/). Stops at the first failed config
# (failed/cNNN, STATUS FAILED). STATUS: "STATE hisq_chain TIMESTAMP pid=P [detail]".
set -u
OUT=$(readlink -f "$1")
source "$HOME/gauge_fields/scripts/milc_env.sh"
HISQ_BIN=$HOME/gauge_fields/bin/ks_spectrum_hisq.9275531c
HISQ_COMMIT=9275531c
source "$OUT/hisq.conf"
export HL HT
mkdir -p "$OUT"/{in,out,corr,done,lock,failed}
st() { local s=$1; shift; echo "$s hisq_chain $(ts) pid=$$ $*" > "$OUT/STATUS.tmp"; mv "$OUT/STATUS.tmp" "$OUT/STATUS"; echo "$s hisq_chain $(ts) pid=$$ $*" >> "$OUT/STATUS.history"; }
rl() { echo "$*" >> "$OUT/RUNLOG.md"; }
[ -f "$OUT/progress.tsv" ] || printf "config\tstart\tend\twall_s\tcg_iters\tlinks_s\tprops_s\tcorr_sha256\n" > "$OUT/progress.tsv"
binsha=$(sha256sum "$HISQ_BIN" | cut -c1-16)
rl "## driver start $(ts) pid=$$ host=$(hostname -s) bin=$HISQ_BIN ($binsha..., commit $HISQ_COMMIT) ENS=$ENS L=$HL T=$HT err=$ERR configs $FIRST-$LAST"
st RUNNING "start"
for c in $(seq "$FIRST" "$LAST"); do
  n=$(printf c%03d "$c")
  [ -e "$OUT/done/$n" ] && continue
  [ -e "$OUT/failed/$n" ] && { st FAILED "previous_failure $n"; rl "$(ts) $n has failed/ marker: stop"; exit 1; }
  mkdir "$OUT/lock/$n" 2>/dev/null || continue
  lat=$ENS/gfix/$n.nersc
  s0=$(awk -v c="$c" '$1==c && $2==0 {print $4, $5, $6, $3}' "$OUT/sites.tsv")
  s1=$(awk -v c="$c" '$1==c && $2==1 {print $4, $5, $6, $3}' "$OUT/sites.tsv")
  fail() { touch "$OUT/failed/$n"; rmdir "$OUT/lock/$n"; rl "$(ts) $n FAILED: $*"; st FAILED "$* $n"; exit 1; }
  [ -r "$lat" ] || fail "missing_config"
  [ -n "$s0" ] && [ -n "$s1" ] || fail "missing_sites"
  "$SCRIPTS/hisq_input.sh" "$lat" reload_serial "hisq_${n}" "$s0" "$s1" "$ERR" "$OUT/corr/$n.tmp" > "$OUT/in/$n.in"
  rm -f "$OUT/corr/$n.tmp"
  t0=$(date +%s.%N); start=$(ts)
  nice -n 19 "$HISQ_BIN" < "$OUT/in/$n.in" > "$OUT/out/$n.out" 2>&1; ex=$?
  t1=$(date +%s.%N); end=$(ts)
  o=$OUT/out/$n.out
  [ $ex -eq 0 ] || fail "exit_$ex"
  /usr/bin/grep -q "Archive style checksum = .* OK" "$o" || fail "no_checksum_ok"
  /usr/bin/grep -qi "violation" "$o" && fail "checksum_violation"
  /usr/bin/grep -q "NOT converged" "$o" && fail "not_converged"
  nok=$(/usr/bin/grep -c "OK converged" "$o")
  [ "$nok" -eq 12 ] || fail "converged_lines_$nok"
  nrow=$(awk -F'\t' 'NF==3' "$OUT/corr/$n.tmp" 2>/dev/null | wc -l)
  [ "$nrow" -eq $((16 * HT)) ] || fail "corr_rows_$nrow"
  mv "$OUT/corr/$n.tmp" "$OUT/corr/$n.txt"
  iters=$(/usr/bin/grep -oE "iters= [0-9]+" "$o" | awk '{s+=$2} END{print s}')
  lk=$(awk '/Aggregate time to create fermion links/{print $NF}' "$o")
  pr=$(awk '/Aggregate time to compute propagators/{print $NF}' "$o")
  sha=$(sha256sum "$OUT/corr/$n.txt" | cut -c1-16)
  wall=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.2f", b-a}')
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$n" "$start" "$end" "$wall" "$iters" "$lk" "$pr" "$sha" >> "$OUT/progress.tsv"
  rl "$(ts) $n OK src0 [$s0] src1 [$s1] (MILC x y z t) wall ${wall}s iters $iters corr sha256 $sha..."
  touch "$OUT/done/$n"; rmdir "$OUT/lock/$n"
  st RUNNING "last $n"
done
ndone=$(ls "$OUT/done" | wc -l)
rl "## driver end $(ts): $ndone configs done"
st DONE "ndone=$ndone"
