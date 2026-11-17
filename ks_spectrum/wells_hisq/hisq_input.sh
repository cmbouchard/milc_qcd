#!/bin/bash
# hisq_input.sh LATFILE RELOADCMD JOBID "x0 y0 z0 t0" "x1 y1 z1 t1" [ERR] [CORRFILE]
# ks_spectrum_hisq input for E1 (12^3x24): Goldstone pion (pion5 x pion5) from two point
# sources, am = 0.2, Naik eps 0, 8 momenta (units 2pi/L) projected with full complex phase
# exp(+i p.x), x relative to coordinate_origin 0 0 0 0 (rephase by exp(-i p.x_src) offline
# to get the B-spline convention exp(-i p.(x - x_src)) up to complex conjugation).
lat=$1; reload=$2; jobid=$3; s0=($4); s1=($5); err=${6:-1e-10}; corr=${7:-corr.txt}
moms="000 100 010 001 011 101 110 111"
cat <<EOI
prompt 0
nx ${HL:-12}
ny ${HL:-12}
nz ${HL:-12}
nt ${HT:-24}
iseed 1224757
job_id $jobid

$reload $lat
u0 1.0
no_gauge_fix
forget
staple_weight 0
ape_iter 0
coordinate_origin 0 0 0 0
time_bc antiperiodic

max_number_of_eigenpairs 0

number_of_pbp_masses 0

number_of_base_sources 2

point
field_type KS
subset full
origin ${s0[*]}
source_label p0
forget_source

point
field_type KS
subset full
origin ${s1[*]}
source_label p1
forget_source

number_of_modified_sources 0

number_of_sets 2
EOI
for i in 0 1; do cat <<EOI

set_type single
inv_type UML
max_cg_iterations 1000
max_cg_restarts 10
check yes
momentum_twist 0 0 0
precision 2

source $i

number_of_propagators 1

mass 0.2
naik_term_epsilon 0
error_for_propagator $err
rel_error_for_propagator 0

fresh_ksprop
forget_ksprop
EOI
done
cat <<EOI

number_of_quarks 2

propagator 0
identity
op_label d
forget_ksprop

propagator 1
identity
op_label d
forget_ksprop

number_of_mesons 2
EOI
for i in 0 1; do
  if [ $i = 0 ]; then s=(${s0[*]}); else s=(${s1[*]}); fi
  cat <<EOI

pair $i $i
spectrum_request meson

save_corr_fnal $corr
r_offset ${s[*]}

number_of_correlators 8

EOI
  for m in $moms; do
    echo "correlator PION_5 p$m 1 * 1 pion5 ${m:0:1} ${m:1:1} ${m:2:1} EO EO EO"
  done
done
cat <<EOI

number_of_baryons 0

EOI
