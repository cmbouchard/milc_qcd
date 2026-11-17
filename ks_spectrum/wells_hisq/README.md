# HISQ valence scripts (wells)

Snapshot copies of the working scripts in ~/gauge_fields/scripts on wells (the working copies
there are the ones that run). Used for HISQ point-to-all propagators and Goldstone-pion 2pt
functions on the quenched Luescher-Weisz ensembles, for comparison with B-spline quarks.

- hisq_input.sh: ks_spectrum_hisq input (two point sources, one mass, Naik eps 0, 8 momenta
  000 100 010 001 011 101 110 111 with EO EO EO parities, coordinate_origin 0 0 0 0).
- hisq_chain.sh OUTDIR: one ks_spectrum_hisq job per config (hisq.conf: ENS HL HT ERR FIRST
  LAST; sites.tsv "config src t x y z"); checks exit status, NERSC checksum, CG convergence
  and correlator row count; restartable; stops at the first failure. Needs
  ~/gauge_fields/scripts/milc_env.sh and bin/ks_spectrum_hisq.9275531c.
- hisq_to_bspline.py: converts the MILC correlator files to the B-spline pion_correlator.py
  format: C_B = exp(+i p.x_src) conj(C_MILC) / (Nx Ny Nz). Sign and staggered taste
  normalisation are not converted.
- build_gpl_hisq.py: copy of the B-spline build_gpl.py for these files (keys pion_hisq_m...;
  SIGN +1 where the original writes -Re C).

Build: ks_spectrum/Makefile, `make ks_spectrum_hisq`.
