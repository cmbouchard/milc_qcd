#!/usr/bin/env python
"""
hisq_to_bspline.py OUTDIR MILC_CORR... [--sites sites.tsv] [--tag hisq_m0p2]

Convert ks_spectrum_hisq correlator files (save_corr_fnal; written by hisq_chain.sh with
coordinate_origin 0 0 0 0 and parities EO EO EO) into the B-spline pion_correlator.py text
format, one file per (source, momentum):

    {ensemble}_c{NNN}_{tag}_src{i}_t{t}x{x}y{y}z{z}_momk{k1}_{k2}_{k3}_pion_corr.txt
    columns  t  Re(C(t))  Im(C(t))   (t = 0 at the source timeslice), '#' header

Conventions (see runs/hisq_e1_test/README.md):
  MILC      C_M(p,t) = sum_x exp(+i p.x) C(x,t), x relative to the lattice origin
  B-spline  C_B(p,t) = (1/V_s) sum_x exp(-i p.(x - x_src)) C(x,t),  V_s = Nx Ny Nz
  C(x,t) real  =>  C_B = exp(+i p.x_src) conj(C_M) / V_s.
NOT converted: the overall sign (B-spline's density has a leading minus, HISQ pion5 x pion5
is positive) and the staggered taste/normalisation factors. Precision is that of the MILC
file (%e, 7 significant digits). The source index i follows the order of the pairs in the
MILC file (pair 0 = src0 = the B-spline site, pair 1 = src1 = opposite site); with --sites
each source origin is checked against sites.tsv ("config src t x y z").

The filename tag replaces B-spline's "n{n}_m{mass}[_u{u0}]" so build_gpl.py's regex (which
requires n<digits>) does not match these files by design.
"""
import argparse
import os
import re
import sys

import numpy as np


def parse_milc(path):
    """Return dict(L=(nx,ny,nz,nt), blocks=[(origin_xyzt, mom(k1,k2,k3), C complex[nt])])."""
    L = None
    origin = None
    mom = None
    blocks = []
    cur = None
    with open(path) as f:
        for line in f:
            s = line.rstrip('\n')
            if s.startswith('lattice_size:'):
                L = tuple(int(v) for v in s.split(':', 1)[1].split(','))
            elif s.startswith('quark_source_origin:'):
                origin = tuple(int(v) for v in re.findall(r'-?\d+', s.split(':', 1)[1]))
            elif s.startswith('momentum:'):
                lab = s.split(':', 1)[1].strip()
                m = re.fullmatch(r'p(\d)(\d)(\d)', lab)
                if m is None:
                    raise ValueError(f'{path}: unexpected momentum label {lab!r}')
                mom = tuple(int(v) for v in m.groups())
            elif s.startswith('correlator:'):
                name = s.split(':', 1)[1].strip()
                if name != 'PION_5':
                    raise ValueError(f'{path}: unexpected correlator {name!r}')
            elif s == '---':
                cur = None
            elif s == '...':
                if mom is not None and origin is not None:
                    cur = []
                    blocks.append((origin, mom, cur))
                    mom = None
            else:
                parts = s.split('\t')
                if cur is not None and len(parts) == 3:
                    t = int(parts[0])
                    if t != len(cur):
                        raise ValueError(f'{path}: t out of order ({t})')
                    cur.append(complex(float(parts[1]), float(parts[2])))
    if L is None:
        raise ValueError(f'{path}: no lattice_size')
    return dict(L=L, blocks=[(o, k, np.array(c)) for o, k, c in blocks])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('outdir')
    ap.add_argument('files', nargs='+')
    ap.add_argument('--sites', help='sites.tsv to cross-check the source origins')
    ap.add_argument('--tag', default='hisq_m0p2')
    ap.add_argument('--ensemble', default=None,
                    help='ensemble id (default: the directory name above hisq_*/corr/)')
    a = ap.parse_args()

    sites = None
    if a.sites:
        sites = {}
        for line in open(a.sites):
            c, i, t, x, y, z = map(int, line.split())
            sites[(c, i)] = (t, x, y, z)

    os.makedirs(a.outdir, exist_ok=True)
    nwritten = 0
    for path in a.files:
        m = re.fullmatch(r'c(\d+)\.txt', os.path.basename(path))
        if m is None:
            sys.exit(f'{path}: expected cNNN.txt')
        cfg = int(m.group(1))
        ens = a.ensemble or os.path.basename(
            os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(path)))))
        d = parse_milc(path)
        nx, ny, nz, nt = d['L']
        Vs = nx * ny * nz
        # source order = order of first appearance of each origin
        order = []
        for o, _, _ in d['blocks']:
            if o not in order:
                order.append(o)
        if len(d['blocks']) != 8 * len(order):
            sys.exit(f'{path}: {len(d["blocks"])} blocks for {len(order)} sources (want 8 each)')
        for o, k, C in d['blocks']:
            if len(C) != nt:
                sys.exit(f'{path}: {len(C)} timeslices, want {nt}')
            i = order.index(o)
            x, y, z, t = o
            if sites is not None and sites.get((cfg, i)) != (t, x, y, z):
                sys.exit(f'{path}: source {i} origin (t,x,y,z)={(t, x, y, z)} '
                         f'!= sites.tsv {sites.get((cfg, i))}')
            ph = np.exp(2j * np.pi * (k[0] * x / nx + k[1] * y / ny + k[2] * z / nz))
            CB = ph * np.conj(C) / Vs
            name = (f'{ens}_c{cfg:03d}_{a.tag}_src{i}_t{t}x{x}y{y}z{z}'
                    f'_momk{k[0]}_{k[1]}_{k[2]}_pion_corr.txt')
            header = (f"pion two-point correlator, momentum k=(k1,k2,k3)={k}\n"
                      f"HISQ (MILC ks_spectrum_hisq) pion5 x pion5, converted by hisq_to_bspline.py "
                      f"from {os.path.abspath(path)}: exp(+i p.x_src) conj(C_MILC) / {Vs}; "
                      f"sign and taste normalisation NOT converted\n"
                      f"{'t':>4s} {'Re(C(t))':>22s} {'Im(C(t))':>22s}")
            data = np.column_stack([np.arange(nt), CB.real, CB.imag])
            np.savetxt(os.path.join(a.outdir, name), data,
                       fmt=["%4d", "%22.14e", "%22.14e"], header=header)
            nwritten += 1
    print(f'{len(a.files)} MILC files -> {nwritten} correlator files in {a.outdir}')


if __name__ == '__main__':
    main()
