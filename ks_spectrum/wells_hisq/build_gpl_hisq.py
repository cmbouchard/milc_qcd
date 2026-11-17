"""
build_gpl_hisq.py
==================
HISQ copy of ~/spline_fermions/build_gpl.py (copied 2026-11-17 from the version
there at that date; sha256 of the original in the comment below SOURCE_SHA256).
Reads the B-spline-format files written by scripts/hisq_to_bspline.py,
    {ensemble_id}_c{config}_hisq_m{mass_str}_src{idx}_t{t}x{x}y{y}z{z}_momk{k1}_{k2}_{k3}_pion_corr.txt
and writes the same .gpl layout with keys pion_hisq_m{mass_str}[_src{i}]_momk{k1}_{k2}_{k3}
(or ..._momavg{a}_{b}_{c}). Differences from the original, nothing else changed:
  * filename regex: 'hisq' replaces 'n{n}' and there is no u0 field;
  * key: pion_hisq_m... instead of pion_n{n}_m..._u{u0};
  * SIGN = +1: the original writes -Re C(t) because the B-spline density has a leading minus
    (its correlators are negative); HISQ pion5 x pion5 files are positive and are written as
    +Re C(t), so both gpl files are positive.
SOURCE_SHA256 5e9f6b8fdae0cbbb4149db33b318d5cfa6f03593b60e6053fcdc1804ec96487a
Original docstring follows (it refers to the B-spline names).

Collect pion two-point correlator .txt files (produced by
pion_correlator.py) into a single .gpl file: one row per correlator, 
with columns 

    key  Re(C(0))  Re(C(1))  ...  Re(C(Nt-1))

where `key` identifies which physical quantity a row belongs to
(spline order, mass, momentum).

By default, src0 and src1 for a given configuration are AVERAGED
together into a single row per configuration (101 rows per key, for
the 20x10^3, n=1 ensemble). A configuration missing one of the two
sources is reported and skipped (can't average what isn't there)
rather than silently written as a single-source row.

Pass --keep-src-separate to revert to the old behavior: src0/src1
kept as distinct rows, with src folded into the key (202 rows per key
in that case).

Pass --average-momenta to also average, per configuration (and per
source, with --keep-src-separate), over momenta related by the cubic
group, i.e. with the same sorted |k_i|. Such rows get the key suffix
momavg{a}_{b}_{c} (|k| sorted descending, e.g. momavg1_1_0 for the
average of (1,1,0), (1,0,1), (0,1,1)) instead of momk{k1}_{k2}_{k3}.
Every momentum of a class that appears anywhere in --corr-dir must be
present for a configuration, otherwise it is reported and skipped.

Usage
-----
    python build_gpl.py [--corr-dir ./correlators] [--out correlators.gpl]
                         [--keep-src-separate] [--average-momenta]
                         [--no-header]

Input filenames are expected to match what pion_correlator.py produces:
    {ensemble_id}_{config}_n{n}_m{mass_str}_src{idx}_t{t}x{x}y{y}z{z}_momk{k1}_{k2}_{k3}_pion_corr.txt
e.g.
    link_W_20x10x10x10_b570_10000_n1_m0p2_src0_t5x0y2z7_momk0_0_0_pion_corr.txt
or, for the NERSC ensembles (config written c<config>, E3 replica chain
letter appended to the ensemble_id),
    l2040b8145A_c005_n1_m0p2_src0_t5x0y2z7_momk0_0_0_pion_corr.txt

Rows are grouped by (ensemble_id, config), so E3's chains A/B/C, which
share config numbers, never overwrite each other; within a key they
are written chain by chain, each in config (Monte Carlo) order. Keys
do not include the ensemble, so keep one ensemble per --corr-dir.
"""

import argparse
import glob
import os
import re
from collections import defaultdict

import numpy as np


# Sign applied to Re C(t) when writing rows (original build_gpl.py: -1).
SIGN = +1

_CORR_FILENAME_RE = re.compile(
    r'^(?P<ensemble_id>.+)_'
    r'c?(?P<config>\d+)_'
    r'hisq_'
    r'm(?P<mass_str>[A-Za-z0-9]+)_'
    r'src(?P<src>\d+)_'
    r't(?P<t0>\d+)x(?P<x0>\d+)y(?P<y0>\d+)z(?P<z0>\d+)_'
    r'momk(?P<k1>-?\d+)_(?P<k2>-?\d+)_(?P<k3>-?\d+)_'
    r'pion_corr\.txt$'
)


def parse_corr_filename(filename):
    """
    Parse a pion_correlator.py output filename into its components.
    Returns a dict, or None if the filename doesn't match the expected
    pattern (callers should skip and warn about unparseable files
    rather than crash the whole batch on one oddly-named file).
    """
    stem = os.path.basename(filename)
    m = _CORR_FILENAME_RE.match(stem)
    if m is None:
        return None
    d = m.groupdict()
    return dict(
        ensemble_id=d['ensemble_id'],
        config=int(d['config']),
        mass_str=d['mass_str'],
        src=int(d['src']),
        site=(int(d['t0']), int(d['x0']), int(d['y0']), int(d['z0'])),
        mom=(int(d['k1']), int(d['k2']), int(d['k3'])),
    )

def read_re_c_t(filename):
    """
    Read a pion_correlator.py output file (columns: t, Re(C(t)),
    Im(C(t)), '#'-prefixed header lines) and return Re(C(t)) as an
    array ordered by increasing t (t=0 first).
    """
    data = np.loadtxt(filename, comments='#')
    if data.ndim == 1:
        data = data.reshape(1, -1)   # guard against a single-row file
    t = data[:, 0].astype(int)
    re_c = data[:, 1]
    order = np.argsort(t)
    if not np.array_equal(t[order], np.arange(len(t))):
        raise ValueError(
            f"{filename}: time column isn't a clean 0..Nt-1 sequence "
            f"(got {sorted(t)}) -- can't safely order Re(C(t))"
        )
    return re_c[order]

def build_key(info, include_src):
    """
    Build the .gpl row key from a parsed-filename info dict.
    Deliberately omits config number and site (they vary per row by
    design). Includes src in the key only when include_src=True (ie.
    --keep-src-separate); when averaging sources (the default), src
    is never part of the key since src0/src1 are combined into one
    row per configuration.
    """
    k1, k2, k3 = info['mom']
    parts = [f"pion_hisq_m{info['mass_str']}"]
    if include_src:
        parts.append(f"src{info['src']}")
    parts.append(f"momk{k1}_{k2}_{k3}")
    return "_".join(parts)

def collect_rows(corr_dir, average_sources=True):
    """
    Scan corr_dir for pion_correlator.py output files and group them
    by .gpl key.

    average_sources=True (default): src0 and src1 for a given
        configuration are averaged into a single row. `key` never
        includes src. A configuration that only has one of the two
        sources present is reported via `errors` and left out of
        `rows` -- averaging silently over a missing source would
        quietly change the effective statistics, so it's treated as
        a reportable problem rather than papered over.
    average_sources=False (--keep-src-separate): old behavior, src0/
        src1 kept as distinct rows with src folded into the key.

    Returns:
        rows: dict key -> list of ((ensemble_id, config), src, re_c_array), sorted
              by config (and by src too, when sources are kept
              separate) for deterministic output order. `src` is
              None for averaged rows.
        skipped: list of filenames that didn't match the expected
              naming pattern (reported, not silently dropped)
        errors: list of (filename_or_description, exception) for
              files that matched the naming pattern but failed to
              read/parse cleanly, or (when averaging) configurations
              missing one of the two sources needed to average.
    """
    files = sorted(glob.glob(os.path.join(corr_dir, "*_pion_corr.txt")))
    skipped = []
    errors = []
    if not average_sources:
        rows = defaultdict(list)
        for f in files:
            info = parse_corr_filename(f)
            if info is None:
                skipped.append(f)
                continue
            try:
                re_c = read_re_c_t(f)
            except Exception as e:
                errors.append((f, e))
                continue
            key = build_key(info, include_src=True)
            rows[key].append(((info['ensemble_id'], info['config']),
                              info['src'], SIGN * re_c))
        for key in rows:
            rows[key].sort(key=lambda r: (r[0], r[1]))
        return dict(rows), skipped, errors
    # average_sources=True: group by (key, (ensemble_id, config)), then
    # combine src0/src1
    by_key_config = defaultdict(dict)  # by_key_config[key][(ens, config)] = {src: re_c}
    for f in files:
        info = parse_corr_filename(f)
        if info is None:
            skipped.append(f)
            continue
        try:
            re_c = read_re_c_t(f)
        except Exception as e:
            errors.append((f, e))
            continue
        key = build_key(info, include_src=False)
        by_key_config[key].setdefault((info['ensemble_id'], info['config']),
                                      {})[info['src']] = SIGN * re_c
    rows = {}
    for key, configs in by_key_config.items():
        entries = []
        for config in sorted(configs):
            srcs = configs[config]
            if 0 in srcs and 1 in srcs:
                avg = 0.5 * (srcs[0] + srcs[1])
                entries.append((config, None, avg))
            else:
                have = sorted(srcs.keys())
                errors.append((
                    f"key={key!r} config={config}",
                    ValueError(
                        f"only source(s) {have} present, need both 0 "
                        f"and 1 to average -- configuration skipped"
                    )
                ))
        rows[key] = entries
    return rows, skipped, errors

def momentum_class(mom):
    """Cubic-group class label of a momentum: |k_i| sorted descending."""
    return tuple(sorted((abs(k) for k in mom), reverse=True))

def average_momenta(rows, infos):
    """
    Average rows over cubic-equivalent momenta.

    rows: output of collect_rows(); infos: dict key -> parsed-filename
    info of one file with that key (for its momentum). Rows are
    averaged per (configuration, src). A (configuration, src) missing
    any member of the class (as seen anywhere in the input) is
    reported in `errors` and skipped.

    Returns (rows, errors) in the same format as collect_rows().
    """
    members = defaultdict(list)  # new key -> [old keys]
    for key, info in infos.items():
        a, b, c = momentum_class(info['mom'])
        new_key = re.sub(r'momk-?\d+_-?\d+_-?\d+$',
                         f"momavg{a}_{b}_{c}", key)
        members[new_key].append(key)
    out = {}
    errors = []
    for new_key, old_keys in members.items():
        by_cfg = defaultdict(dict)  # (config, src) -> {old key: re_c}
        for old in old_keys:
            for config, src, re_c in rows[old]:
                by_cfg[(config, src)][old] = re_c
        entries = []
        for (config, src), vals in sorted(by_cfg.items(),
                                          key=lambda kv: (kv[0][0], kv[0][1] or 0)):
            missing = sorted(set(old_keys) - set(vals))
            if missing:
                errors.append((
                    f"key={new_key!r} config={config} src={src}",
                    ValueError(f"missing {missing} -- configuration skipped")
                ))
                continue
            entries.append((config, src,
                            np.mean([vals[k] for k in old_keys], axis=0)))
        out[new_key] = entries
    return out, errors

def write_gpl(rows, out_path, header=True, corr_dir=None, fmt="%.14e"):
    """
    Write the collected rows to a .gpl-style ascii file: one line per
    row, "key val0 val1 ... valN" whitespace-separated. Rows are
    grouped by key (all rows for one key are contiguous), keys written
    in sorted order for a stable, diffable file across regenerations.

    Checks that every row sharing a key has the same number of time
    slices (a real, if unlikely, sign of trouble -- e.g. a mismatched
    lattice shape -- rather than something to silently paper over).
    """
    with open(out_path, "w") as f:
        if header:
            f.write(f"# HISQ pion5 x pion5 correlator data (build_gpl_hisq.py, SIGN={SIGN:+d}), one row per measurement\n")
            f.write(f"# columns: key  Re(C(t=0))  Re(C(t=1))  ...\n")
            if corr_dir:
                f.write(f"# source directory: {os.path.abspath(corr_dir)}\n")
            f.write(f"# {sum(len(v) for v in rows.values())} rows, "
                     f"{len(rows)} distinct keys\n")
        for key in sorted(rows):
            entries = rows[key]
            lengths = {len(re_c) for _, _, re_c in entries}
            if len(lengths) > 1:
                raise ValueError(
                    f"key {key!r}: rows have inconsistent lengths {lengths} "
                    f"-- refusing to write ragged data"
                )
            for config, src, re_c in entries:
                values = " ".join(fmt % v for v in re_c)
                f.write(f"{key} {values}\n")

def main():
    parser = argparse.ArgumentParser(
        description="Collect pion_correlator.py .txt output into a "
                     "single keyed .gpl ascii data file."
    )
    parser.add_argument("--corr-dir", default="./correlators",
                         help="directory containing *_pion_corr.txt "
                              "files (default: ./correlators)")
    parser.add_argument("--out", default="correlators.gpl",
                         help="output .gpl file path (default: "
                              "correlators.gpl)")
    parser.add_argument("--keep-src-separate", action="store_true",
                         help="don't average src0/src1; keep them as "
                              "distinct rows with src folded into the "
                              "key (202 rows per key) -- the default "
                              "instead averages src0/src1 per "
                              "configuration into a single row (101 "
                              "rows per key)")
    parser.add_argument("--average-momenta", action="store_true",
                         help="also average over cubic-equivalent "
                              "momenta (same sorted |k_i|), per "
                              "configuration; keys become "
                              "...momavg{a}_{b}_{c}")
    parser.add_argument("--no-header", action="store_true",
                         help="omit the '#'-prefixed provenance/column "
                              "header lines, for strict .gpl readers "
                              "that don't skip comment lines")
    args = parser.parse_args()

    rows, skipped, errors = collect_rows(args.corr_dir,
                                          average_sources=not args.keep_src_separate)

    if args.average_momenta and rows:
        infos = {}
        for f in glob.glob(os.path.join(args.corr_dir, "*_pion_corr.txt")):
            info = parse_corr_filename(f)
            if info is not None:
                infos.setdefault(build_key(
                    info, include_src=args.keep_src_separate), info)
        rows, mom_errors = average_momenta(rows, infos)
        errors.extend(mom_errors)

    if not rows:
        print(f"No parseable *_pion_corr.txt files found in {args.corr_dir}")
        return

    write_gpl(rows, args.out, header=not args.no_header, corr_dir=args.corr_dir)

    n_rows = sum(len(v) for v in rows.values())
    print(f"Wrote {n_rows} rows across {len(rows)} keys -> {args.out}")

    # Flag any key with fewer rows than the most common count -- likely
    # missing configs for that momentum/source combination, worth
    # knowing about rather than silently having ragged statistics.
    counts = {key: len(v) for key, v in rows.items()}
    if counts:
        expected = max(counts.values())
        incomplete = {k: c for k, c in counts.items() if c != expected}
        if incomplete:
            print(f"\nWARNING: {len(incomplete)} key(s) have fewer rows "
                  f"than the most complete key ({expected} rows):")
            for k, c in sorted(incomplete.items()):
                print(f"  {k}: {c} rows ({expected - c} missing)")
    if skipped:
        print(f"\n{len(skipped)} file(s) in {args.corr_dir} didn't match "
              f"the expected naming pattern and were skipped:")
        for f in skipped[:10]:
            print(f"  {f}")
        if len(skipped) > 10:
            print(f"  ... and {len(skipped) - 10} more")
    if errors:
        print(f"\n{len(errors)} problem(s) encountered (failed to "
              f"read/parse a file, or a configuration missing a "
              f"source needed for averaging):")
        for f, e in errors[:10]:
            print(f"  {f}: {e}")
        if len(errors) > 10:
            print(f"  ... and {len(errors) - 10} more")

if __name__ == "__main__":
    main()
