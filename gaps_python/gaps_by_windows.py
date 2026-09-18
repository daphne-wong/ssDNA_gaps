#!/usr/bin/env python3
"""
dnascent_gaps.py -- call candidate ssDNA gaps as internal analogue-negative
segments inside otherwise continuously labelled nascent DNA, from DNAscent
detect output (v4.x modbam or legacy human-readable .detect).

IMPORTANT: DNAscent does not detect ssDNA gaps. It reports, per thymidine,
P(BrdU) and P(EdU). This script applies an *operational* definition:

    a candidate gap = a run of analogue-NEGATIVE reference positions that is
    (i)  internal to a single read,
    (ii) flanked on BOTH sides by analogue-POSITIVE track of >= --min-flank bp,
    (iii) shorter than --max-gap bp (default 1000),
    (iv) covered by enough thymidine calls that the absence of signal is
         informative rather than a hole in call density.

Everything downstream depends on that definition, so it is exposed as
parameters rather than hard-coded. Validate against an S1-nuclease-minus
control library before interpreting any of it biologically.

Usage
-----
  python dnascent_gaps.py -i detect_output.bam  -o out_prefix          # modbam
  python dnascent_gaps.py -i output.detect      -o out_prefix          # legacy

Outputs
-------
  <prefix>.gaps.bed          BED6 of candidate gaps (chrom, start, end, readID, len, strand)
  <prefix>.gaps.tsv          per-gap detail incl. flank lengths and call density
  <prefix>.per_chrom.tsv     per-chromosome counts and gaps per Mb of interrogated DNA
  <prefix>.per_read.tsv      per-read gap counts and interrogated length
"""

import argparse
import os
import sys
from collections import defaultdict


# ----------------------------------------------------------------------
# Input parsing: yields (read_id, chrom, strand, [(ref_pos, p_analogue), ...])
# ----------------------------------------------------------------------

def parse_detect(path, analogue):
    """Legacy human-readable .detect: header lines '#', read headers
    '>readID contig refStart refEnd strand', then 'pos  pBrdU  pEdU  6mer'."""
    read_id = chrom = strand = None
    calls = []
    with open(path) as fh:
        for line in fh:
            if not line or line[0] == '#':
                continue
            if line[0] == '>':
                if read_id is not None:
                    yield read_id, chrom, strand, calls
                f = line[1:].split()
                raw = f[4] if len(f) > 4 else 'fwd'
                strand = '-' if raw in ('rev', '-') else '+'
                read_id, chrom = f[0], f[1]
                calls = []
                continue
            f = line.split()
            if len(f) < 2:
                continue
            pos = int(f[0])
            p_brdu = float(f[1])
            p_edu = float(f[2]) if len(f) > 3 else 0.0
            calls.append((pos, _combine(p_brdu, p_edu, analogue)))
    if read_id is not None:
        yield read_id, chrom, strand, calls


def parse_modbam(path, analogue, min_mapq, min_read_len):
    """DNAscent v4 modbam. MM codes: 'b' = BrdU (written as 5caU), 'e' = EdU
    (written as 5fU). pysam exposes ML qualities on 0-255; rescaled to 0-1."""
    try:
        import pysam
    except ImportError:
        sys.exit("modbam input needs pysam:  pip install pysam")

    bam = pysam.AlignmentFile(path, 'rb')
    for read in bam.fetch(until_eof=False):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < min_mapq:
            continue
        if read.reference_length is None or read.reference_length < min_read_len:
            continue

        mods = read.modified_bases or {}
        # collect query-position -> probability for each analogue
        q_brdu, q_edu = {}, {}
        for key, vals in mods.items():
            code = str(key[2])
            target = None
            if code in ('b', '5caU', '21983'):
                target = q_brdu
            elif code in ('e', '5fU', '21982'):
                target = q_edu
            if target is None:
                continue
            for qpos, qual in vals:
                target[qpos] = qual / 255.0

        if not q_brdu and not q_edu:
            continue

        # query position -> reference position
        q2r = {q: r for q, r in read.get_aligned_pairs(matches_only=True)}
        calls = []
        for qpos in set(q_brdu) | set(q_edu):
            rpos = q2r.get(qpos)
            if rpos is None:
                continue
            calls.append((rpos, _combine(q_brdu.get(qpos, 0.0),
                                         q_edu.get(qpos, 0.0), analogue)))
        calls.sort()
        yield (read.query_name, read.reference_name,
               '-' if read.is_reverse else '+', calls)
    bam.close()


def _combine(p_brdu, p_edu, analogue):
    if analogue == 'BrdU':
        return p_brdu
    if analogue == 'EdU':
        return p_edu
    return max(p_brdu, p_edu)


# ----------------------------------------------------------------------
# Core: tile the read into windows, binarise, find internal negative runs
# ----------------------------------------------------------------------

def window_read(calls, window, prob_thresh, min_calls, min_frac):
    """Return list of (win_start, win_end, state, n_calls) with
    state in {1 positive, 0 negative, None uninformative}."""
    if not calls:
        return []
    start, end = calls[0][0], calls[-1][0]
    n_win = max(1, (end - start) // window + 1)
    tot = [0] * n_win
    pos = [0] * n_win
    for p, prob in calls:
        i = min((p - start) // window, n_win - 1)
        tot[i] += 1
        if prob >= prob_thresh:
            pos[i] += 1
    out = []
    for i in range(n_win):
        w0 = start + i * window
        w1 = min(w0 + window, end + 1)
        if tot[i] < min_calls:
            state = None
        elif pos[i] / tot[i] >= min_frac:
            state = 1
        else:
            state = 0
        out.append((w0, w1, state, tot[i]))
    return out


def call_gaps(windows, max_gap, min_gap, min_flank, min_calls_in_gap,
              max_uninformative_frac):
    """Find negative runs flanked by positive runs on both sides."""
    if not windows:
        return []
    # collapse consecutive windows into runs of equal state
    runs = []
    cur = [windows[0][0], windows[0][1], windows[0][2], windows[0][3], 1]
    for w0, w1, st, n in windows[1:]:
        if st == cur[2]:
            cur[1] = w1
            cur[3] += n
            cur[4] += 1
        else:
            runs.append(tuple(cur))
            cur = [w0, w1, st, n, 1]
    runs.append(tuple(cur))

    gaps = []
    for i, (s, e, st, n_calls, n_win) in enumerate(runs):
        if st != 0:
            continue
        length = e - s
        if not (min_gap <= length < max_gap):
            continue
        if n_calls < min_calls_in_gap:
            continue
        # flanks: walk outward, allow uninformative windows to be crossed but
        # count them; require a positive run of >= min_flank on each side
        left = _flank(runs, i, -1, min_flank, max_uninformative_frac)
        right = _flank(runs, i, +1, min_flank, max_uninformative_frac)
        if left is None or right is None:
            continue
        gaps.append({'start': s, 'end': e, 'length': length,
                     'n_calls': n_calls, 'n_win': n_win,
                     'left_flank': left, 'right_flank': right,
                     'calls_per_kb': 1000.0 * n_calls / length})
    return gaps


def _flank(runs, i, step, min_flank, max_uninf):
    """Accumulate positive track length adjacent to runs[i] in direction step.
    Returns flank length, or None if the adjacent run is not positive or the
    flank is too short / too uninformative."""
    j = i + step
    acc = 0
    uninf = 0
    seen_positive = False
    while 0 <= j < len(runs):
        s, e, st, n_calls, n_win = runs[j]
        if st == 1:
            acc += e - s
            seen_positive = True
        elif st is None:
            if seen_positive:
                break
            uninf += e - s
            if uninf > max_uninf * min_flank:
                return None
        else:  # another negative run terminates the flank
            break
        if acc >= min_flank:
            return acc
        j += step
    return acc if acc >= min_flank else None


# ----------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-i', '--input', required=True,
                    help='DNAscent detect output: .bam (modbam) or .detect')
    ap.add_argument('-o', '--out-prefix', required=True)
    ap.add_argument('--analogue', choices=['BrdU', 'EdU', 'either'],
                    default='either',
                    help='which analogue defines "labelled" (default: either)')
    ap.add_argument('--max-gap', type=int, default=1000,
                    help='report gaps strictly shorter than this (default 1000)')
    ap.add_argument('--min-gap', type=int, default=100,
                    help='ignore gaps shorter than this (default 100)')
    ap.add_argument('--window', type=int, default=100,
                    help='tiling window in bp (default 100)')
    ap.add_argument('--prob-threshold', type=float, default=0.5,
                    help='per-thymidine analogue probability cutoff (default 0.5)')
    ap.add_argument('--min-frac', type=float, default=0.05,
                    help='fraction of calls above threshold for a window to be '
                         'labelled; set near your measured substitution rate '
                         '(default 0.05)')
    ap.add_argument('--min-calls-per-window', type=int, default=10,
                    help='below this a window is uninformative (default 10)')
    ap.add_argument('--min-calls-in-gap', type=int, default=20,
                    help='minimum thymidine calls inside a gap (default 20)')
    ap.add_argument('--min-flank', type=int, default=1000,
                    help='labelled track required on each side (default 1000)')
    ap.add_argument('--max-uninformative-frac', type=float, default=0.25)
    ap.add_argument('--min-mapq', type=int, default=20, help='modbam only')
    ap.add_argument('--min-read-length', type=int, default=5000,
                    help='minimum aligned reference length (default 5000)')
    ap.add_argument('--chrom-sizes', default=None,
                    help='optional 2-column file for per-chromosome context')
    args = ap.parse_args()

    if args.input.endswith('.bam'):
        source = parse_modbam(args.input, args.analogue,
                              args.min_mapq, args.min_read_length)
    else:
        source = parse_detect(args.input, args.analogue)

    gaps_out = []
    per_read = []
    interrogated = defaultdict(int)   # chrom -> bp of labelled-flanked read span
    chrom_gaps = defaultdict(int)
    len_by_chrom = defaultdict(list)

    for read_id, chrom, strand, calls in source:
        if len(calls) < 100:
            continue
        span = calls[-1][0] - calls[0][0]
        if span < args.min_read_length:
            continue
        wins = window_read(calls, args.window, args.prob_threshold,
                           args.min_calls_per_window, args.min_frac)
        # interrogated denominator: informative windows only
        informative = sum(w1 - w0 for w0, w1, st, _ in wins if st is not None)
        interrogated[chrom] += informative

        g = call_gaps(wins, args.max_gap, args.min_gap, args.min_flank,
                      args.min_calls_in_gap, args.max_uninformative_frac)
        for x in g:
            x['read_id'] = read_id
            x['chrom'] = chrom
            x['strand'] = strand
            gaps_out.append(x)
            chrom_gaps[chrom] += 1
            len_by_chrom[chrom].append(x['length'])
        per_read.append((read_id, chrom, calls[0][0], calls[-1][0],
                         strand, informative, len(g)))

    p = args.out_prefix
    with open(p + '.gaps.bed', 'w') as fh:
        for x in sorted(gaps_out, key=lambda d: (d['chrom'], d['start'])):
            fh.write(f"{x['chrom']}\t{x['start']}\t{x['end']}\t"
                     f"{x['read_id']}\t{x['length']}\t{x['strand']}\n")

    with open(p + '.gaps.tsv', 'w') as fh:
        fh.write('chrom\tstart\tend\tlength\tread_id\tstrand\t'
                 'n_calls\tcalls_per_kb\tleft_flank\tright_flank\n')
        for x in sorted(gaps_out, key=lambda d: (d['chrom'], d['start'])):
            fh.write(f"{x['chrom']}\t{x['start']}\t{x['end']}\t{x['length']}\t"
                     f"{x['read_id']}\t{x['strand']}\t{x['n_calls']}\t"
                     f"{x['calls_per_kb']:.1f}\t{x['left_flank']}\t"
                     f"{x['right_flank']}\n")

    with open(p + '.per_read.tsv', 'w') as fh:
        fh.write('read_id\tchrom\tref_start\tref_end\tstrand\t'
                 'informative_bp\tn_gaps\n')
        for row in per_read:
            fh.write('\t'.join(str(v) for v in row) + '\n')

    with open(p + '.per_chrom.tsv', 'w') as fh:
        fh.write('chrom\tn_gaps\tinterrogated_Mb\tgaps_per_Mb\t'
                 'median_gap_len\tmean_gap_len\n')
        for chrom in sorted(interrogated, key=_chrom_key):
            mb = interrogated[chrom] / 1e6
            n = chrom_gaps.get(chrom, 0)
            lens = sorted(len_by_chrom.get(chrom, []))
            med = lens[len(lens) // 2] if lens else 0
            mean = sum(lens) / len(lens) if lens else 0
            rate = n / mb if mb > 0 else 0
            fh.write(f"{chrom}\t{n}\t{mb:.3f}\t{rate:.2f}\t{med}\t{mean:.0f}\n")

    tot_mb = sum(interrogated.values()) / 1e6
    sys.stderr.write(
        f"reads used: {len(per_read)}\n"
        f"interrogated: {tot_mb:.2f} Mb\n"
        f"candidate gaps <{args.max_gap} bp: {len(gaps_out)}\n"
        f"gap burden: {len(gaps_out)/tot_mb if tot_mb else 0:.2f} per Mb\n"
        f"wrote {p}.gaps.bed / .gaps.tsv / .per_read.tsv / .per_chrom.tsv\n")


def _chrom_key(c):
    s = c[3:] if c.lower().startswith('chr') else c
    return (0, int(s)) if s.isdigit() else (1, s)


if __name__ == '__main__':
    main()
