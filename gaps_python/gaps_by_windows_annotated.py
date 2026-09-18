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
from collections import defaultdict         # Automatically assigns default value to keys that do not exist (no need to manually check for missing keys; avoid KeyError)


# ----------------------------------------------------------------------
# Input parsing: yields (read_id, chrom, strand, [(ref_pos, p_analogue), ...])
# ----------------------------------------------------------------------

def parse_detect(path, analogue):
    """Legacy human-readable .detect: header lines '#', read headers
    '>readID contig refStart refEnd strand', then 'pos  pBrdU  pEdU  6mer'."""
    read_id = chrom = strand = None
    calls = []                                                              # Holds ref_pos + prob pairs for the calls within each read
    with open(path) as fh:
        for line in fh:                                                     # For each read,
            if not line or line[0] == '#':                                  # Skip line if comment '#''
                continue
            if line[0] == '>':                                              # If line is read header,
                if read_id is not None:
                    yield read_id, chrom, strand, calls                     # Emit completed previous read to calling loop
                f = line[1:].split()                                        # Split header into columns
                raw = f[4] if len(f) > 4 else 'fwd'                         # Read 5th value of header as strand. Default = fwd
                strand = '-' if raw in ('rev', '-') else '+'                # Convert fwd/rev to +/-
                read_id, chrom = f[0], f[1]                                 # Extract read_id and chr from header
                calls = []                                                  # Start new calls list for this read
                continue
            f = line.split()                                                # Split non-header line into fields
            if len(f) < 2:                                                  # Check that at least 1 position and BrdU prob exists
                continue
            pos = int(f[0])                                                 # Convert first field into int(reference coord)
            p_brdu = float(f[1])                                            # Convert second field to BrdU probability
            p_edu = float(f[2]) if len(f) > 3 else 0.0                      # Convert third field to EdU, ONLY if ther eis at least four fields. Otherwise EdU prob = 0
            calls.append((pos, _combine(p_brdu, p_edu, analogue)))          # Call _combine() to choose betwween BrdU, EdU, or max, and stores the call position + prob in calls[]
    if read_id is not None:                                                 # Checks whether final read exists
        yield read_id, chrom, strand, calls                                 # Emits last read


def parse_modbam(path, analogue, min_mapq, min_read_len):
    """DNAscent v4 modbam. MM codes: 'b' = BrdU (written as 5caU), 'e' = EdU
    (written as 5fU). pysam exposes ML qualities on 0-255; rescaled to 0-1."""
    try:
        import pysam
    except ImportError:
        sys.exit("modbam input needs pysam:  pip install pysam")

    bam = pysam.AlignmentFile(path, 'rb')                                               # Open BAM in binary read mode
    for read in bam.fetch(until_eof=False):                                             # Iterate over alignments in each BAM until eof. Requires .bam.bai
        if read.is_unmapped or read.is_secondary or read.is_supplementary:              # Filters out unmapped, secondary, supplementary reads
            continue
        if read.mapping_quality < min_mapq:                                             # Filters out reads that don't meet map_quality thresh
            continue
        if read.reference_length is None or read.reference_length < min_read_len:       # Filters out reads that don't meet read_length thresh
            continue

        mods = read.modified_bases or {}                                                # Gets decoded mod-base annotations from pysam. If none, use empty dict {}
        # collect query-position -> probability for each analogue                       # Key -> value format:        dict[(canonical base, strand, mod)] -> dict[(pos,qual)]
        q_brdu, q_edu = {}, {}                                                          # Create dicts for BrdU and EdU calls
        for key, vals in mods.items():                                                  # Iterate over every mod-base category that the read has
            code = str(key[2])                                                          # Extract mod-code from third element of pysam key and conv to text
            target = None
            if code in ('b', '5caU', '21983'):                                          # If the code is BrdU, direct calls to q_BrdU
                target = q_brdu
            elif code in ('e', '5fU', '21982'):                                         # If the code if EdU, direct calls to q_EdU
                target = q_edu
            if target is None:
                continue
            for qpos, qual in vals:                                                     # Using the call position and call prob (qual),
                target[qpos] = qual / 255.0                                             # Convert 0-255 conf scores to 0-1 and stores it at read position

        if not q_brdu and not q_edu:                                                    # If no recognized read calls, continue
            continue

        # query position -> reference position                                          # Now we are converting from read coord to reference coord
        q2r = {q: r for q, r in read.get_aligned_pairs(matches_only=True)}              # Builds dictionary mapping query(read)_pos to ref_pos
        calls = []                                                                      # Create output list for this read
        for qpos in set(q_brdu) | set(q_edu):                                           # Iterate over all positions within q_BrdU and q_EdU. If both, then only processed once (using union '|'')
            rpos = q2r.get(qpos)                                                        # Look up corresponding ref_pos for the query_pos
            if rpos is None:                                                            # If not aligned to reference base, ignore
                continue
            calls.append((rpos, _combine(q_brdu.get(qpos, 0.0),
                                         q_edu.get(qpos, 0.0), analogue)))              # Add ref_pos and select analogue prob using query_pos. If no BrdU/EdU call, use 0.0
        calls.sort()                                                                    # Sort above output by ref_pos
        yield (read.query_name, read.reference_name,
               '-' if read.is_reverse else '+', calls)                                  # Emit read_id, chr, strand, and calls. Conv fwd/rev to +/-
    bam.close()                                                                         # Close BAM when all reads are iterated through


def _combine(p_brdu, p_edu, analogue):                                                  # Chooses which analogue to "label"
    if analogue == 'BrdU':                                                              # If 'BrdU', return prob(BrdU). Same with 'EdU'. If 'either', return whichever is higher
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
    if not calls:                                           # Check if read has analogues
        return []
    start, end = calls[0][0], calls[-1][0]                  # Define start and end points (det. by first and last call position); need calls[] to be sorted
    n_win = max(1, (end - start) // window + 1)             # Calculate number of windows found in span (creates at least 1)
    tot = [0] * n_win                                       # Create one total-call counter per window
    pos = [0] * n_win                                       # Create one pos-call counter per window
    for p, prob in calls:                                   # For ref_pos (p) and analogue_prob (prob) in calls[]
        i = min((p - start) // window, n_win - 1)           # Calculate call[]'s zero-based window index 'i''
        tot[i] += 1                                         # Add 1 to total-call count
        if prob >= prob_thresh:
            pos[i] += 1                                     # Add 1 to pos-call count if analogue prob >= prob_thresholds
    out = []                                                # Create list to store window classification outputs
    for i in range(n_win):                                  # Iterate over each window index, incl. those with no calls
        w0 = start + i * window                             # Calculate inclusive start coord of window
        w1 = min(w0 + window, end + 1)                      # Calculate exclusive end coord of window. Stops just after final call rather than window end
        if tot[i] < min_calls:                              # Marks window uninformative if total calls < min_calls thresh
            state = None
        elif pos[i] / tot[i] >= min_frac:                   # For a window with enough coverage, calculate fraction of calls > prob_thresh and mark analogue_pos(1) or analogue_neg(0)
            state = 1
        else:
            state = 0
        out.append((w0, w1, state, tot[i]))                 # Stores window coordinates, state, and total call count in out[]
    return out                                              # Returns information on all classified windows


def call_gaps(windows, max_gap, min_gap, min_flank, min_calls_in_gap,
              max_uninformative_frac):
    """Find negative runs flanked by positive runs on both sides."""                        # This is where we define a ssDNA gap -- finding neg runs flanked by pos runs on both sides
    if not windows:                                                                         # If windows is empty, return nothing
        return []
    # collapse consecutive windows into runs of equal state
    runs = []                                                                               # Create list for storing collapsed runs (check for consecutive)
    cur = [windows[0][0], windows[0][1], windows[0][2], windows[0][3], 1]                   # Initialise current run using first window[0]. Values are: (start, end, state, total calls, n_windows)
    for w0, w1, st, n in windows[1:]:                                                       # Iterate through every window after first window
        if st == cur[2]:                                                                    # Does this window have the same state as our current run? -- goal: merge windows of same state. IF SO,
            cur[1] = w1                                                                     # Extend cur[1] end_coord to w1 end_coord
            cur[3] += n                                                                     # Add n total calls to cur[3] total calls
            cur[4] += 1                                                                     # Increment the cur[4] n_windows by 1. Window merging will iterate until it encounters a window with diff state
        else:
            runs.append(tuple(cur))                                                         # If not same state, finalise cur[] by storing it as a tuple in the runs[] list
            cur = [w0, w1, st, n, 1]                                                        # Start a new run from this next window (that doesn't match cur[])
    runs.append(tuple(cur))                                                                 # Stores final run as tuple in runs[]

    gaps = []                                                                               # Create list for storing gap informations
    for i, (s, e, st, n_calls, n_win) in enumerate(runs):                                   # Iterate thru each run; i to inspect neighboring runs, s,e,st,n_calls,n_win describe the run outs[] values
        if st != 0:                                                                         # Ignore if state is positive or uninformative; we want to look at analogue_neg (st=0) runs
            continue
        length = e - s                                                                      # Calculate neg run length using end and start coords
        if not (min_gap <= length < max_gap):                                               # Filter for length; ignore if not between min_ and max_gap values
            continue
        if n_calls < min_calls_in_gap:                                                      # Filter for calls; ignore if less than min_calls_in_gap
            continue
        # flanks: walk outward, allow uninformative windows to be crossed but
        # count them; require a positive run of >= min_flank on each side
        left = _flank(runs, i, -1, min_flank, max_uninformative_frac)                       # Searches to left of candidate gap for analogue_pos runs
        right = _flank(runs, i, +1, min_flank, max_uninformative_frac)                      # Searches to right of candidate gap for analogue_pos runs
        if left is None or right is None:
            continue                                                                        # If either flank returns None (no analogue_pos flank length, too short, or too uninformative), reject gap
        gaps.append({'start': s, 'end': e, 'length': length,                                # Otherwise, append information to gaps[] list
                     'n_calls': n_calls, 'n_win': n_win,
                     'left_flank': left, 'right_flank': right,
                     'calls_per_kb': 1000.0 * n_calls / length})                            # Normalise call count to 1kb gap length
    return gaps                                                                             # Outputs info of all accepted gaps in gaps[]


def _flank(runs, i, step, min_flank, max_uninf):
    """Accumulate positive track length adjacent to runs[i] in direction step.
    Returns flank length, or None if the adjacent run is not positive or the
    flank is too short / too uninformative."""
    j = i + step                                                        # Move to the run immediately next to the candidate gap
    acc = 0                                                             # Initialise accumulated pos-flank length
    uninf = 0                                                           # Initialise accumulated uninformative length btwn gap and pos-signal
    seen_positive = False                                               # Has pos. signal been encountered yet?
    while 0 <= j < len(runs):                                           # Only continue while search remains within the read's run list runs[]
        s, e, st, n_calls, n_win = runs[j]                              # Takes info of neighboring run
        if st == 1:
            acc += e - s                                                # If run is analogue_pos, add length to pos-flank total
            seen_positive = True                                        # Mark pos. signal seen
        elif st is None:                                                
            if seen_positive:                                           # If run is uninformative and has already seen analogue_pos, then stop flank (break).
                break
            uninf += e - s                                              # If run is uninformative and has not seen analogue_pos, add uninformative run length
            if uninf > max_uninf * min_flank:                           # If uninformative run length is too long, return None (adj run is not pos. or flank is too short/uninformative)
                return None
        else:  # another negative run terminates the flank              # If negative run, then stop flank
            break
        if acc >= min_flank:                                            # If sufficient analogue_pos seq length, then return the accepted_flank_length (acc). It only contains pos. base pairs. Any uninformative sequence is not counted towards req. flank length
            return acc
        j += step                                                       # Then move to next run
    return acc if acc >= min_flank else None                            # At end of loop, return acc only if it meets min_flank requirement, otherwise, output None




# If window is:         ++++ ++++ ---- ---- ++++ ++++
# After collapse:       positive run | negative run | positive run

# The negative run only becomes a gap if its LENGTH and CALL COUNT qualify, and BOTH FLANKING POS runs are LONG enough



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

    gaps_out = []                     # Container for every detected gap across all reads
    per_read = []                     # Container for summary tuple (one per accepted-read)
    interrogated = defaultdict(int)   # Store informative bps by chr. Missing chr entries default to 0
    chrom_gaps = defaultdict(int)     # Gap counts by chr
    len_by_chrom = defaultdict(list)  # Gap lengths by chr

    for read_id, chrom, strand, calls in source:                                    # Parsing info per read from BAM or .detect file:
        if len(calls) < 100:                                                        # Filter for reads > 100 analogue calls
            continue
        span = calls[-1][0] - calls[0][0]
        if span < args.min_read_length:                                             # Filter for min_read_length thresh (span between first and last call)
            continue
        wins = window_read(calls, args.window, args.prob_threshold,                 # If pass, then window the read. Outputs window_coords(w0, w1), state(st), and total call count
                           args.min_calls_per_window, args.min_frac)
        positive_bp = sum(w1 - w0 for w0, w1, st, _ in wins if st == 1)             # Note pos_bp count
        negative_bp = sum(w1 - w0 for w0, w1, st, _ in winsif st == 0)              # Note neg_bp count

        informative = positive_bp + negative_bp
        # interrogated denominator: informative windows only
        informative = positive_bp + negative_bp                                     # Add lengths of all windows whose state is 1 or 0 (pos/neg); informative length
        interrogated[chrom] += informative                                          # Add read's informative length to its chr total

        g = call_gaps(wins, args.max_gap, args.min_gap, args.min_flank,             # Detect candidate gaps
                      args.min_calls_in_gap, args.max_uninformative_frac)
        for x in g:                                                                 # For each gap in candidate gaps g[]
            x['read_id'] = read_id
            x['chrom'] = chrom
            x['strand'] = strand
            gaps_out.append(x)                                                      # Add information to gaps_out[]
            chrom_gaps[chrom] += 1                                                  # Increment chr_gap_count by 1
            len_by_chrom[chrom].append(x['length'])                                 # Save gap length for chr-level statistics
        per_read.append((read_id, chrom, calls[0][0], calls[-1][0],                 # Save read_id, chr, first/last call coords, strand, informative length, and n_gaps for each read. ADDED TO PER-READ.TSV FILE
                         strand, informative, len(g)))

    p = args.out_prefix                                                             # Set output_prefix (usually sample_name) as variable
    with open(p + '.gaps.bed', 'w') as fh:                                          # Opens BED file output with filename
        for x in sorted(gaps_out, key=lambda d: (d['chrom'], d['start'])):          # Sort gaps lexographically by chr name, then numerically by start
            fh.write(f"{x['chrom']}\t{x['start']}\t{x['end']}\t"                    # Write BED-like fields (chr, start, end, read_id, gap_length, strand)
                     f"{x['read_id']}\t{x['length']}\t{x['strand']}\n")

    with open(p + '.gaps.tsv', 'w') as fh:                                          # Opens detailed .gaps.tsv output
        fh.write('chrom\tstart\tend\tlength\tread_id\tstrand\t'                     # Header
                 'n_calls\tcalls_per_kb\tleft_flank\tright_flank\n')
        for x in sorted(gaps_out, key=lambda d: (d['chrom'], d['start'])):          # Write gaps info (chr, start, end, gap_length, read_id, strand, n_calls, calls_per_kb, left_flank, right_flank)
            fh.write(f"{x['chrom']}\t{x['start']}\t{x['end']}\t{x['length']}\t"
                     f"{x['read_id']}\t{x['strand']}\t{x['n_calls']}\t"
                     f"{x['calls_per_kb']:.1f}\t{x['left_flank']}\t"
                     f"{x['right_flank']}\n")

    with open(p + '.per_read.tsv', 'w') as fh:                                      # Opens detailed .per_read.tsv output
        fh.write('read_id\tchrom\tref_start\tref_end\tstrand\t'                     # Header
                 'informative_bp\tn_gaps\n')
        for row in per_read:                                                        # For each read, convert every value to text, join with \t, add \n
            fh.write('\t'.join(str(v) for v in row) + '\n')

    with open(p + '.per_chrom.tsv', 'w') as fh:                                     # Opens .per_chr.tsv output
        fh.write('chrom\tn_gaps\tinterrogated_Mb\tgaps_per_Mb\t'                    # Heaer
                 'median_gap_len\tmean_gap_len\n')
        for chrom in sorted(interrogated, key=_chrom_key):                          # Sort chromosomes with _chrom_key (see next function)
            mb = interrogated[chrom] / 1e6                                          # Converts informative bp length to megabases
            n = chrom_gaps.get(chrom, 0)                                            # Get chr gap count
            lens = sorted(len_by_chrom.get(chrom, []))                              # Gets and sorts gap lengths
            med = lens[len(lens) // 2] if lens else 0                               # Get median gap_length (if middle two, get upper)
            mean = sum(lens) / len(lens) if lens else 0                             # Calculate mean_gap_length
            rate = n / mb if mb > 0 else 0                                          # Calculate gaps per informative Mb
            fh.write(f"{chrom}\t{n}\t{mb:.3f}\t{rate:.2f}\t{med}\t{mean:.0f}\n")    # Write chr summary with decimal formatting

    tot_mb = sum(interrogated.values()) / 1e6                                       # Calculate total informative sequence in Mb, across all chr
    sys.stderr.write(                                                               # Write summary to stderr:
        f"reads used: {len(per_read)}\n"                                            # How many reads passed all read-level filters
        f"interrogated: {tot_mb:.2f} Mb\n"                                          # Total informative sequence length in Mb
        f"candidate gaps <{args.max_gap} bp: {len(gaps_out)}\n"                     # Total number of accepted candidate gaps
        f"gap burden: {len(gaps_out)/tot_mb if tot_mb else 0:.2f} per Mb\n"         # Divide all gaps by informative Mb
        f"wrote {p}.gaps.bed / .gaps.tsv / .per_read.tsv / .per_chrom.tsv\n")       # List four types of outputs created


def _chrom_key(c):                                                                  # Defines chr sorting, done numerically then alphabetically (1,2,3...10,11,12,...,X,Y)
    s = c[3:] if c.lower().startswith('chr') else c
    return (0, int(s)) if s.isdigit() else (1, s)


if __name__ == '__main__':                                                          # Controls whether main() runs. Sets __name__ to "__main__" if file is executed directly
    main()                                                                          # Then calls main(). If antoher Python script imports this file, the functions become available without automaticlaly running the analysis.
