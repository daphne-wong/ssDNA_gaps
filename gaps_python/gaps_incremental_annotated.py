#!/usr/bin/env python3
"""
Call short internal analogue patterns from DNAscent detect output without
using fixed, non-overlapping windows to set event boundaries.

The three reported event types are:

    EdU-negative-EdU
    BrdU-negative-BrdU
    EdU-BrdU-EdU

For each supported left flank, the middle interval is tested from just above
--min-gap to just below --max-gap in --gap-step increments (default 1 bp).
The first length with the requested middle state and matching right flank is
reported. If no matching right flank is found before --max-gap, that seed is
not reported as a short event.

Flanks have no bp-length cutoff. Instead, the nearest --min-flank-calls calls
on each side must support the requested analogue. Thus flank lengths can vary
with local thymidine density.

Fixed windows are used only to retain the original definition of informative
bp for per-read/per-chromosome denominators. They do not determine event
coordinates or lengths.

Outputs
-------
  <prefix>.gaps.bed
  <prefix>.gaps.tsv
  <prefix>.per_read.tsv
  <prefix>.per_chrom.tsv
  <prefix>.log
"""

import argparse
import bisect
import logging
import os
import sys
import time
from collections import Counter, defaultdict


EVENT_TYPES = ("EdU-negative-EdU", "BrdU-negative-BrdU", "EdU-BrdU-EdU")


# ----------------------------------------------------------------------
# Input parsing: yields (read_id, chrom, strand,
#                        [(ref_pos, p_brdu, p_edu), ...])
# ----------------------------------------------------------------------

def parse_detect(path):
    """Legacy human-readable .detect: header lines '#', read headers
    '>readID contig refStart refEnd strand', then 'pos pBrdU pEdU 6mer'."""
    read_id = chrom = strand = None
    calls = []                                                              # Holds position and probability values for calls within each read
    with open(path) as fh:
        for line in fh:                                                     # Process the input one line at a time
            if not line or line[0] == '#':                                  # Skip comment lines
                continue
            if line[0] == '>':                                              # A read header starts a new read
                if read_id is not None:
                    calls.sort()
                    yield read_id, chrom, strand, calls                     # Emit the completed previous read
                fields = line[1:].split()                                   # Split the header into columns
                if len(fields) < 2:
                    read_id = None
                    calls = []
                    continue
                raw_strand = fields[4] if len(fields) > 4 else 'fwd'        # Use fwd when the strand field is absent
                strand = '-' if raw_strand in ('rev', '-') else '+'         # Convert fwd/rev to +/-
                read_id, chrom = fields[0], fields[1]                       # Extract read ID and chromosome
                calls = []                                                  # Start the call list for this read
                continue
            fields = line.split()                                           # Split a call line into fields
            if read_id is None or len(fields) < 2:
                continue
            pos = int(fields[0])                                            # Reference coordinate
            p_brdu = float(fields[1])                                       # BrdU probability
            p_edu = float(fields[2]) if len(fields) >= 3 else 0.0           # EdU probability, or 0 when absent
            calls.append((pos, p_brdu, p_edu))                              # Retain both analogue probabilities
    if read_id is not None:                                                 # Emit the final read
        calls.sort()
        yield read_id, chrom, strand, calls


def parse_modbam(path, min_mapq, min_read_len):
    """DNAscent v4 modBAM. MM codes: 'b' = BrdU (5caU), 'e' = EdU
    (5fU). pysam exposes ML qualities on 0-255; these are rescaled to 0-1."""
    try:
        import pysam
    except ImportError:
        sys.exit("modBAM input needs pysam: python -m pip install 'pysam>=0.19'")

    bam = pysam.AlignmentFile(path, 'rb')                                               # Open BAM in binary read mode
    try:
        for read in bam.fetch(until_eof=True):                                          # Iterate through all BAM records
            if read.is_unmapped or read.is_secondary or read.is_supplementary:          # Keep primary mapped reads
                continue
            if read.mapping_quality < min_mapq:                                         # Apply mapping-quality threshold
                continue
            if read.reference_length is None or read.reference_length < min_read_len:   # Apply aligned-length threshold
                continue

            q_brdu, q_edu = {}, {}                                                      # Query-position to probability dictionaries
            for key, values in (read.modified_bases or {}).items():                     # Iterate over decoded modified-base categories
                code = str(key[2])                                                      # Extract the modification code
                target = None
                if code in ('b', '5caU', '21983'):
                    target = q_brdu
                elif code in ('e', '5fU', '21982'):
                    target = q_edu
                if target is not None:
                    for qpos, qual in values:
                        target[qpos] = qual / 255.0                                     # Rescale confidence from 0-255 to 0-1

            if not q_brdu and not q_edu:                                                # Ignore reads without recognized analogue calls
                continue
            q2r = {q: r for q, r in read.get_aligned_pairs(matches_only=True)}          # Map query to reference coordinates
            calls = []                                                                  # Store reference position and both probabilities
            for qpos in set(q_brdu) | set(q_edu):                                       # Process the union of BrdU and EdU positions once
                rpos = q2r.get(qpos)                                                    # Look up the aligned reference coordinate
                if rpos is not None:
                    calls.append((rpos, q_brdu.get(qpos, 0.0), q_edu.get(qpos, 0.0)))
            calls.sort()                                                                # Sort calls by reference coordinate
            yield read.query_name, read.reference_name, '-' if read.is_reverse else '+', calls
    finally:
        bam.close()                                                                     # Close the BAM after iteration


def make_prefix(calls, prob_threshold):
    positions = [x[0] for x in calls]
    pre_brdu = [0]
    pre_edu = [0]
    pre_either = [0]
    for _, p_brdu, p_edu in calls:
        brdu = int(p_brdu >= prob_threshold)
        edu = int(p_edu >= prob_threshold)
        pre_brdu.append(pre_brdu[-1] + brdu)
        pre_edu.append(pre_edu[-1] + edu)
        pre_either.append(pre_either[-1] + int(brdu or edu))
    return positions, pre_brdu, pre_edu, pre_either


def slice_counts(prefix, lo, hi):
    return prefix[hi] - prefix[lo]


def interval_stats(positions, prefixes, start, end):
    lo = bisect.bisect_left(positions, start)
    hi = bisect.bisect_left(positions, end)
    pre_brdu, pre_edu, pre_either = prefixes
    return {
        'lo': lo,
        'hi': hi,
        'n': hi - lo,
        'brdu': slice_counts(pre_brdu, lo, hi),
        'edu': slice_counts(pre_edu, lo, hi),
        'either': slice_counts(pre_either, lo, hi)
    }


def labelled(stats, analogue, min_fraction):
    if stats['n'] == 0:
        return False
    wanted = stats['edu'] if analogue == 'EdU' else stats['brdu']
    other = stats['brdu'] if analogue == 'EdU' else stats['edu']
    return wanted / stats['n'] >= min_fraction and wanted > other


def middle_matches(stats, middle_type, min_calls, min_labelled_fraction,
                   max_negative_fraction, max_other_fraction):
    if stats['n'] < min_calls:
        return False
    if middle_type == 'negative':
        return stats['either'] / stats['n'] < max_negative_fraction
    if middle_type == 'BrdU':
        return (stats['brdu'] / stats['n'] >= min_labelled_fraction and
                stats['edu'] / stats['n'] <= max_other_fraction and
                stats['brdu'] > stats['edu'])
    raise ValueError(f"unsupported middle type: {middle_type}")


def flank_stats(positions, prefixes, boundary, side, min_calls):
    """Summarize the nearest min_calls calls on one side of boundary."""
    pivot = bisect.bisect_left(positions, boundary)
    if side == 'left':
        lo, hi = pivot - min_calls, pivot
    else:
        lo, hi = pivot, pivot + min_calls
    if lo < 0 or hi > len(positions):
        return None
    pre_brdu, pre_edu, pre_either = prefixes
    start = positions[lo]
    end = positions[hi - 1] + 1
    return {
        'lo': lo,
        'hi': hi,
        'n': hi - lo,
        'brdu': slice_counts(pre_brdu, lo, hi),
        'edu': slice_counts(pre_edu, lo, hi),
        'either': slice_counts(pre_either, lo, hi),
        'start': start,
        'end': end,
        'length': boundary - start if side == 'left' else end - boundary
    }


def call_events(calls, min_gap, max_gap, gap_step, min_calls_in_middle,
                min_flank_calls, min_labelled_fraction,
                max_negative_fraction, max_other_fraction,
                prob_threshold):
    """Call non-overlapping short events without tiled gap boundaries."""
    if len(calls) < 2 * min_flank_calls + min_calls_in_middle:
        return []
    positions, pre_brdu, pre_edu, pre_either = make_prefix(calls, prob_threshold)
    prefixes = (pre_brdu, pre_edu, pre_either)
    events = []
    i = min_flank_calls - 1

    while i < len(positions) - min_flank_calls:
        start = positions[i] + 1
        left = flank_stats(positions, prefixes, start, 'left', min_flank_calls)
        left_type = None
        if left and labelled(left, 'EdU', min_labelled_fraction):
            left_type = 'EdU'
        elif left and labelled(left, 'BrdU', min_labelled_fraction):
            left_type = 'BrdU'
        if left_type is None:
            i += 1
            continue

        found = None
        first_length = min_gap + gap_step
        for length in range(first_length, max_gap, gap_step):
            end = start + length
            if end > positions[-1]:
                break
            right = flank_stats(positions, prefixes, end, 'right', min_flank_calls)
            if right is None or not labelled(right, left_type, min_labelled_fraction):
                continue
            middle = interval_stats(positions, prefixes, start, end)

            if middle_matches(middle, 'negative', min_calls_in_middle,
                              min_labelled_fraction, max_negative_fraction,
                              max_other_fraction):
                event_type = f'{left_type}-negative-{left_type}'
                found = (end, event_type, middle, left, right)
                break
            if left_type == 'EdU' and middle_matches(
                    middle, 'BrdU', min_calls_in_middle,
                    min_labelled_fraction, max_negative_fraction,
                    max_other_fraction):
                found = (end, 'EdU-BrdU-EdU', middle, left, right)
                break

        if found is None:
            i += 1
            continue

        end, event_type, middle, left, right = found
        length = end - start
        events.append({
            'start': start,
            'end': end,
            'length': length,
            'event_type': event_type,
            'n_calls': middle['n'],
            'calls_per_kb': 1000.0 * middle['n'] / length,
            'middle_brdu_fraction': middle['brdu'] / middle['n'],
            'middle_edu_fraction': middle['edu'] / middle['n'],
            'middle_labelled_fraction': middle['either'] / middle['n'],
            'left_flank': left['length'],
            'right_flank': right['length'],
            'left_flank_calls': left['n'],
            'right_flank_calls': right['n']
        })
        i = max(i + 1, bisect.bisect_left(positions, end))
    return events


def informative_bp(calls, window, prob_threshold, min_calls, min_fraction):
    """Original tiled denominator only; this does not set event boundaries."""
    if not calls:
        return 0
    start, end = calls[0][0], calls[-1][0]                          # Span between first and last calls
    n_windows = max(1, (end - start) // window + 1)                 # Create at least one denominator window
    total = [0] * n_windows                                         # Total calls per window
    positive = [0] * n_windows                                      # Analogue-positive calls per window
    for pos, p_brdu, p_edu in calls:
        idx = min((pos - start) // window, n_windows - 1)           # Zero-based window index
        total[idx] += 1
        if max(p_brdu, p_edu) >= prob_threshold:
            positive[idx] += 1
    bp = 0
    for idx in range(n_windows):
        if total[idx] < min_calls:                                  # Low-density windows are uninformative
            continue
        window_start = start + idx * window                         # Inclusive window start
        window_end = min(window_start + window, end + 1)            # Exclusive window end
        # Both labelled and analogue-negative classified windows 
        # are considered informative.
        _ = positive[idx] / total[idx] >= min_fraction
        bp += window_end - window_start
    return bp


def chrom_key(chrom):
    """Sort chromosomes numerically first, then alphabetically (1...22, X, Y)."""
    value = chrom[3:] if chrom.lower().startswith('chr') else chrom
    return (0, int(value)) if value.isdigit() else (1, value)


def configure_logger(path):
    logger = logging.getLogger('gaps_incremental')
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter('%(asctime)s\t%(levelname)s\t%(message)s')
    file_handler = logging.FileHandler(path, mode='w')
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler(sys.stderr)
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input', required=True, help='DNAscent .bam/modBAM or legacy .detect file')
    parser.add_argument('-o', '--out-prefix', required=True)
    parser.add_argument('--min-gap', type=int, default=100, help='exclusive lower gap-length bound (default: 100)')
    parser.add_argument('--max-gap', type=int, default=1000, help='exclusive upper gap-length bound (default: 1000)')
    parser.add_argument('--gap-step', type=int, default=1, help='bp increment when growing a candidate (default: 1)')
    parser.add_argument('--prob-threshold', type=float, default=0.5, help='per-call analogue probability threshold (default: 0.5)')
    parser.add_argument('--min-labelled-frac', type=float, default=0.05, help='minimum desired-analogue fraction in a labelled region (default: 0.05)')
    parser.add_argument('--max-negative-frac', type=float, default=0.05, help='maximum either-analogue fraction in a negative middle (default: 0.05)')
    parser.add_argument('--max-other-analogue-frac', type=float, default=0.05, help='maximum EdU fraction in a BrdU middle (default: 0.05)')
    parser.add_argument('--min-calls-in-gap', type=int, default=20, help='minimum calls in the middle interval (default: 20)')
    parser.add_argument('--min-flank-calls', type=int, default=20, help='nearest calls used to support each flank; no bp-length limit (default: 20)')
    parser.add_argument('--window', type=int, default=100, help='window used only for informative-bp denominators (default: 100)')
    parser.add_argument('--min-calls-per-window', type=int, default=10, help='minimum calls in an informative denominator window (default: 10)')
    parser.add_argument('--min-frac', type=float, default=0.05, help='retained for denominator compatibility (default: 0.05)')
    parser.add_argument('--min-mapq', type=int, default=20, help='modBAM only (default: 20)')
    parser.add_argument('--min-read-length', type=int, default=5000, help='minimum called reference span (default: 5000)')
    parser.add_argument('--progress-every', type=int, default=10000, help='log progress every N parsed reads; 0 disables (default: 10000)')
    return parser.parse_args()


def validate_args(args):
    if args.min_gap < 0 or args.max_gap <= args.min_gap + 1:
        raise ValueError('--max-gap must leave at least one integer length strictly above --min-gap')
    if args.gap_step <= 0:
        raise ValueError('--gap-step must be positive')
    if args.min_calls_in_gap <= 0 or args.min_flank_calls <= 0:
        raise ValueError('call-count thresholds must be positive')
    if args.window <= 0 or args.min_calls_per_window <= 0:
        raise ValueError('denominator window settings must be positive')
    for name in ('prob_threshold', 'min_labelled_frac', 'max_negative_frac', 'max_other_analogue_frac', 'min_frac'):
        value = getattr(args, name)
        if not 0 <= value <= 1:
            raise ValueError(f'--{name.replace("_", "-")} must be between 0 and 1')


def main():
    args = parse_args()
    validate_args(args)
    output_dir = os.path.dirname(os.path.abspath(args.out_prefix))
    os.makedirs(output_dir, exist_ok=True)
    logger = configure_logger(args.out_prefix + '.log')
    started = time.time()
    logger.info('started incremental event calling')
    logger.info('input=%s out_prefix=%s', args.input, args.out_prefix)
    logger.info('parameters=%s', vars(args))

    if args.input.lower().endswith('.bam'):
        source = parse_modbam(args.input, args.min_mapq, args.min_read_length)
        logger.info('input_type=modBAM')
    else:
        source = parse_detect(args.input)
        logger.info('input_type=legacy_detect')

    all_events = []                                                 # Container for every detected event across all reads
    per_read = []                                                   # One summary tuple per accepted read
    interrogated = defaultdict(int)                                 # Informative bp by chromosome; missing keys default to 0
    chrom_events = defaultdict(Counter)                             # Event counts by chromosome and event type
    lengths_by_chrom = defaultdict(list)                            # Event lengths used for chromosome-level statistics
    counters = Counter()

    for read_id, chrom, strand, calls in source:                    # Process one parsed BAM or .detect read
        counters['reads_parsed'] += 1
        if args.progress_every and counters['reads_parsed'] % args.progress_every == 0:
            logger.info('progress reads_parsed=%d reads_used=%d events=%d', counters['reads_parsed'], counters['reads_used'], len(all_events))
        if len(calls) < 100:                                        # Require at least 100 analogue calls
            counters['reads_too_few_calls'] += 1
            continue
        span = calls[-1][0] - calls[0][0]                                           # Reference span between first and last calls
        if span < args.min_read_length:                                             # Apply minimum read-span threshold
            counters['reads_too_short'] += 1
            continue

        counters['reads_used'] += 1
        informative = informative_bp(calls, args.window, args.prob_threshold,
                                     args.min_calls_per_window, args.min_frac)
        interrogated[chrom] += informative                                          # Add this read's informative length to its chromosome
        events = call_events(
            calls, args.min_gap, args.max_gap, args.gap_step,
            args.min_calls_in_gap, args.min_flank_calls,
            args.min_labelled_frac, args.max_negative_frac,
            args.max_other_analogue_frac, args.prob_threshold)
        read_counts = Counter(event['event_type'] for event in events)
        for event in events:                                                        # Add read and alignment information to each accepted event
            event.update(read_id=read_id, chrom=chrom, strand=strand)
            all_events.append(event)
            chrom_events[chrom][event['event_type']] += 1
            lengths_by_chrom[chrom].append(event['length'])
        per_read.append((read_id, chrom, calls[0][0], calls[-1][0], strand,         # Save the per-read summary row
                         informative, len(events),
                         read_counts['EdU-negative-EdU'],
                         read_counts['BrdU-negative-BrdU'],
                         read_counts['EdU-BrdU-EdU']))

    events_sorted = sorted(all_events, key=lambda event: (chrom_key(event['chrom']), event['start'], event['read_id'])) # Sort by chromosome and start
    prefix = args.out_prefix
    with open(prefix + '.gaps.bed', 'w') as handle:                # Open the BED6 event output
        for event in events_sorted:
            name = f"{event['read_id']}|{event['event_type']}"
            handle.write(f"{event['chrom']}\t{event['start']}\t{event['end']}\t{name}\t{event['length']}\t{event['strand']}\n")

    with open(prefix + '.gaps.tsv', 'w') as handle:                # Open the detailed event output
        handle.write('chrom\tstart\tend\tlength\tread_id\tstrand\tevent_type\t'
                     'n_calls\tcalls_per_kb\tmiddle_edu_fraction\t'
                     'middle_brdu_fraction\tmiddle_labelled_fraction\t'
                     'left_flank\tright_flank\tleft_flank_calls\tright_flank_calls\n')
        for event in events_sorted:
            handle.write(
                f"{event['chrom']}\t{event['start']}\t{event['end']}\t{event['length']}\t"
                f"{event['read_id']}\t{event['strand']}\t{event['event_type']}\t"
                f"{event['n_calls']}\t{event['calls_per_kb']:.1f}\t"
                f"{event['middle_edu_fraction']:.4f}\t{event['middle_brdu_fraction']:.4f}\t"
                f"{event['middle_labelled_fraction']:.4f}\t{event['left_flank']}\t"
                f"{event['right_flank']}\t{event['left_flank_calls']}\t"
                f"{event['right_flank_calls']}\n")

    with open(prefix + '.per_read.tsv', 'w') as handle:             # Open the per-read summary output
        handle.write('read_id\tchrom\tref_start\tref_end\tstrand\tinformative_bp\t'
                     'n_gaps\tn_EdU_negative_EdU\tn_BrdU_negative_BrdU\tn_EdU_BrdU_EdU\n')
        for row in per_read:                                                
            handle.write('\t'.join(str(value) for value in row) + '\n')     # Convert values to text and join them with tabs

    with open(prefix + '.per_chrom.tsv', 'w') as handle:            # Open the per-chromosome summary output
        handle.write('chrom\tn_gaps\tinterrogated_Mb\tgaps_per_Mb\tmedian_gap_len\t'
                     'mean_gap_len\tn_EdU_negative_EdU\tn_BrdU_negative_BrdU\t'
                     'n_EdU_BrdU_EdU\n')
        for chrom in sorted(interrogated, key=chrom_key):           # Sort chromosomes numerically, then alphabetically
            mb = interrogated[chrom] / 1e6                          # Convert informative bp to megabases
            event_counts = chrom_events[chrom]                      # Get counts for each event type
            n_events = sum(event_counts.values())                   # Count all accepted event types
            lengths = sorted(lengths_by_chrom.get(chrom, []))       # Sort lengths for summary statistics
            median = lengths[len(lengths) // 2] if lengths else 0   # Use the upper middle value for an even count
            mean = sum(lengths) / len(lengths) if lengths else 0    # Calculate mean event length
            rate = n_events / mb if mb else 0                       # Events per informative Mb
            handle.write(
                f"{chrom}\t{n_events}\t{mb:.3f}\t{rate:.2f}\t{median}\t{mean:.0f}\t"
                f"{event_counts['EdU-negative-EdU']}\t"
                f"{event_counts['BrdU-negative-BrdU']}\t"
                f"{event_counts['EdU-BrdU-EdU']}\n")

    total_mb = sum(interrogated.values()) / 1e6                     # Total informative sequence across chromosomes
    event_counts = Counter(event['event_type'] for event in all_events)
    logger.info('reads_parsed=%d reads_used=%d reads_too_few_calls=%d reads_too_short=%d',
                counters['reads_parsed'], counters['reads_used'],
                counters['reads_too_few_calls'], counters['reads_too_short'])
    logger.info('interrogated_Mb=%.3f total_events=%d gap_burden_per_Mb=%.3f',
                total_mb, len(all_events), len(all_events) / total_mb if total_mb else 0)
    for event_type in EVENT_TYPES:
        logger.info('%s=%d', event_type, event_counts[event_type])
    logger.info('wrote %s.gaps.bed, %s.gaps.tsv, %s.per_read.tsv, %s.per_chrom.tsv, %s.log',
                prefix, prefix, prefix, prefix, prefix)
    logger.info('finished elapsed_seconds=%.2f', time.time() - started)


if __name__ == '__main__':                                         # Run main only when executed as a script
    main()
