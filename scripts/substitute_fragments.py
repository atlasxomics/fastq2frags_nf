#!/usr/bin/env python3
"""
substitute_fragments.py

For each affected barcode in masked_positions.csv (mask column == 0), removes its
fragments from the input fragment file and adds copies of an *unaffected* substitute
barcode's fragments, relabeled with the affected barcode. The substitute barcode's own
fragments are kept unchanged.

Substitute selection (build_lookups): the substitute MUST be an unaffected position.
Starting from the affected spot (r, c), the first in-bounds, unaffected position in this
priority order is used (column shifts preferred, distance 1 before 2, + before -, no
diagonals, never beyond distance 2):

    col+1, col-1, row+1, row-1, col+2, col-2, row+2, row-2

If none of those eight is unaffected, the position is left completely untouched
(pass-through): its fragments are kept and no replacement is inserted.

Sort order is preserved: relabeled copies are written immediately after the substitute
barcode's own line at the same (chrom, start, stop), so the output stays tabix-sorted.

Performance: the OneDrive /mnt/c (9p) filesystem is the dominant cost (~11-65 MB/s, no
page cache) vs WSL-native ext4 (~430 MB/s) -- keep inputs/outputs on ext4. The pipeline
decouples I/O from compute and uses a compiled C inner filter (3.7x faster than the
pure-Python loop; the bottleneck then becomes bgzip compression):

    bgzip -d -c -@ D  INPUT  ->  substitute_filter (C)  ->  bgzip -@ C -c  ->  OUTPUT

The C engine is auto-compiled from substitute_filter.c on first run (needs gcc). If gcc
is unavailable, the script falls back to an equivalent (slower) pure-Python filter.

Usage (WSL, conda env latch1 active; files on ext4):
    python substitute_fragments.py \\
        masked_positions.csv fragments.tsv.gz output.tsv.gz [--mask-col 4] [--threads 11]

Runtime deps (must be on PATH): bgzip (htslib), tabix; gcc for the C engine (optional).
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path
from typing import Optional

PROGRESS_EVERY = 10_000_000        # stderr update interval (Python fallback only)

# Substitute search: (d_row, d_col) offsets in priority order, with labels for the summary.
# Column shifts preferred over row shifts; distance 1 before distance 2; + before -.
SUB_OFFSETS = [(0, 1), (0, -1), (1, 0), (-1, 0), (0, 2), (0, -2), (2, 0), (-2, 0)]
SUB_LABELS  = ['col+1', 'col-1', 'row+1', 'row-1', 'col+2', 'col-2', 'row+2', 'row-2']

SCRIPT_DIR = Path(__file__).resolve().parent
C_SOURCE   = SCRIPT_DIR / 'substitute_filter.c'
C_BINARY   = SCRIPT_DIR / 'substitute_filter'


# ---------------------------------------------------------------------------
# Lookup construction
# ---------------------------------------------------------------------------

def build_lookups(masked_csv_path: str, mask_col: int, fold_col: int = 5,
                  threshold: Optional[float] = None,
                  barcode_col: int = 0, row_col: int = 2, col_col: int = 3,
                  skip_header: bool = False):
    """
    Read masked_positions.csv once and return:
      affected_set : set[bytes]               16-mers to discard (have a valid substitute)
      sub_to_aff   : dict[bytes, list[bytes]] substitute_16mer -> [affected_16mer, ...]
      stats        : dict                     per-rule counts for the build-time summary

    The substitute for each affected (r, c) is the first in-bounds, unaffected position in
    SUB_OFFSETS order. If no such position exists within distance 2, the affected position
    is left untouched (NOT added to affected_set) and counted under 'no_substitute'.

    When threshold is set, affected positions whose fold_col value is strictly less than
    threshold are also left untouched (counted under 'below_threshold').
    """
    all_rows = []
    with open(masked_csv_path, newline='') as fh:
        for row in csv.reader(fh):
            if row:
                all_rows.append(row)
    if skip_header:
        all_rows = all_rows[1:]

    # (row, col) -> 16-mer barcode (bytes), fold value, and the set of affected positions.
    pos_to_bc: dict   = {}
    pos_to_fold: dict = {}
    affected_pos: set = set()
    for row in all_rows:
        bc = row[barcode_col].strip().encode()
        r, c = int(row[row_col]), int(row[col_col])
        pos_to_bc[(r, c)] = bc
        try:
            pos_to_fold[(r, c)] = float(row[fold_col])
        except (IndexError, ValueError):
            pos_to_fold[(r, c)] = None
        try:
            if int(row[mask_col]) == 0:
                affected_pos.add((r, c))
        except (IndexError, ValueError):
            pass

    affected_set: set = set()
    sub_to_aff: dict  = defaultdict(list)

    stats = {label: 0 for label in SUB_LABELS}
    stats['no_substitute']   = 0
    stats['below_threshold'] = 0
    stats['total_affected']  = len(affected_pos)

    for (r, c) in affected_pos:
        if threshold is not None:
            fold_val = pos_to_fold.get((r, c))
            if fold_val is None or fold_val < threshold:
                stats['below_threshold'] += 1
                continue
        aff_bc = pos_to_bc[(r, c)]
        for (dr, dc), label in zip(SUB_OFFSETS, SUB_LABELS):
            cand = (r + dr, c + dc)
            if cand in pos_to_bc and cand not in affected_pos:
                affected_set.add(aff_bc)
                sub_to_aff[pos_to_bc[cand]].append(aff_bc)
                stats[label] += 1
                break
        else:
            # No unaffected neighbor within distance 2 -> leave untouched (pass-through).
            stats['no_substitute'] += 1

    return affected_set, sub_to_aff, stats


def print_build_summary(stats: dict, n_sub_barcodes: int) -> None:
    print('Substitution list built:', file=sys.stderr)
    print(f'  Total affected positions : {stats["total_affected"]}', file=sys.stderr)
    for label in SUB_LABELS:
        if stats[label]:
            print(f'    resolved by {label:<6} : {stats[label]}', file=sys.stderr)
    if stats.get('below_threshold'):
        print(f'  Below threshold (left as-is): {stats["below_threshold"]}', file=sys.stderr)
    if stats['no_substitute']:
        print(f'  No substitute (left as-is) : {stats["no_substitute"]}', file=sys.stderr)
    resolved = stats['total_affected'] - stats['no_substitute'] - stats.get('below_threshold', 0)
    print(f'  Affected barcodes removed  : {resolved} '
          f'(-> {n_sub_barcodes} substitute barcodes)', file=sys.stderr)


def write_map_file(affected_set, sub_to_aff, path) -> None:
    """Emit the C engine's map file: 'A\\t<16mer>' and 'S\\t<sub>\\t<aff>' lines."""
    with open(path, 'w') as fh:
        for bc in affected_set:
            fh.write('A\t' + bc.decode() + '\n')
        for sub, affs in sub_to_aff.items():
            s = sub.decode()
            for aff in affs:
                fh.write('S\t' + s + '\t' + aff.decode() + '\n')


# ---------------------------------------------------------------------------
# C engine build
# ---------------------------------------------------------------------------

def ensure_c_binary():
    """Compile substitute_filter.c if needed. Return the binary Path, or None if gcc is
    unavailable / compilation fails."""
    gcc = shutil.which('gcc') or shutil.which('cc')
    if not gcc:
        return None
    if not C_SOURCE.exists():
        print(f'  WARNING: {C_SOURCE.name} not found next to script', file=sys.stderr)
        return None
    needs_build = (not C_BINARY.exists()
                   or C_BINARY.stat().st_mtime < C_SOURCE.stat().st_mtime)
    if needs_build:
        print(f'Compiling C engine: {gcc} -O3 -o {C_BINARY.name} {C_SOURCE.name}', file=sys.stderr)
        try:
            subprocess.run([gcc, '-O3', '-o', str(C_BINARY), str(C_SOURCE)], check=True)
        except subprocess.CalledProcessError as e:
            print(f'  WARNING: C compile failed ({e}); falling back to Python', file=sys.stderr)
            return None
    return C_BINARY


# ---------------------------------------------------------------------------
# Pipelines
# ---------------------------------------------------------------------------

def _decompress_cmd(frag, d_threads):
    return ['bgzip', '-d', '-c', '-@', str(d_threads), frag]


def _compress_cmd(c_threads, level):
    cmd = ['bgzip', '-@', str(c_threads), '-c']
    if level is not None:
        cmd += ['-l', str(level)]
    return cmd


def run_c_pipeline(binary, map_path, frag, out_path, d_threads, c_threads, level):
    """bgzip -d -@D | substitute_filter map | bgzip -@C -c > out"""
    out_fh = open(out_path, 'wb')
    zin = subprocess.Popen(_decompress_cmd(frag, d_threads),
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    cf = subprocess.Popen([str(binary), str(map_path)],
                          stdin=zin.stdout, stdout=subprocess.PIPE)
    zin.stdout.close()                      # cf owns the read end now
    comp = subprocess.Popen(_compress_cmd(c_threads, level),
                            stdin=cf.stdout, stdout=out_fh)
    cf.stdout.close()                       # comp owns the read end now

    comp.wait()
    cf.wait()
    dec_err = zin.stderr.read() if zin.stderr else b''
    zin.wait()
    out_fh.close()

    for ln in dec_err.splitlines():
        if ln and b'EOF marker is absent' not in ln:
            sys.stderr.buffer.write(b'  [bgzip] ' + ln + b'\n')

    for name, rc in (('bgzip(decompress)', zin.returncode),
                     ('substitute_filter', cf.returncode),
                     ('bgzip(compress)', comp.returncode)):
        if rc != 0:
            print(f'\nERROR: {name} exited with code {rc}', file=sys.stderr)
            sys.exit(1)


def run_python_pipeline(affected_set, sub_to_aff, frag, out_path,
                        d_threads, c_threads, level, block_size, flush_size):
    """Pure-Python fallback filter (equivalent output, ~3.7x slower than the C engine)."""
    out_fh = open(out_path, 'wb')
    zin = subprocess.Popen(_decompress_cmd(frag, d_threads),
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=block_size)
    comp = subprocess.Popen(_compress_cmd(c_threads, level),
                            stdin=subprocess.PIPE, stdout=out_fh, bufsize=0)
    comp_in = comp.stdin

    lines_read = lines_written = 0
    tail = b''
    out_buf = bytearray()
    extend, append = out_buf.extend, out_buf.append
    is_affected = affected_set.__contains__
    sub_get = sub_to_aff.get
    t0 = time.monotonic()

    def flush():
        if out_buf:
            comp_in.write(out_buf)
            out_buf.clear()

    def handle(line):
        nonlocal lines_written
        try:
            a = line.index(b'\t'); b = line.index(b'\t', a + 1)
            c = line.index(b'\t', b + 1); d = line.index(b'\t', c + 1)
        except ValueError:
            extend(line); append(10); lines_written += 1; return
        dash = line.rfind(b'-', c + 1, d)
        bc16 = line[c + 1:dash] if dash >= 0 else line[c + 1:d]
        if is_affected(bc16):
            return
        extend(line); append(10); lines_written += 1
        aff_list = sub_get(bc16)
        if aff_list:
            suffix = line[dash:d] if dash >= 0 else b''
            rest = line[d:]
            prefix = line[:c + 1]
            for aff in aff_list:
                extend(prefix); extend(aff); extend(suffix); extend(rest); append(10)
                lines_written += 1

    try:
        while True:
            block = zin.stdout.read(block_size)
            if not block:
                break
            parts = (tail + block).split(b'\n')
            tail = parts.pop()
            for line in parts:
                if not line:
                    continue
                lines_read += 1
                handle(line)
                if len(out_buf) >= flush_size:
                    flush()
                if lines_read % PROGRESS_EVERY == 0:
                    el = time.monotonic() - t0
                    print(f'  {lines_read/1e6:,.0f}M read  {lines_written/1e6:,.0f}M written  '
                          f'{lines_read/el/1e6:.1f}M lines/s', file=sys.stderr, end='\r')
        if tail:
            lines_read += 1
            handle(tail)
        flush()
    finally:
        comp.stdin.close(); comp.wait(); out_fh.close()
        dec_err = zin.stderr.read() if zin.stderr else b''
        zin.stdout.close(); zin.wait()

    for ln in dec_err.splitlines():
        if ln and b'EOF marker is absent' not in ln:
            sys.stderr.buffer.write(b'  [bgzip] ' + ln + b'\n')
    if comp.returncode != 0 or zin.returncode != 0:
        print(f'\nERROR: bgzip exited (decompress={zin.returncode}, compress={comp.returncode})',
              file=sys.stderr)
        sys.exit(1)
    print(f'\n  {lines_read:,} read, {lines_written:,} written', file=sys.stderr)


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def run(masked_csv, fragments_gz, output_gz, mask_col, fold_col, threshold, engine,
        d_threads, c_threads, level, block_size, flush_size,
        barcode_col=0, row_col=2, col_col=3, skip_header=False) -> None:
    t_start = time.monotonic()
    output_path = Path(output_gz)

    print('Building barcode lookup tables ...', file=sys.stderr)
    affected_set, sub_to_aff, stats = build_lookups(
        masked_csv, mask_col, fold_col, threshold,
        barcode_col=barcode_col, row_col=row_col, col_col=col_col, skip_header=skip_header,
    )
    print_build_summary(stats, len(sub_to_aff))

    # Decide engine.
    binary = None
    if engine in ('auto', 'c'):
        binary = ensure_c_binary()
        if binary is None and engine == 'c':
            print('ERROR: --engine c requested but gcc/C engine unavailable', file=sys.stderr)
            sys.exit(1)
    use_c = binary is not None and engine != 'python'

    print(f'Engine: {"C (substitute_filter)" if use_c else "Python fallback"}  |  '
          f'bgzip decompress -@{d_threads}  compress -@{c_threads}'
          f'{f" -l{level}" if level is not None else ""}', file=sys.stderr)
    print(f'Streaming {fragments_gz} -> {output_path} ...', file=sys.stderr)

    map_path = None
    try:
        if use_c:
            fd, map_path = tempfile.mkstemp(suffix='.submap')
            os.close(fd)
            write_map_file(affected_set, sub_to_aff, map_path)
            run_c_pipeline(binary, map_path, fragments_gz, output_path,
                           d_threads, c_threads, level)
        else:
            run_python_pipeline(affected_set, sub_to_aff, fragments_gz, output_path,
                                d_threads, c_threads, level, block_size, flush_size)
    finally:
        if map_path and os.path.exists(map_path):
            os.unlink(map_path)

    elapsed = time.monotonic() - t_start
    print(f'Streaming done in {elapsed:.1f}s', file=sys.stderr)

    print('Running tabix index ...', file=sys.stderr)
    subprocess.run(['tabix', '-p', 'bed', str(output_path)], check=True)
    print(f'Output : {output_path}', file=sys.stderr)
    print(f'Index  : {output_path}.tbi', file=sys.stderr)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def resolve_threads(total, c_threads, d_threads):
    """Split a total thread budget: reserve 1 core for the filter, give the remainder
    mostly to compression (the heavier stage), the rest to decompression."""
    if c_threads is not None and d_threads is not None:
        return c_threads, d_threads
    avail = max(2, total - 1)
    if d_threads is None and c_threads is None:
        d = max(1, avail // 5)
        c = max(1, avail - d)
    elif d_threads is None:
        c = c_threads
        d = max(1, avail - c)
    else:
        d = d_threads
        c = max(1, avail - d)
    return c, d


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            'Substitute affected-barcode fragments with proxy reads from the nearest '
            'UNAFFECTED neighbor barcode (col+/-1, row+/-1, col/row +/-2; pass-through if '
            'none). Input: sorted, BGZF-gzipped TSV (chrom, start, stop, barcode, count). '
            'Output is bgzipped and tabix-indexed. Keep files on WSL-native ext4 for speed.'
        )
    )
    parser.add_argument('masked_positions',
                        help='CSV with per-barcode mask info. Default layout (no header): barcode(0), on_tissue(1), row(2), col(3), mask(4+). Use --barcode-col/--row-col/--col-col to override; use --header if the file has a header row.')
    parser.add_argument('fragments', help='BGZF-gzipped input fragment file (.tsv.gz)')
    parser.add_argument('output', help='Output path for bgzipped fragment file (.tsv.gz)')
    parser.add_argument('--barcode-col', type=int, default=0,
                        help='0-based column index for the barcode/16-mer sequence (default: 0)')
    parser.add_argument('--row-col', type=int, default=2,
                        help='0-based column index for the spatial row coordinate (default: 2)')
    parser.add_argument('--col-col', type=int, default=3,
                        help='0-based column index for the spatial col coordinate (default: 3)')
    parser.add_argument('--mask-col', type=int, default=4,
                        help='0-based mask column index (0=affected/replace, 1=unaffected). Default: 4')
    parser.add_argument('--fold-col', type=int, default=5,
                        help='0-based column index for the fold/magnitude value (default: 5)')
    parser.add_argument('--header', action='store_true',
                        help='Skip the first row of masked_positions as a header (default: no header)')
    parser.add_argument('--threshold', type=float, default=None,
                        help='Minimum fold/magnitude value required for substitution. Affected '
                             'positions with fold < threshold are left unchanged (pass-through). '
                             'Default: substitute all affected positions regardless of fold value.')

    parser.add_argument('--engine', choices=['auto', 'c', 'python'], default='auto',
                        help='Hot-loop engine: auto (C if gcc available, else Python), '
                             'c (force C), python (force pure-Python). Default: auto')
    parser.add_argument('--threads', type=int, default=min(11, os.cpu_count() or 11),
                        help='Total thread budget; auto-split into compress/decompress '
                             '(1 core reserved for the filter). Default: min(11, ncpu)')
    parser.add_argument('--compress-threads', type=int, default=None,
                        help='Override bgzip compression threads')
    parser.add_argument('--decompress-threads', type=int, default=None,
                        help='Override bgzip decompression threads')
    parser.add_argument('--compress-level', type=int, default=None,
                        help='bgzip compression level 1-9 (default: bgzip default 6). '
                             'Lower = faster + larger; compression is the pipeline bottleneck.')

    parser.add_argument('--block-mb', type=int, default=32,
                        help='Python-fallback decompress read-block size in MB (default: 32)')
    parser.add_argument('--flush-mb', type=int, default=16,
                        help='Python-fallback output buffer flush threshold in MB (default: 16)')
    args = parser.parse_args()

    c_threads, d_threads = resolve_threads(args.threads, args.compress_threads, args.decompress_threads)
    block_size = args.block_mb * 1024 * 1024
    flush_size = args.flush_mb * 1024 * 1024

    run(args.masked_positions, args.fragments, args.output,
        args.mask_col, args.fold_col, args.threshold, args.engine,
        d_threads, c_threads, args.compress_level, block_size, flush_size,
        barcode_col=args.barcode_col, row_col=args.row_col, col_col=args.col_col,
        skip_header=args.header)


if __name__ == '__main__':
    main()
