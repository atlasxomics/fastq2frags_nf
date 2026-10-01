#!/usr/bin/env python3
"""Cross-talk correction: substitute fragments of affected barcodes with proxy
reads from the nearest unaffected neighbour (scripts/substitute_fragments.py).

Column indices are resolved from the xtalkTable header, falling back to the
known xtalkTable layout when a name is absent.
"""
import argparse
import csv
import logging

from substitute_fragments import run, resolve_threads

logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s", level=logging.INFO
)


def resolve_xtalk_columns(path):
    name_to_idx = {}
    try:
        with open(path, newline="") as fh:
            header = next(csv.reader(fh))
        for i, name in enumerate(header):
            name_to_idx[name.strip()] = i
    except (OSError, StopIteration):
        pass
    return {
        "barcode": name_to_idx.get("mer", 1),
        "row":     name_to_idx.get("row", 3),
        "col":     name_to_idx.get("col", 4),
        "mask":    name_to_idx.get("xtalk_mask", 5),
        "fold":    name_to_idx.get("fold_above_background", 7),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("table", help="xtalkTable CSV (with header)")
    parser.add_argument("fragments", help="BGZF fragment file (.tsv.gz)")
    parser.add_argument("output", help="Output BGZF fragment file (.tsv.gz); a .tbi is written alongside")
    parser.add_argument("--threshold", type=float, default=2.0)
    parser.add_argument("--threads", type=int, default=40)
    args = parser.parse_args()

    cols = resolve_xtalk_columns(args.table)
    logging.info(f"xtalk column map: {cols}; threshold={args.threshold}")
    c_threads, d_threads = resolve_threads(args.threads, None, None)
    run(
        args.table,
        args.fragments,
        args.output,
        cols["mask"],
        cols["fold"],
        args.threshold,
        "auto",
        d_threads,
        c_threads,
        None,
        32 * 1024 * 1024,
        16 * 1024 * 1024,
        barcode_col=cols["barcode"],
        row_col=cols["row"],
        col_col=cols["col"],
        skip_header=True,
    )
    logging.info("Completed cross-talk correction.")


if __name__ == "__main__":
    main()
