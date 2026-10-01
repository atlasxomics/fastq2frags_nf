#!/usr/bin/env python3
"""Remove out-of-bound and (optionally) blacklisted fragments from a
space-delimited Chromap fragment file; writes a tab-delimited fragment file.

Wraps lib/oob.py (installed in the container on PYTHONPATH).
"""
import argparse
import logging
import os

import oob

logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s", level=logging.INFO
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Space-delimited fragment file (barcodes already suffixed)")
    parser.add_argument("output", help="Tab-delimited output fragment file")
    parser.add_argument("--genome-id", required=True)
    parser.add_argument("--asset-root", default="/root",
                        help="Directory containing chrom_sizes/ and blacklist/ (default: /root)")
    parser.add_argument("--no-blacklist", action="store_true",
                        help="Skip blacklist filtering")
    args = parser.parse_args()

    if args.genome_id not in oob.chromsize_paths:
        raise SystemExit(
            f"No chrom sizes configured for genome '{args.genome_id}'. "
            f"Known: {', '.join(sorted(oob.chromsize_paths))}"
        )

    chromsizes = oob.load_chromsizes(
        os.path.join(args.asset_root, oob.chromsize_paths[args.genome_id])
    )
    logging.info("Filtering out-of-bound fragments...")
    oob.filter_oob(args.input, args.output, chromsizes)

    blacklist = oob.blacklist_paths.get(args.genome_id)
    if args.no_blacklist:
        logging.info("Blacklist filtering disabled by input parameter.")
    elif blacklist is None:
        logging.info(f"No blacklist configured for {args.genome_id}; skipping blacklist filtering.")
    else:
        blacklist = os.path.join(args.asset_root, blacklist)
        logging.info(f"Filtering blacklist regions with: {blacklist}")
        tmp = f"{args.output}.blacklist.tmp"
        oob.filter_blacklist(args.output, tmp, blacklist)
        os.replace(tmp, args.output)
        logging.info("Completed blacklist filtering.")


if __name__ == "__main__":
    main()
