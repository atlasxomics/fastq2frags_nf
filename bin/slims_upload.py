#!/usr/bin/env python3
"""Push run metrics from summary.csv to SLIMS.

Credentials are read from the SLIMS_USERNAME and SLIMS_PASSWORD environment
variables (provided as Nextflow secrets).
"""
import argparse
import csv
import logging
import os
import sys

from slims.criteria import equals
from slims.slims import Slims

logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s", level=logging.INFO
)

SLIMS_URL = "https://slims.atlasxomics.com/slimsrest"

NA = ["na", "#N/A", "=NA()", "NA", "nan", None]

# summary.csv metric name -> SLIMS Result field
MAPPING = {
    "Genome Build": "rslt_cf_refGenome",
    "Atlas Pipeline Version": "rslt_cf_pipelineVersion",
    "Fraction Uniquely Aligned Reads": "rslt_cf_confidentlyMappedReadPairs",
    "Fragment-weighted FRiP": "rslt_cf_fractionOfHighQualityFragmentsOrlapPe",
    "Total Peaks Called": "rslt_cf_numberOfPeaks",
    "Fragment-weighted Duplication": "rslt_cf_percentDuplicates",
    "Number of input read pairs": "rslt_cf_sequencedReadPairs1",
    "TSS Enrichment Score": "rslt_cf_tssEnrichmentScore",
    "Fraction Unaligned Reads": "rslt_cf_unmappedReadPairs",
    "Fraction Reads with Valid Barcode": "rslt_cf_validBarcodes",
}


def csv_to_dict(path):
    with open(path) as f:
        reader = csv.reader(f)
        return dict(zip(next(reader), next(reader)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", help="summary.csv from fragReport.py")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--ng-id", default="NG00000")
    args = parser.parse_args()

    username = os.environ.get("SLIMS_USERNAME")
    password = os.environ.get("SLIMS_PASSWORD")
    if not username or not password:
        sys.exit("SLIMS_USERNAME / SLIMS_PASSWORD are not set; "
                 "define them with `nextflow secrets set`.")

    results = csv_to_dict(args.summary)
    payload = {MAPPING[k]: v for k, v in results.items()
               if k in MAPPING and v not in NA}

    ng_id = args.ng_id
    if ng_id in (None, "", "NG00000"):
        ng_id = args.run_id.split("_")[-1]
    logging.info(f"Looking up SLIMS content {ng_id}")

    slims = Slims("slims", SLIMS_URL, username, password)
    content = slims.fetch("Content", equals("cntn_id", ng_id))
    if not content:
        sys.exit(f"No SLIMS content found for ng_id '{ng_id}'.")

    payload["rslt_fk_content"] = content[0].pk()
    payload["rslt_fk_test"] = 39
    payload["rslt_value"] = "upload"

    logging.info(f"upload succeeded: {slims.add('Result', payload)}")


if __name__ == "__main__":
    main()
