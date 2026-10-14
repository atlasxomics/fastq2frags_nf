#!/usr/bin/env python3
"""Compare run_metrics/summary.csv from two pipeline runs, column by column.

    tests/compare_summaries.py A B [--csv results.csv]

A and B are either two summary.csv files or two results directories (the
pipeline's --outdir). For directories, samples are matched by name using
<dir>/<sample>/run_metrics/summary.csv. Local paths and s3:// URIs both work
(S3 needs the AWS CLI and credentials).

Each column gets one of these statuses:
  PASS    within its rule
  FAIL    outside its rule
  REVIEW  can legitimately flip between runs, or has no rule; check by hand
  skip    timing and run identifiers (shown only with --all)

The rules allow for runs made before FILTER_LINKERS used bbduk ordered=t. In
those runs the read order varied with thread timing, and Chromap's
tie-breaking then moved a few fragments (about 1 in 10,000 in the test data).
Metrics computed before Chromap's tie-breaking must match exactly. Metrics
computed from the final fragments get a small tolerance.

With --csv FILE, every column's result (including skipped and identical ones)
is also written to FILE: sample, column, status, a, b, abs_diff, tolerance, rule.
The two runs are "a" and "b", in the order given.

Exit status: 0 if nothing failed, 1 if any column or sample failed.
"""
import argparse
import csv
import io
import subprocess
import sys
from pathlib import Path

# Upstream of Chromap tie-breaking: read order can't change these.
EXACT = [
    "Genome Build",
    "Cell Barcode Whitelist File",
    "Profile Length",
    "Flank Size",
    "Number of input read pairs",
    "Fraction Reads with Valid Barcode",
    "Number of reads",
    "Number of mapped reads",
    "Number of mappings",
    "Number of uni-mappings",
    "Number of corrected barcodes",
    "Number of barcodes in whitelist",
    "Fraction Uniquely Aligned Reads",
    "Fraction Unaligned Reads",
]
# Computed from the final fragments: relative tolerance (--rel), or one unit
# in the last reported decimal, whichever is larger.
RELATIVE = [
    "Mean Fragments per pixel",
    "Mean Duplication Rate",
    "Median Duplication Rate",
    "Fragment-weighted Duplication",
    "Mean Mitochondrial Fraction",
    "Median Mitochondrial Fraction",
    "Mean FRiP",
    "Median FRiP",
    "Fragment-weighted FRiP",
    "Total Peaks Called",
    "FRiP (calibrated)",
    "Library TSS Enrichment",
    "TSS Overlap Fraction",
    "Background Mean Signal",
    "TSS Enrichment Score",
    "Max Fold Enrichment",
    "Min Fold Enrichment",
]
# Absolute tolerances, for medians of counts and counts near a threshold.
ABSOLUTE = {
    "Median Fragments per pixel": 1,
    "Total Barcodes Analyzed": 2,
}
# A fraction of barcodes: allow --barcodes barcodes' worth of change.
BARCODE_FRACTION = ["Fraction Barcodes >= Calc Floor"]
# Yes/no and flags derived from values near cut-offs: may flip.
REVIEW_IF_DIFFERENT = ["FRiP Calibration Reliable", "FRiP Reliability Flags"]
IGNORE = [
    "Fragment Import Time",
    "Analysis Time",
    "Bigwig Export Time",
    "FRiP Calc Time",
    "Total Runtime",
    "Sample Run ID",
    "Atlas Pipeline Version",
]


def read_text(path):
    if path.startswith("s3://"):
        return subprocess.run(["aws", "s3", "cp", path, "-"], check=True,
                              capture_output=True, text=True).stdout
    return Path(path).read_text()


def read_summary(path):
    """Values row of summary.csv (header, values, then a row of descriptions)."""
    rows = list(csv.DictReader(io.StringIO(read_text(path))))
    if not rows:
        sys.exit(f"{path}: no data row")
    return rows[0]


def find_summaries(root):
    """Map sample name -> summary.csv path under a results directory."""
    suffix = "/run_metrics/summary.csv"
    if root.startswith("s3://"):
        prefix = root.rstrip("/") + "/"
        bucket = prefix[5:].split("/", 1)[0]
        out = subprocess.run(["aws", "s3", "ls", "--recursive", prefix], check=True,
                             capture_output=True, text=True).stdout
        keys = [line.split(None, 3)[3] for line in out.splitlines() if line.strip()]
        paths = [f"s3://{bucket}/{k}" for k in keys if k.endswith(suffix)]
    else:
        paths = [str(p) for p in Path(root).glob("*" + suffix)]
    found = {}
    for p in paths:
        sample = p[: -len(suffix)].rsplit("/", 1)[-1]
        found[sample] = p
    return found


def decimals(s):
    return len(s.split(".", 1)[1]) if "." in s else 0


def as_float(s):
    try:
        return float(s)
    except ValueError:
        return None


def check(col, a, b, a_row, b_row, args):
    """Return (status, rule, tolerance) for one column. tolerance is None
    when the column isn't compared numerically."""
    if col in IGNORE:
        return "skip", "ignored", None
    if a is None or b is None:
        return "REVIEW", "column missing in one file", None
    if a == b:
        return "PASS", "identical", None
    if col in EXACT:
        return "FAIL", "must match exactly", None
    if col in REVIEW_IF_DIFFERENT:
        return "REVIEW", "may flip near a cut-off", None

    x, y = as_float(a), as_float(b)
    if x is None or y is None:
        return "FAIL" if col in RELATIVE + list(ABSOLUTE) + BARCODE_FRACTION else "REVIEW", "non-numeric difference", None
    diff = abs(x - y)

    if col in RELATIVE:
        tol = max(args.rel * max(abs(x), abs(y)), 10 ** -max(decimals(a), decimals(b)))
        rule = f"±{args.rel:.2%} or last decimal (tol {tol:.4g})"
    elif col in ABSOLUTE:
        tol = ABSOLUTE[col]
        rule = f"±{tol} absolute"
    elif col in BARCODE_FRACTION:
        total = max(as_float(a_row.get("Total Barcodes Analyzed", "")) or 0,
                    as_float(b_row.get("Total Barcodes Analyzed", "")) or 0)
        if not total:
            return "REVIEW", "no Total Barcodes Analyzed to scale by", None
        tol = max(args.barcodes / total, 10 ** -max(decimals(a), decimals(b)))
        rule = f"±{args.barcodes} barcodes of {total:.0f} (tol {tol:.4g})"
    else:
        return "REVIEW", "no rule for this column", None
    return ("PASS" if diff <= tol + 1e-12 else "FAIL"), rule, tol


def compare(sample, a_path, b_path, args, records):
    """Print one sample's comparison, append a record per column, return pass/fail."""
    a_row, b_row = read_summary(a_path), read_summary(b_path)
    cols = list(a_row) + [c for c in b_row if c not in a_row]
    counts = {"PASS": 0, "FAIL": 0, "REVIEW": 0, "skip": 0}
    lines = []
    for col in cols:
        a, b = a_row.get(col), b_row.get(col)
        status, rule, tol = check(col, a, b, a_row, b_row, args)
        counts[status] += 1
        x, y = as_float(a or ""), as_float(b or "")
        records.append({
            "sample": sample, "column": col, "status": status, "a": a, "b": b,
            "abs_diff": "" if x is None or y is None else f"{abs(x - y):.6g}",
            "tolerance": "" if tol is None else f"{tol:.6g}", "rule": rule,
        })
        if args.all or status in ("FAIL", "REVIEW") or (status == "PASS" and a != b):
            lines.append(f"  {status:6s} {col:36s} {str(a):>14s} {str(b):>14s}  {rule}")
    print(f"== {sample}: {counts['PASS']} pass, {counts['FAIL']} fail, "
          f"{counts['REVIEW']} review, {counts['skip']} skipped")
    print("\n".join(lines) if lines else "  all compared columns identical")
    return counts["FAIL"] == 0


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("a", help="summary.csv or results directory (local or s3://)")
    parser.add_argument("b", help="summary.csv or results directory (local or s3://)")
    parser.add_argument("--rel", type=float, default=0.001,
                        help="relative tolerance for fragment-based metrics (default 0.001 = 0.1%%)")
    parser.add_argument("--barcodes", type=int, default=2,
                        help="barcodes allowed to cross the fragment floor (default 2)")
    parser.add_argument("--all", action="store_true",
                        help="list every column, not just differences")
    parser.add_argument("--csv", metavar="FILE",
                        help="also write every column's result to FILE (one row per sample and column)")
    args = parser.parse_args()

    records = []

    if args.a.endswith(".csv") and args.b.endswith(".csv"):
        ok = compare(Path(args.a).parent.parent.name or "sample", args.a, args.b, args, records)
    else:
        a, b = find_summaries(args.a), find_summaries(args.b)
        if not a or not b:
            sys.exit(f"no */run_metrics/summary.csv found under {args.a if not a else args.b}")
        ok = True
        for sample in sorted(set(a) | set(b)):
            if sample not in a or sample not in b:
                where = args.a if sample in a else args.b
                print(f"== {sample}: FAIL, only in {where}")
                records.append({"sample": sample, "column": "", "status": "FAIL",
                                "rule": f"sample only in {where}"})
                ok = False
                continue
            ok = compare(sample, a[sample], b[sample], args, records) and ok
    if args.csv:
        fields = ["sample", "column", "status", "a", "b", "abs_diff", "tolerance", "rule"]
        with open(args.csv, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, restval="")
            writer.writeheader()
            writer.writerows(records)
        print(f"wrote {len(records)} rows to {args.csv}")
    print("RESULT:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
