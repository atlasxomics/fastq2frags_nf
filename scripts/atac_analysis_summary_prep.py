import snapatac2 as snap
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import argparse
import time
import sys
import json

from datetime import datetime
from typing import Optional, Dict

import tss_norm

# globals Defs
# change to accomodate pipeline norms
frag_genome_dict = {
    "GRCh38": {
        "annotation_path": "/root/references/gencode_v41_GRCh38.gff3.gz",
        "sizes_path": "/root/chrom_sizes/hg38_chrom_sizes.txt",
        "blacklist_path": "/root/blacklist/hg38-blacklist.v2.bed"
    },
    "GRCm38": {
        "annotation_path": "/root/references/gencode_vM25_GRCm38.gff3.gz",
        "sizes_path": "/root/chrom_sizes/mm10_chrom_sizes.txt",
        "blacklist_path": "/root/blacklist/mm10-blacklist.v2.bed"
    },
    "GRCm39": {
        "annotation_path": "/root/references/gencode_vM33_GRCm39.gff3.gz",
        "sizes_path": "/root/chrom_sizes/mm39_chrom_sizes.txt",
        "blacklist_path": "/root/blacklist/mm39-blacklist.v2.bed"
    },
    "mm39": {
        "annotation_path": "/root/references/gencode_vM33_GRCm39.gff3.gz",
        "sizes_path": "/root/chrom_sizes/mm39_chrom_sizes.txt",
        "blacklist_path": "/root/blacklist/mm39-blacklist.v2.bed"
    },
    "Rnor6": {
        "annotation_path": "/root/references/ensemble.Rnor_6.0.104.chr.gtf.gz",
        "sizes_path": "/root/chrom_sizes/rn6_chrom_sizes.txt",
        "blacklist_path": None  # no rat blacklist available
    }
}

# Calibration mapping the variant-C snapatac2 FRiP (max(mean, median)) onto the
# historical pycisTopic FRIP reference: pycisTopic ~= a * snapC + b. Fit over 15 recent
# spatial (D0####) runs (a=1.097, b=-0.069, R^2=0.896, LOO-MAE=0.027); an independent
# high-barcode bulk (B0####) cohort gives a nearly identical fit (a=1.096, b=-0.055,
# R^2~0.77). The relationship holds for samples with >~12k barcodes REGARDLESS of
# protocol; it does NOT hold for sparse chips (e.g. 2500-pixel B0, R^2~0), where the
# calibrated value is unreliable. Override per-run via --frip-calib-slope/-intercept.
FRIP_CALIB_SLOPE = 1.097
FRIP_CALIB_INTERCEPT = -0.069
# Two-step (variant E: peaks/bigwig on all barcodes, FRiP on >=100) uses a different
# fit -- cross-run D0 variant-E fit (a=0.900, b=-0.114, R^2=0.920, n=15).
FRIP_CALIB_SLOPE_E = 0.900
FRIP_CALIB_INTERCEPT_E = -0.114

# FRiP reliability gate (coverage-collapse). The snapatac2<->pycisTopic calibration is
# fit on populations where most barcodes clear the calc floor (>=100 fragments). On
# sparse/background libraries (e.g. IgG controls) MACS3 overfits tight peaks onto the few
# reads, so FRiP/TSSe/fold-enrichment all look spuriously HIGH -- peak-specificity metrics
# cannot flag it. The reliable discriminator is coverage: when the calc-floor set is
# (near) empty the calibration is out of scope, so the calibrated FRiP is withheld (NA).
# See the B02745 (IgG) vs B02746 (real) investigation and bin/frip_coverage_gate.py.
FRIP_GATE_FRAC_PASS_MIN = 0.5

# Function Defs


def frip_reliability(calc_floor_fallback, frac_pass, median_nfrag,
                     min_fragments, frac_pass_min):
    """Coverage-collapse reliability gate for the calibrated FRiP.

    The calibration is only valid for populations with adequate per-barcode coverage.
    Trips (calibrated FRiP unreliable -> withheld) when any coverage condition fails:
      - calc_floor_fallback: the >=min_fragments import returned 0 barcodes, so the
        pipeline fell back to a floor of 1 (the calc set is empty at the intended floor).
      - low_frac_pass: fraction of barcodes clearing min_fragments < frac_pass_min.
      - low_median_nfrag: median fragments/barcode < min_fragments.

    Args:
        calc_floor_fallback (bool): whether the calc-floor import fell back to floor 1.
        frac_pass (float): n(barcodes >= min_fragments) / n(all input barcodes).
        median_nfrag (float): median fragments per barcode (input-floor set).
        min_fragments (int): the calc floor.
        frac_pass_min (float): minimum acceptable frac_pass (<=0 disables that condition).

    Returns:
        (reliable: bool, flags: list of str). flags is empty when reliable.
    """
    flags = []
    if calc_floor_fallback:
        flags.append("calc_floor_fallback")
    if frac_pass < frac_pass_min:
        flags.append("low_frac_pass")
    if median_nfrag < min_fragments:
        flags.append("low_median_nfrag")
    return (len(flags) == 0), flags


def _convert_value(value, value_type, line_num):
    """
    Convert string value to specified type.
    
    Args:
        value (str): The value to convert
        value_type (str): Target type for conversion
        line_num (int): Line number for error reporting
        
    Returns:
        Converted value
        
    Raises:
        ValueError: If conversion fails
    """
    try:
        if value_type == 'numeric':
            # Try int first, then float, fallback to string
            try:
                return int(value)
            except ValueError:
                try:
                    return float(value)
                except ValueError:
                    return value
        elif value_type == 'int':
            return int(value)
        elif value_type == 'float':
            return float(value)
        elif value_type == 'str':
            return value
        elif value_type == 'bool':
            if value.lower() in ('true', '1', 'yes', 'on'):
                return True
            elif value.lower() in ('false', '0', 'no', 'off'):
                return False
            else:
                raise ValueError(f"Cannot convert '{value}' to boolean")
        else:
            raise ValueError(f"Unsupported value_type: {value_type}")
    except ValueError as e:
        raise ValueError(f"Line {line_num}: {e}")


def read_file_as_dict(filepath, value_type='numeric'):
    """
    Try to open a JSON file first. If that fails, open as tab-delimited file
    where each line is a label-value pair.
    
    Args:
        filepath (str): Path to the file
        value_type (str): Type conversion for tab-delimited values. Options:
            - 'numeric': Try int first, then float, fallback to string
            - 'int': Convert to integer
            - 'float': Convert to float
            - 'str': Keep as string
            - 'bool': Convert to boolean (True/False, case-insensitive)
        
    Returns:
        dict: Dictionary containing the parsed data
        
    Raises:
        FileNotFoundError: If the file doesn't exist
        ValueError: If tab-delimited format is invalid or type conversion fails
    """
    try:
        # First, try to open as JSON
        with open(filepath, 'r', encoding='utf-8') as f:
            return json.load(f)
    
    except (json.JSONDecodeError, UnicodeDecodeError):
        # If JSON parsing fails, try tab-delimited format
        try:
            result = {}
            with open(filepath, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line:  # Skip empty lines
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) != 2:
                        raise ValueError(f"Line {line_num} doesn't have exactly 2 tab-separated values: '{line}'")
                    
                    label, value = parts
                    result[label] = _convert_value(value.strip(), value_type, line_num)
            
            return result
            
        except Exception as e:
            raise ValueError(f"Failed to parse as both JSON and tab-delimited format: {e}")


def get_unique_log_file(output_dir, prefix):
    """
    Create a unique log file name. If the base name exists, append timestamp.
    """
    base_log_file = os.path.join(output_dir, f"{prefix}_log.txt")
    
    if not os.path.exists(base_log_file):
        return base_log_file
    else:
        # Create timestamp (you can adjust format as needed)
        timestamp = datetime.now().strftime("%H%M%S")
        timestamped_log_file = os.path.join(output_dir, f"{prefix}_log_{timestamp}.txt")
        return timestamped_log_file


def parse_arguments():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(description='ATAC-seq fragment analysis script')
    parser.add_argument('fragment_file', 
                        help='Path to the fragment file (e.g., fragments.tsv.gz)')
    parser.add_argument('--skip-bigwig', 
                        action='store_true',
                        help='Skip the final bigwig generation step')
    parser.add_argument('--genome', 
                        choices=['GRCh38', 'GRCm38', 'GRCm39', 'mm39', 'Rnor6'],
                        default='GRCm38',
                        help='Genome to use (default: GRCm38)')
    parser.add_argument('--min-fragments',
                        type=int,
                        default=100,
                        help='Calc floor: min fragments per cell for FRiP + per-cell '
                             'metrics (default: 100).')
    parser.add_argument('--input-fragment-floor',
                        type=int, default=1,
                        help='Input floor: min fragments per barcode for the import used for '
                             'PEAK CALLING and BIGWIG export (default: 1 = all barcodes). When '
                             'this is < --min-fragments, the task runs the two-step "variant E" '
                             '(peaks/bigwig on all barcodes, FRiP/metrics on the calc-floor set). '
                             'Set equal to --min-fragments (e.g. 100) for single-set "variant C".')
    parser.add_argument('--n-jobs', 
                        type=int, 
                        default=12,
                        help='Number of parallel jobs (default: 12)')
    parser.add_argument('--skip-frip',
                        action='store_true',
                        help='Skip the FRiP calculation')
    # FRiP peak-set construction. Defaults produce a fixed-width (500 bp),
    # blacklist-filtered consensus peak set at MACS3 q<=0.05 -- a standardized FRiP.
    # Pass --frip-merge-half-width 0 and --frip-no-blacklist to recover the older
    # raw-MACS3-pseudobulk behaviour.
    parser.add_argument('--frip-merge-half-width',
                        type=int, default=250,
                        help='Merge MACS3 peaks into fixed-width consensus via '
                             'snap.tl.merge_peaks; peaks become 2*hw bp wide '
                             '(default: 250 -> 500 bp). Pass 0 to disable merging.')
    parser.add_argument('--frip-groupby',
                        default=None,
                        help='obs column to group cells by for MACS3 peak calling '
                             '(default: None = single pseudobulk).')
    parser.add_argument('--frip-blacklist',
                        default=None,
                        help='BED of blacklist regions subtracted from peaks before FRiP. '
                             'Default: the genome-appropriate blacklist from frag_genome_dict.')
    parser.add_argument('--frip-no-blacklist',
                        action='store_true',
                        help='Disable blacklist subtraction (default: blacklist applied).')
    parser.add_argument('--frip-macs3-qvalue',
                        type=float, default=0.05,
                        help='q-value cutoff forwarded to snap.tl.macs3 (default: 0.05).')
    # Calibration of FRiP onto the pycisTopic reference. Default (None) selects by mode:
    # two-step/variant-E -> (0.900, -0.114); single-set/variant-C -> (1.097, -0.069).
    parser.add_argument('--frip-calib-slope',
                        type=float, default=None,
                        help='Slope a for FRiP (calibrated) = a*max(mean,median)+b. '
                             'Default: mode-selected (E: 0.900, C: 1.097).')
    parser.add_argument('--frip-calib-intercept',
                        type=float, default=None,
                        help='Intercept b for FRiP (calibrated). '
                             'Default: mode-selected (E: -0.114, C: -0.069).')
    parser.add_argument('--frip-gate-frac-pass-min',
                        type=float, default=FRIP_GATE_FRAC_PASS_MIN,
                        help='Coverage-collapse gate: withhold the calibrated FRiP (NA) '
                             'when the fraction of barcodes clearing --min-fragments is '
                             'below this (default: %.2f). The calc-floor fallback and a '
                             'median depth below --min-fragments also withhold it. '
                             'Set <=0 to disable the frac_pass condition.'
                             % FRIP_GATE_FRAC_PASS_MIN)

    return parser.parse_args()


def format_time(seconds):
    """
    Format time in seconds to human readable format
    """
    if seconds < 60:
        return f"{seconds:.2f} seconds"
    elif seconds < 3600:
        minutes = seconds / 60
        return f"{minutes:.2f} minutes"
    else:
        hours = seconds / 3600
        return f"{hours:.2f} hours"


def _peaks_to_bed_df(peaks_obj):
    """
    Normalise a snapatac2 peak object (pandas or polars DataFrame) to a pandas
    DataFrame with the first three columns being Chromosome/Start/End. snapatac2
    2.8.0 stores macs3 output with a 'Chromosome','Start','End' schema; merge_peaks
    returns a frame with a single 'Peaks' string column like 'chr1:100-600'.
    """
    import pandas as pd

    # polars -> pandas
    if hasattr(peaks_obj, "to_pandas"):
        df = peaks_obj.to_pandas()
    else:
        df = peaks_obj.copy()

    cols = list(df.columns)
    if "Peaks" in cols and not {"Chromosome", "Start", "End"}.issubset(cols):
        # merge_peaks style: 'chr:start-end'
        split = df["Peaks"].str.split(":", n=1, expand=True)
        chrom = split[0]
        se = split[1].str.split("-", n=1, expand=True)
        out = pd.DataFrame({
            "Chromosome": chrom,
            "Start": se[0].astype(int),
            "End": se[1].astype(int),
        })
        return out
    # macs3 style: assume first three columns are chrom/start/end
    out = df.iloc[:, :3].copy()
    out.columns = ["Chromosome", "Start", "End"]
    out["Start"] = out["Start"].astype(int)
    out["End"] = out["End"].astype(int)
    return out


def _subtract_blacklist(bed_df, blacklist_path):
    """Remove peaks overlapping any blacklist region using pyranges."""
    import pyranges as pr

    peaks_pr = pr.PyRanges(bed_df.rename(columns={
        "Chromosome": "Chromosome", "Start": "Start", "End": "End"}))
    black_pr = pr.read_bed(blacklist_path)
    kept = peaks_pr.overlap(black_pr, invert=True)
    return kept.df[["Chromosome", "Start", "End"]].reset_index(drop=True)


def compute_frip(adata, output_dir, prefix,
                 groupby=None,               # type: Optional[str]
                 merge_half_width=None,      # type: Optional[int]
                 blacklist=None,             # type: Optional[str]
                 chrom_sizes=None,           # type: Optional[Dict[str, int]]
                 macs3_kwargs=None,          # type: Optional[Dict]
                 run_macs3=True,             # type: bool
                 score_adata=None):          # type: Optional[object]
    """
    Compute Fraction of Reads in Peaks with snapatac2.

    Default call (groupby=None, merge_half_width=None, blacklist=None,
    macs3_kwargs=None) reproduces the historical production behaviour exactly:
    ``snap.tl.macs3(adata, groupby=None)`` on a single pseudobulk, peaks written to
    ``{prefix}_peakList.bed.gz``, then ``snap.metrics.frip`` -> obs['pseudobulk_frip'].

    Optional convergence levers toward the pycisTopic reference:
      - merge_half_width: run ``snap.tl.merge_peaks(half_width=...)`` to build a
        fixed-width consensus (250 -> 500 bp), matching pycisTopic get_consensus_peaks.
      - blacklist: BED of regions subtracted from the peak set before FRiP.
      - macs3_kwargs: extra kwargs (e.g. qvalue, nolambda, shift, extsize) forwarded
        to snap.tl.macs3.
      - score_adata: if given, peaks are called on `adata` but FRiP is scored on
        `score_adata` (two-step: peaks from all barcodes, FRiP on the >=100 subset).

    Returns a dict with mean_frip, median_frip, frag_wt_frip, total_npeak,
    per_cell (pandas Series indexed by barcode) and peak_bed_path.
    """
    if macs3_kwargs is None:
        macs3_kwargs = {}

    # --- call peaks (skip to reuse peaks already in adata.uns from a prior call) ---
    if run_macs3:
        snap.tl.macs3(adata, groupby=groupby, **macs3_kwargs)
    if groupby is None:
        peaks_obj = adata.uns['macs3_pseudobulk']
    else:
        peaks_obj = adata.uns['macs3']

    # --- optional consensus merge to fixed width (matches pycisTopic) ---
    if merge_half_width is not None:
        if chrom_sizes is None:
            raise ValueError("merge_half_width requires chrom_sizes")
        if groupby is None:
            merge_input = {'pseudobulk': peaks_obj}
        else:
            merge_input = peaks_obj
        merged = snap.tl.merge_peaks(merge_input, chrom_sizes,
                                     half_width=merge_half_width)
        bed_df = _peaks_to_bed_df(merged)
    else:
        bed_df = _peaks_to_bed_df(peaks_obj)

    # --- optional blacklist subtraction ---
    if blacklist is not None:
        bed_df = _subtract_blacklist(bed_df, blacklist)

    # --- write peak BED and score FRiP ---
    peakFilePath = os.path.join(output_dir, f"{prefix}_peakList.bed.gz")
    bed_df.to_csv(peakFilePath, sep="\t", index=False, header=False)
    total_npeak = len(bed_df)

    # Peaks come from `adata` (the input-floor set); FRiP may be scored on a
    # different set (score_adata, e.g. the >=100 calc-floor set) -> two-step mode.
    frip_target = score_adata if score_adata is not None else adata
    snap.metrics.frip(frip_target, {'pseudobulk_frip': peakFilePath})
    per_cell = frip_target.obs['pseudobulk_frip']

    return {
        "mean_frip": float(per_cell.mean()),
        "median_frip": float(per_cell.median()),
        "total_npeak": int(total_npeak),
        "per_cell": per_cell,
        "peak_bed_path": peakFilePath,
    }


def main():
    # Record start time
    start_time = time.time()

    # Parse command line arguments
    args = parse_arguments()

    # Validate fragment file exists
    if not os.path.exists(args.fragment_file):
        print(f"Error: Fragment file '{args.fragment_file}' not found!")
        sys.exit(1)

    fragFile_localPath = args.fragment_file

    # load chsizes from file and set path to annotations:
    chSizes = read_file_as_dict(frag_genome_dict[args.genome]['sizes_path'])

    annot_path = frag_genome_dict[args.genome]['annotation_path']

    outputDir = os.path.dirname(fragFile_localPath)
    preFix = fragFile_localPath.split("/")[-1].split('.')[0]
    logFile = get_unique_log_file(outputDir, preFix)

    logging.basicConfig(
        force=True,
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s]\t%(message)s",
        handlers=[
            logging.FileHandler(logFile),
            logging.StreamHandler()
        ]
    )
    
    # Log analysis parameters
    logging.info("Starting analysis with parameters:")
    logging.info(f"  Fragment file: {fragFile_localPath}")
    logging.info(f"  Genome: {args.genome}")
    logging.info(f"  Min fragments: {args.min_fragments}")
    logging.info(f"  Number of jobs: {args.n_jobs}")
    logging.info(f"  Skip bigwig: {args.skip_bigwig}")
    logging.info(f"  Skip FRiP: {args.skip_frip}")

    # Two-step ("variant E") when the input floor (peaks + bigwig) is below the calc
    # floor (FRiP + per-cell metrics). Equal floors -> single-set "variant C".
    two_step = args.input_fragment_floor < args.min_fragments
    logging.info(f"  Input floor (peaks/bigwig): {args.input_fragment_floor}")
    logging.info(f"  Calc floor (FRiP/metrics):  {args.min_fragments}")
    logging.info(f"  Mode: {'two-step (variant E)' if two_step else 'single-set (variant C)'}")

    # Record import start time
    import_start = time.time()
    logging.info("Starting fragment import...")

    # load a fragment file from chromap
    # For shallow sequencing, may need to reduce min_frags
    min_frag = args.min_fragments
    adataDone = False
    retry_count = 0
    max_retries = 2
    while not adataDone and retry_count < max_retries:
        logging.info(f"Attempting to import fragments with {min_frag=}")
        adata = snap.pp.import_fragments(
            fragFile_localPath,
            chrom_sizes=chSizes,
            min_num_fragments=min_frag,
            sorted_by_barcode=False,
            n_jobs=args.n_jobs
        )
        if adata.obs.shape[0] != 0:
            adataDone = True
            logging.info(f"Success: Loaded {adata.obs.shape[0]} barcodes")
        else:
            min_frag = 1
            logging.warning(f"Attempting with {min_frag=}, check whether low coverage is expected")
            retry_count += 1
    if max_retries == retry_count and not adataDone:
        logging.error("Failed to load any barcodes, exiting")
        sys.exit(1)

    # Fallback fired when the calc-floor import returned 0 barcodes and min_frag was
    # dropped to 1 -> the intended calc set (>=min_fragments) is empty (coverage collapse).
    calc_floor_fallback = adataDone and (min_frag < args.min_fragments)
    if calc_floor_fallback:
        logging.warning(
            f"calc-floor fallback: no barcode reached {args.min_fragments} fragments; "
            f"imported at floor {min_frag}. Calibrated FRiP will be withheld (out of scope).")

    import_time = time.time() - import_start
    logging.info(f"Fragment import completed in {format_time(import_time)}")

    # Two-step: a separate all-(input-floor)-barcode import drives peak calling and
    # bigwig export, while `adata` (calc floor) drives FRiP and every per-cell metric.
    if two_step:
        logging.info(f"Importing input-floor set for peaks/bigwig "
                     f"(min_num_fragments={args.input_fragment_floor})...")
        adata_input = snap.pp.import_fragments(
            fragFile_localPath,
            chrom_sizes=chSizes,
            min_num_fragments=args.input_fragment_floor,
            sorted_by_barcode=False,
            n_jobs=args.n_jobs
        )
        logging.info(f"Input-floor set: {adata_input.obs.shape[0]} barcodes")
    else:
        adata_input = adata

    # Coverage metrics for the FRiP reliability gate, computed on the input-floor set
    # (all barcodes; robust to the calc-floor fallback since it counts >=min_fragments
    # directly rather than trusting the imported floor).
    n_input = adata_input.obs.shape[0]
    n_pass = int((adata_input.obs['n_fragment'] >= args.min_fragments).sum())
    frac_pass = (n_pass / n_input) if n_input else 0.0
    median_nfrag_input = float(adata_input.obs['n_fragment'].median())
    logging.info(f"frac_pass:\t{frac_pass}")
    logging.info(f"n_pass_calc_floor:\t{n_pass}")
    logging.info(f"calc_floor_fallback:\t{calc_floor_fallback}")

    # create a groupby column for the bigwig-export set (the input-floor set)
    gb = ['cell'] * adata_input.obs.shape[0]

    # Record plotting start time
    plot_start = time.time()
    logging.info("Generating plots and metrics...")

    # plot fragment size distribution
    frag_size_plot = os.path.join(outputDir, f"{preFix}_frag_size.png")
    snap.pl.frag_size_distr(adata, interactive=False, out_file=frag_size_plot, show=False)
    logging.info(f"saving plot\t{frag_size_plot}")

    # compute tss and plot
    # Use the genome specified by user
    snap.metrics.tsse(adata, annot_path, inplace=True)
    tsse_plot = os.path.join(outputDir, f"{preFix}_tsse.png")
    try:
        snap.pl.tsse(adata, interactive=False, out_file=tsse_plot, show=False)
        logging.info(f"saving plot\t{tsse_plot}")
    except ValueError as e:
        logging.warning(f"ValueError occurred plotting snapatac2:plot tsse: {e}")

    # create obs data and output summary metrics
    obs_df = adata.obs[['n_fragment', 'frac_dup', 'tsse', 'frac_mito']]

    total_fragments = obs_df.n_fragment.sum()
    per_cell_frag_fraction = obs_df['n_fragment'] / total_fragments
    per_cell_frag_wt_dups = per_cell_frag_fraction * obs_df['frac_dup']
    frag_wt_dups = per_cell_frag_wt_dups.sum()

    # fragment metric
    logging.info(f"mean nfrag:\t{obs_df['n_fragment'].mean()}")
    logging.info(f"median nfrag:\t{obs_df['n_fragment'].median()}")

    # duplicate metric
    logging.info(f"mean dup:\t{obs_df['frac_dup'].mean()}")
    logging.info(f"median dup:\t{obs_df['frac_dup'].median()}")
    logging.info(f"frag_wt_dup:\t{frag_wt_dups}")

    # mito fraction
    logging.info(f"mito_mean:\t{obs_df['frac_mito'].mean()}")
    logging.info(f"mito_median:\t{obs_df['frac_mito'].median()}")

    # do Frip calculations if not skipping
    if not args.skip_frip:
        logging.info("starting FRiP calculation")
        frip_start = time.time()

        # Peak-set construction levers. Defaults: 500 bp consensus + genome
        # blacklist + MACS3 q<=0.05. See parse_arguments for how to opt out.
        frip_groupby = getattr(args, "frip_groupby", None)
        frip_qvalue = getattr(args, "frip_macs3_qvalue", None)
        macs3_kwargs = {} if frip_qvalue is None else {"qvalue": frip_qvalue}
        if two_step:
            # pyc-like peak calling on the all-barcode input set (variant E)
            macs3_kwargs.update({"nolambda": True, "shift": 73, "extsize": 146})

        frip_merge_hw = getattr(args, "frip_merge_half_width", 250)
        if not frip_merge_hw:          # 0 or None disables merging
            frip_merge_hw = None

        if getattr(args, "frip_no_blacklist", False):
            frip_blacklist = None
        else:
            frip_blacklist = getattr(args, "frip_blacklist", None)
            if frip_blacklist is None:
                frip_blacklist = frag_genome_dict[args.genome].get("blacklist_path")
        if frip_blacklist and not os.path.exists(frip_blacklist):
            logging.warning(f"blacklist not found, skipping subtraction:\t{frip_blacklist}")
            frip_blacklist = None
        logging.info(f"FRiP peaks:\tinput_barcodes={adata_input.obs.shape[0]} "
                     f"merge_half_width={frip_merge_hw} macs3_kwargs={macs3_kwargs} "
                     f"blacklist={frip_blacklist}")

        # Peaks called on the input-floor set; FRiP scored on the calc-floor set
        # (adata) when two-step, else on the single set.
        frip_res = compute_frip(
            adata_input, outputDir, preFix,
            groupby=frip_groupby,
            merge_half_width=frip_merge_hw,
            blacklist=frip_blacklist,
            chrom_sizes=chSizes,
            macs3_kwargs=macs3_kwargs,
            score_adata=(adata if two_step else None),
        )
        peakFilePath = frip_res["peak_bed_path"]

        # frip metric
        obs_df['pseudobulk_frip'] = frip_res["per_cell"]
        per_cell_frag_wt_frip = per_cell_frag_fraction * obs_df['pseudobulk_frip']
        frag_wt_frip = per_cell_frag_wt_frip.sum()
        logging.info(f"mean frip:\t{obs_df['pseudobulk_frip'].mean()}")
        logging.info(f"median frip:\t{obs_df['pseudobulk_frip'].median()}")
        logging.info(f"frag_wt_frip:\t{frag_wt_frip}")
        logging.info(f"total_npeak:\t{frip_res['total_npeak']}")

        # FRiP (calibrated): map the summary scalar (max(mean, median)) onto the
        # pycisTopic reference. Constants are mode-selected by default (variant E vs C)
        # unless explicitly overridden on the CLI. Clamped to [0, 1].
        frip_maxmm = max(obs_df['pseudobulk_frip'].mean(),
                         obs_df['pseudobulk_frip'].median())
        calib_slope = args.frip_calib_slope
        if calib_slope is None:
            calib_slope = FRIP_CALIB_SLOPE_E if two_step else FRIP_CALIB_SLOPE
        calib_intercept = args.frip_calib_intercept
        if calib_intercept is None:
            calib_intercept = FRIP_CALIB_INTERCEPT_E if two_step else FRIP_CALIB_INTERCEPT
        frip_calibrated = min(1.0, max(0.0, calib_slope * frip_maxmm + calib_intercept))
        logging.info(f"frip_calib_params:\ta={calib_slope} b={calib_intercept}")

        # Coverage-collapse gate: withhold (NA) the calibrated FRiP when the calibration
        # is out of scope. Raw FRiP mean/median/frag_wt and total_npeak are unaffected.
        reliable, frip_flags = frip_reliability(
            calc_floor_fallback, frac_pass, median_nfrag_input,
            args.min_fragments, args.frip_gate_frac_pass_min)
        frip_calibrated_out = f"{frip_calibrated}" if reliable else "NA"
        logging.info(f"frip_reliable:\t{reliable}")
        # ';'-separated (fragReport strips ',' as a thousands separator).
        logging.info(f"frip_flags:\t{';'.join(frip_flags) if frip_flags else 'none'}")
        # frip_calibrated_raw always carries the numeric value for debugging/harness;
        # frip_calibrated is NA when withheld (what the report/summary consume).
        logging.info(f"frip_calibrated_raw:\t{frip_calibrated}")
        logging.info(f"frip_calibrated:\t{frip_calibrated_out}")

        # Plot Distribution for frip fraction
        frip_hist_plot = os.path.join(outputDir, f"{preFix}_frip_hist.png")
        plt.hist(adata.obs['pseudobulk_frip'], density=False, bins=100)
        plt.savefig(frip_hist_plot)
        plt.close()
        logging.info(f"saving plot\t{frip_hist_plot}")
        frip_time = time.time() - frip_start

    # Saving obs file
    obs_csv = os.path.join(outputDir, f"{preFix}_adata_obs.csv")
    obs_df.to_csv(obs_csv, index=True)
    logging.info(f"saving file\t{obs_csv}")
    
    # Plot Distribution for mito fraction
    mito_hist_plot = os.path.join(outputDir, f"{preFix}_mito_hist.png")
    plt.hist(adata.obs['frac_mito'], density=False, bins=100)
    plt.savefig(mito_hist_plot)
    plt.close()
    logging.info(f"saving plot\t{mito_hist_plot}")

    data = adata.obs['frac_dup']
    # Create the figure and axis
    fig, ax = plt.subplots(figsize=(8, 6))

    # Create histogram with range from 0 to 1
    n, bins, patches = ax.hist(data, bins=100, range=(0, 1), alpha=0.7, color='skyblue', edgecolor='black')

    # Set x-axis range and ticks at 0.05 intervals
    ax.set_xlim(0, 1)
    x_ticks = np.arange(0, 1.05, 0.05)  # 0.05 intervals from 0 to 1
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([f'{tick:.2f}' for tick in x_ticks], rotation=45)

    # Add vertical red dotted line at x=1
    ax.axvline(x=1, color='red', linestyle='--', linewidth=2, label='x = 1')

    # Add labels and title
    ax.set_xlabel('Cell Duplicate Rate')
    ax.set_ylabel('Counts')
    ax.set_title(f'{preFix}: Duplicate Rate Histogram')

    # Add grid for better readability
    ax.grid(True, alpha=0.3)

    # Add legend
    ax.legend()

    # Adjust layout to prevent label cutoff
    plt.tight_layout()

    dup_hist_plot = os.path.join(outputDir, f"{preFix}_dup_hist.png")
    plt.savefig(dup_hist_plot)
    logging.info(f"saving plot\t{dup_hist_plot}")

    # Display the plot
    plt.close()

    # Log tss metrics
    logging.info(f"library_tsse:\t{adata.uns['library_tsse']}")
    logging.info(f"frac_overlap_TSS:\t{adata.uns['frac_overlap_TSS']}")

    # PLot profile for TSS
    raw_tss_plot = os.path.join(outputDir, f"{preFix}_tss_profile.png")
    raw_tss_profile = adata.uns['TSS_profile']
    plt.plot(raw_tss_profile)
    plt.savefig(raw_tss_plot)
    plt.close()
    logging.info(f"saving plot\t{raw_tss_plot}")

    # Plot normalized TSS profile
    norm_tss_plot = os.path.join(outputDir, f"{preFix}_norm_tss.png")
    tss_fig, _, _, _ = tss_norm.demo_tss_normalization(raw_tss_profile, flank_size=50, window_size=1000)
    tss_fig.savefig(norm_tss_plot)

    # Plot frag count versus duplication rate
    sns.set_theme(style="dark")

    x, y = np.log10(adata.obs['n_fragment']), adata.obs['frac_dup']

    # Draw a combo histogram and scatterplot with density contours
    f, ax = plt.subplots(figsize=(6, 6))
    sns.scatterplot(x=x, y=y, s=5, color=".15")
    sns.histplot(x=x, y=y, bins=50, pthresh=.1, cmap="mako")
    sns.kdeplot(x=x, y=y, levels=5, color="w", linewidths=1)
    frag_v_dup_plot = os.path.join(outputDir, f"{preFix}_frag_v_dup.png")
    plt.savefig(frag_v_dup_plot)
    plt.close()
    logging.info(f"saving plot\t{frag_v_dup_plot}")

    plot_time = time.time() - plot_start
    logging.info(f"Plotting and metrics completed in {format_time(plot_time)}")

    # Export big wig coverage (conditional)
    if not args.skip_bigwig:
        bigwig_start = time.time()
        logging.info("Starting bigwig export...")

        bigwig_file = os.path.join(outputDir, f"{preFix}_cell.bw")
        # Bigwig from the input-floor set (all barcodes in two-step mode); gb is sized
        # to adata_input above.
        snap.ex.export_coverage(
            adata_input,
            prefix=preFix,
            groupby=gb,
            out_dir=outputDir,
            suffix='.bw',
            normalization=None,
            n_jobs=args.n_jobs
        )
        logging.info(f"saving file\t{bigwig_file}")
        
        bigwig_time = time.time() - bigwig_start
        logging.info(f"Bigwig export completed in {format_time(bigwig_time)}")
    else:
        logging.info("Skipping bigwig generation (--skip-bigwig flag specified)")
        bigwig_time = 0

    # Calculate and log total time
    total_time = time.time() - start_time
    
    logging.info("="*60)
    logging.info("ANALYSIS SUMMARY")
    logging.info("="*60)
    logging.info(f"Fragment file: {fragFile_localPath}")
    logging.info(f"Total cells analyzed: {adata.obs.shape[0]}")
    logging.info(f"Genome: {args.genome}")
    logging.info(f"Output directory: {outputDir}")
    logging.info(f"Output prefix: {preFix}")
    logging.info("")
    logging.info("TIMING BREAKDOWN:")
    logging.info(f"  Fragment import: {format_time(import_time)}")
    logging.info(f"  Plotting/metrics: {format_time(plot_time)}")
    if not args.skip_bigwig:
        logging.info(f"  Bigwig export: {format_time(bigwig_time)}")
    else:
        logging.info("  Bigwig export: skipped")
    if not args.skip_frip:
        logging.info(f"  FRiP calculation: {format_time(frip_time)}")
    else:
        logging.info("  FRiP calculation: skipped")
    logging.info(f"  Total runtime: {format_time(total_time)}")
    logging.info("="*60)


if __name__ == "__main__":
    main()
