#!/usr/bin/env python3
"""Write lane_qc.png/.pdf (row/column median fragment counts with outlier lanes
highlighted) next to the per-barcode QC table. Wraps lib/outliers.py."""
import argparse
from pathlib import Path

from outliers import plotting_task


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("singlecell", help="fragments_adata_obs.csv")
    parser.add_argument("positions", help="tissue_positions_list.csv for the barcode set")
    args = parser.parse_args()
    plotting_task(Path(args.singlecell).resolve(), Path(args.positions).resolve())


if __name__ == "__main__":
    main()
