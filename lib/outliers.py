import os
import matplotlib.pyplot as plt
import numpy as np
#import pandas as pd

from pathlib import Path
from typing import List, Tuple, Union


def get_axis_avgs(
   singlecell_path: Path,
   positions_path: Path
):
    """
    Given singlecell.csv and barcode file, return pd.DataFrame with median
    row and column counts for each lane.
    """
    import pandas as pd
    sc = pd.read_csv(singlecell_path)

    # Normalize barcode column name (some files use 'barcode' vs 'barcodes')
    if 'barcodes' not in sc.columns and 'barcode' in sc.columns:
        sc = sc.rename(columns={'barcode': 'barcodes'})

    # Normalize fragment count column name: accept 'passed_filters' or 'n_fragments'/'n_fragment'
    if 'passed_filters' not in sc.columns:
        if 'n_fragments' in sc.columns:
            sc = sc.rename(columns={'n_fragments': 'passed_filters'})
        elif 'n_fragment' in sc.columns:
            sc = sc.rename(columns={'n_fragment': 'passed_filters'})

    # Drop a leading NO_BARCODE row if present in the first column
    if sc.shape[0] > 0 and str(sc.iloc[0, 0]) == "NO_BARCODE":
        sc = sc.drop(0, axis=0)

    # Ensure barcodes column exists, find it if needed, then strip trailing '-1'
    if 'barcodes' not in sc.columns:
        barcode_col = next((col for col in sc.columns 
                           if sc[col].astype(str).str.contains(r'^[AGCT]+-1$', regex=True).any()), 
                          None)
        if barcode_col is None:
            raise ValueError("No barcode column found in singlecell CSV.")
        sc = sc.rename(columns={barcode_col: 'barcodes'})
    
    sc['barcodes'] = sc['barcodes'].astype(str).str.replace(r'-1$', '', regex=True)

    positions = pd.read_csv(positions_path, header=None)
    positions.columns = ["barcodes", "on_tissue", "row", "column"]

    # Normalize positions barcodes to match the fragments/singlecell format
    positions['barcodes'] = positions['barcodes'].astype(str).str.replace(r'-1$', '', regex=True)

    merged = pd.merge(sc, positions, how='outer')

    avgs = [merged.groupby([axis]).median(numeric_only=True)["passed_filters"]
            for axis in ["row", "column"]]

    df = pd.merge(
        avgs[0], avgs[1], right_index=True, left_index=True, how="outer"
    )
    df.columns = ["row_avg", "col_avg"]
    df.index.name = None

    # If there are less than min4Median data points, change value to NA
    min4Median = 3
    df.loc[merged.groupby(['row']).count()['passed_filters'] < min4Median,
           'row_avg'] = pd.NA
    df.loc[merged.groupby(['column']).count()['passed_filters'] < min4Median,
           'col_avg'] = pd.NA

    return df


def get_upper_bounds(
    avgs,
    sigma = 1
):
    """
    Given pd.DataFrame with median row and column counts for each lane,
    return values above with a lane average is considered an outlier for
    """
    import pandas as pd

    mean = np.mean(avgs)
    std = np.std(avgs)

    row, col = mean + sigma * std

    return row, col


def get_outliers(
    avgs,
    row_bound,
    col_bound
):
    """Add boolean column identifying outliers to lane averages;
    merger with average to return a dataframe."""
    import pandas as pd
    row_outliers = avgs["row_avg"] > row_bound
    col_outliers = avgs["col_avg"] > col_bound

    avgs = avgs.merge(row_outliers, left_index=True, right_index=True)
    avgs = avgs.merge(col_outliers, left_index=True, right_index=True)

    avgs.columns = ["row_avg", "col_avg", "row_outlier", "col_outlier"]

    return avgs


def plotting_task(
    singlecell_path: Path,
    positions_path: Path
)-> str:
    """Save a barplot.pdf containing row/column medians with outlier
    highlighted at sigma 1 and 2"""
    import pandas as pd
    output_dir = singlecell_path.parent
    dfs = []
    for sigma in [1, 2]:
        avgs = get_axis_avgs(singlecell_path, positions_path)
        row_bound, col_bound = get_upper_bounds(avgs, sigma=sigma)
        df = get_outliers(avgs, row_bound, col_bound)
        dfs.append((sigma, df))

    plt.rc("figure", figsize=(15, 10))

    fig, ax = plt.subplots(4, 1)
    plt.subplots_adjust(wspace=0, hspace=0.7)

    i = 0
    for sigma, df in dfs:

        row_x = df.index
        row_y = df["row_avg"]
        row_colors = df["row_outlier"]

        ax[i].bar(
            [x + 1 for x in row_x[~row_colors]],
            row_y[~row_colors],
            color="blue",
            edgecolor=(0, 0, 0),
        )

        ax[i].bar(
            [x + 1 for x in row_x[row_colors]],
            row_y[row_colors],
            color="red",
            edgecolor=(0, 0, 0),
            label="outliers"
        )

        ax[i].set_title(f"row medians (sigma={sigma})")
        ax[i].set_ylabel("median frag counts")
        ax[i].legend()

        i += 1

        col_x = df.index
        col_y = df["col_avg"]
        col_colors = df["col_outlier"]

        ax[i].bar(
            [x + 1 for x in col_x[~col_colors]],
            col_y[~col_colors],
            color="blue",
            edgecolor=(0, 0, 0),
        )

        ax[i].bar(
            [x + 1 for x in col_x[col_colors]],
            col_y[col_colors],
            color="red",
            edgecolor=(0, 0, 0),
            label="outlier"
        )

        ax[i].set_title(f"column medians (sigma={sigma})")
        ax[i].set_ylabel("median frag counts")
        ax[i].legend()

        xticks = [x + 1 for x in range(len(df.index))]
        font_size = 5 if len(df.index) == 96 else 9
        ax[i].set_xticks(xticks)
        ax[i].set_xticklabels(xticks, fontsize=font_size)

        i += 1
    
    plt.savefig(os.path.join(output_dir, "lane_qc.pdf"))
    plt.savefig(os.path.join(output_dir, "lane_qc.png"))

    return os.path.join(output_dir, "lane_qc.png")