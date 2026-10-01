import logging


logging.basicConfig(
    format="%(levelname)s - %(asctime)s - %(message)s", level=logging.INFO
)

chromsize_paths = {
    "GRCh38":  "chrom_sizes/hg38_chrom_sizes.txt",
    "GRCm38": "chrom_sizes/mm10_chrom_sizes.txt",
    "mm10": "chrom_sizes/mm10_chrom_sizes.txt",
    "GRCm39": "chrom_sizes/mm39_chrom_sizes.txt",
    "mm39": "chrom_sizes/mm39_chrom_sizes.txt",
    "Rnor6": "chrom_sizes/rn6_chrom_sizes.txt"
}

blacklist_paths = {
    "GRCh38": "blacklist/hg38-blacklist.v2.bed",
    "hg38": "blacklist/hg38-blacklist.v2.bed",
    "hg19": "blacklist/hg19-blacklist.v2.bed",
    "GRCm38": "blacklist/mm10-blacklist.v2.bed",
    "mm10": "blacklist/mm10-blacklist.v2.bed",
    "GRCm39": "blacklist/mm39-blacklist.v2.bed",
    "mm39": "blacklist/mm39-blacklist.v2.bed",
    "dm6": "blacklist/dm6-blacklist.v2.bed",
}


def call_oob(row, chromsizes):
    import pandas as pd
    chrom, start, end = row["chrom"], row["start"], row["end"]
    if chrom in chromsizes:
        chrom_length = chromsizes[chrom]
        if start > chrom_length or end > chrom_length:
            logging.warning(
                f"Found one beyond {chrom=}, {chrom_length=}: \
                {(chrom, start, end)}"
            )
            return False
        else:
            return True
    else:
        logging.warning(f"Wrong chromosome found: {(chrom, start, end)}")
        return False


def filter_oob(
    input_file: str,
    output_file: str,
    chromsizes: dict,
    chunksize: int = 10000000
):
    """Reads a large input file in chunks and applies call_oob() filtering.
    Instead of keeping everything in memory, writes results directly to a file.
    """
    import pandas as pd
    columns = ["chrom", "start", "end", "cellBarcode", "duplicates"]
    first_write = True  # Track whether it's the first chunk

    for chunk in pd.read_csv(
        input_file,
        sep=" ",
        header=None,
        names=columns,
        chunksize=chunksize
    ):

        # Apply the filtering condition to each chunk
        condition = chunk.apply(lambda row: call_oob(row, chromsizes), axis=1)
        filtered_chunk = chunk[condition]

        # Write filtered data directly to a file
        filtered_chunk.to_csv(
            output_file,
            sep="\t",
            index=False,
            header=False,
            mode="w" if first_write else "a"
        )
        first_write = False
        logging.info(f"Checked {chunksize} lines for out-of-bound fragments...")

    logging.info(f"Out-of-bound cleaned data saved to {output_file}")


def load_chromsizes(path: str) -> dict:
    sizes = {}
    with open(path, 'r') as f:
        for line in f:
            chrom, length = line.strip().split()
            sizes[chrom] = int(length)

    return sizes


def filter_blacklist(
    input_file: str,
    output_file: str,
    blacklist_file: str,
    chunksize: int = 2_000_000
):
    """Remove fragments overlapping ENCODE blacklist regions.

    Operates in chunks to avoid loading the full fragment file in memory.
    """
    import shutil
    import pandas as pd
    import pyranges as pr

    blacklist_df = pd.read_csv(
        blacklist_file,
        sep="\t",
        comment="#",
        header=None,
        usecols=[0, 1, 2],
        names=["Chromosome", "Start", "End"]
    )

    blacklist_df = blacklist_df.dropna()
    blacklist_df["Start"] = blacklist_df["Start"].astype(int)
    blacklist_df["End"] = blacklist_df["End"].astype(int)

    if blacklist_df.empty:
        logging.warning(
            f"Blacklist file is empty ({blacklist_file}). Skipping blacklist filtering."
        )
        shutil.copyfile(input_file, output_file)
        return

    blacklist_pr = pr.PyRanges(blacklist_df)

    columns = ["Chromosome", "Start", "End", "cellBarcode", "duplicates"]
    row_id_col = "__row_id"
    row_offset = 0
    first_write = True

    for chunk in pd.read_csv(
        input_file,
        sep="\t",
        header=None,
        names=columns,
        chunksize=chunksize
    ):
        if chunk.empty:
            continue

        chunk["Start"] = chunk["Start"].astype(int)
        chunk["End"] = chunk["End"].astype(int)
        chunk = chunk[chunk["Start"] < chunk["End"]]

        if chunk.empty:
            continue

        # Preserve global row order after interval subtraction; tabix expects
        # chromosome blocks to remain continuous in the original sorted stream.
        chunk[row_id_col] = range(row_offset, row_offset + len(chunk))
        chunk_pr = pr.PyRanges(chunk)
        filtered = chunk_pr.overlap(blacklist_pr, invert=True).as_df()

        if filtered.empty:
            row_offset += len(chunk)
            continue

        filtered = filtered.sort_values(row_id_col)
        filtered = filtered[columns]
        filtered.to_csv(
            output_file,
            sep="\t",
            index=False,
            header=False,
            mode="w" if first_write else "a"
        )
        first_write = False
        row_offset += len(chunk)

    if first_write:
        open(output_file, "w").close()

    logging.info(f"Blacklist-cleaned data saved to {output_file}")
