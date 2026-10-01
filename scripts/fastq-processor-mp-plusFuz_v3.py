#!/usr/bin/env python3
import argparse
import csv
import gzip
import logging
import numpy as np
import os
import random
import re
import threading
import time

from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from functools import partial
from pathlib import Path
from sklearn.metrics import pairwise_distances
from typing import Dict, List, Optional, Tuple


@dataclass
class FastqRead:
    name: str
    sequence: str
    plus: str
    quality: str

    def to_string(self) -> str:
        return f"{self.name}\n{self.sequence}\n{self.plus}\n{self.quality}"


class BarcodeProcessor:

    def dna_to_numeric_vector(self, dna_sequence: str):
        mapping = {'A': [0], 'C': [1], 'G': [2], 'T': [3]}
        list_of_lists = [mapping[base] for base in dna_sequence]
        return [item for sublist in list_of_lists for item in sublist]

    def extract_barcode(self, read2_seq: str) -> Optional[str]:
        """Extract 16-mer barcode from two 8-base positions in read2."""
        if self.bc_extract == 'position':
            try:
                first_part = read2_seq[22:30]
                second_part = read2_seq[60:68]
                return first_part + second_part
            except IndexError:
                return None

        elif self.bc_extract == 'regex':
            try:
                matchString = r"(?P<cell_B>.{8})ATCCACGTGCTTGAGAGGCCAGAGCATTCG(?P<cell_A>.{8})"
                pattern = re.compile(matchString)
                bc_match = pattern.search(read2_seq)

                first_part = bc_match.group('cell_B')
                second_part = bc_match.group('cell_A')
                return first_part + second_part
            except AttributeError:
                return None

    def should_keep_read(self, barcode: str):
        """Determine if read should be kept based on sampling fraction."""
        seq_barcode = barcode  # preserve raw bc from read
        if barcode not in self.sampling_fractions:
            if not self.minDist:
                return None, seq_barcode, barcode
            else:
                # Convert query sequence to numerical vector
                query_vector = np.array(self.dna_to_numeric_vector(barcode)).reshape(1, -1)
                # Calculate pairwise distances between query and reference sequences using scikit-learn
                distances = pairwise_distances(self.reference_array, query_vector, metric='hamming')
                minIdx = np.argmin(distances)
                minErrors = distances[minIdx][0]*16
                logging.debug(f"barcode: {barcode}\tfrom matcher: {self.bc_list[minIdx]=}\tminErrors: {minErrors}")  # Output: ('apple', 90.0)
                numNeighboors = len(np.where(distances<=distances[minIdx]))-1
                if (numNeighboors == 1) and (minErrors <= self.minDist):  # a single match found
                    logging.debug(f"rescued: {barcode}->{self.bc_list[minIdx]}")
                    barcode = self.bc_list[minIdx]
                else:
                    return None, seq_barcode, barcode
        logging.debug(f"final barcode: {barcode}")
        sampling_fraction = self.sampling_fractions[barcode]
        return random.random() < sampling_fraction, seq_barcode, barcode

    def update_stats(self, local_stats: dict):
        """Thread-safe update of global statistics."""
        with self.stats_lock:
            for key, value in local_stats.items():
                self.stats[key] += value

    def __init__(self, barcode_file: str, bc_fraction_col: str, bc_extract: str, allowed_error: int):
        """Initialize with barcode sampling fractions from CSV file."""
        self.sampling_fractions: Dict[str, float] = {}
        self.stats = defaultdict(int)
        self.stats_lock = threading.Lock()
        self.bc_extract = bc_extract
        # calculate score cuttoff
        # 1 - (max edit distance allowed / query len) * 100
        # So for 16bp barcode, a single change is a 93.75 score
        self.bc_cutoff = (1 - (int(allowed_error)/16)) * 100
        self.minDist = int(allowed_error)

        # Read barcode table
        with open(barcode_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                self.sampling_fractions[row['mer']] = float(
                    row[bc_fraction_col]
                )
        self.bc_list = list(self.sampling_fractions.keys())
        # Convert reference sequences to numerical vectors and store in numpy array
        self.reference_array = np.array([self.dna_to_numeric_vector(seq) for seq in self.bc_list])
        logging.info(f"Loaded {len(self.sampling_fractions)} barcodes")
        logging.info(f"Extracting Barcode using {self.bc_extract}")


class ChunkWriter:
    def __init__(self, filename: str):
        self.filename = filename
        self.lock = threading.Lock()

    def write_chunk(self, chunk: str):
        with self.lock:
            with open(self.filename, 'a') as f:
                f.write(chunk)


def process_chunk(chunk_data: Tuple[List[str], List[str]],
                  processor: BarcodeProcessor,
                  kept_writers: Tuple[ChunkWriter, ChunkWriter],
                  filtered_writers: Tuple[ChunkWriter, ChunkWriter]) -> dict:
    """Process a chunk of reads in parallel."""
    local_stats = defaultdict(int)
    kept_chunk1 = []
    kept_chunk2 = []
    filtered_chunk1 = []
    filtered_chunk2 = []

    r1_chunk, r2_chunk = chunk_data

    for i in range(0, len(r1_chunk), 4):
        local_stats['total_reads'] += 1

        # Extract reads
        try:
            r1 = FastqRead(
                name=r1_chunk[i].strip(),
                sequence=r1_chunk[i+1].strip(),
                plus=r1_chunk[i+2].strip(),
                quality=r1_chunk[i+3].strip()
            )
            r2 = FastqRead(
                name=r2_chunk[i].strip(),
                sequence=r2_chunk[i+1].strip(),
                plus=r2_chunk[i+2].strip(),
                quality=r2_chunk[i+3].strip()
            )
        except IndexError:
            continue

        # Extract barcode
        barcode = processor.extract_barcode(r2.sequence)
        if not barcode:
            local_stats['invalid_barcode'] += 1
            continue

        # Sort into appropriate output streams
        keep_it, seqBC, foundBC = processor.should_keep_read(barcode)
        # Append barcode to read names
        r1.name = f"{r1.name}_SB::{seqBC}::CB{foundBC}"
        r2.name = f"{r2.name}_SB::{seqBC}::CB{foundBC}"
        if keep_it is None:
            local_stats['unknown_barcode'] += 1
            continue
        if keep_it:
            local_stats['kept_reads'] += 1
            kept_chunk1.append(r1.to_string() + '\n')
            kept_chunk2.append(r2.to_string() + '\n')
        else:
            local_stats['filtered_reads'] += 1
            # filtered_chunk1.append(r1.to_string() + '\n')
            # filtered_chunk2.append(r2.to_string() + '\n')

    # Write accumulated chunks
    if kept_chunk1:
        kept_writers[0].write_chunk(''.join(kept_chunk1))
        kept_writers[1].write_chunk(''.join(kept_chunk2))
    if filtered_chunk1:
        filtered_writers[0].write_chunk(''.join(filtered_chunk1))
        filtered_writers[1].write_chunk(''.join(filtered_chunk2))

    return local_stats


def chunk_reader(file_path: str, chunk_size: int = 100000):
    """Generate chunks of lines from a file."""
    open_fn = gzip.open if file_path.endswith('.gz') else open
    mode = 'rt' if file_path.endswith('.gz') else 'r'

    chunk = []
    with open_fn(file_path, mode) as f:
        for line in f:
            chunk.append(line)
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []
    if chunk:
        yield chunk


def process_fastq_files(read1_path: str, read2_path: str,
                        out1_path: str, out2_path: str,
                        filtered1_path: str, filtered2_path: str,
                        processor: BarcodeProcessor,
                        num_threads: int,
                        reads_per_chunk: int):
    """Process paired FASTQ files using multiple threads."""

    # Initialize chunk writers
    kept_writers = (ChunkWriter(out1_path), ChunkWriter(out2_path))
    filtered_writers = (
        ChunkWriter(filtered1_path), ChunkWriter(filtered2_path)
    )

    # Clear output files
    for writer in kept_writers + filtered_writers:
        open(writer.filename, 'w').close()

    # Create thread pool
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        # Process chunks in parallel
        chunk_size = int(reads_per_chunk) * 4  # lines per chunk is 4 * reads per chunk
        r1_chunks = chunk_reader(read1_path, chunk_size)
        r2_chunks = chunk_reader(read2_path, chunk_size)

        process_func = partial(
            process_chunk,
            processor=processor,
            kept_writers=kept_writers,
            filtered_writers=filtered_writers,
        )

        # Submit chunks for processing
        futures = []
        total_reads = 0
        total_chunks = 0
        for r1_chunk, r2_chunk in zip(r1_chunks, r2_chunks):
            futures.append(executor.submit(process_func, (r1_chunk, r2_chunk)))

            # Update total reads count (each chunk is 100,000 reads)
            total_reads += chunk_size // 4
            total_chunks += 1
            completed = sum(1 for f in futures if f.done())
            logging.info(f"submitted {total_chunks} chunks\tcompleted {completed}")

            # Log progress every 5 million reads
            if total_reads % 5000000 == 0:
                completed = sum(1 for f in futures if f.done())
                completed_reads = (completed * chunk_size) // 4
                logging.info(f"Processed {completed_reads:,} reads...")

        # Collect results
        for future in futures:
            local_stats = future.result()
            processor.update_stats(local_stats)


def main():
    parser = argparse.ArgumentParser(
        description='Process paired FASTQ files with barcode-based sampling')

    parser.add_argument('--read1', required=True,
                        help='Path to read1 FASTQ file')
    parser.add_argument('--read2', required=True,
                        help='Path to read2 FASTQ file')
    parser.add_argument('--barcode-table', required=True,
                        help='Path to barcode sampling fraction CSV')
    parser.add_argument('--output-prefix', required=True,
                        help='Prefix for output files')
    parser.add_argument('--threads', type=int, default=min(16, os.cpu_count()),
                        help='Number of processing threads')
    parser.add_argument('--log-level', default='INFO',
                        help='Logging level')
    parser.add_argument('--final-fraction', default='finalDownFract',
                        help='Column indicating per bc fraction')
    parser.add_argument('--max-error', default=0,
                        help='Maximum error in barcode')
    parser.add_argument('--bc-extract', default='position',
                        help='How bc bases are extracted: position | regex')
    parser.add_argument('--reads-per-chunk', default=100000,
                        help='How many reads in a chunk')

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=args.log_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler("./fastq_clean.log"),
            logging.StreamHandler()
        ]
    )

    logging.info(f"Starting processing with {args.threads} threads and {args.reads_per_chunk} reads per chunk")
    start_time = time.time()

    # Initialize processor
    processor = BarcodeProcessor(args.barcode_table, args.final_fraction, args.bc_extract, args.max_error)

    # Setup output paths
    out_dir = Path(args.output_prefix).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    out1 = f"{args.output_prefix}_cleaned_R1.fastq"
    out2 = f"{args.output_prefix}_cleaned_R2.fastq"
    filtered1 = f"{args.output_prefix}_filtered_out_R1.fastq"
    filtered2 = f"{args.output_prefix}_filtered_out_R2.fastq"

    # Process files
    process_fastq_files(
        args.read1, args.read2,
        out1, out2,
        filtered1, filtered2,
        processor,
        args.threads,
        reads_per_chunk=args.reads_per_chunk
    )

    # Log final statistics
    elapsed = time.time() - start_time
    logging.info(f"Processing complete in {elapsed/60:.1f} minutes!")
    logging.info(f"Total reads processed: {processor.stats['total_reads']:,}")
    logging.info(f"Reads kept: {processor.stats['kept_reads']:,}")
    logging.info(f"Reads filtered: {processor.stats['filtered_reads']:,}")
    logging.info(f"Unknown barcodes: {processor.stats['unknown_barcode']:,}")
    logging.info(f"Invalid barcodes: {processor.stats['invalid_barcode']:,}")


if __name__ == '__main__':
    main()
