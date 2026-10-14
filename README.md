# fastq2frags (Nextflow)

AtlasXomics epigenomic preprocessing for DBiT-seq ATAC-style data, as a
Nextflow (DSL2) pipeline. It filters reads by ligation linker, aligns them with
Chromap, builds a tabix-indexed fragment file, and generates fragment-level QC
metrics and an HTML report.

This is the Nextflow version of the Latch workflow in
[atlasxomics/atlasXomics_pipeline_wf](https://github.com/atlasxomics/atlasXomics_pipeline_wf).
It runs the same scripts and produces the same output layout.

Read layout ([Zhang et al. 2023](https://www.nature.com/articles/s41586-023-05795-1)):

- read 1: genomic sequence
- read 2: linker1 | barcodeA | linker2 | barcodeB | genomic sequence

## Steps

| Process | What it does |
| --- | --- |
| `CLEAN_FASTQ` | Optional. Cleans reads against a correction table (`--cleaning`). |
| `FILTER_LINKERS` | Keeps read pairs whose read 2 matches ligation linker 1 and linker 2 with bbduk (at most 3 mismatches). In bulk modes, it then assigns random barcodes. |
| `CHROMAP_ALIGN` | Aligns with Chromap and extracts barcodes from read 2. |
| `PROCESS_FRAGMENTS` | Appends `-1` to barcodes, removes out-of-bound fragments and, when enabled, blacklisted ones, then bgzips and tabix-indexes the fragment file. |
| `XTALK_CORRECT` | Optional. Cross-talk correction (`--xtalk_correction`). |
| `FRAGMENT_METRICS` | Runs SnapATAC2 QC, calls pseudobulk MACS3 peaks, computes FRiP, exports the bigwig, and builds the lane QC plot and the HTML report. |
| `SLIMS_UPLOAD` | Optional. Pushes run metrics to SLIMS (`--upload_to_slims`). |

## Build the container

All processes run in one image, built from the repository root:

```bash
docker build -f containers/Dockerfile -t atx-fastq2frags:$(cat version) .
```

The pipeline uses `atx-fastq2frags:<version>` by default. After pushing the
image to a registry, point the pipeline at it with
`--container <registry>/atx-fastq2frags:<version>`.

### Published images (ECR)

The "Publish container" GitHub Actions workflow
(`.github/workflows/publish-container.yml`) builds the image on GitHub's x86
runners and pushes it to Amazon ECR. Pushing a tag that matches the `version`
file publishes two images:

| Tag | Image |
|---|---|
| `<version>` | Release image, for real runs. |
| `<version>-local` | Chromap sort buffer cut to 1 GiB, for local tests only. |

```bash
git tag "v$(cat version)"
git push origin "v$(cat version)"
```

Running the workflow by hand from the Actions tab publishes
`<version>-<short sha>` and `<version>-<short sha>-local` instead, so release
tags are never overwritten. One-time AWS setup is in `.github/aws/README.md`.

To pull from ECR:

```bash
aws ecr get-login-password --region us-east-1 | \
  docker login --username AWS --password-stdin <account>.dkr.ecr.us-east-1.amazonaws.com
docker pull <account>.dkr.ecr.us-east-1.amazonaws.com/atx-fastq2frags:<version>
```

## Run

Single sample:

```bash
nextflow run main.nf -profile docker \
    --run_id D01234_NG01234 \
    --r1 sample_R1_001.fastq.gz \
    --r2 sample_R2_001.fastq.gz \
    --genome_dir /refs/Chromap_references/Human \
    --barcode_file bc220_25-APR.txt \
    --outdir results
```

For several samples, pass a samplesheet with `--input` in place of `--run_id`, `--r1` and `--r2`:

```bash
nextflow run main.nf -profile docker --input samplesheet.csv --genome_dir ... --outdir results
```

```csv
run_id,r1,r2,ng_id,ds_table
D01234_NG01234,s3://bucket/D01234_R1.fastq.gz,s3://bucket/D01234_R2.fastq.gz,,
D01235,s3://bucket/D01235_R1.fastq.gz,s3://bucket/D01235_R2.fastq.gz,NG01235,s3://bucket/xtalkTable.csv
```

The `ng_id` and `ds_table` columns are optional. When a cell is empty, the
pipeline uses the value of `--ng_id` or `--ds_table`.

Inputs can be local paths or `s3://`, `gs://`, `az://` or `https://` URLs.
Choose an executor (AWS Batch, Slurm, Kubernetes, and so on) in your own config
or profile.

### Resources

Resource requests are set per label in `nextflow.config`:

| Label | Processes | Default |
| --- | --- | --- |
| `process_high` | CLEAN_FASTQ, FILTER_LINKERS, CHROMAP_ALIGN, XTALK_CORRECT | 42 CPU / 192 GB |
| `process_metrics` | FRAGMENT_METRICS | 16 CPU / 64 GB |
| `process_medium` | PROCESS_FRAGMENTS | 8 CPU / 32 GB |
| `process_single` | SLIMS_UPLOAD | 1 CPU / 2 GB |

On a small machine, add `-profile docker,test` to cap every process at 4 CPU / 12 GB. For laptop test runs, see [Local testing](#local-testing).

Chromap needs more than 12 GB at any input size. Its sort/dedupe step
allocates a fixed 10 GiB buffer (`src/mapping_writer.h`) on top of the genome
index, which is several GB for a full human or mouse genome. Give
`CHROMAP_ALIGN` at least 16 GB for small references and about 32 GB for full
genomes.

## Local testing

A small test runs on a laptop, including a 16 GB Apple Silicon Mac with
Docker Desktop memory set to about 12 GB. It covers the pipeline's code paths,
not real analysis. The test reference is only hg38 chr21 + chrM, so just
about 3% of reads align.

1. Build an image with a 1 GiB Chromap sort buffer instead of 10 GiB. On Apple
   Silicon, keep `--platform linux/amd64`: Chromap needs x86 instructions, so
   the image runs under emulation and is slower.

   ```bash
   docker build --platform linux/amd64 \
       --build-arg CHROMAP_SORT_BUFFER_GIB=1 \
       -f containers/Dockerfile -t atx-fastq2frags:local .
   ```

2. Download the test data and build its Chromap index. This takes about 650 MB
   of disk and downloads about 350 MB.

   ```bash
   tests/get_test_data.sh test_data
   ```

3. Run the tests below. The `local` profile caps every process at 4 CPU / 10 GB.
   With the 1 GiB buffer, Chromap peaks at about 4.7 GB and the other steps
   stay under 1 GB.

   ```bash
   # Single sample
   nextflow run main.nf -profile docker,local --container atx-fastq2frags:local \
       --run_id demo --r1 test_data/demo_R1.fastq.gz --r2 test_data/demo_R2.fastq.gz \
       --genome_dir test_data/ref --barcode_file bc50.txt --outdir test_results

   # Samplesheet: two samples, cross-talk correction on the second only
   nextflow run main.nf -profile docker,local --container atx-fastq2frags:local \
       --input test_data/samplesheet.csv --xtalk_correction \
       --genome_dir test_data/ref --barcode_file bc50.txt --outdir test_results

   # Bulk mode, linker 2 filtering skipped, no blacklist filtering
   nextflow run main.nf -profile docker,local --container atx-fastq2frags:local \
       --run_id demobulk --r1 test_data/demo_R1.fastq.gz --r2 test_data/demo_R2.fastq.gz \
       --genome_dir test_data/ref --barcode_file bc50.txt \
       --bulk --skip_linker2 --blacklist_filtering false --outdir test_results
   ```

Use the small-buffer image only for tests. Real runs should use the default
build (10 GiB buffer) on a machine with at least 32 GB of RAM.

### Test data

| Data | Source | Details |
| --- | --- | --- |
| `demo_R1.fastq.gz`, `demo_R2.fastq.gz` | First 2,000,000 read pairs of `s3://latch-public/test-data/13502/atx_demo_R{1,2}_001.fastq.gz` (public) | ATX demo run with 50x50 barcodes (`bc50.txt`). Read 2 is 151 bp. |
| `ref/GRCh38_chr21_chrM.fa` | hg38 `chr21.fa` + `chrM.fa` from AWS iGenomes (`s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/`) | 48 MB. The `GRCh38_` prefix makes the pipeline use the GRCh38 chrom sizes, blacklist and annotation. |
| `ref/GRCh38_chr21_chrM.index` | Built with `chromap -i` | 291 MB. |
| `xtalkTable.csv` | Generated from the public `x50_all_tissue_positions_list.csv` | Synthetic. Every 10th spot is marked affected (`xtalk_mask` 0, `fold_above_background` 3.5), so 250 of 2,500 barcodes. It only exercises the code path. |
| Positions list | `x50_all_tissue_positions_list.csv`, which the pipeline downloads | Used for the lane QC plot. |

### Expected results

These are the results from the runs used to validate the pipeline. Counts
can shift slightly between Chromap builds.

| Test | Checks | Result |
| --- | --- | --- |
| Single sample | All steps end to end and the full output layout | 2M pairs → 1,796,345 pass linker 1 → 1,764,259 pass linker 2. Chromap: 105,204 mapped reads, 96.4% valid barcodes. 31,037 raw fragments → about 29,180 after out-of-bound and blacklist filtering, across 2,395 barcodes. About 12,350 peaks, TSS enrichment 8.39. Calibrated FRiP is `NA` because the data is too shallow (`low_frac_pass;low_median_nfrag`). |
| Samplesheet + `--xtalk_correction` | Per-sample branching; cross-talk correction with the compiled C engine | `demo` logs a warning and skips correction (no table). `demoxt`: 250 affected barcodes substituted (200 by row+1, 50 by row-1). `fragments.final.tsv.gz` is published at the top level. |
| `--bulk --skip_linker2 --blacklist_filtering false`, as UID 1000 | Bulk barcodes, linker skip, no blacklist, non-root container user | Random barcodes written into read 2, empty linker 2 stats file, all outputs present. |

Not tested: `--cleaning` (no table with `finalDownFract` was available),
`--no_ligation_bulk`, `--upload_to_slims`, and full-genome references.

## Parameters

| Parameter | Default | Description |
| --- | --- | --- |
| `--input` | | Samplesheet CSV (`run_id,r1,r2[,ng_id][,ds_table]`). |
| `--run_id`, `--r1`, `--r2` | | Single-sample input, used when `--input` is not given. The run ID cannot start with or contain `/`, `-`, `.` or a space. |
| `--genome_dir` | required | Directory containing the Chromap `*.fa` and `*.index`. |
| `--genome_id` | inferred | Overrides the genome ID inferred from the `*.fa` name (`GRCh38`, `GRCm38`, `GRCm39`, `Rnor6`, `mm10`, `mm39`). |
| `--barcode_file` | `bc220_25-APR.txt` | Barcode whitelist: a file name in `barcodes/`, or a path. |
| `--positions_file` | matches barcode file | Spatial positions list for lane QC. Needed only with custom barcode files. |
| `--skip_linker1`, `--skip_linker2` | `false` | Skip bbduk linker filtering. |
| `--bulk`, `--no_ligation_bulk` | `false` | Assign random barcodes, for bulk runs. |
| `--blacklist_filtering` | `true` | Remove fragments that overlap ENCODE blacklist regions. |
| `--ds_table` | | Correction table for cleaning and cross-talk correction. |
| `--cleaning` | `false` | Clean FASTQs against `ds_table`. Needs columns `mer` and `finalDownFract`. |
| `--xtalk_correction` | `false` | Cross-talk correction. Needs columns `mer`, `row`, `col`, `xtalk_mask` and `fold_above_background`. |
| `--xtalk_threshold` | `2.0` | Minimum `fold_above_background` for a barcode to be substituted. |
| `--min_fragments` | `100` | Fragment floor for FRiP and per-cell metrics. |
| `--input_fragment_floor` | `1` | Fragment floor for peak calling and the bigwig. |
| `--upload_to_slims` | `false` | Push `summary.csv` metrics to SLIMS. |
| `--ng_id` | `NG00000` | SLIMS content ID. When unset, it is taken from the `_NGxxxxx` suffix of `run_id`. |
| `--outdir` | `results` | Output root. |
| `--container` | `atx-fastq2frags:<version>` | Container image for all processes. |

If you turn on cleaning or cross-talk correction without a `ds_table`, the
pipeline logs a warning and skips that step. A `ds_table` that is missing
required columns stops the run.

`nextflow_schema.json` describes these parameters for Seqera Platform's launch
form. When you add or change a parameter, update it there too.

### SLIMS credentials

```bash
nextflow secrets set SLIMS_USERNAME <user>
nextflow secrets set SLIMS_PASSWORD <password>
nextflow run main.nf ... --upload_to_slims
```

## Outputs

```
<outdir>/<run_id>/
├── fragments.tsv.gz                  # top-level copy (fragments.final.tsv.gz when xtalk-corrected)
├── fragment_analysis_report.html     # top-level copy
├── fragmentscell.bw                  # top-level copy
├── filtered_fastqs/                  # <run_id>_linker2_R{1,2}.fastq.gz, <run_id>_l{1,2}_stats.txt
├── chromap_output/                   # aln.bed, chromap_log.txt, fragments(.final).tsv.gz(.tbi)
└── run_metrics/                      # see below
<outdir>/pipeline_info/               # Nextflow timeline, report, trace
```

| `chromap_output/` file | Description |
| --- | --- |
| `aln.bed` | Raw Chromap output, before barcode suffixing and filtering. |
| `chromap_log.txt` | Chromap log with alignment, mapping and barcode statistics. |
| `fragments.tsv.gz` (+ `.tbi`) | Final fragment file: `-1` barcode suffixes, out-of-bound fragments removed and, when enabled, blacklisted fragments removed. |
| `fragments.final.tsv.gz` (+ `.tbi`) | Cross-talk-corrected fragment file. Present only when `--xtalk_correction` runs. |

| `run_metrics/` file | Description |
| --- | --- |
| `fragment_analysis_report.html` | Self-contained HTML report summarizing fragment metrics, with the key QC plots embedded. |
| `summary.csv` | Run-level metrics, the same values the report shows and SLIMS receives. |
| `run_id.txt` | The `run_id` used to label the report. |
| `chromap_log.txt` | Copy of the Chromap log, used by the report. |
| `fragments_log.txt` | SnapATAC2 analysis log: parameters, metric summaries and timings. |
| `fragments_adata_obs.csv` | Per-barcode QC table: fragment count, duplicate fraction, TSS enrichment, mitochondrial fraction and FRiP. |
| `fragments_peakList.bed.gz` | Pseudobulk MACS3 peak calls, used for FRiP. |
| `fragmentscell.bw` | BigWig coverage track. |
| `fragments_frag_size.png` | Fragment size distribution. |
| `fragments_frip_hist.png` | Histogram of per-barcode FRiP. |
| `fragments_mito_hist.png` | Histogram of mitochondrial fraction. |
| `fragments_dup_hist.png` | Histogram of duplicate fraction. |
| `fragments_tss_profile.png` | Raw aggregate TSS profile. |
| `fragments_norm_tss.png` | Normalized TSS profile. |
| `fragments_frag_v_dup.png` | Fragment count against duplicate fraction. |
| `lane_qc.png`, `lane_qc.pdf` | Row and column median fragment counts per lane, with outlier lanes highlighted. |

## Repository layout

| Path | Contents |
| --- | --- |
| `main.nf`, `nextflow.config` | Pipeline entry point and configuration. |
| `modules/local/` | One file per process. |
| `bin/` | Small CLIs that the processes call. Nextflow puts them on `PATH`. |
| `lib/` | Fragment filtering (`oob.py`) and lane QC (`outliers.py`) helpers, installed in the image. |
| `scripts/` | Analysis scripts: cleaning, bulk barcodes, cross-talk substitution, SnapATAC2 metrics and the report. |
| `barcodes/`, `blacklist/`, `chrom_sizes/`, `references/` | Reference assets, staged or baked into the image. |
| `containers/Dockerfile` | Image for all processes. |
| `tests/get_test_data.sh` | Downloads the local test dataset (see Local testing). |

## Differences from the Latch workflow

- The Latch Registry upsert and execution renaming are not included, because
  both depend on the Latch SDK.
- A failed Chromap run fails the pipeline. The Latch task only logged the
  error and continued.

Questions? Contact support@atlasxomics.com.
