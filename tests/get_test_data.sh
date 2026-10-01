#!/usr/bin/env bash
# Download a small test dataset and build a matching Chromap index.
#
#   tests/get_test_data.sh [outdir] [read_pairs]
#
# outdir      defaults to test_data/
# read_pairs  defaults to 2000000 (first N pairs of the public ATX demo run)
#
# Produces:
#   <outdir>/demo_R1.fastq.gz, demo_R2.fastq.gz   demo reads (bc50 barcodes)
#   <outdir>/ref/GRCh38_chr21_chrM.{fa,index}    hg38 chr21 + chrM reference
#   <outdir>/xtalkTable.csv                       synthetic cross-talk table
#   <outdir>/samplesheet.csv                      2 samples, one with xtalkTable
#
# Needs curl, gzip and docker. The index is built with the pipeline image
# ($IMAGE, default atx-fastq2frags:local).
set -euo pipefail

OUT=${1:-test_data}
PAIRS=${2:-2000000}
IMAGE=${IMAGE:-atx-fastq2frags:local}

DEMO=https://latch-public.s3.amazonaws.com/test-data/13502
IGENOMES=https://ngi-igenomes.s3.amazonaws.com/igenomes/Homo_sapiens/UCSC/hg38/Sequence/Chromosomes

mkdir -p "$OUT/ref"
OUT=$(cd "$OUT" && pwd)

for r in R1 R2; do
    f="$OUT/demo_${r}.fastq.gz"
    if [ -s "$f" ]; then echo "exists: $f"; continue; fi
    echo "Downloading first ${PAIRS} reads of demo ${r}..."
    # head closes the pipe early, so curl/gzip exiting non-zero is expected here.
    (set +o pipefail; curl -sS "$DEMO/atx_demo_${r}_001.fastq.gz" 2>/dev/null \
        | gzip -dc 2>/dev/null | head -n $((PAIRS * 4)) | gzip -1 > "$f")
    n=$(gzip -dc "$f" | wc -l | tr -d ' ')
    [ "$n" -eq $((PAIRS * 4)) ] || { echo "ERROR: $f has $n lines, expected $((PAIRS * 4))" >&2; exit 1; }
done

FA="$OUT/ref/GRCh38_chr21_chrM.fa"
if [ ! -s "$FA" ]; then
    echo "Downloading hg38 chr21 + chrM..."
    { curl -fsS "$IGENOMES/chr21.fa"; curl -fsS "$IGENOMES/chrM.fa"; } > "$FA.tmp"
    mv "$FA.tmp" "$FA"
fi

IDX="$OUT/ref/GRCh38_chr21_chrM.index"
if [ ! -s "$IDX" ]; then
    echo "Building Chromap index with $IMAGE..."
    docker run --rm -v "$OUT/ref:/ref" "$IMAGE" \
        chromap -i -r /ref/GRCh38_chr21_chrM.fa -o /ref/GRCh38_chr21_chrM.index
fi

# Synthetic cross-talk table on the 50x50 grid: every 10th spot is "affected"
# (xtalk_mask 0, fold 3.5). Only meant to exercise the XTALK_CORRECT path.
echo "Writing xtalkTable.csv and samplesheet.csv..."
curl -fsS "$DEMO/x50_all_tissue_positions_list.csv" | tr -d '\r' | awk -F, '
    BEGIN { OFS = ","; print "idx,mer,on_tissue,row,col,xtalk_mask,note,fold_above_background" }
    { m = (NR % 10 == 0) ? 0 : 1; print NR - 1, $1, $2, $3, $4, m, "synthetic", (m == 0 ? 3.5 : 1.0) }
' > "$OUT/xtalkTable.csv"

cat > "$OUT/samplesheet.csv" <<EOF
run_id,r1,r2,ng_id,ds_table
demo,$OUT/demo_R1.fastq.gz,$OUT/demo_R2.fastq.gz,,
demoxt,$OUT/demo_R1.fastq.gz,$OUT/demo_R2.fastq.gz,NG00001,$OUT/xtalkTable.csv
EOF

echo "Done: $OUT"
