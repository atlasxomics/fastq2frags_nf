process CLEAN_FASTQ {
    tag "${meta.id}"
    label 'process_high'

    input:
    tuple val(meta), path(r1), path(r2), path(ds_table)

    output:
    tuple val(meta), path("${meta.id}_cleaned_R1.fastq.gz"), path("${meta.id}_cleaned_R2.fastq.gz"), emit: reads

    script:
    """
    python /root/scripts/fastq-processor-mp-plusFuz_v3.py \\
        --read1 ${r1} \\
        --read2 ${r2} \\
        --barcode-table ${ds_table} \\
        --output-prefix ${meta.id} \\
        --threads ${task.cpus} \\
        --max-error 0 \\
        --bc-extract regex \\
        --reads-per-chunk 800000

    pigz -p ${task.cpus} ${meta.id}_cleaned_R1.fastq ${meta.id}_cleaned_R2.fastq
    """
}
