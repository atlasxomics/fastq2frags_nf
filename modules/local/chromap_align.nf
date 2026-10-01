process CHROMAP_ALIGN {
    tag "${meta.id}"
    label 'process_high'

    publishDir path: { "${params.outdir}/${meta.id}/chromap_output" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(r1), path(r2)
    path fasta
    path index
    path barcodes

    output:
    tuple val(meta), path("aln.bed"),         emit: bed
    tuple val(meta), path("chromap_log.txt"), emit: log

    script:
    """
    set -o pipefail

    chromap \\
        -t ${task.cpus} \\
        --preset atac \\
        -x ${index} \\
        -r ${fasta} \\
        -1 ${r1} \\
        -2 ${r2} \\
        -o aln.bed \\
        -b ${r2} \\
        --barcode-whitelist ${barcodes} \\
        --read-format bc:22:29,bc:60:67,r1:0:-1,r2:117:-1 \\
        2>&1 | tee chromap_log.txt
    """
}
