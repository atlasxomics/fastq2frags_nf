process XTALK_CORRECT {
    tag "${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}/${meta.id}/chromap_output", mode: params.publish_dir_mode
    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'fragments.final.tsv.gz' ? fn : null }

    input:
    tuple val(meta), path(fragments), path(index), path(ds_table)

    output:
    tuple val(meta), path("fragments.final.tsv.gz"), path("fragments.final.tsv.gz.tbi"), emit: fragments

    script:
    """
    xtalk_correct.py \\
        ${ds_table} \\
        ${fragments} \\
        fragments.final.tsv.gz \\
        --threshold ${params.xtalk_threshold} \\
        --threads ${task.cpus}
    """
}
