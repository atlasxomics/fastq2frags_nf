/*
 * Turn Chromap output into the final fragment file: suffix barcodes with -1,
 * drop out-of-bound fragments, optionally drop blacklisted fragments, then
 * bgzip + tabix.
 */
process PROCESS_FRAGMENTS {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/${meta.id}/chromap_output", mode: params.publish_dir_mode
    // Top-level copy of the final fragment file, unless cross-talk correction replaces it.
    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode,
        saveAs: { fn -> (!meta.xtalk && fn == 'fragments.tsv.gz') ? fn : null }

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("fragments.tsv.gz"), path("fragments.tsv.gz.tbi"), emit: fragments

    script:
    def blacklist_flag = params.blacklist_filtering ? '' : '--no-blacklist'
    """
    # Append -1 to barcodes (column 4)
    awk 'BEGIN{FS=OFS=" "}{\$4=\$4"-1"}4' ${bed} > suffixed.bed

    filter_fragments.py suffixed.bed fragments.tsv \\
        --genome-id ${meta.genome_id} \\
        ${blacklist_flag}
    rm suffixed.bed

    bgzip -@ ${task.cpus} -c fragments.tsv > fragments.tsv.gz
    tabix -p bed -f fragments.tsv.gz
    rm fragments.tsv
    """
}
