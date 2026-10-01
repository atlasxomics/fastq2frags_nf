/*
 * SnapATAC2 fragment QC, pseudobulk peaks, bigwig, lane QC plot, and the HTML
 * report. Everything is written to run_metrics/.
 */
process FRAGMENT_METRICS {
    tag "${meta.id}"
    label 'process_metrics'

    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode,
        saveAs: { fn -> fn == 'run_metrics' ? fn : null }
    // Top-level copies of the primary deliverables.
    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode,
        saveAs: { fn -> fn.startsWith('top_level/') ? fn.tokenize('/')[-1] : null }

    input:
    tuple val(meta), path(fragments), path(index), path(chromap_log)
    path positions

    output:
    tuple val(meta), path("run_metrics"),             emit: metrics_dir
    tuple val(meta), path("run_metrics/summary.csv"), emit: summary
    path "top_level/*",                               emit: top_level

    script:
    // The summary script writes its outputs next to the fragment file.
    """
    mkdir run_metrics
    ln -s ../${fragments} run_metrics/${fragments}
    cp -L ${chromap_log} run_metrics/chromap_log.txt
    echo "${meta.id}" > run_metrics/run_id.txt

    python /root/scripts/atac_analysis_summary_prep.py \\
        run_metrics/${fragments} \\
        --genome ${meta.genome_id} \\
        --min-fragments ${params.min_fragments} \\
        --input-fragment-floor ${params.input_fragment_floor} \\
        --n-jobs ${task.cpus}

    lane_qc.py run_metrics/fragments_adata_obs.csv ${positions}

    python /root/scripts/fragReport.py run_metrics

    rm run_metrics/${fragments}

    # Nextflow publishes files nested in an output directory only as part of
    # that directory, so stage the top-level copies separately.
    mkdir top_level
    cp run_metrics/fragment_analysis_report.html top_level/
    bw=\$(ls run_metrics/*cell.bw 2>/dev/null | head -n 1 || true)
    if [ -n "\$bw" ]; then
        cp "\$bw" top_level/fragmentscell.bw
    else
        echo "WARNING: no '*cell.bw' bigwig found to copy" >&2
    fi
    """
}
