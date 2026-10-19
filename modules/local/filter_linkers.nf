/*
 * Keep read pairs whose read 2 carries ligation linker 1 and linker 2 (bbduk,
 * <=3 mismatches). In bulk / no-ligation bulk mode, read 2 barcodes are then
 * replaced with random barcodes from the whitelist.
 */
process FILTER_LINKERS {
    tag "${meta.id}"
    label 'process_high'

    publishDir path: { "${params.outdir}/${meta.id}/filtered_fastqs" }, mode: params.publish_dir_mode

    input:
    tuple val(meta), path(r1), path(r2)
    path barcodes

    output:
    tuple val(meta), path("${meta.id}_linker2_R1.fastq.gz"), path("${meta.id}_linker2_R2.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}_l1_stats.txt"), path("${meta.id}_l2_stats.txt"),                emit: stats

    script:
    def prefix = meta.id
    // Leave headroom for the JVM on top of the heap.
    def xmx = "-Xmx${Math.max(1, (task.memory.toGiga() * 0.85) as int)}g"
    // ordered=t keeps input read order. With several threads bbduk otherwise
    // writes reads in a varying order, and Chromap's output depends on it.
    def bbduk_common = "skipr1=t threads=${task.cpus} ordered=t ${xmx} k=30 mm=f rcomp=f hdist=3"
    // Bulk modes write the barcode-replaced read 2 under the published R2 name.
    def l2_r2 = (params.bulk || params.no_ligation_bulk) ? "linker2_R2.prebulk.fastq.gz" : "${prefix}_linker2_R2.fastq.gz"
    def bulk_flags = [params.bulk ? '-b' : '', params.no_ligation_bulk ? '-nl' : ''].join(' ').trim()
    """
    # If one of bbduk's threads fails (e.g. a writer when the disk fills),
    # bbduk doesn't exit: its main thread waits for that thread forever. Pass
    # its output through to stderr and kill it at the first Java exception,
    # so the task fails instead of hanging.
    run_bbduk() {
        bbduk.sh "\$@" 2>&1 | while IFS= read -r line; do
            printf '%s\\n' "\$line" >&2
            case "\$line" in
                "Exception in thread"*)
                    echo "bbduk: a thread failed, stopping bbduk" >&2
                    pkill -f jgi.BBDuk || true ;;
            esac
        done
        return "\${PIPESTATUS[0]}"
    }

    if [ "${params.skip_linker1}" = "true" ]; then
        echo "linker1 filtering skipped"
        ln -s ${r1} linker1_R1.fastq.gz
        ln -s ${r2} linker1_R2.fastq.gz
        touch ${prefix}_l1_stats.txt
    else
        run_bbduk \\
            in1=${r1} in2=${r2} \\
            outm1=linker1_R1.fastq.gz outm2=linker1_R2.fastq.gz \\
            ${bbduk_common} restrictleft=103 \\
            stats=${prefix}_l1_stats.txt \\
            literal=GTGGCCGATGTTTCGCATCGGCGTACGACT
    fi

    if [ "${params.skip_linker2}" = "true" ]; then
        echo "linker2 filtering skipped"
        cp -L linker1_R1.fastq.gz ${prefix}_linker2_R1.fastq.gz
        cp -L linker1_R2.fastq.gz ${l2_r2}
        touch ${prefix}_l2_stats.txt
    else
        run_bbduk \\
            in1=linker1_R1.fastq.gz in2=linker1_R2.fastq.gz \\
            outm1=${prefix}_linker2_R1.fastq.gz outm2=${l2_r2} \\
            ${bbduk_common} restrictleft=65 \\
            stats=${prefix}_l2_stats.txt \\
            literal=ATCCACGTGCTTGAGAGGCCAGAGCATTCG
    fi

    if [ -n "${bulk_flags}" ]; then
        python /root/scripts/bc_process_newbulk.py \\
            -i ${l2_r2} \\
            -o2 ${prefix}_linker2_R2.fastq.gz \\
            -bcf ${barcodes} \\
            ${bulk_flags}
        rm ${l2_r2}
    fi

    rm -f linker1_R1.fastq.gz linker1_R2.fastq.gz
    """
}
