process SLIMS_UPLOAD {
    tag "${meta.id}"
    label 'process_single'

    secret 'SLIMS_USERNAME'
    secret 'SLIMS_PASSWORD'

    input:
    tuple val(meta), path(summary)

    output:
    tuple val(meta), path("slims_upload.log"), emit: log

    script:
    """
    set -o pipefail
    slims_upload.py ${summary} \\
        --run-id ${meta.id} \\
        --ng-id ${meta.ng_id} \\
        2>&1 | tee slims_upload.log
    """
}
