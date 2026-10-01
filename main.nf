#!/usr/bin/env nextflow
/*
 * AtlasXomics epigenomic preprocessing (fastq2frags) for DBiT-seq ATAC-style data.
 *
 *   read1: genomic sequence
 *   read2: linker1 | barcodeA | linker2 | barcodeB | genomic sequence
 *
 * 1. Optional fastq cleaning against a correction table.
 * 2. Linker 1 / linker 2 filtering with bbduk (+ optional bulk barcode assignment).
 * 3. Chromap alignment -> fragment file (-1 suffix, out-of-bound and blacklist
 *    filtering, bgzip + tabix), optional cross-talk correction.
 * 4. SnapATAC2 fragment metrics, lane QC, and HTML report.
 * 5. Optional SLIMS upload of run metrics.
 */
nextflow.enable.dsl = 2

include { CLEAN_FASTQ       } from './modules/local/clean_fastq'
include { FILTER_LINKERS    } from './modules/local/filter_linkers'
include { CHROMAP_ALIGN     } from './modules/local/chromap_align'
include { PROCESS_FRAGMENTS } from './modules/local/process_fragments'
include { XTALK_CORRECT     } from './modules/local/xtalk_correct'
include { FRAGMENT_METRICS  } from './modules/local/fragment_metrics'
include { SLIMS_UPLOAD      } from './modules/local/slims_upload'


// Infer the genome ID from the reference FASTA name
def inferGenomeId(String name) {
    def known = ['GRCh38', 'GRCm38', 'GRCm39', 'Rnor6', 'mm10', 'mm39']
    def prefix = name.tokenize('_')[0]
    if (prefix in known) return prefix
    def hit = known.find { g ->
        name == g || name.startsWith("${g}.") || name.contains(".${g}.") ||
        name.contains("_${g}_") || name.contains("-${g}-")
    }
    return hit ?: prefix
}

// Spatial positions list for each barcode whitelist
def defaultPositionsFile(String barcodeName) {
    def base = 'https://latch-public.s3.amazonaws.com/test-data/13502'
    def positions = [
        'bc50.txt'        : 'x50_all_tissue_positions_list.csv',
        'bc50_old.txt'    : 'x50-old_tissue_positions_list.csv',
        'bc96.txt'        : 'x96_all_tissue_positions_list.csv',
        'bc96_fg.txt'     : 'xfg96_11DEC_all_tissue_positions_list.csv',
        'bc220-20-MAY.txt': 'xbc220-20-MAY_alltissue_positions_list.csv',
        'bc220-18-SEP.txt': 'xbc220-18-SEP_alltissue_positions_list.csv',
        'bc220_05-OCT.txt': 'xbc220_05-OCT_alltissue_positions_list.csv',
        'bc220_25-APR.txt': 'xbc220_25-APR_alltissue_positions_list.csv',
    ]
    return positions.containsKey(barcodeName) ? "${base}/${positions[barcodeName]}" : null
}

def missingColumns(tablePath, List required) {
    def header = []
    try {
        file(tablePath).withReader { r -> header = (r.readLine() ?: '').split(',')*.trim() }
    } catch (Exception _e) {
        return required
    }
    return required.findAll { !(it in header) }
}

def validateRunId(String id) {
    if (!(id ==~ /[^\/][^-.\s\/]+/)) {
        error "Invalid run_id '${id}': it cannot start with or contain a '/', dash (-), period (.) or space."
    }
}

// Decide per sample whether cleaning / cross-talk correction run: a toggle
// without a table is skipped with a warning, a table missing required columns
// is an error.
def resolveCorrections(String id, ds) {
    def clean = false
    def xtalk = false
    if (params.cleaning) {
        if (!ds) {
            log.warn "[${id}] Cleaning enabled but no ds_table supplied; skipping fastq cleaning."
        } else {
            def miss = missingColumns(ds, ['mer', 'finalDownFract'])
            if (miss) error "[${id}] ds_table missing cleaning columns: ${miss.join(', ')}"
            clean = true
        }
    }
    if (params.xtalk_correction) {
        if (!ds) {
            log.warn "[${id}] Cross-talk correction enabled but no ds_table supplied; skipping."
        } else {
            def miss = missingColumns(ds, ['mer', 'row', 'col', 'xtalk_mask', 'fold_above_background'])
            if (miss) error "[${id}] ds_table missing xtalk columns: ${miss.join(', ')}"
            xtalk = true
        }
    }
    return [clean: clean, xtalk: xtalk]
}


workflow {

    // ---------------------------------------------------------------- references
    if (!params.genome_dir) error "--genome_dir is required (directory with the Chromap *.fa and *.index)"
    def fastas  = file("${params.genome_dir}/*.fa")
    def indexes = file("${params.genome_dir}/*.index")
    if (!fastas)  error "No *.fa file found in genome_dir: ${params.genome_dir}"
    if (!indexes) error "No *.index file found in genome_dir: ${params.genome_dir}"
    def fasta = fastas.sort()[0]
    def index = indexes.sort()[0]
    def genome_id = params.genome_id ?: inferGenomeId(fasta.name)

    def barcode_file = params.barcode_file.contains('/')
        ? file(params.barcode_file, checkIfExists: true)
        : file("${projectDir}/barcodes/${params.barcode_file}", checkIfExists: true)

    def positions_src = params.positions_file ?: defaultPositionsFile(barcode_file.name)
    if (!positions_src) {
        error "No default positions file for barcode file '${barcode_file.name}'; provide --positions_file"
    }
    def positions_file = file(positions_src, checkIfExists: true)

    log.info """
    ATX fastq2frags v${workflow.manifest.version}
    ------------------------------------------------
    input        : ${params.input ?: "${params.r1}, ${params.r2}"}
    genome_dir   : ${params.genome_dir}
    genome_id    : ${genome_id}
    barcode_file : ${barcode_file}
    positions    : ${positions_file.toUriString()}
    outdir       : ${params.outdir}
    """.stripIndent()

    // ------------------------------------------------------------------ samples
    def rows
    if (params.input) {
        rows = Channel.fromPath(params.input, checkIfExists: true).splitCsv(header: true)
    } else {
        if (!(params.run_id && params.r1 && params.r2)) {
            error "Provide either --input <samplesheet.csv> or --run_id, --r1 and --r2"
        }
        rows = Channel.of([run_id: params.run_id, r1: params.r1, r2: params.r2])
    }

    ch_samples = rows.map { row ->
        def id = row.run_id?.trim()
        if (!id) error "Samplesheet row missing run_id: ${row}"
        validateRunId(id)
        def ds = row.ds_table?.trim() ?: params.ds_table
        def fix = resolveCorrections(id, ds)
        def meta = [
            id       : id,
            ng_id    : row.ng_id?.trim() ?: params.ng_id,
            genome_id: genome_id,
            clean    : fix.clean,
            xtalk    : fix.xtalk,
        ]
        [meta, file(row.r1, checkIfExists: true), file(row.r2, checkIfExists: true),
         ds ? file(ds, checkIfExists: true) : []]
    }

    ch_ds_table = ch_samples.map { meta, _r1, _r2, ds -> [meta.id, ds] }

    // ---------------------------------------------------------------- filtering
    ch_to_clean = ch_samples.branch { meta, _r1, _r2, _ds ->
        clean: meta.clean
        skip : true
    }

    CLEAN_FASTQ(ch_to_clean.clean)

    ch_reads = CLEAN_FASTQ.out.reads.mix(
        ch_to_clean.skip.map { meta, r1, r2, _ds -> [meta, r1, r2] }
    )

    FILTER_LINKERS(ch_reads, barcode_file)

    // ---------------------------------------------------------------- alignment
    CHROMAP_ALIGN(FILTER_LINKERS.out.reads, fasta, index, barcode_file)

    PROCESS_FRAGMENTS(CHROMAP_ALIGN.out.bed)

    ch_frags = PROCESS_FRAGMENTS.out.fragments.branch { meta, _frag, _tbi ->
        xtalk: meta.xtalk
        done : true
    }

    XTALK_CORRECT(
        ch_frags.xtalk
            .map { meta, frag, tbi -> [meta.id, meta, frag, tbi] }
            .join(ch_ds_table)
            .map { _id, meta, frag, tbi, ds -> [meta, frag, tbi, ds] }
    )

    ch_final_frags = ch_frags.done.mix(XTALK_CORRECT.out.fragments)

    // ------------------------------------------------------------------ metrics
    FRAGMENT_METRICS(
        ch_final_frags
            .map { meta, frag, tbi -> [meta.id, meta, frag, tbi] }
            .join(CHROMAP_ALIGN.out.log.map { meta, log -> [meta.id, log] })
            .map { _id, meta, frag, tbi, log -> [meta, frag, tbi, log] },
        positions_file
    )

    if (params.upload_to_slims) {
        SLIMS_UPLOAD(FRAGMENT_METRICS.out.summary)
    }
}
