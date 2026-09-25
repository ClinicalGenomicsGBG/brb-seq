/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BARCODESWHITELIST      } from '../modules/local/barcodeswhitelist/main'
include { BCLCONVERT             } from '../modules/nf-core/bclconvert/main'
include { CONVERTMATRIX          } from '../modules/local/convertmatrix/main'
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FQTK                   } from '../modules/nf-core/fqtk/main'
include { GUNZIP as GUNZIP_FASTA } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF   } from '../modules/nf-core/gunzip/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { STAGEFASTQDIR          } from '../modules/local/stagefastqdir/main'
include { STAR_GENOMEGENERATE    } from '../modules/nf-core/star/genomegenerate/main'
include { STARSOLO               } from '../modules/nf-core/star/starsolo/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_brb-seq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BRB_SEQ {

    take:
    ch_samplesheet            // channel: samplesheet read in from --input
    ch_bclconvert_samplesheet // channel: samplesheet read in from --bclconvert_samplesheet (optional)
    ch_rundir                 // channel: run directory path from --rundir
    ch_fasta                  // channel: FASTA file path from --fasta
    ch_gtf                    // channel: GTF file path from --gtf
    ch_star_index             // channel: pre-built STAR index directory from --star_index (optional)
    unzip_fasta               // boolean parameter: whether to unzip FASTA file (if gzipped) for STAR genome generation
    unzip_gtf                 // boolean parameter: whether to unzip GTF file (if gzipped) for STAR genome generation

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    
    ch_bclconvert_in = ch_bclconvert_samplesheet
        .combine(ch_rundir)
        .map { samplesheet, rundir -> tuple([id:rundir.name], samplesheet, rundir) }

    BCLCONVERT(ch_bclconvert_in)
    ch_multiqc_files = ch_multiqc_files.mix(BCLCONVERT.out.reports.map { _meta, file -> file })

    ch_demuxed = BCLCONVERT.out.fastq
        .map { _meta, fq ->
            def match = (fq =~ /^(.+)_S\d+/)
            def udi = match[0][1]
            tuple([id: udi], fq)
        }
        .branch { _meta, fq ->
            reads1: fq =~ /_R1_/
            reads2: fq =~ /_R2_/
        }

    ch_whitelist = ch_samplesheet
        .collectFile(newLine: true) {meta, barcode ->
            ["${meta.uid}.whitelist.txt", barcode]
        }
        .map { file ->
            def udi = file =~ (/^(.+)\.whitelist\.txt/)[0][1]
            tuple([id: udi], file)
        }

    ch_fqtk_samplesheet = ch_samplesheet
        .collectFile(newLine: true) {meta, barcode ->
            ["${meta.uid}.fqtk.txt", meta.id + "\t" + barcode]
        }
        .map { file ->
            def udi = file =~ (/^(.+)\.fqtk\.txt/)[0][1]
            tuple([id: udi], file)
        }
    
    ch_demuxed.reads1
        .join( ch_demuxed.reads2 )
        .join( ch_whitelist )
        .multiMap { meta, reads1, reads2, whitelist ->
            star_fq: [meta, "CB_UMI_Simple", [reads1, reads2]]
            star_barcodes: whitelist
        }
        .set { ch_star_input }

    ch_demuxed.reads1
        .concat( ch_demuxed.reads2 )
        .groupTuple()
        .set { ch_grouped_reads }
      

    FASTQC ( ch_grouped_reads )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{_meta, file -> file})

    //
    // Demultiplex the multiplexed FASTQs with FQTK for QC/archival purposes only.
    // STARsolo further below still consumes the original multiplexed FASTQs.
    //
    STAGEFASTQDIR ( ch_grouped_reads)

    ch_fqtk_samplesheet
        .join( STAGEFASTQDIR.out.dir )
        .join( ch_demuxed.reads1 )
        .join( ch_demuxed.reads2 )
        .map { meta, samplesheet, fastq_dir, reads1, reads2 ->
            def read_structure_pairs = [[reads1.name, '14B14T'], [reads2.name, '90T']]
            [meta, samplesheet, fastq_dir, read_structure_pairs]
        }
        .set { ch_fqtk_input }

    FQTK ( ch_fqtk_input )
    ch_multiqc_files = ch_multiqc_files.mix(FQTK.out.metrics.map { _meta, file -> file })

    ch_star_index_outputs = channel.empty()

    if (params.star_index) {
        ch_star_index_final = ch_star_index
    } else {
        if (unzip_fasta) {
            GUNZIP_FASTA ( ch_fasta )
            ch_fasta = GUNZIP_FASTA.out.gunzip.collect()
        }

        if (unzip_gtf) {
            GUNZIP_GTF ( ch_gtf )
            ch_gtf = GUNZIP_GTF.out.gunzip.collect()
        }

        STAR_GENOMEGENERATE (
            ch_fasta,
            ch_gtf
        )
        ch_star_index_final = STAR_GENOMEGENERATE.out.index.collect()

        // Only publish the index when the pipeline generated it itself
        // (a user-supplied --star_index is never re-published).
        if (params.save_star_index) {
            ch_star_index_outputs = STAR_GENOMEGENERATE.out.index
        }
    }

    STARSOLO (
        ch_star_input.star_fq,
        ch_star_input.star_barcodes,
        ch_star_index_final,
    )
    ch_multiqc_files = ch_multiqc_files.mix(STARSOLO.out.log_final.map { _meta, file -> file } )

    // All STARsolo outputs, published together under a single "starsolo" directory
    // (replicates the process-wide publishDir behavior previously configured).
    ch_starsolo_outputs = STARSOLO.out.counts
        .mix(STARSOLO.out.log_final)
        .mix(STARSOLO.out.log_out)
        .mix(STARSOLO.out.log_progress)
        .mix(STARSOLO.out.summary)

    CONVERTMATRIX (
        STARSOLO.out.counts.join( ch_fqtk_samplesheet )
    )

    // All FQTK outputs (demuxed FASTQs, metrics, unmatched reads), published together.
    ch_fqtk_outputs = FQTK.out.sample_fastq
        .mix(FQTK.out.metrics)
        .mix(FQTK.out.most_frequent_unmatched)
        .transpose()

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'brb-seq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    // All MultiQC outputs (report, data, plots), published together.
    ch_multiqc_outputs = MULTIQC.out.report
        .mix(MULTIQC.out.data)
        .mix(MULTIQC.out.plots)

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    multiqc        = ch_multiqc_outputs          // channel: multiqc report + data + plots, for publishing
    starsolo       = ch_starsolo_outputs         // channel: STARsolo alignment + count outputs, for publishing
    umi_counts     = CONVERTMATRIX.out.tsv       // channel: per-sample UMI/read count matrices, for publishing
    fqtk           = ch_fqtk_outputs             // channel: FQTK demultiplexed FASTQs + metrics, for publishing
    star_index     = ch_star_index_outputs       // channel: generated STAR index, for publishing (only when --save_star_index is set)
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
