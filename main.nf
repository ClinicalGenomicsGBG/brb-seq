#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ClinicalGenomicsGBG/brb-seq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/ClinicalGenomicsGBG/brb-seq
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl=2

include { BRB_SEQ  } from './workflows/brb-seq'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_brb-seq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_brb-seq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow CLINICALGENOMICSGBG_BRB_SEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    val_bclconvert_samplesheet
    val_fasta
    val_gtf
    val_rundir
    val_star_index



    main:

    ch_fasta = val_fasta
        ? channel.fromPath(params.fasta).map { file -> [ [ id: file.simpleName], file] }.collect()
        : channel.empty()
    ch_gtf = val_gtf
        ? channel.fromPath(params.gtf).map { file -> [ [ id: file.simpleName], file] }.collect()
        : channel.empty()
    ch_star_index = val_star_index
        ? channel.fromPath(params.star_index).map { file -> [ [ id: file.simpleName], file] }.collect()
        : channel.empty()
    ch_rundir = channel.fromPath(val_rundir)
    ch_bclconvert_samplesheet = channel.fromPath(val_bclconvert_samplesheet)

    BRB_SEQ (
        ch_samplesheet,
        ch_bclconvert_samplesheet,
        ch_fasta,
        ch_gtf,
        ch_rundir,
        ch_star_index,
        val_fasta ? val_fasta : false,
        val_gtf ? val_gtf.endsWith('.gz') : false
    )

    emit:
    multiqc_report = BRB_SEQ.out.multiqc_report // channel: /path/to/multiqc_report.html
    multiqc        = BRB_SEQ.out.multiqc        // channel: multiqc report + data + plots, for publishing
    starsolo       = BRB_SEQ.out.starsolo       // channel: STARsolo alignment + count outputs, for publishing
    umi_counts     = BRB_SEQ.out.umi_counts     // channel: per-sample UMI/read count matrices, for publishing
    fqtk           = BRB_SEQ.out.fqtk           // channel: FQTK demultiplexed FASTQs + metrics, for publishing
    star_index     = BRB_SEQ.out.star_index     // channel: generated STAR index, for publishing (only when --save_star_index is set)
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        workflow.outputDir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    CLINICALGENOMICSGBG_BRB_SEQ (
        PIPELINE_INITIALISATION.out.samplesheet,
        params.bclconvert_samplesheet,
        params.fasta,
        params.gtf,
        params.rundir,
        params.star_index
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        workflow.outputDir,
        params.monochrome_logs,
        params.hook_url,
        CLINICALGENOMICSGBG_BRB_SEQ.out.multiqc_report
    )

    publish:
    multiqc    = CLINICALGENOMICSGBG_BRB_SEQ.out.multiqc
    starsolo   = CLINICALGENOMICSGBG_BRB_SEQ.out.starsolo
    umi_counts = CLINICALGENOMICSGBG_BRB_SEQ.out.umi_counts
    fqtk       = CLINICALGENOMICSGBG_BRB_SEQ.out.fqtk
    star_index = CLINICALGENOMICSGBG_BRB_SEQ.out.star_index
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

output {
    multiqc {
        path 'multiqc'
    }
    starsolo {
        path 'starsolo'
    }
    umi_counts {
        path 'umi_counts'
    }
    fqtk {
        path { meta, file -> 
            file >> "fastq/${meta.id}/${file.name}"
        }
    }
    star_index {
        path 'star_index'
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
