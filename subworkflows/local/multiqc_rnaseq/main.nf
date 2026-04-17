//
// MultiQC report assembly for nf-core/rnaseq.
//

include { MULTIQC                 } from '../../../modules/nf-core/multiqc'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc    } from '../../nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../utils_nfcore_rnaseq_pipeline'
include { multiqcNameReplacements } from '../utils_nfcore_rnaseq_pipeline'

workflow MULTIQC_RNASEQ {

    take:
    ch_multiqc_files           // channel: [ val(meta), path(file) ]
    ch_fastq                   // channel: [ val(meta), [ reads ] ]
    ch_collated_versions       // channel: path(versions yaml)
    mqc_default_config         // path: pipeline-bundled MultiQC config
    mqc_custom_config          // path (or []): optional user MultiQC config
    mqc_logo                   // path (or []): optional custom logo
    methods_description_yml    // path: methods-description YAML template
    skip_quantification_merge  // boolean

    main:

    // Workflow summary and methods description rendered as MultiQC sections.
    ch_workflow_summary = channel.value(
        paramsSummaryMultiqc(
            paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        )
    ).collectFile(name: 'workflow_summary_mqc.yaml')

    ch_methods_description = channel.value(
        methodsDescriptionText(methods_description_yml)
    ).collectFile(name: 'methods_description_mqc.yaml')

    // Everything MultiQC will ingest: the collected per-sample QC outputs
    // plus the global context files, tagged with an empty meta.
    ch_multiqc_all = ch_multiqc_files.mix(
        ch_workflow_summary
            .mix(ch_collated_versions)
            .mix(ch_methods_description)
            .map { f -> [[:], f] }
    )

    // --replace-names TSV so MultiQC uses sample IDs rather than FASTQ basenames.
    ch_name_replacements = multiqcNameReplacements(ch_fastq)

    def mqc_configs = [mqc_default_config, mqc_custom_config].findAll { it }

    if (skip_quantification_merge) {
        // One MultiQC report per sample. Split incoming files into per-sample
        // and global buckets, then attach the global bucket to every sample.
        ch_branched = ch_multiqc_all
            .branch { meta, _file ->
                per_sample: meta.id != null
                global: true
            }

        ch_global_files = ch_branched.global
            .map { _meta, f -> f }
            .collect()

        ch_multiqc_input = ch_branched.per_sample
            .map { meta, f -> [meta.id, f] }
            .groupTuple()
            .combine(ch_global_files.toList())
            .map { id, sample_files, global_files ->
                [
                    [id: id],
                    sample_files + (global_files ?: []),
                    mqc_configs,
                    mqc_logo,
                    [],  // no replace_names — each report contains one sample's files
                    [],
                ]
            }
    } else {
        // One merged MultiQC report. 'multiqc_report' is a sentinel meta.id
        // used by conf/modules/multiqc.config to pick the merged output
        // path/prefix. Wrap the collected file list in a 1-tuple so
        // .combine() doesn't spread it across the downstream closure args.
        ch_all_files = ch_multiqc_all
            .map { _meta, f -> f }
            .collect()
            .map { files -> [files] }

        ch_multiqc_input = ch_all_files
            .combine(ch_name_replacements.ifEmpty([]).toList())
            .map { files, replace_names ->
                [
                    [id: 'multiqc_report'],
                    files,
                    mqc_configs,
                    mqc_logo,
                    replace_names ?: [],
                    [],
                ]
            }
    }

    MULTIQC(ch_multiqc_input)

    emit:
    report = MULTIQC.out.report.map { _meta, report -> report }
}
