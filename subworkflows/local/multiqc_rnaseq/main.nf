//
// MultiQC report assembly for nf-core/rnaseq.
//

include { MULTIQC                 } from '../../../modules/nf-core/multiqc'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc    } from '../../nf-core/utils_nfcore_pipeline'
include { workflowVersionToYAML   } from '../../nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../utils_nfcore_rnaseq_pipeline'
include { multiqcNameReplacements } from '../utils_nfcore_rnaseq_pipeline'
include { multiqcSampleMergeYaml  } from '../utils_nfcore_rnaseq_pipeline'

workflow MULTIQC_RNASEQ {

    take:
    ch_multiqc_files           // channel: [ val(meta), path(file) ]       - flat, contributor outputs
    ch_per_sample_bundle_raw   // channel: [ id, meta, f1, f2, ... ]       - per-sample, grown by `.join(..., remainder: true)` at each subworkflow aggregation site
    ch_trim_read_count         // channel: [ val(meta), val(num_reads) ]   - for fail_trimmed section
    ch_percent_mapped_pass     // channel: [ id, percent_mapped, pass ]    - for fail_mapped section
    ch_strand_comparison       // channel: [ val(meta), status, lines ]    - for fail_strand section
    ch_fastq                   // channel: [ val(meta), [ reads ] ]
    ch_collated_versions       // channel: path(versions yaml)
    samplesheet_path           // path: pipeline input samplesheet
    samplesheet_schema         // path: samplesheet JSON schema
    mqc_default_config         // path: pipeline-bundled MultiQC config
    mqc_custom_config          // path (or []): optional user MultiQC config
    mqc_logo                   // path (or []): optional custom logo
    methods_description_yml    // path: methods-description YAML template
    sample_status_header       // path: MultiQC custom content header for fail_* tables
    min_trimmed_reads          // integer: threshold for fail_trimmed classification
    skip_quantification_merge  // boolean

    main:

    // Per-sample fail_* TSVs and their merged-mode aggregates. The anchor
    // is derived from the bundle itself so every bundle sample has a match
    // on every fail_* stream, keeping the downstream `.join(..., remainder: true)`
    // progressive (unmatched samples would otherwise wait for channel close).
    // `skip:` tracks the header length so editing `sample_status_header.txt`
    // doesn't silently mis-skip the merged aggregate's concatenation.
    def status_header_lines = sample_status_header.readLines().size() + 1  // parent header + one column row
    ch_sample_anchor_by_id = ch_per_sample_bundle_raw.map { row -> [row[0], row[1]] }

    ch_fail_trimmed_fail_by_id = ch_trim_read_count
        .filter { _meta, n -> n <= min_trimmed_reads.toFloat() }
        .collectFile { meta, n ->
            ["${meta.id}_fail_trimmed_samples_mqc.tsv",
             "Sample\tReads after trimming\n${meta.id}\t${n}\n"]
        }
        .map { f -> [f.baseName.replace('_fail_trimmed_samples_mqc', ''), f] }
    ch_fail_trimmed_all = ch_sample_anchor_by_id
        .join(ch_fail_trimmed_fail_by_id, remainder: true)
        .map { _id, meta, f -> [meta, f ?: []] }
    ch_fail_trimmed_merged = ch_fail_trimmed_all
        .map { _meta, f -> f }
        .flatten()
        .collectFile(name: 'fail_trimmed_samples_mqc.tsv', keepHeader: true)
        .map { f -> [[:], f] }

    ch_fail_mapped_fail_by_id = ch_percent_mapped_pass
        .filter { _id, _pm, pass -> pass != null && !pass }
        .collectFile { id, percent_mapped, _pass ->
            ["${id}_fail_mapped_samples_mqc.tsv",
             sample_status_header.text +
             "Sample\tSTAR uniquely mapped reads (%)\n${id}\t${percent_mapped}\n"]
        }
        .map { f -> [f.baseName.replace('_fail_mapped_samples_mqc', ''), f] }
    ch_fail_mapped_all = ch_sample_anchor_by_id
        .join(ch_fail_mapped_fail_by_id, remainder: true)
        .map { _id, meta, f -> [meta, f ?: []] }
    ch_fail_mapped_merged = ch_fail_mapped_all
        .map { _meta, f -> f }
        .flatten()
        .collectFile(name: 'fail_mapped_samples_mqc.tsv', keepHeader: true, skip: status_header_lines)
        .map { f -> [[:], f] }

    // fail_strand shows status for every sample — pass or fail — that ran
    // through RSeQC infer_experiment. Samples without a comparison still
    // emit an empty placeholder so the bundle join stays progressive.
    def strand_header = "Sample\tStrand inference method\tStatus\tProvided strandedness\tInferred strandedness\tSense (%)\tAntisense (%)\tUnstranded (%)\n"
    ch_fail_strand_by_id = ch_strand_comparison
        .collectFile { meta, _s, lines ->
            ["${meta.id}_fail_strand_check_mqc.tsv",
             sample_status_header.text + strand_header + lines.join('\n') + '\n']
        }
        .map { f -> [f.baseName.replace('_fail_strand_check_mqc', ''), f] }
    ch_fail_strand_all = ch_sample_anchor_by_id
        .join(ch_fail_strand_by_id, remainder: true)
        .map { _id, meta, f -> [meta, f ?: []] }
    ch_fail_strand_merged = ch_fail_strand_all
        .map { _meta, f -> f }
        .flatten()
        .collectFile(name: 'fail_strand_check_mqc.tsv', keepHeader: true, skip: status_header_lines)
        .map { f -> [[:], f] }

    // Extend the raw bundle with fail_* rows, then collapse to `[meta, files_list]`.
    ch_per_sample_bundle = ch_per_sample_bundle_raw
        .join(ch_fail_trimmed_all.map { meta, f -> [meta.id, f] }, remainder: true)
        .join(ch_fail_mapped_all.map { meta, f -> [meta.id, f] },  remainder: true)
        .join(ch_fail_strand_all.map { meta, f -> [meta.id, f] },  remainder: true)
        .map { row ->
            [row[1], row.drop(2).findAll { it != null }.collectMany { entry -> (entry instanceof List) ? entry : [entry] }]
        }

    // Per-run table_sample_merge config: only PE samples from the samplesheet
    // get their _1/_2 rows grouped in the General Stats table.
    ch_mqc_dynamic_config = channel.of(multiqcSampleMergeYaml(samplesheet_path, samplesheet_schema))
        .collectFile(name: 'multiqc_sample_merge.yml')

    // Workflow summary and methods description rendered as MultiQC sections.
    ch_workflow_summary = channel.value(
        paramsSummaryMultiqc(
            paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        )
    ).collectFile(name: 'workflow_summary_mqc.yaml')

    ch_methods_description = channel.value(
        methodsDescriptionText(methods_description_yml)
    ).collectFile(name: 'methods_description_mqc.yaml')

    // --replace-names TSV so MultiQC uses sample IDs rather than FASTQ basenames.
    ch_name_replacements = multiqcNameReplacements(ch_fastq)

    // Merged-mode file set: contributor outputs + run-level context + fail_* aggregates.
    ch_multiqc_files_merged = ch_multiqc_files
        .mix(ch_fail_trimmed_merged)
        .mix(ch_fail_mapped_merged)
        .mix(ch_fail_strand_merged)
        .mix(ch_workflow_summary.mix(ch_collated_versions).mix(ch_methods_description).map { f -> [[:], f] })

    if (skip_quantification_merge) {
        // One report per sample. Per-sample reports carry the pipeline-
        // identity YAML (name + version + commit SHA + Nextflow version)
        // instead of the full collated tool versions — the latter closes
        // only once the workflow-global `versions` topic closes, which
        // would reintroduce a full-run barrier. The full YAML is still
        // published unchanged to `pipeline_info/`.
        ch_manifest_versions = channel.value(workflowVersionToYAML())
            .collectFile(name: 'nf_core_rnaseq_software_mqc_versions.yml')

        ch_static_globals = ch_workflow_summary
            .mix(ch_methods_description)
            .mix(ch_manifest_versions)
            .collect()

        ch_global_files = ch_fail_trimmed_merged
            .mix(ch_fail_mapped_merged)
            .mix(ch_fail_strand_merged)
            .map { _meta, f -> f }
            .collect()
            .ifEmpty([])

        ch_multiqc_input = ch_per_sample_bundle
            .combine(ch_static_globals.toList())
            .combine(ch_global_files.toList())
            .combine(ch_mqc_dynamic_config)
            .map { meta, sample_files, static_globals, run_globals, dyn ->
                [
                    [id: meta.id],
                    sample_files + (static_globals ?: []) + (run_globals ?: []),
                    [mqc_default_config, dyn, mqc_custom_config].findAll { it },
                    mqc_logo,
                    [],  // no replace_names — each report contains one sample
                    [],
                ]
            }
    } else {
        // Merged mode: `multiqc_report` is a sentinel meta.id used by
        // conf/modules/multiqc.config to pick the merged output path.
        // Wrap files in a 1-tuple so `.combine()` keeps them grouped.
        ch_multiqc_input = ch_multiqc_files_merged
            .map { _meta, f -> f }
            .collect()
            .map { files -> [files] }
            .combine(ch_name_replacements.ifEmpty([]).toList())
            .combine(ch_mqc_dynamic_config)
            .map { files, replace_names, dyn ->
                [
                    [id: 'multiqc_report'],
                    files,
                    [mqc_default_config, dyn, mqc_custom_config].findAll { it },
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
