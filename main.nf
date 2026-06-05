#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.worms_db = params.worms_db ?: null
params.metadata = params.metadata ?: null
params.run_ampliseq = params.run_ampliseq == null ? true : params.run_ampliseq
params.outdir = params.outdir ?: 'results'
params.ampliseq_outdir = params.ampliseq_outdir ?: "${params.outdir}/ampliseq"
params.darwincore_outdir = params.darwincore_outdir ?: "${params.outdir}/darwincore"
params.ampliseq_profile = params.ampliseq_profile ?: 'standard'
params.clean_prefix = params.clean_prefix ?: false

workflow {
    def ampliseq_done_ch
    def ampliseq_results_path
    def launch_dir = new File(workflow.launchDir.toString())
    def root_outdir = params.outdir ?: 'results'
    def root_outdir_abs = new File(root_outdir.toString()).isAbsolute()
        ? new File(root_outdir.toString()).absolutePath
        : new File(launch_dir, root_outdir.toString()).absolutePath
    def ampliseq_outdir_abs = params.ampliseq_outdir
        ? (new File(params.ampliseq_outdir.toString()).isAbsolute()
            ? new File(params.ampliseq_outdir.toString()).absolutePath
            : new File(launch_dir, params.ampliseq_outdir.toString()).absolutePath)
        : "${root_outdir_abs}/ampliseq"
    def darwincore_outdir_abs = params.darwincore_outdir
        ? (new File(params.darwincore_outdir.toString()).isAbsolute()
            ? new File(params.darwincore_outdir.toString()).absolutePath
            : new File(launch_dir, params.darwincore_outdir.toString()).absolutePath)
        : "${root_outdir_abs}/darwincore"

    params.ampliseq_outdir = ampliseq_outdir_abs
    params.darwincore_outdir = darwincore_outdir_abs

    if (!params.worms_db) error('Set --worms_db to your WoRMS SQLite database (aphiasync)')
    if (!params.metadata) error('Set --metadata to your sample metadata TSV')

    file(params.worms_db, checkIfExists: true)
    file(params.metadata, checkIfExists: true)

    if (params.run_ampliseq) {
        if (!params.input) error('Set --input samplesheet when --run_ampliseq true')
        RUN_AMPLISEQ(channel.value(ampliseq_outdir_abs))
        ampliseq_done_ch = RUN_AMPLISEQ.out
        ampliseq_results_path = ampliseq_outdir_abs
    } else if (!params.ampliseq_results) {
        error('Set --ampliseq_results when --run_ampliseq false')
    } else {
        ampliseq_done_ch = channel.value('existing_ampliseq')
        ampliseq_results_path = params.ampliseq_results
        file(ampliseq_results_path, checkIfExists: true)
    }

    params.ampliseq_results = ampliseq_results_path

    def darwincore_outdir_ch = channel.value(darwincore_outdir_abs)
    WORMS_MATCH(channel.of('sintax', 'vsearch'), ampliseq_done_ch, darwincore_outdir_ch)
    BUILD_DARWIN_CORE(WORMS_MATCH.out.worms.collect(), darwincore_outdir_ch)
    WORMS_MATCH_OCCURRENCE(BUILD_DARWIN_CORE.out.publishing, darwincore_outdir_ch)
}

process RUN_AMPLISEQ {
    tag 'ampliseq'
    publishDir { ampliseq_outdir_abs }, mode: 'copy', pattern: 'ampliseq.done'

    input:
    val ampliseq_outdir_abs

    output:
    path 'ampliseq.done'

    script:
    def ampliseq_profile = params.ampliseq_profile ?: 'standard'
    def ampliseq_config_flag = ampliseq_profile == 'docker'
        ? " -c ${projectDir}/conf/ampliseq_docker.config"
        : ''
    """
    set -euo pipefail
    nextflow run "${projectDir}/third_party/nf-core-ampliseq" \\
      -profile ${ampliseq_profile}${ampliseq_config_flag} \\
      --input ${params.input} \\
      --outdir ${ampliseq_outdir_abs} \\
      ${params.single_end ? '--single_end' : ''} \\
      ${params.skip_cutadapt ? '--skip_cutadapt' : ''} \\
      ${params.skip_dada_taxonomy ? '--skip_dada_taxonomy' : ''} \\
      ${params.sintax_ref_tax_custom ? "--sintax_ref_tax_custom ${params.sintax_ref_tax_custom}" : ''} \\
      ${params.sintax_assign_taxlevels ? "--sintax_assign_taxlevels ${params.sintax_assign_taxlevels}" : ''} \\
      ${params.vsearch_lca_ref_tax_custom ? "--vsearch_lca_ref_tax_custom ${params.vsearch_lca_ref_tax_custom}" : ''} \\
      ${params.vsearch_lca_assign_taxlevels ? "--vsearch_lca_assign_taxlevels ${params.vsearch_lca_assign_taxlevels}" : ''} \\
      ${params.vsearch_lca_id != null ? "--vsearch_lca_id ${params.vsearch_lca_id}" : ''} \\
      ${params.vsearch_lca_maxaccepts != null ? "--vsearch_lca_maxaccepts ${params.vsearch_lca_maxaccepts}" : ''} \\
      ${params.vsearch_lca_maxrejects != null ? "--vsearch_lca_maxrejects ${params.vsearch_lca_maxrejects}" : ''} \\
      ${params.vsearch_lca_lca_cutoff != null ? "--vsearch_lca_lca_cutoff ${params.vsearch_lca_lca_cutoff}" : ''}
    touch ampliseq.done
    """
}

process WORMS_MATCH {
    tag "$method"
    publishDir { "${darwincore_outdir}/worms/${method}" }, mode: 'copy'

    input:
    val method
    val _amp_done
    val darwincore_outdir

    output:
    tuple val(method), path("worms_matched.${method}.tsv"), emit: worms

    script:
    def out = "worms_matched.${method}.tsv"
    """
    set -euo pipefail
    python3 "${projectDir}/bin/worms_match.py" \\
      --method ${method} \\
      --input ${params.ampliseq_results} \\
      --worms-db ${params.worms_db} \\
      --output ${out}
    """
}

process BUILD_DARWIN_CORE {
    tag 'dwc'
    publishDir { darwincore_outdir }, mode: 'copy', overwrite: true

    input:
    val _worms_done
    val darwincore_outdir

    output:
    path 'publishing', emit: publishing

    script:
    """
    set -euo pipefail
    python3 "${projectDir}/bin/build_darwin_core.py" \\
      --ampliseq-results ${params.ampliseq_results} \\
      --metadata ${params.metadata} \\
      --output ${darwincore_outdir} \\
      ${params.clean_prefix ? '--clean-prefix' : ''}
    mkdir -p publishing
    cp ${darwincore_outdir}/publishing/*.tsv publishing/
    """
}

process WORMS_MATCH_OCCURRENCE {
    tag 'occurrence'
    publishDir { darwincore_outdir }, mode: 'copy', overwrite: true

    input:
    path publishing_in
    val darwincore_outdir

    output:
    path 'publishing', emit: publishing

    script:
    """
    set -euo pipefail
    python3 "${projectDir}/bin/worms_match.py" \\
      --input ${publishing_in}/occurrence.tsv \\
      --output ${publishing_in}/occurrence.tsv \\
      --worms-db ${params.worms_db}
    if [ "${publishing_in}" != "publishing" ]; then
      ln -sfn ${publishing_in} publishing
    fi
    """
}
