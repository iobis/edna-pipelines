#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.worms_db = params.worms_db ?: null
params.metadata = params.metadata ?: null
params.outdir = params.outdir ?: 'results'
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
    // Fixed layout under --outdir (not user-overridable).
    def ampliseq_outdir_abs = "${root_outdir_abs}/ampliseq"
    def darwincore_outdir_abs = "${root_outdir_abs}/darwincore"

    if (!params.worms_db) error('Set --worms_db to your WoRMS SQLite database (aphiasync)')
    if (!params.metadata) error('Set --metadata to your sample metadata TSV')

    file(params.worms_db, checkIfExists: true)
    file(params.metadata, checkIfExists: true)

    // Skip nested ampliseq when --ampliseq_results points at an existing outdir.
    if (params.ampliseq_results) {
        ampliseq_done_ch = channel.value('existing_ampliseq')
        ampliseq_results_path = params.ampliseq_results
        file(ampliseq_results_path, checkIfExists: true)
    } else {
        if (!params.input) error('Set --input samplesheet, or --ampliseq_results to reuse a prior run')
        RUN_AMPLISEQ(channel.value(ampliseq_outdir_abs))
        ampliseq_done_ch = RUN_AMPLISEQ.out.done
        ampliseq_results_path = ampliseq_outdir_abs
    }

    params.ampliseq_results = ampliseq_results_path

    def darwincore_outdir_ch = channel.value(darwincore_outdir_abs)

    // Report image, build from containers/Dockerfile when missing or Dockerfile has changed.
    def report_image_ch
    if (workflow.profile.tokenize(',').contains('docker')) {
        ENSURE_REPORT_CONTAINER()
        report_image_ch = ENSURE_REPORT_CONTAINER.out.image.map { it.text.trim() }
    } else {
        report_image_ch = channel.value('')
    }

    WORMS_MATCH(channel.of('sintax', 'vsearch'), ampliseq_done_ch, darwincore_outdir_ch)
    BUILD_DARWIN_CORE(WORMS_MATCH.out.worms.collect(), darwincore_outdir_ch)
    WORMS_MATCH_OCCURRENCE(BUILD_DARWIN_CORE.out.publishing, darwincore_outdir_ch)
    SUMMARY_REPORT(
        WORMS_MATCH.out.worms.filter { method, _path -> method == 'sintax' }.map { _method, path -> path },
        WORMS_MATCH.out.worms.filter { method, _path -> method == 'vsearch' }.map { _method, path -> path },
        darwincore_outdir_ch,
        report_image_ch
    )
}

process RUN_AMPLISEQ {
    tag 'ampliseq'
    publishDir { ampliseq_outdir_abs }, mode: 'copy', pattern: 'ampliseq.{done,revision}'

    input:
    val ampliseq_outdir_abs

    output:
    path 'ampliseq.done', emit: done
    path 'ampliseq.revision', emit: revision

    script:
    def ampliseq_profile = params.ampliseq_profile ?: 'standard'
    def ampliseq_config_flag = ampliseq_profile == 'docker'
        ? " -c ${projectDir}/conf/ampliseq_docker.config"
        : ''
    def ampliseq_dir = "${projectDir}/third_party/nf-core-ampliseq"
    def ampliseq_update = params.ampliseq_update ? 'true' : 'false'
    def ampliseq_sequencing_type = params.sequencing_type ?: 'illumina_pe'
    """
    set -euo pipefail
    {
      echo "ampliseq_repo=${params.ampliseq_repo}"
      echo "ampliseq_ref=${params.ampliseq_ref}"
      bash "${projectDir}/scripts/ensure_ampliseq.sh" \\
        --repo "${params.ampliseq_repo}" \\
        --ref "${params.ampliseq_ref}" \\
        --update ${ampliseq_update} \\
        --target-dir "${ampliseq_dir}"
    } | tee ampliseq.revision
    echo "Running nf-core/ampliseq from ${ampliseq_dir} (${params.ampliseq_repo} @ ${params.ampliseq_ref})" >&2
    nextflow run "${ampliseq_dir}" \\
      -profile ${ampliseq_profile}${ampliseq_config_flag} \\
      --input ${params.input} \\
      --outdir ${ampliseq_outdir_abs} \\
      --sequencing_type ${ampliseq_sequencing_type} \\
      ${params.binned_quality ? "--binned_quality ${params.binned_quality}" : ''} \\
      ${params.ignore_binned_quality ? '--ignore_binned_quality' : ''} \\
      ${params.skip_cutadapt ? '--skip_cutadapt' : ''} \\
      ${params.skip_dada_taxonomy ? '--skip_dada_taxonomy' : ''} \\
      ${params.skip_phyloseq ? '--skip_phyloseq' : ''} \\
      ${params.skip_tse ? '--skip_tse' : ''} \\
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
    publishDir { "${darwincore_outdir}/worms/${method}" }, mode: 'copy', overwrite: true

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

process ENSURE_REPORT_CONTAINER {
    tag 'report-image'
    // Always check, rebuild only when the hash-tagged image is missing.
    cache false

    output:
    path 'report_image.txt', emit: image

    script:
    """
    set -euo pipefail
    "${projectDir}/scripts/ensure_report_container.sh" \\
      "${projectDir}/containers/Dockerfile" > report_image.txt
    """
}

process SUMMARY_REPORT {
    tag 'summary'
    publishDir { "${darwincore_outdir}/report" }, mode: 'copy', overwrite: true
    container { report_image ? report_image.toString() : null }
    containerOptions '-u $(id -u):$(id -g)'

    input:
    path sintax_tsv
    path vsearch_tsv
    val darwincore_outdir
    val report_image

    output:
    path 'summary_report.html', emit: report
    path 'report_params.tsv', emit: params

    script:
    """
    set -euo pipefail
    cat > report_params.tsv <<'EOF'
parameter	value
pipeline_version	${workflow.manifest.version}
run_date	${workflow.start.toString().take(19).replace('T', ' ')}
outdir	${params.outdir}
ampliseq_results	${params.ampliseq_results}
darwincore_outdir	${darwincore_outdir}
metadata	${params.metadata}
worms_db	${params.worms_db}
run_ampliseq	${params.ampliseq_results ? false : true}
ampliseq_repo	${params.ampliseq_repo ?: ''}
ampliseq_ref	${params.ampliseq_ref ?: ''}
ampliseq_update	${params.ampliseq_update ?: false}
sequencing_type	${params.sequencing_type ?: ''}
skip_cutadapt	${params.skip_cutadapt ?: false}
skip_dada_taxonomy	${params.skip_dada_taxonomy ?: false}
skip_phyloseq	${params.skip_phyloseq ?: false}
skip_tse	${params.skip_tse ?: false}
sintax_ref_tax_custom	${params.sintax_ref_tax_custom ?: ''}
sintax_assign_taxlevels	${params.sintax_assign_taxlevels ?: ''}
vsearch_lca_ref_tax_custom	${params.vsearch_lca_ref_tax_custom ?: ''}
vsearch_lca_assign_taxlevels	${params.vsearch_lca_assign_taxlevels ?: ''}
vsearch_lca_id	${params.vsearch_lca_id != null ? params.vsearch_lca_id : ''}
vsearch_lca_maxaccepts	${params.vsearch_lca_maxaccepts != null ? params.vsearch_lca_maxaccepts : ''}
vsearch_lca_maxrejects	${params.vsearch_lca_maxrejects != null ? params.vsearch_lca_maxrejects : ''}
vsearch_lca_lca_cutoff	${params.vsearch_lca_lca_cutoff != null ? params.vsearch_lca_lca_cutoff : ''}
clean_prefix	${params.clean_prefix != null ? params.clean_prefix : false}
EOF
    cp "${projectDir}/bin/summary_report.Rmd" summary_report.Rmd
    cp "${projectDir}/bin/summary_report.R" summary_report.R
    cp "${projectDir}/bin/render_summary_report.R" render_summary_report.R
    Rscript render_summary_report.R \\
      --sintax "\${PWD}/${sintax_tsv}" \\
      --vsearch "\${PWD}/${vsearch_tsv}" \\
      --params "\${PWD}/report_params.tsv" \\
      --output "\${PWD}/summary_report.html" \\
      --rmd "\${PWD}/summary_report.Rmd"
    """
}
