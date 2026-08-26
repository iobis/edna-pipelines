# edna-pipelines

This is an eDNA metabarcoding pipeline with WoRMS aligned Darwin Core output building on nf-core/ampliseq.

## Parameters

Default parameters values are set in `nextflow.config`.

| Parameter | Description | Forwarded to ampliseq | Default |
| --- | --- | --- | --- |
| `--input` | Sample sheet (required unless `--ampliseq_results`) | yes | |
| `--metadata` | Metadata sheet | no | |
| `--worms_db` | WoRMS SQLite for taxon matching | no | |
| `--outdir` | Root output path | ampliseq subfolder as `--outdir` | `results` |
| `--ampliseq_results` | Reuse existing ampliseq outdir (skips nested ampliseq) | no | |
| `--clean_prefix` | Strip `s_` sample prefix in DwC builder | no | |
| `--ampliseq_results` | Reuse an existing ampliseq outdir (skips nested ampliseq) | no | |
| `--ampliseq_repo` | Git URL for ampliseq checkout | no | `https://github.com/nf-core/ampliseq.git` |
| `--ampliseq_ref` | Git ref / branch for ampliseq | no | `dev` |
| `--ampliseq_update` | Re-fetch ampliseq on each run | no | `false` |
| `--ampliseq_profile` | Nextflow profile for nested ampliseq | As `-profile` | `standard`, or `docker` when `-profile docker` is set |
| `--sequencing_type` | Sequencing type | yes | `illumina_pe` |
| `--primer_fwd` | Forward primer (required unless `--skip_cutadapt`) | yes | |
| `--primer_rev` | Reverse primer (required unless `--skip_cutadapt`) | yes | |
| `--binned_quality` | Comma separated quality bins | yes | |
| `--ignore_binned_quality` | Ignore binned quality warnings | yes | `false` |
| `--skip_cutadapt` | | yes | `false` |
| `--skip_dada_taxonomy` | | yes | `false` |
| `--skip_phyloseq` | | yes | `false` |
| `--skip_tse` | | yes | `false` |
| `--sintax_ref_tax_custom` | | yes | |
| `--sintax_assign_taxlevels` | | yes | |
| `--vsearch_lca_ref_tax_custom` | | yes | |
| `--vsearch_lca_assign_taxlevels` | | yes | |
| `--vsearch_lca_id` | | yes | `0.9` |
| `--vsearch_lca_maxaccepts` | | yes | `0` |
| `--vsearch_lca_maxrejects` | | yes | `0` |
| `--vsearch_lca_lca_cutoff` | | yes | `0.9` |

The following key ampliseq parameters are currently not configurable:

| Parameter | Default |
| --- | --- |
| `--truncq` | `2` |
| `--max_ee` | `2` |
| `--min_len` | `50` |
| `--sample_inference` | `pooled` |
| `--trunclenf` / `--trunclenr` | |
| `--cutadapt_min_overlap` | `3` |
| `--cutadapt_max_error_rate` | `0.1` |
| `--vsearch_lca_query_cov` | `1.0` |

In addition, `--sintax_cutoff` (default `0.8`) is [currently not configurable](https://github.com/nf-core/ampliseq/blob/master/conf/modules.config#L581) in ampliseq.

## Steps
### Read QC

| Nextflow process | Description |
| --- | --- |
| `FASTQC` | Per-sample quality plots for raw reads |
| `MULTIQC` | |

### Primer removal

The subworkflow `CUTADAPT_WORKFLOW` can be disabled with `--skip_cutadapt`.

| Nextflow process | Description |
| --- | --- |
| `CUTADAPT_BASIC` | Trim primer sequences from reads |
| `CUTADAPT_SUMMARY_STD` |  |
| `CUTADAPT_SUMMARY_MERGE` |  |
| `CUTADAPT_READTHROUGH` | optional |
| `CUTADAPT_DOUBLEPRIMER` | optional |
| `CUTADAPT_SUMMARY_DOUBLEPRIMER` | optional |

### Quality trimming

| Nextflow process | Description |
| --- | --- |
| `DADA2_QUALITY1` | Quality profiles of primer-trimmed reads |
| `TRUNCLEN` | Auto-pick truncation lengths from quality plots |
| `DADA2_FILTNTRIM` | Filter and truncate reads by quality and length |
| `DADA2_QUALITY2` | Quality profiles after filtering/truncation |

### ASV inference

| Nextflow process | Description |
| --- | --- |
| `DADA2_ERR` | Learn DADA2 sequencing error models |
| `DADA2_DENOISING` | Infer ASVs from filtered reads |
| `DADA2_RMCHIMERA` | Remove chimeric ASVs |
| `DADA2_STATS` | Track per-sample read counts through DADA2 |
| `DADA2_MERGE` | Combine ASV tables across samples |

### Taxonomic classification
#### SINTAX

| Nextflow process | Description |
| --- | --- |
| `FORMAT_TAXONOMY_SINTAX` |  |
| `VSEARCH_SINTAX` | Assign taxonomy with SINTAX |
| `FORMAT_TAXRESULTS_SINTAX` | Turn SINTAX output into a taxonomy table |

#### VSEARCH + LCA

| Nextflow process | Description |
| --- | --- |
| `FORMAT_TAXONOMY_VSEARCH_LCA` | Convert the reference DB for VSEARCH + LCA |
| `VSEARCH_USEARCHGLOBAL_LCA` | Assign taxonomy with VSEARCH + LCA |
| `FORMAT_TAXRESULTS_VSEARCH_LCA` | Turn VSEARCH + LCA output into a taxonomy table |

### WoRMS matching

| Nextflow process | Description |
| --- | --- |
| `WORMS_MATCH` |  |
| `WORMS_MATCH_OCCURRENCE` |  |

### Darwin core packaging

| Nextflow process | Description |
| --- | --- |
| `BUILD_DARWIN_CORE` |  |

### Summary report

| Nextflow process | Description |
| --- | --- |
| `SUMMARY_REPORT` (ampliseq) | HTML summary of the ampliseq run |
| `SUMMARY_REPORT` (wrapper) | HTML summary of WoRMS matching and Darwin Core output |

## Tests

```bash
pip install pytest && pytest
```
