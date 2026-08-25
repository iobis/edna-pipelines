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
| `--primer_fwd` / `primer_rev` | |

In addition, `--sintax_cutoff` (default `0.8`) is [currently not configurable](https://github.com/nf-core/ampliseq/blob/master/conf/modules.config#L581) in ampliseq.

## Steps

### FastQC and MultiQC (ampliseq)

### ASV inference with DADA2 (ampliseq)

### Taxonomic classification: SINTAX (ampliseq)

### Taxonomic classification: VSEARCH + LCA (ampliseq)

### WoRMS matching

### Darwin core packaging

### Summary report

## Tests

```bash
pip install pytest && pytest
```
