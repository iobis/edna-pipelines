# edna-pipelines

This is an eDNA metabarcoding pipeline with WoRMS aligned Darwin Core output building on nf-core/ampliseq.

## Run

Single ended:

```bash
nextflow run prototyping/edna-platform-pipelines/main.nf -profile docker \
  --input "$(pwd)/prototyping/ampliseq-wilderlab-poc/samplesheet_local_two_samples_ci.tsv" \
  --single_end \
  --skip_cutadapt \
  --skip_dada_taxonomy \
  --sintax_ref_tax_custom /Volumes/acasis/reference_databases/MIDORI2_UNIQ_NUC_GB269_CO1_SINTAX.fasta.gz \
  --sintax_assign_taxlevels Kingdom,Phylum,Class,Order,Family,Genus,Species \
  --vsearch_lca_ref_tax_custom /Volumes/acasis/reference_databases/MIDORI2_UNIQ_NUC_GB269_CO1_SINTAX.fasta.gz \
  --vsearch_lca_assign_taxlevels Kingdom,Phylum,Class,Order,Family,Genus,Species \
  --vsearch_lca_id 1 \
  --vsearch_lca_maxaccepts 0 \
  --vsearch_lca_maxrejects 0 \
  --vsearch_lca_lca_cutoff 1 \
  --metadata "$(pwd)/prototyping/ampliseq-wilderlab-poc/metadata_wilderlab.tsv" \
  --worms_db /Volumes/acasis/worms/worms_draft_20260522.db \
  --outdir results
```

## Output

After a full run with `--outdir results`:

```
results/
├── ampliseq/
└── darwincore/
    ├── worms/
    └── publishing/
```
