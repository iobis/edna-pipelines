# edna-pipelines

This is an eDNA metabarcoding pipeline with WoRMS aligned Darwin Core output building on nf-core/ampliseq.

## Run

Single ended:

```bash
nextflow run main.nf -profile docker \
  --input "$(pwd)/data/samplesheet.tsv" \
  --single_end \
  --skip_cutadapt \
  --skip_dada_taxonomy \
  --sintax_ref_tax_custom "$(pwd)/data/ncbi_18s_bu_pga_derep_filtered_sintax_20260702.fasta.gz" \
  --sintax_assign_taxlevels Kingdom,Phylum,Class,Order,Family,Genus,Species \
  --vsearch_lca_ref_tax_custom "$(pwd)/data/ncbi_18s_bu_pga_derep_filtered_sintax_20260702.fasta.gz" \
  --vsearch_lca_assign_taxlevels Kingdom,Phylum,Class,Order,Family,Genus,Species \
  --vsearch_lca_id 1 \
  --vsearch_lca_maxaccepts 0 \
  --vsearch_lca_maxrejects 0 \
  --vsearch_lca_lca_cutoff 1 \
  --metadata "$(pwd)/data/metadata.tsv" \
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

## Tests

```bash
pip install pytest && pytest
```

## Debug

After a pipeline run, use VS Code launch configs in `.vscode/launch.json`:

- **Debug build_darwin_core (local results)** — `bin/build_darwin_core.py`
- **Debug worms_match (vsearch)** — `bin/worms_match.py` on `results/ampliseq/vsearch_lca/ASV_tax_vsearch_lca.user.tsv`
