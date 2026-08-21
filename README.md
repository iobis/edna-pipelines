# edna-pipelines

This is an eDNA metabarcoding pipeline with WoRMS aligned Darwin Core output building on nf-core/ampliseq.

## Run

Single ended:

```bash
nextflow run main.nf -profile docker \
  --input "$(pwd)/data/samplesheet.tsv" \
  --sequencing_type illumina_se \
  --binned_quality 2,12,24,40 \
  --skip_cutadapt \
  --skip_dada_taxonomy \
  --skip_phyloseq \
  --skip_tse \
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

`-profile docker` runs nested nf-core/ampliseq in Docker and runs `SUMMARY_REPORT` in a local report image. WoRMS matching and Darwin Core still use host Python (from the pipeline `.venv` when launched via the platform).

### Report container

With `-profile docker`, the pipeline ensures a report image exists before `SUMMARY_REPORT`:

- Image tag is `edna-pipelines-report:<dockerfile-sha256-12>`
- Rebuilds only when `containers/Dockerfile` changes (or the tagged image is missing)
- Report R/Rmd scripts are copied into the task at runtime, so code changes do not require a rebuild

### nf-core/ampliseq source

Defaults: `--ampliseq_repo https://github.com/nf-core/ampliseq.git`, `--ampliseq_ref dev`, `--ampliseq_update false`.

On the first run, the pipeline clones ampliseq into `third_party/nf-core-ampliseq` and checks out `ampliseq_ref`. With `--ampliseq_update true`, each run fetches and checks out that ref again (use when tracking `dev` or testing a PR branch). The checked-out commit is logged and written to `results/ampliseq/ampliseq.revision`.

Example: test a PR branch from a fork:

```bash
nextflow run main.nf -profile docker ... \
  --ampliseq_repo https://github.com/pieterprovoost/ampliseq.git \
  --ampliseq_ref fix-taxonomy-hash-comment \
  --ampliseq_update true
```

## Output

After a full run with `--outdir results`:

```
results/
├── ampliseq/
│   ├── ampliseq.revision   # repo, ref, and commit used for this run
│   └── ampliseq.done
└── darwincore/
    ├── worms/
    ├── publishing/
    └── report/
        ├── summary_report.html
        └── report_params.tsv
```

The HTML summary covers run parameters, unmatched WoRMS names, and SINTAX vs VSEARCH name agreement. It links to the AmpliSeq summary report under `../ampliseq/summary_report/`.

## Tests

```bash
pip install pytest && pytest
```

## Debug

After a pipeline run, use VS Code launch configs in `.vscode/launch.json`:

- **Debug build_darwin_core (local results)** — `bin/build_darwin_core.py`
- **Debug worms_match (vsearch)** — `bin/worms_match.py` on `results/ampliseq/vsearch_lca/ASV_tax_vsearch_lca.user.tsv`
- **Render summary report (local results)** — R Debugger on `bin/debug_render_summary_report.R`; set breakpoints in `bin/summary_report.R` (needs [R Debugger](https://marketplace.visualstudio.com/items?itemName=RDebugger.r-debugger) + `vscDebugger`)
