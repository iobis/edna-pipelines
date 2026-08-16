# Local Cursor/VS Code R Debugger entrypoint (no CLI flags).
# Set breakpoints in bin/summary_report.R (analysis) or bin/render_summary_report.R.
# Working directory must be the edna-pipelines repo root.

options(edna.skip_render_cli = TRUE)
source(file.path("bin", "summary_report.R"), local = FALSE)
source(file.path("bin", "render_summary_report.R"), local = FALSE)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

render_summary_report(
  sintax = file.path(root, "results/darwincore/worms/sintax/worms_matched.sintax.tsv"),
  vsearch = file.path(root, "results/darwincore/worms/vsearch/worms_matched.vsearch.tsv"),
  params_tsv = file.path(root, "results/darwincore/report/report_params.tsv"),
  output = file.path(root, "results/darwincore/report/summary_report.html"),
  rmd = file.path(root, "bin/summary_report.Rmd"),
  revision = file.path(root, "results/ampliseq/ampliseq.revision")
)
