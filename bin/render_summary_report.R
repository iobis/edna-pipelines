#!/usr/bin/env Rscript
# Render bin/summary_report.Rmd to HTML.

abs_path <- function(path, label) {
  if (is.null(path) || !nzchar(path)) {
    stop(label, " path is empty", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop(label, " not found: ", path, " (cwd=", getwd(), ")", call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

ensure_summary_report_helpers <- function() {
  if (exists("build_summary_report", mode = "function", inherits = TRUE)) {
    return(invisible(NULL))
  }
  candidates <- c(
    file.path("bin", "summary_report.R"),
    "summary_report.R"
  )
  ca <- commandArgs(trailingOnly = FALSE)
  file_arg <- ca[startsWith(ca, "--file=")]
  if (length(file_arg) == 1) {
    candidates <- c(
      candidates,
      file.path(dirname(normalizePath(sub("^--file=", "", file_arg))), "summary_report.R")
    )
  }
  helper <- candidates[file.exists(candidates)][1]
  if (is.na(helper)) {
    stop("Could not find summary_report.R", call. = FALSE)
  }
  sys.source(helper, envir = .GlobalEnv)
}

#' Render the eDNA summary report.
#'
#' Builds the report object in R (debuggable), then knits a thin Rmd.
render_summary_report <- function(
    sintax,
    vsearch,
    params_tsv,
    output,
    rmd,
    revision = NULL
) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package 'rmarkdown' is required. Install with install.packages('rmarkdown').", call. = FALSE)
  }

  ensure_summary_report_helpers()

  sintax <- abs_path(sintax, "sintax")
  vsearch <- abs_path(vsearch, "vsearch")
  params_tsv <- abs_path(params_tsv, "params")
  rmd <- abs_path(rmd, "rmd")
  revision_path <- if (!is.null(revision) && nzchar(revision) && file.exists(revision)) {
    abs_path(revision, "revision")
  } else {
    NULL
  }

  # Compute everything here so breakpoints in summary_report.R work.
  report <- build_summary_report(
    sintax_tsv = sintax,
    vsearch_tsv = vsearch,
    params_tsv = params_tsv,
    revision_file = revision_path
  )

  if (dirname(output) == ".") {
    output_abs <- file.path(getwd(), basename(output))
  } else {
    output_abs <- output
  }
  outdir <- dirname(output_abs)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  outdir <- abs_path(outdir, "output dir")

  rmarkdown::render(
    input = rmd,
    output_file = basename(output_abs),
    output_dir = outdir,
    knit_root_dir = getwd(),
    params = list(report = report),
    quiet = TRUE
  )

  out <- file.path(outdir, basename(output_abs))
  message("Wrote ", out)
  invisible(out)
}

parse_render_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  usage <- function() {
    message(
      "Usage: render_summary_report.R ",
      "--sintax PATH --vsearch PATH --params PATH --output PATH [--revision PATH] [--rmd PATH]"
    )
  }

  get_opt <- function(flag, required = TRUE) {
    i <- match(flag, args)
    if (is.na(i) || i == length(args)) {
      if (required) {
        usage()
        stop("Missing ", flag, call. = FALSE)
      }
      return(NULL)
    }
    args[[i + 1]]
  }

  rmd <- get_opt("--rmd", required = FALSE)
  if (is.null(rmd)) {
    ca <- commandArgs(trailingOnly = FALSE)
    file_arg <- ca[startsWith(ca, "--file=")]
    if (length(file_arg) == 1) {
      rmd <- file.path(dirname(normalizePath(sub("^--file=", "", file_arg))), "summary_report.Rmd")
    }
  }
  if (is.null(rmd) || !file.exists(rmd)) {
    stop("summary_report.Rmd not found; pass --rmd", call. = FALSE)
  }

  list(
    sintax = get_opt("--sintax"),
    vsearch = get_opt("--vsearch"),
    params_tsv = get_opt("--params"),
    output = get_opt("--output"),
    revision = get_opt("--revision", required = FALSE),
    rmd = rmd
  )
}

# CLI entry when run via Rscript (skipped when sourced for debugging).
if (!isTRUE(getOption("edna.skip_render_cli"))) {
  opts <- parse_render_args()
  render_summary_report(
    sintax = opts$sintax,
    vsearch = opts$vsearch,
    params_tsv = opts$params_tsv,
    output = opts$output,
    rmd = opts$rmd,
    revision = opts$revision
  )
}
