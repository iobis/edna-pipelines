# Summary report analysis helpers (debuggable from R Debugger).

TAX_RANKS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
TAX_PLACEHOLDERS <- c("", "na", "n/a", "unassigned")

is_assigned <- function(x) {
  !is.na(x) & !(tolower(trimws(as.character(x))) %in% TAX_PLACEHOLDERS)
}

finest_name <- function(row, ranks = TAX_RANKS) {
  for (rank in rev(ranks)) {
    val <- as.character(row[[rank]])
    if (nzchar(val)) return(val)
  }
  ""
}

read_tax <- function(path, ranks = TAX_RANKS) {
  df <- read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, quote = "")
  for (rank in ranks) {
    if (!rank %in% names(df)) df[[rank]] <- ""
    df[[rank]] <- ifelse(is_assigned(df[[rank]]), trimws(df[[rank]]), "")
  }
  if (!"ASV_ID" %in% names(df)) stop("ASV_ID column missing in ", path)
  if (!"scientificNameID" %in% names(df)) df$scientificNameID <- ""
  if (!"scientificName" %in% names(df)) df$scientificName <- ""
  df$scientificNameID <- trimws(as.character(df$scientificNameID))
  df$scientificName <- trimws(as.character(df$scientificName))
  df
}

has_taxonomy <- function(df, ranks = TAX_RANKS) {
  apply(df[ranks], 1, function(r) any(nzchar(r)))
}

read_parameters_table <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  read.delim(path, sep = "\t", stringsAsFactors = FALSE, quote = "")
}

read_revision_text <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

unmatched_names_table <- function(df, method, ranks = TAX_RANKS) {
  unmatched <- has_taxonomy(df, ranks = ranks) & !nzchar(df$scientificNameID)
  names_u <- vapply(
    which(unmatched),
    function(i) finest_name(df[i, , drop = FALSE], ranks = ranks),
    character(1)
  )
  names_u <- sort(unique(names_u[nzchar(names_u)]))
  if (length(names_u) == 0) {
    return(data.frame(method = character(), name = character(), stringsAsFactors = FALSE))
  }
  data.frame(method = method, name = names_u, stringsAsFactors = FALSE)
}

unmatched_summary_text <- function(sintax, vsearch, unmatched, ranks = TAX_RANKS) {
  sprintf(
    "SINTAX: %d ASVs unmatched (%d unique names). VSEARCH: %d ASVs unmatched (%d unique names).",
    sum(has_taxonomy(sintax, ranks = ranks) & !nzchar(sintax$scientificNameID)),
    sum(unmatched$method == "sintax"),
    sum(has_taxonomy(vsearch, ranks = ranks) & !nzchar(vsearch$scientificNameID)),
    sum(unmatched$method == "vsearch")
  )
}

agreement_by_rank <- function(sintax, vsearch, ranks = TAX_RANKS) {
  shared_ids <- intersect(sintax$ASV_ID, vsearch$ASV_ID)
  s <- sintax[match(shared_ids, sintax$ASV_ID), , drop = FALSE]
  v <- vsearch[match(shared_ids, vsearch$ASV_ID), , drop = FALSE]

  agree_rows <- lapply(ranks, function(rank) {
    a <- tolower(s[[rank]])
    b <- tolower(v[[rank]])
    both <- nzchar(a) & nzchar(b)
    data.frame(
      rank = rank,
      both_assigned = sum(both),
      same = sum(both & a == b),
      different = sum(both & a != b),
      only_sintax = sum(nzchar(a) & !nzchar(b)),
      only_vsearch = sum(!nzchar(a) & nzchar(b)),
      neither = sum(!nzchar(a) & !nzchar(b)),
      stringsAsFactors = FALSE
    )
  })
  list(
    shared_ids = shared_ids,
    sintax = s,
    vsearch = v,
    by_rank = do.call(rbind, agree_rows)
  )
}

finest_agreement_text <- function(sintax_aligned, vsearch_aligned, ranks = TAX_RANKS) {
  s_fine <- vapply(
    seq_len(nrow(sintax_aligned)),
    function(i) tolower(finest_name(sintax_aligned[i, , drop = FALSE], ranks = ranks)),
    character(1)
  )
  v_fine <- vapply(
    seq_len(nrow(vsearch_aligned)),
    function(i) tolower(finest_name(vsearch_aligned[i, , drop = FALSE], ranks = ranks)),
    character(1)
  )
  both_fine <- nzchar(s_fine) & nzchar(v_fine)
  sprintf(
    "Finest name (any rank): %d same, %d different, %d only SINTAX, %d only VSEARCH, %d neither (of %d shared ASVs).",
    sum(both_fine & s_fine == v_fine),
    sum(both_fine & s_fine != v_fine),
    sum(nzchar(s_fine) & !nzchar(v_fine)),
    sum(!nzchar(s_fine) & nzchar(v_fine)),
    sum(!nzchar(s_fine) & !nzchar(v_fine)),
    nrow(sintax_aligned)
  )
}

#' Build all summary-report tables/text (call this under the debugger).
build_summary_report <- function(
    sintax_tsv,
    vsearch_tsv,
    params_tsv = NULL,
    revision_file = NULL,
    ranks = TAX_RANKS
) {
  sintax <- read_tax(sintax_tsv, ranks = ranks)
  vsearch <- read_tax(vsearch_tsv, ranks = ranks)

  unmatched <- rbind(
    unmatched_names_table(sintax, "sintax", ranks = ranks),
    unmatched_names_table(vsearch, "vsearch", ranks = ranks)
  )
  agreement <- agreement_by_rank(sintax, vsearch, ranks = ranks)

  list(
    parameters = read_parameters_table(params_tsv),
    revision_text = read_revision_text(revision_file),
    unmatched = unmatched,
    unmatched_summary = unmatched_summary_text(sintax, vsearch, unmatched, ranks = ranks),
    agreement = agreement$by_rank,
    finest_summary = finest_agreement_text(agreement$sintax, agreement$vsearch, ranks = ranks),
    n_sintax = nrow(sintax),
    n_vsearch = nrow(vsearch),
    n_shared = length(agreement$shared_ids)
  )
}
