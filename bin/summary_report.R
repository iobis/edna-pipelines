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

param_value <- function(params, name) {
  if (is.null(params) || !"parameter" %in% names(params)) {
    return(NULL)
  }
  rows <- params$parameter == name
  if (!any(rows)) {
    return(NULL)
  }
  value <- params$value[which(rows)[1]]
  if (is.na(value) || !nzchar(value)) {
    return(NULL)
  }
  value
}

ampliseq_revision_path <- function(params) {
  ampliseq_results <- param_value(params, "ampliseq_results")
  if (is.null(ampliseq_results)) {
    return(NULL)
  }
  file.path(ampliseq_results, "ampliseq.revision")
}

read_revision_table <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines[nzchar(trimws(lines))])
  if (length(lines) == 0) {
    return(NULL)
  }
  parts <- strsplit(lines, "=", fixed = TRUE)
  keys <- vapply(parts, function(x) x[[1]], character(1))
  values <- vapply(parts, function(x) paste(x[-1], collapse = "="), character(1))
  data.frame(property = keys, value = values, stringsAsFactors = FALSE)
}

revision_section <- function(params) {
  path <- ampliseq_revision_path(params)
  revision <- read_revision_table(path)
  if (!is.null(revision)) {
    return(list(revision = revision, revision_message = NULL))
  }
  message <- if (is.null(path)) {
    "ampliseq_results was not recorded in report parameters; ampliseq version is unavailable."
  } else {
    sprintf("ampliseq.revision not found at %s.", path)
  }
  list(revision = NULL, revision_message = message)
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

unmatched_summary_table <- function(sintax, vsearch, unmatched, ranks = TAX_RANKS) {
  data.frame(
    method = c("sintax", "vsearch"),
    asvs = c(
      sum(has_taxonomy(sintax, ranks = ranks) & !nzchar(sintax$scientificNameID)),
      sum(has_taxonomy(vsearch, ranks = ranks) & !nzchar(vsearch$scientificNameID))
    ),
    names = c(
      sum(unmatched$method == "sintax"),
      sum(unmatched$method == "vsearch")
    ),
    stringsAsFactors = FALSE
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

sintax_species_change_text <- function(sintax, vsearch) {
  v_idx <- match(sintax$ASV_ID, vsearch$ASV_ID)
  s_species <- sintax$Species
  s_genus <- sintax$Genus
  v_species <- ifelse(is.na(v_idx), "", vsearch$Species[v_idx])
  v_genus <- ifelse(is.na(v_idx), "", vsearch$Genus[v_idx])

  genera_agree <- nzchar(s_genus) & nzchar(v_genus) & tolower(s_genus) == tolower(v_genus)
  v_accepted <- ifelse(nzchar(v_species) & genera_agree, v_species, "")

  has_sintax_species <- nzchar(s_species)
  n_removed <- sum(has_sintax_species & !nzchar(v_accepted))
  n_updated <- sum(
    has_sintax_species & nzchar(v_accepted) & tolower(s_species) != tolower(v_accepted)
  )

  sprintf(
    "%d SINTAX species names were removed and %d were updated because VSEARCH did not agree.",
    n_removed,
    n_updated
  )
}

#' Build all summary-report tables/text (call this under the debugger).
build_summary_report <- function(
    sintax_tsv,
    vsearch_tsv,
    params_tsv = NULL,
    ranks = TAX_RANKS
) {
  parameters <- read_parameters_table(params_tsv)
  revision_info <- revision_section(parameters)

  sintax <- read_tax(sintax_tsv, ranks = ranks)
  vsearch <- read_tax(vsearch_tsv, ranks = ranks)

  unmatched <- rbind(
    unmatched_names_table(sintax, "sintax", ranks = ranks),
    unmatched_names_table(vsearch, "vsearch", ranks = ranks)
  )
  agreement <- agreement_by_rank(sintax, vsearch, ranks = ranks)

  list(
    parameters = parameters,
    revision = revision_info$revision,
    revision_message = revision_info$revision_message,
    unmatched = unmatched,
    unmatched_summary = unmatched_summary_table(sintax, vsearch, unmatched, ranks = ranks),
    agreement = agreement$by_rank,
    sintax_species_changes = sintax_species_change_text(sintax, vsearch),
    n_sintax = nrow(sintax),
    n_vsearch = nrow(vsearch),
    n_shared = length(agreement$shared_ids)
  )
}
