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

unmatched_mask <- function(df, ranks = TAX_RANKS) {
  has_taxonomy(df, ranks = ranks) & !nzchar(df$scientificNameID)
}

row_finest_names <- function(df, ranks = TAX_RANKS) {
  if (nrow(df) == 0) {
    return(character())
  }
  vapply(
    seq_len(nrow(df)),
    function(i) finest_name(df[i, , drop = FALSE], ranks = ranks),
    character(1)
  )
}

unmatched_names_table <- function(sintax, vsearch, ranks = TAX_RANKS) {
  s <- sintax[unmatched_mask(sintax, ranks), , drop = FALSE]
  v <- vsearch[unmatched_mask(vsearch, ranks), , drop = FALSE]
  pairs <- rbind(
    data.frame(
      asv = s$ASV_ID,
      name = row_finest_names(s, ranks = ranks),
      stringsAsFactors = FALSE
    ),
    data.frame(
      asv = v$ASV_ID,
      name = row_finest_names(v, ranks = ranks),
      stringsAsFactors = FALSE
    )
  )
  pairs <- pairs[nzchar(pairs$name), , drop = FALSE]
  if (nrow(pairs) == 0) {
    return(data.frame(name = character(), asvs = integer(), stringsAsFactors = FALSE))
  }
  pairs <- unique(pairs)
  counts <- as.data.frame(table(pairs$name), stringsAsFactors = FALSE)
  names(counts) <- c("name", "asvs")
  counts[order(-counts$asvs, counts$name), ]
}

fmt_pct <- function(n, d) {
  if (d == 0) {
    return("0%")
  }
  sprintf("%.1f%%", 100 * n / d)
}

unmatched_summary_text <- function(sintax, vsearch, ranks = TAX_RANKS) {
  unmatched_asvs <- unique(c(
    sintax$ASV_ID[unmatched_mask(sintax, ranks)],
    vsearch$ASV_ID[unmatched_mask(vsearch, ranks)]
  ))
  tax_asvs <- unique(c(
    sintax$ASV_ID[has_taxonomy(sintax, ranks)],
    vsearch$ASV_ID[has_taxonomy(vsearch, ranks)]
  ))
  unmatched_names <- unique(c(
    row_finest_names(sintax[unmatched_mask(sintax, ranks), , drop = FALSE], ranks = ranks),
    row_finest_names(vsearch[unmatched_mask(vsearch, ranks), , drop = FALSE], ranks = ranks)
  ))
  unmatched_names <- unmatched_names[nzchar(unmatched_names)]
  tax_names <- unique(c(
    row_finest_names(sintax[has_taxonomy(sintax, ranks), , drop = FALSE], ranks = ranks),
    row_finest_names(vsearch[has_taxonomy(vsearch, ranks), , drop = FALSE], ranks = ranks)
  ))
  tax_names <- tax_names[nzchar(tax_names)]

  n_asv <- length(unmatched_asvs)
  n_asv_total <- length(tax_asvs)
  n_name <- length(unmatched_names)
  n_name_total <- length(tax_names)

  sprintf(
    "%d of %d ASVs (%s) and %d of %d names (%s) could not be matched to WoRMS.",
    n_asv,
    n_asv_total,
    fmt_pct(n_asv, n_asv_total),
    n_name,
    n_name_total,
    fmt_pct(n_name, n_name_total)
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

classification_rate_by_rank <- function(occurrence_path, ranks = tolower(TAX_RANKS)) {
  if (is.null(occurrence_path) || !file.exists(occurrence_path)) {
    return(NULL)
  }
  occ <- read.delim(occurrence_path, sep = "\t", stringsAsFactors = FALSE,
                     check.names = FALSE, quote = "")
  # Extract ASV ID from occurrenceID (format: sampleid_asvhash)
  occ$asv_id <- sub(".*_", "", occ$occurrenceID)
  occ_dedup <- occ[!duplicated(occ$asv_id), , drop = FALSE]
  n_total <- nrow(occ_dedup)
  if (n_total == 0) {
    return(NULL)
  }
  pcts <- vapply(ranks, function(rank) {
    if (!rank %in% names(occ_dedup)) return(0)
    sum(nzchar(trimws(occ_dedup[[rank]]))) / n_total * 100
  }, numeric(1))
  data.frame(
    rank = factor(rev(names(pcts)), levels = rev(names(pcts))),
    classified = round(rev(pcts), 1),
    stringsAsFactors = FALSE
  )
}

taxonomy_sunburst_data <- function(occurrence_path, ranks = tolower(TAX_RANKS)) {
  if (is.null(occurrence_path) || !file.exists(occurrence_path)) {
    return(NULL)
  }
  occ <- read.delim(occurrence_path, sep = "\t", stringsAsFactors = FALSE,
                     check.names = FALSE, quote = "")
  occ$asv_id <- sub(".*_", "", occ$occurrenceID)
  occ_dedup <- occ[!duplicated(occ$asv_id), , drop = FALSE]
  ranks <- intersect(ranks, names(occ_dedup))
  if (length(ranks) == 0 || nrow(occ_dedup) == 0) return(NULL)

  for (r in ranks) {
    occ_dedup[[r]] <- trimws(occ_dedup[[r]])
    occ_dedup[[r]][!nzchar(occ_dedup[[r]])] <- NA
  }

  rows <- list()
  for (depth in seq_along(ranks)) {
    r <- ranks[depth]
    parent_r <- if (depth == 1) NULL else ranks[depth - 1]
    cols <- ranks[seq_len(depth)]
    sub <- occ_dedup[!is.na(occ_dedup[[r]]), cols, drop = FALSE]
    if (nrow(sub) == 0) next
    agg <- aggregate(rep(1, nrow(sub)), by = lapply(cols, function(c) sub[[c]]), FUN = sum)
    names(agg) <- c(cols, "n")
    agg$label <- agg[[r]]
    if (is.null(parent_r)) {
      agg$id <- agg[[r]]
      agg$parent <- ""
    } else {
      agg$id <- apply(agg[cols], 1, paste, collapse = " - ")
      agg$parent <- apply(agg[, cols[-length(cols)], drop = FALSE], 1, paste, collapse = " - ")
    }
    rows[[depth]] <- agg[, c("id", "label", "parent", "n"), drop = FALSE]
  }
  do.call(rbind, rows)
}

detected_species_table <- function(occurrence_path) {
  if (is.null(occurrence_path) || !file.exists(occurrence_path)) {
    return(NULL)
  }
  occ <- read.delim(occurrence_path, sep = "\t", stringsAsFactors = FALSE,
                    check.names = FALSE, quote = "")
  needed <- c("species", "phylum", "class", "sample_id", "organismQuantity")
  missing <- setdiff(needed, names(occ))
  if (length(missing) > 0) {
    return(NULL)
  }
  sp <- occ[
    occ$taxonRank == "species" & nzchar(trimws(occ$species)),
    needed,
    drop = FALSE
  ]
  if (nrow(sp) == 0) {
    return(data.frame(
      species = character(),
      phylum = character(),
      class = character(),
      asvs = integer(),
      samples = integer(),
      stringsAsFactors = FALSE
    ))
  }
  sp$species <- trimws(as.character(sp$species))
  sp$phylum <- trimws(as.character(sp$phylum))
  sp$class <- trimws(as.character(sp$class))
  sp$organismQuantity <- suppressWarnings(as.integer(sp$organismQuantity))
  sp$organismQuantity[is.na(sp$organismQuantity)] <- 0L

  qty_counts <- aggregate(
    organismQuantity ~ species,
    data = sp,
    FUN = sum
  )
  sample_counts <- aggregate(
    sample_id ~ species,
    data = sp,
    FUN = function(x) length(unique(x))
  )
  # One taxonomy label per species (first non-empty if multiple)
  tax <- aggregate(
    cbind(phylum, class) ~ species,
    data = sp,
    FUN = function(x) {
      x <- unique(x[nzchar(x)])
      if (length(x) == 0) "" else x[[1]]
    }
  )

  out <- data.frame(
    species = qty_counts$species,
    phylum = tax$phylum[match(qty_counts$species, tax$species)],
    class = tax$class[match(qty_counts$species, tax$species)],
    asvs = as.integer(qty_counts$organismQuantity),
    samples = as.integer(sample_counts$sample_id[match(qty_counts$species, sample_counts$species)]),
    stringsAsFactors = FALSE
  )
  out[order(-out$asvs, out$species), , drop = FALSE]
}

#' Build all summary-report tables/text (call this under the debugger).
build_summary_report <- function(
    sintax_tsv,
    vsearch_tsv,
    params_tsv = NULL,
    occurrence_tsv = NULL,
    ranks = TAX_RANKS
) {
  parameters <- read_parameters_table(params_tsv)
  revision_info <- revision_section(parameters)

  if (is.null(occurrence_tsv)) {
    darwincore_outdir <- param_value(parameters, "darwincore_outdir")
    if (!is.null(darwincore_outdir)) {
      candidate <- file.path(darwincore_outdir, "publishing", "occurrence.tsv")
      if (file.exists(candidate)) {
        occurrence_tsv <- candidate
      }
    }
  }

  sintax <- read_tax(sintax_tsv, ranks = ranks)
  vsearch <- read_tax(vsearch_tsv, ranks = ranks)

  unmatched <- unmatched_names_table(sintax, vsearch, ranks = ranks)
  agreement <- agreement_by_rank(sintax, vsearch, ranks = ranks)

  pipeline_version <- param_value(parameters, "pipeline_version")
  if (is.null(pipeline_version)) {
    pipeline_version <- "unknown"
  }
  run_date <- param_value(parameters, "run_date")
  if (is.null(run_date)) {
    run_date <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  }

  list(
    pipeline_version = pipeline_version,
    run_date = run_date,
    parameters = parameters,
    revision = revision_info$revision,
    revision_message = revision_info$revision_message,
    unmatched = unmatched,
    unmatched_summary = unmatched_summary_text(sintax, vsearch, ranks = ranks),
    agreement = agreement$by_rank,
    sintax_species_changes = sintax_species_change_text(sintax, vsearch),
    classification_rate = classification_rate_by_rank(occurrence_tsv),
    sunburst = taxonomy_sunburst_data(occurrence_tsv),
    detected_species = detected_species_table(occurrence_tsv),
    n_sintax = nrow(sintax),
    n_vsearch = nrow(vsearch),
    n_shared = length(agreement$shared_ids)
  )
}
