library(dplyr)
library(tidyr)
library(tibble)
library(purrr)

TAX_RANKS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
TAX_PLACEHOLDERS <- c("", "na", "n/a", "unassigned")

is_assigned <- function(x) {
  !is.na(x) & !(tolower(trimws(as.character(x))) %in% TAX_PLACEHOLDERS)
}

blank_to_na <- function(x) {
  x <- trimws(as.character(x))
  x[!nzchar(x) | is.na(x)] <- NA_character_
  x
}

finest_name <- function(df, ranks = TAX_RANKS) {
  if (nrow(df) == 0) {
    return(character())
  }
  vals <- lapply(ranks, function(rank) blank_to_na(df[[rank]]))
  names(vals) <- ranks
  coalesce(!!!rev(vals))
}

read_tsv_safe <- function(path) {
  read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, quote = "")
}

read_tax <- function(path, ranks = TAX_RANKS) {
  df <- as_tibble(read_tsv_safe(path))
  for (rank in ranks) {
    if (!rank %in% names(df)) {
      df[[rank]] <- ""
    }
    df[[rank]] <- ifelse(is_assigned(df[[rank]]), trimws(df[[rank]]), "")
  }
  if (!"ASV_ID" %in% names(df)) {
    stop("ASV_ID column missing in ", path, call. = FALSE)
  }
  if (!"scientificNameID" %in% names(df)) {
    df$scientificNameID <- ""
  }
  if (!"scientificName" %in% names(df)) {
    df$scientificName <- ""
  }
  df %>%
    mutate(
      scientificNameID = trimws(as.character(.data$scientificNameID)),
      scientificName = trimws(as.character(.data$scientificName))
    )
}

has_taxonomy <- function(df, ranks = TAX_RANKS) {
  Reduce(`|`, lapply(ranks, function(rank) nzchar(df[[rank]])))
}

read_parameters_table <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  as_tibble(read_tsv_safe(path))
}

param_value <- function(params, name) {
  if (is.null(params) || !"parameter" %in% names(params)) {
    return(NULL)
  }
  value <- params %>%
    filter(.data$parameter == !!name) %>%
    pull(.data$value)
  if (length(value) == 0 || is.na(value[[1]]) || !nzchar(value[[1]])) {
    return(NULL)
  }
  value[[1]]
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
  lines <- trimws(readLines(path, warn = FALSE))
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) {
    return(NULL)
  }
  tibble(line = lines) %>%
    separate(
      .data$line,
      into = c("property", "value"),
      sep = "=",
      extra = "merge",
      fill = "right"
    ) %>%
    mutate(value = coalesce(.data$value, ""))
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

empty_unmatched <- function() {
  tibble(name = character(), asvs = integer())
}

unmatched_names_table <- function(sintax, vsearch, ranks = TAX_RANKS) {
  pairs <- bind_rows(
    sintax %>%
      filter(unmatched_mask(., ranks = ranks)) %>%
      transmute(asv = .data$ASV_ID, name = finest_name(., ranks = ranks)),
    vsearch %>%
      filter(unmatched_mask(., ranks = ranks)) %>%
      transmute(asv = .data$ASV_ID, name = finest_name(., ranks = ranks))
  ) %>%
    filter(nzchar(.data$name)) %>%
    distinct()

  if (nrow(pairs) == 0) {
    return(empty_unmatched())
  }

  pairs %>%
    count(.data$name, name = "asvs") %>%
    arrange(desc(.data$asvs), .data$name)
}

fmt_pct <- function(n, d) {
  if (d == 0) {
    return("0%")
  }
  sprintf("%.1f%%", 100 * n / d)
}

unmatched_summary_text <- function(sintax, vsearch, ranks = TAX_RANKS) {
  unmatched_sintax <- sintax %>% filter(unmatched_mask(., ranks = ranks))
  unmatched_vsearch <- vsearch %>% filter(unmatched_mask(., ranks = ranks))
  tax_sintax <- sintax %>% filter(has_taxonomy(., ranks = ranks))
  tax_vsearch <- vsearch %>% filter(has_taxonomy(., ranks = ranks))

  unmatched_asvs <- unique(c(unmatched_sintax$ASV_ID, unmatched_vsearch$ASV_ID))
  tax_asvs <- unique(c(tax_sintax$ASV_ID, tax_vsearch$ASV_ID))
  unmatched_names <- unique(c(
    finest_name(unmatched_sintax, ranks = ranks),
    finest_name(unmatched_vsearch, ranks = ranks)
  ))
  unmatched_names <- unmatched_names[!is.na(unmatched_names) & nzchar(unmatched_names)]
  tax_names <- unique(c(
    finest_name(tax_sintax, ranks = ranks),
    finest_name(tax_vsearch, ranks = ranks)
  ))
  tax_names <- tax_names[!is.na(tax_names) & nzchar(tax_names)]

  sprintf(
    "%d of %d ASVs (%s) and %d of %d names (%s) could not be matched to WoRMS.",
    length(unmatched_asvs),
    length(tax_asvs),
    fmt_pct(length(unmatched_asvs), length(tax_asvs)),
    length(unmatched_names),
    length(tax_names),
    fmt_pct(length(unmatched_names), length(tax_names))
  )
}

agreement_by_rank <- function(sintax, vsearch, ranks = TAX_RANKS) {
  shared <- inner_join(
    sintax %>% select(ASV_ID, all_of(ranks)),
    vsearch %>% select(ASV_ID, all_of(ranks)),
    by = "ASV_ID",
    suffix = c("_sintax", "_vsearch")
  )

  by_rank <- map_dfr(ranks, function(rank) {
    a <- tolower(shared[[paste0(rank, "_sintax")]])
    b <- tolower(shared[[paste0(rank, "_vsearch")]])
    both <- nzchar(a) & nzchar(b)
    tibble(
      rank = rank,
      both_assigned = sum(both),
      same = sum(both & a == b),
      different = sum(both & a != b),
      only_sintax = sum(nzchar(a) & !nzchar(b)),
      only_vsearch = sum(!nzchar(a) & nzchar(b)),
      neither = sum(!nzchar(a) & !nzchar(b))
    )
  })

  list(
    shared_ids = shared$ASV_ID,
    sintax = sintax %>% filter(.data$ASV_ID %in% shared$ASV_ID),
    vsearch = vsearch %>% filter(.data$ASV_ID %in% shared$ASV_ID),
    by_rank = by_rank
  )
}

sintax_species_change_text <- function(sintax, vsearch) {
  compared <- sintax %>%
    select(ASV_ID, Species, Genus) %>%
    left_join(
      vsearch %>% select(ASV_ID, Species, Genus),
      by = "ASV_ID",
      suffix = c("_sintax", "_vsearch")
    ) %>%
    mutate(
      across(c("Species_vsearch", "Genus_vsearch"), ~ coalesce(.x, "")),
      genera_agree = nzchar(.data$Genus_sintax) &
        nzchar(.data$Genus_vsearch) &
        tolower(.data$Genus_sintax) == tolower(.data$Genus_vsearch),
      v_accepted = if_else(
        nzchar(.data$Species_vsearch) & .data$genera_agree,
        .data$Species_vsearch,
        ""
      ),
      has_sintax_species = nzchar(.data$Species_sintax)
    )

  n_removed <- sum(compared$has_sintax_species & !nzchar(compared$v_accepted))
  n_updated <- sum(
    compared$has_sintax_species &
      nzchar(compared$v_accepted) &
      tolower(compared$Species_sintax) != tolower(compared$v_accepted)
  )

  sprintf(
    "%d SINTAX species names were removed and %d were updated because VSEARCH did not agree.",
    n_removed,
    n_updated
  )
}

read_occurrence <- function(occurrence_path) {
  if (is.null(occurrence_path) || !file.exists(occurrence_path)) {
    return(NULL)
  }
  as_tibble(read_tsv_safe(occurrence_path)) %>%
    mutate(asv_id = sub(".*_", "", .data$occurrenceID))
}

classification_rate_by_rank <- function(occurrence_path, ranks = tolower(TAX_RANKS)) {
  occ <- read_occurrence(occurrence_path)
  if (is.null(occ)) {
    return(NULL)
  }
  occ_dedup <- occ %>% distinct(.data$asv_id, .keep_all = TRUE)
  n_total <- nrow(occ_dedup)
  if (n_total == 0) {
    return(NULL)
  }

  pcts <- vapply(ranks, function(rank) {
    if (!rank %in% names(occ_dedup)) {
      return(0)
    }
    sum(nzchar(trimws(occ_dedup[[rank]]))) / n_total * 100
  }, numeric(1))

  tibble(
    rank = factor(rev(names(pcts)), levels = rev(names(pcts))),
    classified = round(rev(pcts), 1)
  )
}

taxonomy_sunburst_data <- function(occurrence_path, ranks = tolower(TAX_RANKS)) {
  occ <- read_occurrence(occurrence_path)
  if (is.null(occ)) {
    return(NULL)
  }
  ranks <- intersect(ranks, names(occ))
  occ_dedup <- occ %>% distinct(.data$asv_id, .keep_all = TRUE)
  if (length(ranks) == 0 || nrow(occ_dedup) == 0) {
    return(NULL)
  }

  occ_dedup <- occ_dedup %>%
    mutate(across(all_of(ranks), blank_to_na))

  rows <- lapply(seq_along(ranks), function(depth) {
    r <- ranks[[depth]]
    cols <- ranks[seq_len(depth)]
    sub <- occ_dedup %>%
      filter(!is.na(.data[[r]])) %>%
      select(all_of(cols))
    if (nrow(sub) == 0) {
      return(NULL)
    }
    sub %>%
      count(across(all_of(cols)), name = "n") %>%
      mutate(
        label = .data[[r]],
        id = if (depth == 1) {
          .data[[r]]
        } else {
          do.call(paste, c(across(all_of(cols)), sep = " - "))
        },
        parent = if (depth == 1) {
          ""
        } else {
          do.call(paste, c(across(all_of(cols[-length(cols)])), sep = " - "))
        }
      ) %>%
      select("id", "label", "parent", "n")
  })

  bind_rows(rows)
}

empty_detected_species <- function() {
  tibble(
    species = character(),
    phylum = character(),
    class = character(),
    asvs = integer(),
    reads = integer(),
    samples = integer()
  )
}

empty_sample_summary <- function() {
  tibble(
    sample_id = character(),
    asvs = integer(),
    species = integer()
  )
}

sample_summary_table <- function(occurrence_path) {
  occ <- read_occurrence(occurrence_path)
  if (is.null(occ)) {
    return(NULL)
  }
  needed <- c("sample_id", "asv_id", "species", "taxonRank")
  if (length(setdiff(needed, names(occ))) > 0) {
    return(NULL)
  }

  occ <- occ %>%
    mutate(
      sample_id = trimws(as.character(.data$sample_id)),
      species = trimws(as.character(.data$species))
    ) %>%
    filter(nzchar(.data$sample_id))

  if (nrow(occ) == 0) {
    return(empty_sample_summary())
  }

  asv_counts <- occ %>%
    group_by(.data$sample_id) %>%
    summarise(asvs = n_distinct(.data$asv_id), .groups = "drop")

  species_counts <- occ %>%
    filter(.data$taxonRank == "species", nzchar(.data$species)) %>%
    group_by(.data$sample_id) %>%
    summarise(species = n_distinct(.data$species), .groups = "drop")

  asv_counts %>%
    left_join(species_counts, by = "sample_id") %>%
    mutate(species = coalesce(.data$species, 0L)) %>%
    arrange(.data$sample_id)
}

first_non_empty <- function(x) {
  x <- unique(x[nzchar(x) & !is.na(x)])
  if (length(x) == 0) "" else x[[1]]
}

detected_species_table <- function(occurrence_path) {
  occ <- read_occurrence(occurrence_path)
  if (is.null(occ)) {
    return(NULL)
  }
  needed <- c("species", "phylum", "class", "sample_id", "organismQuantity", "asv_id")
  if (length(setdiff(needed, names(occ))) > 0 || !"taxonRank" %in% names(occ)) {
    return(NULL)
  }

  sp <- occ %>%
    filter(.data$taxonRank == "species", nzchar(trimws(.data$species))) %>%
    mutate(
      species = trimws(as.character(.data$species)),
      phylum = trimws(as.character(.data$phylum)),
      class = trimws(as.character(.data$class)),
      organismQuantity = coalesce(
        suppressWarnings(as.integer(.data$organismQuantity)),
        0L
      )
    )

  if (nrow(sp) == 0) {
    return(empty_detected_species())
  }

  sp %>%
    group_by(.data$species) %>%
    summarise(
      phylum = first_non_empty(.data$phylum),
      class = first_non_empty(.data$class),
      asvs = n_distinct(.data$asv_id),
      reads = sum(.data$organismQuantity),
      samples = n_distinct(.data$sample_id),
      .groups = "drop"
    ) %>%
    arrange(desc(.data$reads), desc(.data$asvs), .data$species)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

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
  agreement <- agreement_by_rank(sintax, vsearch, ranks = ranks)

  list(
    pipeline_version = param_value(parameters, "pipeline_version") %||% "unknown",
    run_date = param_value(parameters, "run_date") %||% format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    parameters = parameters,
    revision = revision_info$revision,
    revision_message = revision_info$revision_message,
    unmatched = unmatched_names_table(sintax, vsearch, ranks = ranks),
    unmatched_summary = unmatched_summary_text(sintax, vsearch, ranks = ranks),
    agreement = agreement$by_rank,
    sintax_species_changes = sintax_species_change_text(sintax, vsearch),
    classification_rate = classification_rate_by_rank(occurrence_tsv),
    sunburst = taxonomy_sunburst_data(occurrence_tsv),
    sample_summary = sample_summary_table(occurrence_tsv),
    detected_species = detected_species_table(occurrence_tsv),
    n_sintax = nrow(sintax),
    n_vsearch = nrow(vsearch),
    n_shared = length(agreement$shared_ids)
  )
}
