# =============================================================================
# Identify leading edge genes
# -----------------------------------------------------------------------------
# For each GSEA result directory, this script:
#   1. Reads the merged GO enrichment table and splits the leading-edge
#      ("core_enrichment") gene lists into one row per gene.
#   2. Counts, per gene, how many GO terms of each ontology (BP/MF/CC, ...)
#      the gene appears in as a leading-edge member.
#   3. Joins those counts with the matching differential gene expression
#      (DGE) results (log2 fold change, p-value, adjusted p-value).
#
# Expected input layout (relative to the working directory):
#   data/gsea/<name>_All/GO_merged.csv   columns: source_ID, ontology,
#                                                 core_enrichment ("/"-separated)
#   data/dge/<name>.csv                  columns: gene, log2FoldChange, pvalue,
#                                                 padj
#
# Requires R >= 4.1 (native pipe `|>` and `\(x)` lambdas) and the tidyverse
# (>= 2.0 / dplyr >= 1.1 for `separate_longer_delim()`).
# =============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
})

# ---- Configuration ----------------------------------------------------------

# Directory containing one sub-directory per GSEA result
path_prefix <- file.path("data", "gsea")

# Directory containing one DGE CSV per comparison
dge_path_prefix <- file.path("data", "dge")

# ---- Helpers ----------------------------------------------------------------

#' Stop with a clear message if a file does not exist.
check_file_exists <- function(path) {
    if (!file.exists(path)) {
        stop("Required input file not found: ", path, call. = FALSE)
    }
    invisible(path)
}

#' Stop with a clear message if a data frame lacks required columns.
check_columns <- function(data, required, source) {
    missing_cols <- setdiff(required, colnames(data))
    if (length(missing_cols) > 0) {
        stop(
            "Missing column(s) [", paste(missing_cols, collapse = ", "),
            "] in ", source,
            call. = FALSE
        )
    }
    invisible(data)
}

# ---- Discover GSEA result directories ---------------------------------------

if (!dir.exists(path_prefix)) {
    stop("GSEA directory not found: ", path_prefix,
        " (is the working directory the project root?)",
        call. = FALSE
    )
}

# Sub-directory names, named by the same name with the "_All" suffix removed.
# Values  -> directory names on disk, with "_All" (passed to leading_edge(),
#            and used for both the GSEA folder and the DGE file name)
# Names   -> comparison names without "_All" (only used to label the results)
# Only directories are kept (stray files in data/gsea are ignored).
gsea_files <- list.files(path_prefix) |>
    keep(\(x) dir.exists(file.path(path_prefix, x))) |>
    set_names(\(x) str_remove(x, "_All"))

# ---- Loaders ----------------------------------------------------------------

#' Load a GSEA result and expand leading-edge gene lists.
#'
#' @param filepath Name of the GSEA sub-directory inside `path_prefix`.
#' @return Tibble with one row per (GO term, leading-edge gene):
#'   ID, Ontology, GeneSymbol.
load_gsea_result <- function(filepath) {
    file <- file.path(path_prefix, filepath, "GO_merged.csv")
    check_file_exists(file)

    raw <- read_csv(file, show_col_types = FALSE)
    check_columns(raw, c("source_ID", "ontology", "core_enrichment"), file)

    raw |>
        select(ID = source_ID, Ontology = ontology, LeadingEdge = core_enrichment) |>
        # "GENE1/GENE2/GENE3" -> one row per gene
        separate_longer_delim(LeadingEdge, "/") |>
        rename(GeneSymbol = LeadingEdge)
}

#' Load the differential expression table for one comparison.
#'
#' @param filepath Comparison name (file name without ".csv") in `dge_path_prefix`.
#' @return Tibble: GeneSymbol, LFC, PValue, AdjPValue.
load_dge_result <- function(filepath) {
    file <- file.path(dge_path_prefix, paste0(filepath, ".csv"))
    check_file_exists(file)

    raw <- read_csv(file, show_col_types = FALSE)
    check_columns(raw, c("gene", "log2FoldChange", "pvalue", "padj"), file)

    raw |>
        select(
            GeneSymbol = gene, LFC = log2FoldChange, PValue = pvalue,
            AdjPValue = padj
        )
}

# ---- Transformations --------------------------------------------------------

#' Count leading-edge occurrences per gene, one column per ontology.
#'
#' Genes absent from an ontology get NA in that column (not 0).
format_leadingedge_result <- function(data) {
    data |>
        count(Ontology, GeneSymbol) |>
        pivot_wider(names_from = Ontology, values_from = n)
}

#' Attach DGE statistics to the leading-edge table (keeps all leading-edge
#' genes; genes missing from the DGE table get NA statistics).
enrich_leadingedge_result <- function(data, dge_table) {
    data |>
        left_join(dge_table, by = "GeneSymbol")
}

# ---- Main function ----------------------------------------------------------

#' Build the enriched leading-edge table for one comparison.
#'
#' @param filepath Name of the GSEA sub-directory as on disk (e.g.
#'   "Treatment_2h_All"). The same string is used to locate the DGE file
#'   (`data/dge/<filepath>.csv`), so the "_All" suffix is kept for both.
#' @return Tibble: GeneSymbol, one count column per ontology, LFC, PValue,
#'   AdjPValue.
leading_edge <- function(filepath) {
    gsea_result <- load_gsea_result(filepath)
    dge_table <- load_dge_result(filepath)

    enriched_result <- gsea_result |>
        format_leadingedge_result() |>
        enrich_leadingedge_result(dge_table)

    # Explicit return (the original returned this invisibly via assignment)
    return(enriched_result)
}

# ---- Run --------------------------------------------------------------------

# Only run when executed as a script (Rscript leading_edge.R) or sourced
# interactively; skipped if you set `RUN_LEADING_EDGE <- FALSE` beforehand
# and just want the function definitions.
if (!isFALSE(get0("RUN_LEADING_EDGE", ifnotfound = TRUE))) {
    if (length(gsea_files) == 0) {
        stop("No GSEA result directories found in ", path_prefix, call. = FALSE)
    }

    # Named list of result tables, one per comparison.
    # Each comparison is processed independently; a failure is reported as a
    # warning and the remaining comparisons still run.
    leading_edge_results <- imap(gsea_files, \(dir_name, comparison) {
        tryCatch(
            leading_edge(dir_name),
            error = \(e) {
                warning("Skipping '", comparison, "': ", conditionMessage(e),
                    call. = FALSE
                )
                NULL
            }
        )
    }) |>
        compact() # drop failed (NULL) comparisons

    message(
        "Processed ", length(leading_edge_results), " of ",
        length(gsea_files), " comparison(s)."
    )
}

dir.create(file.path("data", "leading_edge"))


iwalk(leading_edge_results, \(df, name) write_csv(df, file.path("data", "leading_edge", paste0(name, ".csv"))))
