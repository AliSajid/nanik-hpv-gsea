# 01_run_gsea.R
# Reads data/dge/*.csv (written by 00_export_dge.R), runs GSEA against the three
# GO sub-ontologies (BP, CC, MF) and writes results to data/gsea/<dge name>/.

library(clusterProfiler)
library(msigdbr)
library(dplyr)
library(BiocParallel)
register(SerialParam())

# ---- Config ------------------------------------------------------------------
dge_dir <- file.path("data", "dge")
out_root <- file.path("data", "gsea")

gene_col <- "gene"
rank_col <- "log2FoldChange"

# Name used for file names -> msigdbr collection / subcollection.
# Add further entries here to analyse more MSigDB sets.
ontologies <- list(
    GOBP = list(collection = "C5", subcollection = "GO:BP"),
    GOCC = list(collection = "C5", subcollection = "GO:CC"),
    GOMF = list(collection = "C5", subcollection = "GO:MF")
)
merged_name <- "GO_merged"

# ---- Helpers -----------------------------------------------------------------
# Named, finite, de-duplicated, descending log2FC vector
make_ranked_list <- function(df) {
    ranks <- setNames(df[[rank_col]], df[[gene_col]])
    ranks <- ranks[is.finite(ranks) & !duplicated(names(ranks))]
    sort(ranks, decreasing = TRUE)
}

run_gsea <- function(ranks, gene_sets) {
    set.seed(2020)
    tryCatch(
        GSEA(
            geneList      = ranks,
            TERM2GENE     = dplyr::select(gene_sets, gs_name, gene_symbol),
            minGSSize     = 10,
            maxGSSize     = 1000,
            pvalueCutoff  = 1L,
            pAdjustMethod = "BH",
            verbose       = TRUE,
            seed          = TRUE,
            by            = "fgsea"
        ),
        error = function(e) {
            message("    GSEA failed: ", conditionMessage(e))
            NULL
        }
    )
}

# gseaResult -> tidy data frame.
# `leading_edge` (the "tags=..%, list=..%, signal=..%" string) is dropped;
# `core_enrichment` holds the leading edge genes, separated by "/".
format_results <- function(gsea, ontology, sheet, term_ids) {
    as.data.frame(gsea) %>%
        mutate(
            source_ID = term_ids$gs_exact_source[match(ID, term_ids$gs_name)],
            ontology  = ontology,
            sheet     = sheet
        ) %>%
        dplyr::select(-any_of("leading_edge")) %>%
        relocate(ID, source_ID, Description, ontology, sheet) %>%
        `rownames<-`(NULL)
}

# ---- Input files -------------------------------------------------------------
dge_files <- list.files(dge_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(dge_files) == 0) stop("No CSV files in ", dge_dir, " - run 00_export_dge.R first")

# ---- Load gene sets once -----------------------------------------------------
gene_sets <- lapply(names(ontologies), function(name) {
    o <- ontologies[[name]]
    gs <- msigdbr(
        species = "Homo sapiens",
        collection = o$collection, subcollection = o$subcollection
    )
    if (nrow(gs) == 0) stop("No gene sets found for ", name)
    gs
})
names(gene_sets) <- names(ontologies)

term_ids <- lapply(gene_sets, function(gs) {
    distinct(dplyr::select(gs, gs_name, gs_exact_source), gs_name, .keep_all = TRUE)
})

# ---- Main loop ---------------------------------------------------------------
failed <- character()

for (file in dge_files) {
    sheet <- tools::file_path_sans_ext(basename(file))
    tryCatch(
        {
            message("Sheet: ", sheet)

            df <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
            missing_cols <- setdiff(c(gene_col, rank_col), colnames(df))
            if (length(missing_cols) > 0) {
                stop("Missing column(s): ", paste(missing_cols, collapse = ", "))
            }

            ranks <- make_ranked_list(df)

            out_dir <- file.path(out_root, sheet)
            dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

            results <- list()
            for (ont in names(ontologies)) {
                message("  ", ont)
                gsea <- run_gsea(ranks, gene_sets[[ont]])

                if (is.null(gsea) || nrow(gsea@result) == 0) {
                    message("    no significant gene sets")
                    next
                }

                results[[ont]] <- format_results(gsea, ont, sheet, term_ids[[ont]])
                write.csv(results[[ont]], file.path(out_dir, paste0(ont, ".csv")), row.names = FALSE)
                saveRDS(gsea, file.path(out_dir, paste0(ont, ".rds"))) # needed for enrichplots
            }

            if (length(results) > 0) {
                write.csv(bind_rows(results),
                    file.path(out_dir, paste0(merged_name, ".csv")),
                    row.names = FALSE
                )
            }

            message("  -> saved to ", out_dir)
        },
        error = function(e) {
            message("Sheet ", sheet, " failed: ", conditionMessage(e))
            failed <<- c(failed, sheet)
        }
    )
}

if (length(failed) > 0) {
    warning("Failed sheets: ", paste(failed, collapse = ", "))
}
