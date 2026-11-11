# =============================================================================
# Leading-edge recurrence vs. log2 fold change
# -----------------------------------------------------------------------------
# For every comparison in `leading_edge_results` (produced by leading_edge.R),
# draws a scatter plot:
#   x = log2 fold change (LFC)
#   y = number of leading-edge GO terms the gene appears in (all ontologies
#       summed)
#   colour = significant (adjusted p-value < threshold) or not
# and labels the most recurrent genes. One PNG is saved per comparison.
#
# Genes in the upper left/right corners are strongly regulated AND drive many
# enriched GO terms, making them good candidates for follow-up.
#
# Usage:  Rscript plot_leading_edge.R      (from the project root, next to
#                                           leading_edge.R and data/)
# =============================================================================

suppressPackageStartupMessages({
    library(tidyverse)
})

# ggrepel is needed for non-overlapping labels
if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required. Install it with: ",
        "install.packages(\"ggrepel\")",
        call. = FALSE
    )
}

# ---- Configuration ----------------------------------------------------------

leading_edge_script <- "R/05_leading_edge_analysis.R" # script that builds the result tables
output_dir <- "figures" # PNGs are written here
padj_threshold <- 0.05 # significance cut-off (AdjPValue)
n_labels <- 15 # top recurrent genes to label per plot
plot_width <- 8 # inches
plot_height <- 6 # inches
plot_dpi <- 300

# Columns in each result table that are NOT ontology counts
stat_columns <- c("GeneSymbol", "LFC", "PValue", "AdjPValue")

# ---- Load leading-edge results ----------------------------------------------

# Reuse an existing `leading_edge_results` object if the session already has
# one; otherwise run the pipeline script to create it.
if (!exists("leading_edge_results")) {
    if (!file.exists(leading_edge_script)) {
        stop("Cannot find '", leading_edge_script,
            "'. Run this script from the project root.",
            call. = FALSE
        )
    }
    source(leading_edge_script)
}

if (!is.list(leading_edge_results) || length(leading_edge_results) == 0) {
    stop("`leading_edge_results` is empty; nothing to plot.", call. = FALSE)
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Data preparation -------------------------------------------------------

#' Add TotalCount and Significant columns to one comparison's table.
#'
#' Ontology columns are detected automatically (everything that is not a
#' statistics column), so this works whichever ontologies are present.
#' NA counts mean "gene not in this ontology" and count as 0.
prepare_plot_data <- function(data) {
    ontology_cols <- setdiff(colnames(data), stat_columns)

    data |>
        mutate(
            TotalCount = rowSums(pick(all_of(ontology_cols)), na.rm = TRUE),
            Significance = if_else(
                !is.na(AdjPValue) & AdjPValue < padj_threshold,
                paste0("Adj. p < ", padj_threshold),
                "Not significant"
            ) |>
                factor(levels = c(
                    paste0("Adj. p < ", padj_threshold),
                    "Not significant"
                ))
        ) |>
        # Genes with no DGE match (NA LFC) cannot be placed on the x axis
        filter(!is.na(LFC))
}

# ---- Plotting ---------------------------------------------------------------

#' Draw the recurrence-vs-LFC scatter for one comparison.
plot_leading_edge <- function(data, comparison) {
    plot_data <- prepare_plot_data(data)

    if (nrow(plot_data) == 0) {
        warning("Skipping '", comparison, "': no genes with an LFC value.",
            call. = FALSE
        )
        return(NULL)
    }

    label_data <- slice_max(plot_data, TotalCount,
        n = n_labels,
        with_ties = FALSE
    )

    ggplot(plot_data, aes(x = LFC, y = TotalCount, colour = Significance)) +
        geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
        geom_point(alpha = 0.6, size = 1.8) +
        ggrepel::geom_text_repel(
            data = label_data, aes(label = GeneSymbol),
            max.overlaps = Inf, size = 3, show.legend = FALSE
        ) +
        scale_colour_manual(
            values = setNames(c("#D55E00", "grey55"), levels(plot_data$Significance)),
            drop = FALSE
        ) +
        labs(
            title = comparison,
            subtitle = paste0("Top ", n_labels, " most recurrent genes labelled"),
            x = "log2 fold change",
            y = "Leading-edge GO terms (all ontologies)",
            colour = NULL
        ) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom")
}

# ---- Render and save one figure per comparison ------------------------------

iwalk(leading_edge_results, \(data, comparison) {
    p <- tryCatch(
        plot_leading_edge(data, comparison),
        error = \(e) {
            warning("Plot failed for '", comparison, "': ",
                conditionMessage(e),
                call. = FALSE
            )
            NULL
        }
    )
    if (is.null(p)) {
        return(invisible(NULL))
    }

    # Make the comparison name safe for use as a file name
    file_name <- paste0("leading_edge_", str_replace_all(comparison, "[^A-Za-z0-9_-]", "_"), ".png")
    out_file <- file.path(output_dir, file_name)

    ggsave(out_file,
        plot = p, width = plot_width, height = plot_height,
        dpi = plot_dpi, bg = "white"
    )
    message("Saved ", out_file)
})
