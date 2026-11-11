# 02_make_plots.R
# Reads the GSEA output of 01_run_gsea.R (data/gsea/<sheet>/) and creates all
# figures in figures/<sheet>/.

library(dplyr)
library(ggplot2)
library(stringr)
library(clusterProfiler)
library(enrichplot)
library(GSEAmining)

# ---- Config ------------------------------------------------------------------
data_root        <- file.path("data", "gsea")
fig_root         <- "figures"
top_n            <- 20      # per direction
make_enrichplots <- TRUE    # one figure per significant gene set (can be many files)

# ---- Helpers -----------------------------------------------------------------
save_plot <- function(plot, base, width, height) {
  dir.create(dirname(base), recursive = TRUE, showWarnings = FALSE)
  ggsave(paste0(base, ".pdf"), plot, width = width, height = height)
  ggsave(paste0(base, ".png"), plot, width = width, height = height)
}

# Top `top_n` strongest positive and negative NES, ordered by descending NES
make_dotplot <- function(df, title) {
  df <- df %>%
    filter(!is.na(Description), !is.na(NES), !is.na(p.adjust)) %>%
    mutate(Count = str_count(core_enrichment, "/") + 1)   # leading edge size

  plot_df <- bind_rows(
    df %>% filter(NES > 0) %>% slice_max(NES, n = top_n, with_ties = FALSE),
    df %>% filter(NES < 0) %>% slice_min(NES, n = top_n, with_ties = FALSE)
  ) %>%
    arrange(desc(NES)) %>%
    mutate(
      label = ifelse(!is.na(source_ID) & nzchar(source_ID),
                     paste0(Description, " (", source_ID, ")"), Description),
      label = factor(label, levels = rev(unique(label)))   # highest NES on top
    )

  if (nrow(plot_df) == 0) return(NULL)

  p <- ggplot(plot_df, aes(x = NES, y = label)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_continuous(name = "p.adjust", low = "#d73027", high = "#4575b4",
                           trans = "reverse") +
    scale_size_continuous(name = "Count") +
    labs(title = paste0("Top ", top_n, " positive and negative NES - ", title),
         x = "NES", y = NULL) +
    theme_bw() +
    theme(plot.title = element_text(size = 9, face = "bold"),
          axis.text.y = element_text(size = 6),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank())

  list(plot = p, height = max(6, min(18, 4 + nrow(plot_df) * 0.30)))
}

make_gseamining_plots <- function(res, out_dir) {
  filtered <- tryCatch(gm_filter(res), error = function(e) NULL)
  if (is.null(filtered) || nrow(filtered) == 0) {
    message("    GSEAmining: nothing to plot")
    return(invisible(NULL))
  }
  clustered <- tryCatch(gm_clust(filtered), error = function(e) NULL)
  if (is.null(clustered)) {
    message("    GSEAmining: clustering failed")
    return(invisible(NULL))
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  png(file.path(out_dir, "dendrogram_plot.png")); gm_dendplot(filtered, clustered); dev.off()
  pdf(file.path(out_dir, "dendrogram_plot.pdf")); gm_dendplot(filtered, clustered); dev.off()

  for (fun in c("enrichterms", "enrichcores")) {
    p <- tryCatch(
      if (fun == "enrichterms") gm_enrichterms(filtered, clustered)
      else gm_enrichcores(filtered, clustered),
      error = function(e) NULL
    )
    if (inherits(p, "ggplot")) {
      save_plot(p, file.path(out_dir, paste0(fun, "_plot")), 10, 8)
    }
  }
}

make_enrichplots_for <- function(rds_file, out_dir) {
  gsea <- readRDS(rds_file)
  for (id in gsea@result$ID) {
    p <- tryCatch(
      gseaplot(gsea, geneSetID = id, title = id, color.line = "#0d76ff"),
      error = function(e) { message("    gseaplot failed: ", id); NULL }
    )
    if (!is.null(p)) {
      save_plot(p, file.path(out_dir, gsub("[^[:alnum:]_.-]", "_", id)), 10, 8)
    }
  }
}

# ---- Main loop ---------------------------------------------------------------
sheet_dirs <- list.dirs(data_root, recursive = FALSE)
if (length(sheet_dirs) == 0) stop("No result folders found in ", data_root,
                                  " - run 01_run_gsea.R first")

for (sheet_dir in sheet_dirs) {
  sheet <- basename(sheet_dir)
  tryCatch({
    message("Sheet: ", sheet)
    out_dir <- file.path(fig_root, sheet)

    # Dotplots: each ontology + merged
    for (csv in list.files(sheet_dir, pattern = "\\.csv$", full.names = TRUE)) {
      name <- tools::file_path_sans_ext(basename(csv))
      res  <- read.csv(csv, stringsAsFactors = FALSE)

      dp <- make_dotplot(res, name)
      if (!is.null(dp)) {
        save_plot(dp$plot, file.path(out_dir, "dotplots", paste0("dotplot_", name)),
                  width = 10, height = dp$height)
      }

      # GSEAmining on the merged GO table only
      if (grepl("merged", name)) {
        make_gseamining_plots(res, file.path(out_dir, "gseamining"))
      }
    }

    # Enrichplots from the saved GSEA objects
    if (make_enrichplots) {
      for (rds in list.files(sheet_dir, pattern = "\\.rds$", full.names = TRUE)) {
        make_enrichplots_for(
          rds, file.path(out_dir, "enrichplots", tools::file_path_sans_ext(basename(rds)))
        )
      }
    }
  }, error = function(e) {
    message("Plots for ", sheet, " failed: ", conditionMessage(e))
  })
}
