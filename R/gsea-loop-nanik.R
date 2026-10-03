# --- GSEA_LOOP_HPV_coding (Excel‑Version für mehrere Sheets) ---
# alle MSIG Kategorien& Subkategorien
# geändert zu je top 20 pathways nach NES (stärkste Effekte in beide Richtungen)
# NES auch x-Achse
# 1 Graph
# Sortierung der y-Achse: absteigende NES

# install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "enrichplot", "DOSE", "GSEABase", "GSEAmining"))
# install.packages(c("readxl", "dplyr", "msigdbr", "ggplot2", "openxlsx"))
# if (!require("BiocManager", quietly = TRUE)) {
# install.packages("BiocManager")
# }
# BiocManager::install("GSEABase")

library(readxl)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(org.Hs.eg.db)
library(msigdbr)
library(GSEABase)
library(GSEAmining)    
library(ggplot2)
library(openxlsx)
library(stringr)
library(patchwork)
library(BiocParallel)
register(SerialParam())

setwd("C:/Users/Nanik Wahrhausen/Documents/R Projekte")

# 1. Excel‑Pfad und Sheets
path <- "C:/Users/Nanik Wahrhausen/Documents/R Projekte/04.1_DESeq2_coding_total_optimized/20251204_analysis_files/20251204_coding_total_results.xlsx"

sheet_names <- c(
  "HPV_Main_Effect_All",
  "coding_total_Untreated_All",
  "Treatment_2h_All",
  "Treatment_6h_All",
  "Treatment_24h_All",
  "HPV_Interaction_2h_All",
  "HPV_Interaction_6h_All",
  "HPV_Interaction_24h_All"
)

# Prüfen, ob alle Sheets da sind
available_sheets <- excel_sheets(path)
missing_sheets <- setdiff(sheet_names, available_sheets)

if (length(missing_sheets) > 0) {
  stop(paste("Diese Sheets fehlen in der Excel‑Datei:", paste(missing_sheets, collapse = ", ")))
}

# 2. Gene und Rank spalte
#  Angenommen: erste Spalte hat keinen Namen -> wird von readxl zu "...1"
gene_col  <- "...1"
rank_col  <- "log2FoldChange"

# 3. gewünschte MSigDB‑Collections + Subkategorien für Mensch (Homo sapiens)
target_collections <- list(
  H   = "H",
  C1  = "C1",
  C2  = "C2",
  C3  = "C3",
  C4  = "C4",
  C5  = "C5",
  C6  = "C6",
  C7  = "C7",
  C8  = "C8",
  
  C2_REACTOME      = "C2",
  C2_BIOCARTA      = "C2",
  C2_KEGG_LEGACY   = "C2",
  C3_MIRDB         = "C3",
  C5_GOBP          = "C5",
  C7_IMMUNESIGDB   = "C7"
)

# 4. Funktion um GSEA pro Collection / Sheet durchzuführen
perform_gsea <- function(category_df, input_data, gene_col, rank_col) {
  lfc_vector <- input_data[[rank_col]]
  names(lfc_vector) <- input_data[[gene_col]]
  lfc_vector <- sort(lfc_vector, decreasing = TRUE)
  
  # NA/Inf entfernen
  lfc_vector <- lfc_vector[!is.na(lfc_vector)]
  lfc_vector <- lfc_vector[is.finite(lfc_vector)]
  
  # doppelte Namen entfernen
  lfc_vector <- lfc_vector[!duplicated(names(lfc_vector))]
  
  set.seed(2020)
  
  tryCatch({
    clusterProfiler::GSEA(
      geneList = lfc_vector,
      TERM2GENE = dplyr::select(category_df, gs_name, gene_symbol),
      minGSSize = 10,
      maxGSSize = 1000,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose = FALSE,
      seed = TRUE,
      by = "fgsea"
    )
  }, error = function(e) {
    message("GSEA failed: ", conditionMessage(e))
    return(NULL)
  })
}

# 5. Name des Resultats‑Ordners aus Sheet‑Name erzeugen
clean_sheet_name <- function(sheet) {
  x <- tolower(sheet)
  x <- gsub("_all$", "", x)
  x <- gsub("_main_effect", "", x)
  x <- gsub("^coding_total_", "", x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

# 6. FUNKTION: merged GSEA‑results pro Sheet speichern
save_merged_results <- function(all_results, out_dir, sheet_name) {
  merged_results <- data.frame()
  excel_list <- list()
  
  for (category_name in names(all_results)) {
    result <- all_results[[category_name]]
    
    if (!is.null(result) && nrow(as.data.frame(result)) > 0) {
      result_df <- as.data.frame(result)
      result_df$Category <- category_name
      result_df$Sheet <- sheet_name
      
      merged_results <- rbind(merged_results, result_df)
      excel_list[[category_name]] <- result_df
    }
  }
  
  write.csv(merged_results, file = file.path(out_dir, "merged_GSEA_results.csv"), row.names = FALSE)
  write.xlsx(excel_list, file = file.path(out_dir, "merged_GSEA_results.xlsx"), overwrite = TRUE)
  saveRDS(merged_results, file = file.path(out_dir, "merged_GSEA_results.rds"))
  
  merged_results
}

# 7. GSEAmining‑Plots (dendrogram, enriched terms, cores); optional auskommentierbar
generate_gseamining_plots <- function(merged_results, input_file_directory) {
  
  # Falls GSEAmining nicht geladen oder fehlerhaft
  if (!exists("gm_filter") || is.null(merged_results) || nrow(merged_results) == 0) {
    message("GSEAmining‑Plots übersprungen für: ", input_file_directory)
    return(NULL)
  }
  
  gs_filtered <- tryCatch(gm_filter(merged_results), error = function(e) NULL)
  gs_clustered <- tryCatch(gm_clust(gs_filtered), error = function(e) NULL)
  
  if (is.null(gs_filtered) || is.null(gs_clustered) || nrow(gs_filtered) == 0) {
    message("Keine Daten für GSEAmining‑Plots: ", input_file_directory)
    return(NULL)
  }
  
  # dendrogram
  png(filename = file.path(input_file_directory, "dendrogram_plot.png"))
  gm_dendplot(gs_filtered, gs_clustered)
  dev.off()
  
  pdf(file = file.path(input_file_directory, "dendrogram_plot.pdf"))
  gm_dendplot(gs_filtered, gs_clustered)
  dev.off()
  
  # enriched terms
  enriched_terms_plot <- tryCatch(gm_enrichterms(gs_filtered, gs_clustered), error = function(e) NULL)
  if (!is.null(enriched_terms_plot) && any(c("gg", "ggplot") %in% class(enriched_terms_plot))) {
    ggsave(filename = file.path(input_file_directory, "enriched_terms_plot.png"), plot = enriched_terms_plot, width = 10, height = 8)
    ggsave(filename = file.path(input_file_directory, "enriched_terms_plot.pdf"), plot = enriched_terms_plot, width = 10, height = 8)
  }
  
  # enrichment cores
  enrichment_cores_plot <- tryCatch(gm_enrichcores(gs_filtered, gs_clustered), error = function(e) NULL)
  if (!is.null(enrichment_cores_plot) && any(c("gg", "ggplot") %in% class(enrichment_cores_plot))) {
    ggsave(filename = file.path(input_file_directory, "enrichment_cores_plot.png"), plot = enrichment_cores_plot, width = 10, height = 8)
    ggsave(filename = file.path(input_file_directory, "enrichment_cores_plot.pdf"), plot = enrichment_cores_plot, width = 10, height = 8)
  }
}

# 8. FUNKTION: Plots pro Gene Set / Collection speichern
save_plots <- function(plots_list, base_directory, category_name, prefix = NULL, width = 10, height = 8) {
  dir.create(base_directory, showWarnings = FALSE, recursive = TRUE)
  
  safe_category <- gsub("[^[:alnum:]_\\.-]", "_", category_name)
  
  for (plot_name in names(plots_list)) {
    plot <- plots_list[[plot_name]]
    
    if (is.null(prefix)) {
      file_name <- paste0(plot_name, "_", safe_category)
    } else {
      file_name <- paste0(prefix, "_", safe_category)
    }
    
    file_base <- file.path(base_directory, file_name)
    
    ggsave(filename = paste0(file_base, ".pdf"), plot = plot, width = width, height = height)
    ggsave(filename = paste0(file_base, ".png"), plot = plot, width = width, height = height)
  }
}

# 9. Hilfsfunktion für automatische Graph-Höhe
get_dotplot_height <- function(result_df, max_categories = 20, split_sign = TRUE) {
  if (is.null(result_df) || nrow(result_df) == 0) {
    return(6)
  }
  
  n_terms <- min(nrow(result_df), max_categories)
  
  if (split_sign && ".sign" %in% colnames(result_df)) {
    n_pos <- sum(result_df$.sign == "activated", na.rm = TRUE)
    n_neg <- sum(result_df$.sign == "suppressed", na.rm = TRUE)
    
    n_pos <- min(n_pos, max_categories)
    n_neg <- min(n_neg, max_categories)
    
    tallest_panel <- max(n_pos, n_neg, 1)
    height <- 4 + tallest_panel * 0.35
  } else {
    height <- 4 + n_terms * 0.35
  }
  
  height <- max(5, min(height, 16))
  return(height)
}

# 10. MAIN LOOP über alle Sheets
today <- format(Sys.Date(), "%Y%m%d")

for (sheet in sheet_names) {
  tryCatch({
  message("Verarbeite Sheet: ", sheet)
  
  # Excel einlesen
  df <- read_excel(path, sheet = sheet) %>% as.data.frame()
  
  # Check gen_spalte und rank_spalte
  if (!(gene_col %in% colnames(df))) {
    stop(paste("Spalte fehlt in Sheet '", sheet, "': ", gene_col))
  }
  if (!(rank_col %in% colnames(df))) {
    stop(paste("Spalte fehlt in Sheet '", sheet, "': ", rank_col))
  }
  
  df <- df %>%
    filter(!is.na(.data[[gene_col]]), !is.na(.data[[rank_col]])) %>%
    distinct(.data[[gene_col]], .keep_all = TRUE)
  
  # Ordner pro Sheet
  sheet_clean <- clean_sheet_name(sheet)
  results_dir <- paste0(today, "_GSEA_coding_", sheet_clean)
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 2. Ergebnis‑Liste pro Collection
  all_results_list <- list()
  
  # Alle Ziel‑Collections durchgehen
  for (target_name in names(target_collections)) {
    main_collection <- target_collections[[target_name]]
    
    message("  Starte Collection: ", target_name, " (main collection: ", main_collection, ")")
    
    # Lade Gensets
    msig_df <- msigdbr(species = "Homo sapiens", collection = main_collection)
    
    # Spezielle Subkategorien filtern
    if (target_name == "C2_REACTOME") {
      msig_df <- msig_df %>% filter(gs_subcollection == "CP:REACTOME")
    } else if (target_name == "C2_BIOCARTA") {
      msig_df <- msig_df %>% filter(gs_subcollection == "CP:BIOCARTA")
    } else if (target_name == "C2_KEGG_LEGACY") {
      msig_df <- msig_df %>% filter(gs_subcollection == "CP:KEGG_LEGACY")
    } else if (target_name == "C5_GOBP") {
      msig_df <- msig_df %>% filter(gs_subcollection == "GO:BP")
    } else if (target_name == "C3_MIRDB") {
      msig_df <- msig_df %>% filter(gs_subcollection == "MIR:MIRDB")
    } else if (target_name == "C7_IMMUNESIGDB") {
      msig_df <- msig_df %>% filter(gs_subcollection == "IMMUNESIGDB")
    }
    
    if (is.null(msig_df) || nrow(msig_df) == 0) {
      message("Keine Gensets gefunden für: ", target_name)
      next
    }
    
    # GSEA durchführen
    gsea_result <- perform_gsea(msig_df, df, gene_col, rank_col)
    all_results_list[[target_name]] <- gsea_result
  }
  
  # 3. Ergebnisse pro Collection speichern (gsea_results_<Category>_<Sheet>.csv)
  for (category_name in names(all_results_list)) {
    result <- all_results_list[[category_name]]
    if (!is.null(result) && nrow(as.data.frame(result)) > 0) {
      output_file <- file.path(
        results_dir,
        paste0("gsea_results_", category_name, "_", sheet_clean, ".csv")
      )
      write.csv(as.data.frame(result), file = output_file, row.names = FALSE)
    }
  }
  
  # 4. merged results pro Sheet
  merged_results <- save_merged_results(all_results_list, results_dir, sheet)
  
  # 5. GSEAmining‑Plots, falls Paket vorhanden
  generate_gseamining_plots(merged_results, results_dir)
  
  # 6. Unterordner für enrichplots und doseplots
  enrichplots_dir <- file.path(results_dir, "enrichplots")
  dotplots_dir  <- file.path(results_dir, "dotplots")
  dir.create(enrichplots_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(dotplots_dir,  showWarnings = FALSE, recursive = TRUE)
  
  # 7. Loop über alle Collections, um Enrichplots und Dotplots zu erstellen
  for (category_name in names(all_results_list)) {
    result <- all_results_list[[category_name]]
    if (is.null(result) || nrow(as.data.frame(result)) == 0L) {
      next
    }
    
    # Enrichplots pro Gene Set
    geneSetIDs <- as.data.frame(result)$ID
    enrich_plot_list <- list()
    
    for (geneSetID in geneSetIDs) {
      enrich_plot <- tryCatch({
        gseaplot(
          result,
          geneSetID = geneSetID,
          title = geneSetID,
          color.line = "#0d76ff"
        )
      }, error = function(e) {
        message("gseaplot fehlgeschlagen für: ", geneSetID)
        NULL
      })
      
      if (!is.null(enrich_plot)) {
        safe_name <- gsub("[^[:alnum:]_\\.]|_", "_", geneSetID)
        safe_name <- paste0(substr(safe_name, 1, 80), "_", seq_along(geneSetIDs)
                            [which(geneSetIDs == geneSetID)][1])
        enrich_plot_list[[safe_name]] <- enrich_plot
      }
    }
    
    if (length(enrich_plot_list) > 0) {
      save_plots(
        enrich_plot_list,
        enrichplots_dir,
        category_name,
        prefix = "enrichplot",
        width = 10,
        height = 8
      )
    }
    
    # Dotplot
    r_df <- as.data.frame(result)
    
    if (!"Count" %in% colnames(r_df)) {
      if ("core_enrichment" %in% colnames(r_df)) {
        r_df$Count <- stringr::str_count(r_df$core_enrichment, "/") + 1
      } else {
        r_df$Count <- NA_real_
      }
    }
    
    # Nur sinnvolle Einträge behalten
    r_df <- r_df %>%
      filter(!is.na(Description), !is.na(NES), !is.na(p.adjust))
    
    # Top 20 positive NES
    pos_df <- r_df %>%
      filter(NES > 0) %>%
      arrange(desc(NES), p.adjust) %>%
      slice_head(n = 20)
    
    # Top 20 negative NES
    neg_df <- r_df %>%
      filter(NES < 0) %>%
      arrange(desc(NES), p.adjust) %>%   # z.B. -1 vor -1.5 vor -2.3
      slice_head(n = 20)
    
    # Erst positive, dann negative
    plot_df <- bind_rows(pos_df, neg_df)
    
    if (nrow(plot_df) == 0) {
      next
    }
    
    # Y-Achse: oben nach unten genau wie in plot_df
    plot_df <- plot_df %>%
      mutate(Description = factor(Description, levels = rev(unique(Description))))
    
    gsea_dotplot <- ggplot(plot_df, aes(x = NES, y = Description)) +
      geom_point(aes(size = Count, color = p.adjust)) +
      scale_color_continuous(
        name = "p.adjust",
        low = "#d73027",
        high = "#4575b4",
        trans = "reverse"
      ) +
      scale_size_continuous(name = "Count") +
      labs(
        title = paste0("Top 20 positive and top 20 negative NES pathways - ", category_name),
        x = "NES",
        y = NULL
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(size = 9, face = "bold"),
        axis.text.y = element_text(size = 6),   #size 7 wäre eventuell besser lesbar
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    plot_height <- max(6, min(18, 4 + nrow(plot_df) * 0.30))
    
    gsea_dotplot_list <- list(gsea_dotplot = gsea_dotplot)
    
    save_plots(
      gsea_dotplot_list,
      dotplots_dir,
      category_name,
      prefix = "dotplot",
      width = 10,
      height = plot_height
    )
  }
  
  message("--> Bearbeitung von ", sheet, " ist komplett und wurde in Ordner ", results_dir, " gespeichert.")
  
  }, error = function(e) {
    message("Bearbeitung von ", sheet, " ist fehlgeschlagen: ", conditionMessage(e))
  })
}
