# Script to load individual GSEA results and combine them for PAVER

library(tidyverse)

path_prefix <- file.path("data", "gsea")

gsea_files <- list.files(path_prefix) |>
    set_names(\(x) str_remove(x, "_All"))

load_gsea_result <- function(filepath) {
    read_csv(file.path(path_prefix, filepath, "GO_merged.csv")) |>
        select(ID = source_ID, Score = NES) |>
        rename({{ filepath }} := Score)
}

combined_dataset <- gsea_files |>
    map(load_gsea_result) |>
    reduce(\(df_a, df_b) merge(df_a, df_b, by = "ID"))

combined_dataset |>
    write_csv("data/combined_gsea_results.csv")
