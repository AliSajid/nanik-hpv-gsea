# Generate PAVER Plots for Nanik Data


library(PAVER)
library(tidyverse)

data_files <- "data/combined_gsea_results.csv"

embeddings_url <- "~/experiments/PAVER_embeddings/2025-02-06/anc2vec_02_06_25_embeddings.rds"

term2name_url <- "~/experiments/PAVER_embeddings/2025-02-06/anc2vec_02_06_25_term2name.rds"

embeddings <- readRDS(embeddings_url)
term2name <- readRDS(term2name_url)

paver_data <- read_csv(data_files)

minClusterSize <- 3

maxCoreScatter <- 0.33

minGap <- (1 - maxCoreScatter) * 3 / 4


paver_input <- prepare_data(paver_data, embeddings, term2name)

PAVER_result <- generate_themes(
    paver_input,
    maxCoreScatter = maxCoreScatter,
    minGap = minGap,
    minClusterSize = minClusterSize
)
