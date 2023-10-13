## ----------------------------------------------------------------------------------------------------------------------------------------
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(knitr)
library(msigdbr)


## ----------------------------------------------------------------------------------------------------------------------------------------
INPUT_DIRS <- c("/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_clean_outputs/deseq/noRTNormalized")
OUTPUT_DIR <- "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_explore_outputs/gsea"
PATH_TO_SEURAT_OBJECT <- "/raleighlab/data1/liuj/gbm_perturb/analysis/GBM43_1_malignant_only_annotated_20230817.Rds"
WORKING_DIR <- "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_explore_outputs/gsea/noRTNormalized"
setwd(WORKING_DIR)


## ----------------------------------------------------------------------------------------------------------------------------------------
data.list <- list()
for (dir in INPUT_DIRS) {
  file_names <- list.files(dir)
  for (i in file_names) {
    name <- gsub(".csv", "", i)
    df <- read.table(paste(dir, "/", i, sep = ""))
    data.list[[name]] <- df
  }
}


## ----------------------------------------------------------------------------------------------------------------------------------------
ranked_gene_set_list <- list()
for (name in names(data.list)) {
  deseq_table <- data.list[[name]]
  deseq_table <- deseq_table[!is.na(deseq_table$log_fc), ]
  deseq_table <- deseq_table[order(-deseq_table$log_fc), ]
  gene_list <- deseq_table$log_fc
  gene_names <- deseq_table$feature
  gene_names <- gsub("GRCh38-", "", gene_names)
  names(gene_list) <- gene_names
  ranked_gene_set_list[[name]] <- gene_list
}


## ----------------------------------------------------------------------------------------------------------------------------------------
bp.ontology.list <- list()
for (name in names(ranked_gene_set_list)) {
  cat("Processing", name)
  gse <- gseGO(
    ranked_gene_set_list[[name]],
    ont = "BP",
    keyType = "SYMBOL",
    OrgDb = org.Hs.eg.db,
    eps = 1e-300
  )
  bp.ontology.list[[name]] <- gse
}
saveRDS(bp.ontology.list, file = paste(WORKING_DIR, "bp_ontology_list.rds", sep = "/"))


## ----------------------------------------------------------------------------------------------------------------------------------------
hallmark_term_to_gene <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)


## ----------------------------------------------------------------------------------------------------------------------------------------
hallmark.msigdb.ontology.list <- list()
for (name in names(ranked_gene_set_list)) {
  gsea <- GSEA(ranked_gene_set_list[[name]], 
               TERM2GENE = hallmark_term_to_gene)
  hallmark.msigdb.ontology.list[[name]] <- gsea
}
saveRDS(hallmark.msigdb.ontology.list, file = paste(WORKING_DIR, "hallmark_ontology_list.rds", sep = "/"))


## ----------------------------------------------------------------------------------------------------------------------------------------
bp.ontology.list <- readRDS(paste(WORKING_DIR, "bp_ontology_list.rds", sep = "/"))


## ----------------------------------------------------------------------------------------------------------------------------------------


## ----------------------------------------------------------------------------------------------------------------------------------------

