# Given the CED perturbations, retrieve the Stat3 scores that we need to overlay
# on the barplot

fgsea_output_list <- readRDS("output/gsea/condNormalized/CED/Combined/fgsea_ensembl_ontology_list_nes.rds")

stat3_scores <- list()
for (name in names(fgsea_output_list)) {
    df <- fgsea_output_list[[name]]
    df <- subset(df, pathway == "HALLMARK_IL6_JAK_STAT3_SIGNALING")$NES
    stat3_scores[[name]] <- df
}
stat3_df <- data.frame(
    perturb = names(stat3_scores),
    NES = unlist(stat3_scores)
)
rownames(stat3_df) <- NULL
stat3_df$perturb <- gsub("_non-targeting_noRT|_non-targeting_RT", "", 
                         stat3_df$perturb)
stat3_df$perturb <- gsub("preinf_", "", stat3_df$perturb)
write.csv(stat3_df, "data/chromatin_analysis/preinf_stat3_scores.csv", 
          row.names = FALSE, quote = FALSE)
