# Comparison of downsampled analysis with the regular analysis stream to see:
# 1. Whether or not there are meaningful differences in DESeq output with
#    downsampling.
# 2. Whether projecting the log fold changes from downsampled DESeq onto a
#    heatmap made with non-downsampled gene modules qualitatively changes the
#    heatmap.
# 
# Author: Christopher Zou

##############################################################################
##############################################################################
# Dependencies:

library(Seurat)

##############################################################################
##############################################################################
# Inputs:

SEED = 325
set.seed(SEED)
DESEQ_INPUT_WD = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/deseq/"
OUTPUT_WD = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/downsampled_correlations"

cat("SEED: ", SEED, "\n")
cat("DESeq WORKING DIRECTORY: ", DESEQ_WORKING_DIRECTORY, "\n")

##############################################################################
##############################################################################
# Code:

##############################################################################
# Generate Pearson correlation coefficients per perturbation for each context

# Import all of the data

setwd(DESEQ_INPUT_WD)

downsampled_input_dirs = c("GL261_integrated_20230819_ced_noRTNormalized_downsampled_all",
                           "GL261_integrated_20230819_invitro_noRTNormalized_downsampled_all",
                           "GL261_integrated_20230819_preinf_noRTNormalized_downsampled_sorted")
data_list_downsampled = list()
for (dir in downsampled_input_dirs) {
  file_names = list.files(dir)
  for (i in file_names) {
    name = gsub(".csv", "", i)
    df = read.table(paste(dir, "/", i, sep = ""))
    data_list_downsampled[[name]] = df
  }
}

regular_input_dirs = c("GL261_integrated_20230705_ced_noRTNormalized_all",
                       "GL261_integrated_20230705_preinf_noRTNormalized_sorted",
                       "GL261_integrated_20230705_invitro_noRTNormalized_all")
data_list_regular = list()
for (dir in regular_input_dirs) {
  file_names = list.files(dir)
  for (i in file_names) {
    name = gsub(".csv", "", i)
    df = read.table(paste(dir, "/", i, sep = ""))
    data_list_regular[[name]] = df
  }
}

perturbs_downsampled = names(data_list_downsampled)
perturbs_regular = names(data_list_regular)
perturbs_shared = intersect(perturbs_downsampled, perturbs_regular)

# Find correlation between all genes

pearson_cor_list = list()
spearman_cor_list = list()
for (perturb in perturbs_shared) {
  print(paste("Analyzing correlation for: ", perturb))
  data_downsampled = data_list_downsampled[[perturb]]
  data_regular = data_list_regular[[perturb]]
  merged = merge(data_regular, data_downsampled, by = "feature")
  merged = merged[!is.na(merged$log_fc.x),]
  merged = merged[!is.na(merged$log_fc.y),]
  pearson_cor_list[[perturb]] = cor(merged$log_fc.x, merged$log_fc.y)
  spearman_cor_list[[perturb]] = cor(merged$log_fc.x, merged$log_fc.y, method = "spearman")
}

cor_df = data.frame(Pearson_cor = unlist(pearson_cor_list), 
                    Spearman_cor = unlist(spearman_cor_list))
cor_df$context = sapply(rownames(cor_df), function(x) (strsplit(x, split = "_")[[1]][1]))

average_invitro_pearson_corr <- mean(cor_df[cor_df$context == "invitro", "Pearson_cor"])
average_preinf_pearson_corr <- mean(cor_df[cor_df$context == "preinf", "Pearson_cor"])
average_invitro_spearman_corr <- mean(cor_df[cor_df$context == "invitro", "Spearman_cor"])
average_preinf_spearman_corr <- mean(cor_df[cor_df$context == "preinf", "Spearman_cor"])

# Write the outputs

setwd(OUTPUT_WD)

write.table(cor_df, file = "3Context_pearson_correlations_origDEGenes_RT.csv", sep = "\t")
png("3Context_pearson_correlations_origDEGenes_RT_hist.png")
hist(cor_df$Pearson_cor, main = "Correlations, orig DE genes RT", xlab = "perturbs", ylab = "Pearson correlation coefficient")
dev.off()
png("3Context_spearman_correlations_origDEGenes_RT_hist.png")
hist(cor_df$Spearman_cor, main = "Correlations, orig DE genes RT", xlab = "perturbs", ylab = "Spearman correlation coefficient")
dev.off()

# What if we only look at the DE genes?

# Split up the noRT and RT perturbations, and look at each only with their
# previously determined DE genes. We can also generate DE genes based off
# the downsampled and look at those.

# Get noRT previous DE genes

PATH_TO_DE_GENES = "/raleighlab/data1/liuj/gbm_perturb/analysis/gbm_pdx_perturb_GL261_integrate_xfp/GL261_integrated_20230626_noRT_all/GL261_integrated_20230626_noRT_all_deMatSig_padj05_abs0.1.txt"
deMatSig = read.table(PATH_TO_DE_GENES, sep = "\t", header = TRUE)
deMatSig = na.omit(deMatSig)
rownames(deMatSig) = deMatSig$X
deMatSig = deMatSig[,-1]
de_genes_noRT_previous = rownames(deMatSig)

# Get RT previous DE genes

PATH_TO_DE_GENES = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/de_genes/GL261_integrated_20230729_3Context_noRTNormalized_RTOnly/deMatSig_adjp05_lfc01.txt"
deMatSig = read.table(PATH_TO_DE_GENES, header = TRUE)
deMatSig = na.omit(deMatSig)
de_genes_RT_previous = rownames(deMatSig)

# Screen to only these ad see what correlations look like

perturbs_noRT = perturbs_shared[grep("noRT_non-targeting_noRT", perturbs_shared)]
perturbs_RT = perturbs_shared[-grep("noRT_non-targeting_noRT", perturbs_shared)]
pearson_cor_list = list()
spearman_cor_list = list()
for (perturb in perturbs_RT) {
  print(paste("Analyzing correlation for: ", perturb))
  data_downsampled = data_list_downsampled[[perturb]]
  data_regular = data_list_regular[[perturb]]
  merged = merge(data_regular, data_downsampled, by = "feature")
  merged = merged[!is.na(merged$log_fc.x),]
  merged = merged[!is.na(merged$log_fc.y),]
  merged = merged[merged$feature %in% de_genes_RT_previous, ]
  pearson_cor_list[[perturb]] = cor(merged$log_fc.x, merged$log_fc.y)
  spearman_cor_list[[perturb]] = cor(merged$log_fc.x, merged$log_fc.y, method = "spearman")
}

cor_df = data.frame(Pearson_cor = unlist(pearson_cor_list), 
                    Spearman_cor = unlist(spearman_cor_list))
cor_df$context = sapply(rownames(cor_df), function(x) (strsplit(x, split = "_")[[1]][1]))

average_invitro_pearson_corr <- mean(cor_df[cor_df$context == "invitro", "Pearson_cor"])
average_preinf_pearson_corr <- mean(cor_df[cor_df$context == "preinf", "Pearson_cor"])
average_invitro_spearman_corr <- mean(cor_df[cor_df$context == "invitro", "Spearman_cor"])
average_preinf_spearman_corr <- mean(cor_df[cor_df$context == "preinf", "Spearman_cor"])

# Now take a look at generated differentially expressed genes from the new DESeq
# run.

# Take a look at the tables that we generate for averaged out LFC values and
# look at correlation per module

downsampled_modules = read.table("output_matrix_lognormalized_lfc_rank25_MaxModules_downsampled.txt",
                                 sep = "\t")
full_modules = read.table("output_matrix_lognormalized_lfc_rank25_MaxModules_full.txt",
                                 sep = "\t")
pearson_correlations <- numeric(length(ncol(downsampled_modules)))
spearman_correlations <- numeric(length(ncol(downsampled_modules)))
for (i in 1:ncol(downsampled_modules)) {
  pearson_correlations[i] <- cor(downsampled_modules[[i]], full_modules[[i]], method="pearson")
  spearman_correlations[i] <- cor(downsampled_modules[[i]], full_modules[[i]], method="spearman")
}
annotations = read.table("RT_annotation_table.txt", sep = "\t")
annotations$col_name <- rownames(annotations)
result <- data.frame(
  col_name = colnames(downsampled_modules),
  Pearson_Correlation = pearson_correlations,
  Spearman_Correlation = spearman_correlations
)
final_result <- merge(result, annotations, by="col_name")
final_result <- final_result %>%
  group_by(sgRNA) %>%
  mutate(coverage = min(coverage)) %>%
  ungroup()

write.table(final_result, "by_perturb_correlations.txt", sep = "\t")

# HORIZONTAL (PER MODULE CORRELATIONS)

downsampled_modules_t <- data.frame(t(downsampled_modules))
full_modules_t <- data.frame(t(full_modules))

# Initialize vectors to store correlations
pearson_correlations <- numeric(length(ncol(downsampled_modules_t)))
spearman_correlations <- numeric(length(ncol(full_modules_t)))

# Iterate over rows and calculate correlations
for (i in 1:ncol(downsampled_modules_t)) {
  pearson_correlations[i] <- cor(downsampled_modules_t[[i]], full_modules_t[[i]], method="pearson")
  spearman_correlations[i] <- cor(downsampled_modules_t[[i]], full_modules_t[[i]], method="spearman")
}

# Create a dataframe to store the results
result <- data.frame(
  col_name = colnames(downsampled_modules_t),
  Pearson_Correlation = pearson_correlations,
  Spearman_Correlation = spearman_correlations
)
write.table(result, "by_module_correlations.txt", sep = "\t")





