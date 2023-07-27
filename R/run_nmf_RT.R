# Run non-negative matrix factorization given a list of differentially
# expressed genes.

##############################################################################
##############################################################################
# Dependencies:

library(Seurat)
library(DElegate)
library(plyr)
library(dplyr)
library(NMF)
library(Matrix)
library(parallel)
library(doParallel)
set.seed(NULL)

##############################################################################
##############################################################################
# Inputs:

PATH_TO_DE_GENES = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/de_genes/GL261_integrated_20230712_ced_noRTNormalized_all/deMatSig_adjp05_lfc01.txt"
PATHS_TO_LFC = c("/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/deseq/GL261_integrated_20230705_ced_noRTNormalized_all")
PATH_TO_SEURAT_OBJECT = "/raleighlab/data1/liuj/gbm_perturb/analysis/GL261_integrated_20230619.Rds"
PATH_TO_NORMALIZED_DATA = NULL
CONTEXTS = c("CED")
RADIATION_CONDS = c("noRT", "RT")
PERTURBS_TO_REMOVE = c("NA_RT", "NA_noRT", "non-targeting_B_RT", "non-targeting_B_noRT")
OUTPUT_DIR = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/nmf/RT_noRT/ced_noRTNormalized"
OUTPUT_FILE_NAME = "ced_res_list_log_noRTNormalized_ranks2_40.rds"
NCORES = 10

print(paste("PATH_TO_DE_GENES:", PATH_TO_DE_GENES))
print(paste("PATHS_TO_LFC:", PATHS_TO_LFC))
print(paste("PATH_TO_SEURAT_OBJECT:", PATH_TO_SEURAT_OBJECT))
print(paste("CONTEXTS:", CONTEXTS))
print(paste("RADIATION_CONDS:", RADIATION_CONDS))
print(paste("PERTURBS_TO_REMOVE:", PERTURBS_TO_REMOVE))
print(paste("OUTPUT_DIR:", OUTPUT_DIR))
print(paste("OUTPUT_FILE_NAME:", OUTPUT_FILE_NAME))
print(paste("NCORES:", NCORES))

##############################################################################
##############################################################################
# Code

# Retrieve the expression data and subset it by context, RT and DE genes. Note
# that because we eventually want to visualize expression with SCT, we only
# include genes that pass the SCT filters.

data.main = readRDS(PATH_TO_SEURAT_OBJECT)
DefaultAssay(data.main) = "RNA"
data = NormalizeData(data.main, normalization.method = "LogNormalize", scale.factor = 1e6, verbose = TRUE)

data = subset(data, source %in% CONTEXTS)
data = subset(data, cond %in% RADIATION_CONDS)

sgRNAConds = unique(data$sgRNACond)
to_keep = setdiff(sgRNAConds, PERTURBS_TO_REMOVE)
data = subset(data, sgRNACond %in% to_keep)

deMatSig = read.table(PATH_TO_DE_GENES)
deMatSig = na.omit(deMatSig)

de_genes = rownames(deMatSig)
data = subset(data, features = de_genes)

# Create and then pseudobulk by sgRNAContextTreatment

data$sgRNAContextTreatment = paste(data$sgRNA, data$source, data$cond, sep = "_")

counts_data = as.matrix(GetAssayData(data))
metadata = data@meta.data
expr_df = as.data.frame(counts_data)
expr_df$gene = rownames(expr_df)
expr_long = tidyr::pivot_longer(expr_df, -gene, 
                                names_to = "cell_barcode",
                                values_to = "expression")
expr_long$sgRNAContextTreatment = metadata[expr_long$cell_barcode, "sgRNAContextTreatment"]
avg_data = expr_long %>%
  dplyr::group_by(gene, sgRNAContextTreatment) %>%
  dplyr::summarise(mean_expression = mean(expression, na.rm = TRUE))
avg_matrix = tidyr::pivot_wider(avg_data, names_from = sgRNAContextTreatment, values_from = mean_expression)
avg_matrix = as.data.frame(avg_matrix)
rownames(avg_matrix) = avg_matrix$gene
avg_matrix = as.matrix(avg_matrix[, -1])
avg_matrix = avg_matrix[apply(avg_matrix, 1, function(row) any(row != 0)), ]

# Filter by those perturbs that we have deseq output for. This means they passed
# a coverage filter.

perturb_list = c()
for (directory in PATHS_TO_LFC) {
  file_names = list.files(directory)
  for (i in file_names) {
    split_string = strsplit(i, "_")[[1]]
    perturb = paste(split_string[2], split_string[1], split_string[3], sep = "_")
    perturb_list = c(perturb_list, perturb)
  }
}
perturb_list = c(perturb_list, paste("non-targeting", CONTEXTS[1], "noRT", sep = "_"))
avg_matrix = avg_matrix[,perturb_list]

# Obtain a rank estimate by getting NMF at many different ranks.

rank_range = 2:40
ncores = NCORES

run_nmf = function(r) {
  model = NMF::nmf(avg_matrix, r, method = "brunet", nrun = 30)
}
res_list = mclapply(rank_range, run_nmf, mc.cores = ncores)
saveRDS(res_list, paste(OUTPUT_DIR, OUTPUT_FILE_NAME, sep = "/"))

print("Finished executing")