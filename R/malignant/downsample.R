# Downsampling that forms part of the preprocessing pipeline. May
# unify with preprocess_malignant.R in future refactors. Running this script
# will take in the Seurat object, subject it to a coverage constraint, restrict
# to only perturbations found in all contexts, and downsample each perturbation
# to the minimum number of cells across contexts.
# 
# Author: Christopher Zou

##############################################################################
##############################################################################
# Dependencies:

library(Seurat)

##############################################################################
##############################################################################
# Inputs:

PATH_TO_SEURAT_OBJECT = "/raleighlab/data1/liuj/gbm_perturb/analysis/GL261_integrated_20230619.Rds"
SEED = 325
set.seed(SEED)
COVERAGE_FILTER = 5
OUTPUT_PATH = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/downsampled_objects/GL261_integrated_downsampled_20230819.rds"

cat("PATH TO SEURAT OBJECT: ", PATH_TO_SEURAT_OBJECT, "\n")
cat("SEED: ", SEED, "\n")
cat("COVERAGE FILTER: ", COVERAGE_FILTER, "\n")
cat("OUTPUT PATH: ", OUTPUT_PATH, "\n")

##############################################################################
##############################################################################
# Code:

data = readRDS(PATH_TO_SEURAT_OBJECT)
DefaultAssay(data.main) = "RNA"
data = NormalizeData(data.main, normalization.method = "LogNormalize", scale.factor = 1e6, verbose = TRUE)

# Filter for only those perturbations that have passed a coverage filter.

data$sgRNACondContext = paste(data$sgRNA, data$cond, data$source, sep = "_")
sgRNACondContext_counts = table(data$sgRNACondContext)
data = subset(data, sgRNACondContext %in% 
                names(sgRNACondContext_counts[sgRNACondContext_counts > COVERAGE_FILTER]))

# Filter for only those perturbation x cond that are represented in all three 
# contexts.

all_sgRNACondContexts = unique(data$sgRNACondContext)
data$sgRNACond = paste(data$sgRNA, data$cond, sep = "_")
all_sgRNAConds = unique(data$sgRNACond)
qualified_sgRNAConds = c()
for (sgRNACond in all_sgRNAConds) {
  necessary_sgRNACondContexts = c(
    paste(sgRNACond, "CED", sep = "_"),
    paste(sgRNACond, "invitro", sep = "_"),
    paste(sgRNACond, "preinf", sep = "_")
  )
  if (all(necessary_sgRNACondContexts %in% all_sgRNACondContexts)) {
    qualified_sgRNAConds = unique(c(qualified_sgRNAConds, sgRNACond))
  }
}
data = subset(data, sgRNACond %in% qualified_sgRNAConds)

# Downsample each perturbation x cond to its minimum among all three contexts

metadata_df = data@meta.data
all_sgRNAConds = unique(data$sgRNACond)
all_sgRNACondContext_counts = table(data$sgRNACondContext)
cells_to_include = c()
for (sgRNACond in all_sgRNAConds) {
  sgRNACondContexts = c(
    paste(sgRNACond, "CED", sep = "_"),
    paste(sgRNACond, "invitro", sep = "_"),
    paste(sgRNACond, "preinf", sep = "_")
  )
  min_count <- min(all_sgRNACondContext_counts[sgRNACondContexts])
  for (perturb_cond_context in sgRNACondContexts) {
    df = subset(metadata_df, sgRNACondContext == perturb_cond_context)
    cells = rownames(df)
    cells_to_add = sample(cells, size = min_count)
    cells_to_include = c(cells_to_include, cells_to_add)
  }
}
data = subset(data, cells = cells_to_include)
saveRDS(data, file = OUTPUT_PATH)

print("Finished executing!")