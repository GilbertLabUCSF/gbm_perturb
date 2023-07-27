# Run DESeq2 on a Seurat object according to our filters and configurations
# 
# Users should configure:
# - The intput Seurat object
# - The identities of the sorted cells
# - The output directory
# - The experimental context to generate output for. We default to "iv"
# - A random seed
# - Whether we normalize to non-targeting noRT or normalize to the condition
# 
# Output: deseq output per perturbation. By default, we use 3 cells per
# pseudoreplicate.

##############################################################################
##############################################################################
# Dependencies:

library(Seurat)
library(DElegate)
library(dplyr)
set.seed(5220)

##############################################################################
##############################################################################
# Inputs:

PATH_TO_SEURAT_OBJECT = "/raleighlab/data1/liuj/gbm_perturb/analysis/SB28_micro4_20230619.Rds"
NT_GUIDE = "sgNegCtrl"
OUTPUT_DIR = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_sb28_clean_outputs/deseq/sb28_integrated_micro4_macrophagesOnly"
SEED = 5220
MINIMUM_COVERAGE = 5

print(paste("PATH_TO_SEURAT_OBJECT =", PATH_TO_SEURAT_OBJECT))
print(paste("NT_GUIDES =", toString(NT_GUIDES)))
print(paste("OUTPUT_DIR =", OUTPUT_DIR))
print(paste("SEED =", SEED))
print(paste("MINIMUM_COVERAGE =", MINIMUM_COVERAGE))

##############################################################################
##############################################################################
# Code

data = readRDS(PATH_TO_SEURAT_OBJECT)

# Subset the data by sgRNA guide positive cells

data.context = subset(data, sgRNA_binary == TRUE)

# Screen for macrophages only

data.context = subset(data.context, new.cluster.ids == "Macrophage")

# Screen for having a guide that's not NA

perturbs = unique(data.context$GeneA)
perturbs <- perturbs[!is.na(perturbs)]
data.context = subset(data.context, GeneA %in% perturbs)

# Screen for cells part of a group with coverage of > 5 cells

sgRNACond_counts = table(data.context$GeneA)
data.context = subset(data.context, GeneA %in% 
                        names(sgRNACond_counts[sgRNACond_counts >= MINIMUM_COVERAGE]))

# Define functions that help us to find differential genes

find_deseq_differential_genes = function(data.obj, seed, group1, group2, group_column=NULL) {
  # Takes in a normalized Seurat object with metadata
  # indicating group1 and group2 and a group_column argument
  # indicating the metadata column to be used. Sets a seed, then
  # finds differential expressed genes and outputs a data frame.
  # Note that log fold comparisons will be output as log2(group1/group2)
  set.seed(seed)
  
  # Set futures to greater than max capacity; you may need to tweak this
  options(future.globals.maxSize = 5000 * 1024^2)
  
  df = findDE(object = data.obj, group_column = group_column,
               compare = c(group1, group2), method = 'deseq')
  return(df)
}

build_filename = function(group1, group2) {
  filename = paste(group1, group2, sep = "_")
  directory = paste(OUTPUT_DIR, "/", sep = "")
  extension = ".csv"
  return(paste(directory, filename, extension, sep = ""))
}

# Find differential genes, normalizing to sgNegCtrl

for (perturb in perturbs) {
  print(sprintf("Calculating diff genes for %s", perturb))
  found = FALSE
  tryCatch({
    df = find_deseq_differential_genes(data.context, SEED, perturb, NT_GUIDE, group_column = "GeneA")
    found = TRUE
  }, error = function(err) {
    print(paste("Failed to get model because of ", err, "in", perturb))
  })
  if (found) {
    filename = build_filename(perturb, NT_GUIDE)
    write.table(df, filename)
  }
}

print("Finished running!")
