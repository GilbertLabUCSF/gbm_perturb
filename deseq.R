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
set.seed(5220)

##############################################################################
##############################################################################
# Inputs:

PATH_TO_SEURAT_OBJECT = "/raleighlab/data1/liuj/gbm_perturb/analysis/GL261_integrated_20230619.Rds"
EXP_CONTEXT = "iv"
SORTED_IDENTITIES = c()
OUTPUT_DIR = "/JJJJ/JJJ/JJ"
SEED = 5220
NORMALIZE_TO_NT_NORT = FALSE
MINIMUM_COVERAGE = 5

##############################################################################
##############################################################################
# Code

data = readRDS(PATH_TO_SEURAT_OBJECT)

# Isolate the data that we need for running DElegate

cell_barcodes = colnames(data)
guides = data$sgRNA
radiation = data$cond
origin = data$source
sgRNA_metadata = data.frame(
  origin = origin,
  cell_barcode = cell_barcodes,
  radiation = radiation,
  guide = guides
)

data.context = subset(data, source == EXP_CONTEXT)
sgRNA_metadata.context = subset(sgRNAMetadata, origin == EXP_CONTEXT) 
perturbations.context = unique(filter(sgRNA_metadata.context, Guide != "non-targeting", 
                                   Guide != "non-targeting_B", Guide != "None")$Guide)

# Screen for sorted cells only (or other malignancy indicator)

data.context <- subset(data.context, orig.ident %in% SORTED_IDENTITIES)

# Screen for coverage



# Functions that help us to find differential genes

find_deseq_differential_genes = function(data.obj, seed, group1, group2, group_column=NULL) {
  # Takes in a normalized Seurat object with metadata
  # indicating group1 and group2 and a group_column argument
  # indicating the metadata column to be used. Sets a seed, then
  # finds differential expressed genes and outputs a data frame.
  # Note that log fold comparisons will be output as log2(group1/group2)
  set.seed(seed)
  
  # Set futures to greater than max capacity; you may need to tweak this
  options(future.globals.maxSize = 5000 * 1024^2)
  
  df <- findDE(object = data.obj, group_column = group_column,
               compare = c(group1, group2), method = 'deseq')
  return(df)
}

build_filename = function(group1, group2) {
  filename = paste(EXP_CONTEXT, group1, group2, sep = "_")
  directory = paste(OUTPUT_DIR, "/", sep = "")
  extension = ".csv"
  return(paste(directory, filename, extension, sep = ""))
}

# Find differential genes, starting with the non-targeting case

df.nt = find_deseq_differential_genes(data.context, 5220, "non-targeting_RT", 
                                       "non-targeting_noRT", group_column = "sgRNACond")
write.table(df.ced, build_filename("non-targeting_RT", "non-targeting_noRT"))

for (perturb in perturbations.context) {
  print(sprintf("Calculating diff genes for %s", perturb))
  for (cond in c("RT", "noRT")) {
    perturb_condition = paste(perturb, cond, sep = "_")
    if (NORMALIZE_TO_NT_noRT) {
      nt_condition = paste("non-targeting", "noRT", sep = "_")
    } else {
      nt_condition = paste("non-targeting", cond, sep = "_")
    }
    found = FALSE
    tryCatch({
      df <- find_deseq_differential_genes(data.ced,
                                          SEED,
                                          perturb_condition,
                                          nt_condition,
                                          group_column = "sgRNACond"
      )
      found = TRUE
    }, error = function(err) {
      print(paste("Failed to get model because of ", err, "in", perturb_condition))
    })
    if (found) {
      filename = build_filename(perturb_condition, nt_condition)
      write.table(df, filename)
    }
  }
}
