# Run DESeq2 using DElegate on a Seurat object according to our filters 
# and configurations.
# 
# Author: Christopher Zou
# 
# Users should configure:
# - The input Seurat object, which should have source, sorted, cond and 
#   sgRNACond metadata
# - The experimental context
# - The identities of the sorted cells
# - Whether to screen for sorted cells only
# - The non-targeting guides
# - Guides to ignore
# - The output directory
# - A random seed
# - Whether we normalize to non-targeting noRT or normalize to the condition
# - The bottom coverage threshold. We only consider cells with coverage above
#   this.
# This file can be run as a script once configured.
# 
# Output: DElegate output as a CSV file per perturbation. By default, DElegate
# uses 3 cells per pseudoreplicate.

##############################################################################
##############################################################################
# Dependencies:

library(Seurat)
library(DElegate)
library(dplyr)

##############################################################################
##############################################################################
# Inputs:

PATH_TO_SEURAT_OBJECT = "/raleighlab/data1/liuj/gbm_perturb/analysis/GL261_integrated_20230619.Rds"
SORTED_IDENTITIES = c("FACS")
SORTED_IDENTITES_ONLY = TRUE
NT_GUIDES = c("non-targeting")
IGNORE_GUIDES = c("")
OUTPUT_DIR = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/deseq/condNormalized_invitro"
SEED = 5220
NORMALIZE_TO_NT_NORT = TRUE
MINIMUM_COVERAGE = 5

print(paste("PATH_TO_SEURAT_OBJECT =", PATH_TO_SEURAT_OBJECT))
print(paste("SORTED_IDENTITIES =", toString(SORTED_IDENTITIES)))
print(paste("SORTED_IDENTITES_ONLY =", SORTED_IDENTITES_ONLY))
print(paste("NT_GUIDES =", toString(NT_GUIDES)))
print(paste("IGNORE_GUIDES =", toString(IGNORE_GUIDES)))
print(paste("OUTPUT_DIR =", OUTPUT_DIR))
print(paste("SEED =", SEED))
print(paste("NORMALIZE_TO_NT_NORT =", NORMALIZE_TO_NT_NORT))
print(paste("MINIMUM_COVERAGE =", MINIMUM_COVERAGE))
set.seed(SEED)
print(paste("SEED SET TO ", SEED))

##############################################################################
##############################################################################
# Code

data = readRDS(PATH_TO_SEURAT_OBJECT)

# Screen for only human genes

DefaultAssay(data) = "RNA"
all_genes = rownames(data)
human_genes = all_genes[grep("GRCh38-", all_genes)]
data = subset(data, features = human_genes)

# Screen for sorted cells only (or other malignancy indicator)

if (SORTED_IDENTITES_ONLY) {
  data = subset(data, sorted %in% SORTED_IDENTITIES)
}

# Screen for cells part of a group with coverage of >= MINIMUM_COVERAGE cells

sgRNACond_counts = table(data$sgRNACond)
data = subset(data, sgRNACond %in% 
                        names(sgRNACond_counts[sgRNACond_counts >= MINIMUM_COVERAGE]))

# Isolate the data that we need for running DElegate

cell_barcodes = colnames(data)
guides = data$sgRNA
radiation = data$cond
origin = "CED"
sgRNA_metadata = data.frame(
  origin = origin,
  cell_barcode = cell_barcodes,
  radiation = radiation,
  guide = guides
)
perturbations = unique(subset(sgRNA_metadata, !(guide %in% NT_GUIDES))$guide)

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

# Find differential genes, starting with the non-targeting case(s)

for (guide in NT_GUIDES) {
  RT = paste(guide, "RT", sep = "_")
  noRT = paste(guide, "noRT", sep = "_")
  found = FALSE
  tryCatch({
    df.nt = find_deseq_differential_genes(data, SEED, RT, 
                                          noRT, group_column = "sgRNACond")
    found = TRUE
  }, error = function(err) {
    print(paste("Failed to get model because of ", err, "in", guide))
  })
  if (found) {
    write.table(df.nt, build_filename(RT, noRT))
  }
}

for (perturb in perturbations) {
  print(sprintf("Calculating diff genes for %s", perturb))
  for (cond in c("RT", "noRT")) {
    perturb_condition = paste(perturb, cond, sep = "_")
    if (NORMALIZE_TO_NT_NORT) {
      nt_condition = paste("non-targeting", "noRT", sep = "_")
    } else {
      nt_condition = paste("non-targeting", cond, sep = "_") 
    }
    found = FALSE
    tryCatch({
      df = find_deseq_differential_genes(data,
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

print("Finished running!")
