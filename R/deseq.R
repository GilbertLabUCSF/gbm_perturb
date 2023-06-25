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

PATH_TO_SEURAT_OBJECT = "/raleighlab/data1/liuj/gbm_perturb/analysis/SB28_integrated_20230620.Rds"
EXP_CONTEXT = "invitro"
SORTED_IDENTITIES = c("sort")
NT_GUIDES = c("non-targeting_A","non-targeting_B","non-targeting_C","non-targeting_D")
IGNORE_GUIDES = c("NA_RT", "NA_noRT")
OUTPUT_DIR = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_sb28_clean_outputs/de_genes/SB28_integrated_202306025_invitro_condNormalized"
SEED = 5220
NORMALIZE_TO_NT_NORT = FALSE
MINIMUM_COVERAGE = 5

##############################################################################
##############################################################################
# Code

data = readRDS(PATH_TO_SEURAT_OBJECT)

# Isolate the data by context

data.context = subset(data, source == EXP_CONTEXT)

# Screen for sorted cells only (or other malignancy indicator)

data.context = subset(data.context, sorted %in% SORTED_IDENTITIES)

# Screen for cells part of a group with coverage of >= 5 cells

sgRNACond_counts = table(data.context$sgRNACond)
data.context = subset(data.context, sgRNACond %in% 
                        names(sgRNACond_counts[sgRNACond_counts > 5]))

# Get rid of all guides in the ignore category

guides = unique(data.context$sgRNACond)
guides <- guides[!(guides %in% IGNORE_GUIDES)]
data.context = subset(data.context, sgRNACond %in% guides)

# Isolate the data that we need for running DElegate

cell_barcodes = colnames(data.context)
guides = data.context$sgRNA
radiation = data.context$cond
origin = data.context$source
sgRNA_metadata = data.frame(
  origin = origin,
  cell_barcode = cell_barcodes,
  radiation = radiation,
  guide = guides
)
perturbations.context = unique(subset(sgRNA_metadata, !(guide %in% NT_GUIDES))$guide)

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
  filename = paste(EXP_CONTEXT, group1, group2, sep = "_")
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
    df.nt = find_deseq_differential_genes(data.context, SEED, RT, 
                                          noRT, group_column = "sgRNACond")
    found = TRUE
  }, error = function(err) {
    print(paste("Failed to get model because of ", err, "in", guide))
  })
  if (found) {
    write.table(df.nt, build_filename(RT, noRT))
  }
}

# If there are multiple non-targeting perturbations, combine them into a single
# non-targeting perturbation and analyze that.
for (guide in NT_GUIDES) {
  data.context$sgRNACond = gsub(guide, "non-targeting", data.context$sgRNACond)
}

for (perturb in perturbations.context) {
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
      df = find_deseq_differential_genes(data.context,
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
