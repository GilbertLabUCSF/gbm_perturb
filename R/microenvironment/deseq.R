# Microenvironment DESeq
# 
# This script takes in an example Seurat object called microenvironment_data.rds 
# and outputs DESeq data for cells that are labeled with metadata scMRMA_manual
# %in% DESIRED.CELL.TYPES. Note that this script shows how we used
# only sgNegCtrl_3 and relabels oligodendrocyte progenitor cells as OPCs.
# 
# Inputs:
# - The input Seurat object, which should have scMRMA_manual, sgRNA_full, and
#   sgRNA_gene
# - The output directory
# - Desired cell types
# 
# Outputs: 
# - DElegate output as a CSV file per perturbation.

## -----------------------------------------------------------------------------------------------------------------------------------------------
# Dependencies:

library(Seurat)
library(DElegate)
library(dplyr)

## -----------------------------------------------------------------------------------------------------------------------------------------------
# Inputs:

PATH.TO.SEURAT.OBJECT = "microenvironment_data.rds"
OUTPUT.DIR = "data/microenvironment/deseq"
SEED = 5220
DESIRED.CELL.TYPES = c("Oligodendrocytes", "Astrocytes",
                       "Microglia", "Macrophages", "OPCs")

## -----------------------------------------------------------------------------------------------------------------------------------------------
# Code

data = readRDS(PATH.TO.SEURAT.OBJECT)

# Define a function that helps us to find differential genes

find.deseq.differential.genes = function(data.obj, seed, group1, group2, group.column=NULL) {
  # Takes in a normalized Seurat object with metadata
  # indicating group1 and group2 and a group_column argument
  # indicating the metadata column to be used. Sets a seed, then
  # finds differential expressed genes and outputs a data frame.
  # Note that log fold comparisons will be output as log2(group1/group2)
  set.seed(seed)
  
  # Set futures to greater than max capacity; you may need to tweak this
  options(future.globals.maxSize = 5000 * 1024^2)
  
  df = findDE(object = data.obj, group_column = group.column,
              compare = c(group1, group2), method = 'deseq')
  return(df)
}

build.filename = function(group1, group2) {
  filename = paste(group1, group2, sep = "_")
  directory = paste(OUTPUT.DIR, "/", sep = "")
  extension = ".csv"
  return(paste(directory, filename, extension, sep = ""))
}

# Do some metadata cleaning
data$scMRMA_manual <- gsub("Oligodendrocyte progenitor cells", "OPCs", data$scMRMA_manual)

# Keep only the cells that have guides from sgRNA_full we want to keep
data$guide.cell <- paste(data$sgRNA_full, data$scMRMA_manual, sep = "_")
knockouts <- read.table("/raleighlab/data1/czou/gbm_perturb/gbm_perturb/R/microenvironment/pdx_micro_sb28_3_pseudobulk_KD_greaterthan_030.txt",
                        sep = ",")
colnames(knockouts) <- c("guide", "cell_type")
guide.cells.to.keep <- paste(knockouts$guide, knockouts$cell_type, sep = "_")
guide.cells.to.keep <- c(
  guide.cells.to.keep,
  "sgNegCtrl_3_Astrocytes",
  "sgNegCtrl_3_Macrophages",
  "sgNegCtrl_3_OPCs",
  "sgNegCtrl_3_Oligodendrocytes",
  "sgNegCtrl_3_Microglia",
  "sgNegCtrl_3_SB28"
)
guide.cells.to.keep <- gsub("Oligodendrocyte progenitor cells", "OPCs", guide.cells.to.keep)

data <- subset(data, guide.cell %in% guide.cells.to.keep)

# Since we only want to mark sgNegCtrl_3 as the negative control, find all of the
# cells that are affiliated with it and mark their sgRNA_gene differently
data <- data %>% 
  AddMetaData(
    metadata = data@meta.data %>%
      mutate(
        sgRNA_gene = ifelse(sgRNA_full == "sgNegCtrl_3", "sgNegCtrl_3", sgRNA_gene)
      )
  )

# For each of the desired cell types, run DESeq against the non-targeting control.
# This run will use the sgRNA_gene metadata.
for (cell.type in DESIRED.CELL.TYPES) {
  data.cell.type <- subset(data, scMRMA_manual == cell.type)
  
  # Screen for only those perturbations that meet coverage requirements
  sgRNA_gene.counts <- table(data.cell.type$sgRNA_gene)
  data.cell.type <- subset(data.cell.type, sgRNA_gene %in% names(sgRNA_gene.counts[sgRNA_gene.counts > 5]))

  # Screen out control perturbations
  output.dir <- sprintf("%s/%s/%s", OUTPUT.DIR, "sgRNA_gene_sgNegCtrl3_knockouts", cell.type)
  print(sprintf("Outputting %s sgRNA_gene deseq into %s", cell.type, output.dir))
  perturbations <- unique(data.cell.type$sgRNA_gene)
  perturbations <- perturbations[!perturbations == "sgNegCtrl"]
  perturbations <- perturbations[!perturbations == "sgNegCtrl_3"]
  found = FALSE
  for (perturb in perturbations) {
    tryCatch({
      df.deseq <- find.deseq.differential.genes(
        data.cell.type,
        SEED,
        perturb,
        "sgNegCtrl_3",
        group.column = "sgRNA_gene"
      )
      found = TRUE
    }, error = function(err) {
      print(sprintf("Failed to get model because of %s in %s", err, perturb))
    })
    if (found) {
      write.table(df.deseq, sprintf("%s/%s_%s.csv",
                                    output.dir,
                                    perturb,
                                    "sgNegCtrl_3"))
    }
  }
}
