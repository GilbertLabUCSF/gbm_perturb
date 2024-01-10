## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
.libPaths(c("/raleighlab/data1/czou/gbm_perturb/gbm_perturb/renv/library/R-4.3/x86_64-pc-linux-gnu", .libPaths()))


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(Seurat)
library(DElegate)
library(dplyr)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
data <- readRDS(
  "/raleighlab/data1/liuj/gbm_perturb/analysis/SB28_micro_3_all_20231209.Rds"
)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
data.sorted <- subset(data, sorted == "sort")
data.sorted@meta.data <- data.sorted@meta.data %>% mutate(sgRNA_full = ifelse(is.na(sgRNA_full) & sgRNA_binary == FALSE, "unidentified_guide", sgRNA_full))


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
data.sorted$scMRMA_manual <- gsub("Oligodendrocyte progenitor cells", "OPCs", data.sorted$scMRMA_manual)


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
OUTPUT.DIR = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_sb28_clean_outputs/deseq/SB28_integrated_20231207_ME"
SEED = 5220
DESIRED.CELL.TYPES = c("SB28", "Oligodendrocytes", "Astrocytes",
                       "Microglia", "Macrophages", "OPCs")


## ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
for (cell.type in DESIRED.CELL.TYPES) {
  data.cell.type <- subset(data.sorted, scMRMA_manual == cell.type)
  output.dir <- sprintf("%s/%s/%s", OUTPUT.DIR, "sgRNA_full_unidentified_guides", cell.type)
  print(sprintf("Outputting %s sgRNA_full deseq into %s", cell.type, output.dir))
  perturbations <- unique(data.cell.type$sgRNA_full)
  perturbations <- perturbations[!perturbations == "unidentified_guide"]
  perturbations <- perturbations[!is.na(perturbations)]
  found = FALSE
  for (perturb in perturbations) {
    tryCatch({
      df.deseq <- find.deseq.differential.genes(
        data.cell.type,
        SEED,
        perturb,
        "unidentified_guide",
        group.column = "sgRNA_full"
      )
      found = TRUE
    }, error = function(err) {
      print(sprintf("Failed to get model because of %s in %s", err, perturb))
    })
    if (found) {
      write.table(df.deseq, sprintf("%s/%s_%s.csv",
                                    output.dir,
                                    perturb,
                                    "unidentified_guide"))
    }
  }
}

print("finished running")
