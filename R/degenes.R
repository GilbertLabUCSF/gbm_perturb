# Transform DESeq output into a set of differentially expressed genes per
# perturbation per context per condition.
# 
# Users should configure:
# - An input directory, with DESeq outputs
# - An output directory
# - The percentile of LFC that we consider "differentially expressed"
# 
# Output: a design matrix with log fold changes across differentially expressed
# genes.

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

INPUT_DIR = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/deseq/GL261_integrated_20230626_ced_condNormalized_all"
OUTPUT_DIR = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/de_genes/GL261_integrated_20230626_ced_condNormalized_all"
PERCENTILE = 0.1
SEED = 5220

##############################################################################
##############################################################################
# Code

file_names <- list.files(INPUT_DIR)
data.list <- list()
for (i in file_names) {
  name <- gsub(".csv", "", i)
  df <- read.table(paste("differential_gene_sets/deseq_3cond_malignantOnly_ntNoRTNormalized/", i, sep = ""))
  # DE <- filter(df, log_fc != 0 & abs(log_fc) >= 0.1 & padj < 0.05)
  data.list[[name]] <- df
}

# get all DE genes 

deGenesAll <- unique(unlist(lapply(data.list,function(x) x[x$padj < 0.05 & abs(x$log_fc) > 0.1,"feature"])))

# build log2FC matrix 

deMat <- ldply(lapply(data.list,function(x) x[,c("feature","log_fc")]),data.frame)
deMat <- reshape(deMat, idvar = "feature", v.names = "log_fc",timevar = ".id", direction = "wide")
deMat <- deMat[!is.na(deMat$feature),]
row.names(deMat) <- deMat$feature
colnames(deMat) <- gsub("log_fc.","",colnames(deMat))
deMat <- deMat[,-1]
deMat[is.na(deMat)] <- 0
deMat <- deMat[, !colnames(deMat) %in% c("ced_ntRT_ntNoRT", "iv_ntRT_ntNoRT", "preinf_ntRT_ntNoRT")]
deMatsig <- deMat[intersect(deGenesAll,row.names(deMat)),]

write.table(deMat, "linear_modeling/deseq_3cond_malignantOnly_wgcna/datasets/3Cond_deMat.txt")
write.table(deMatsig, "linear_modeling/deseq_3cond_malignantOnly_wgcna/datasets/3cond_deMatsig.txt")