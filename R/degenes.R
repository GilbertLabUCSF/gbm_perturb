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
library(plyr)
library(dplyr)
set.seed(5220)

##############################################################################
##############################################################################
# Inputs:

INPUT_DIR = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/deseq/GL261_integrated_20230626_ced_condNormalized_all"
OUTPUT_DIR = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/de_genes/GL261_integrated_20230626_ced_condNormalized_all"
NT_COLS = c("CED_non-targeting_RT_non-targeting_noRT")
PERCENTILE = 0.1
SEED = 5220

##############################################################################
##############################################################################
# Code

file_names = list.files(INPUT_DIR)
data.list = list()
for (i in file_names) {
  name = gsub(".csv", "", i)
  df = read.table(paste(INPUT_DIR, "/", i, sep = ""))
  data.list[[name]] = df
}

# Filter the dataframes for the top 2.5% upregulated and downregulated

deGenesAll = unique(unlist(lapply(data.list, function(x) {
  x = x[!is.na(x$log_fc), ]
  x = x[order(x$log_fc),]  # Sort data based on log_fc values
  
  n = nrow(x)
  cutoff_n = round(n * 0.025)  # We set a 2.5% cutoff on either end
  top_genes = tail(x, cutoff_n)
  bottom_genes = head(x, cutoff_n)
  
  # Print minimum log_fc of top genes and maximum log_fc of bottom genes
  print(paste("Minimum log_fc of top genes: ", min(top_genes$log_fc)))
  print(paste("Maximum log_fc of bottom genes: ", max(bottom_genes$log_fc)))
  
  print(length(c(top_genes$feature, bottom_genes$feature)))
  c(top_genes$feature, bottom_genes$feature)  # Return the combined vector
})))

# Filter the dataframes based on only the absolute values - but top 5% of those.

deGenesAll = unique(unlist(lapply(data.list, function(x) {
  x = x[!is.na(x$log_fc), ]
  x = x[order(abs(x$log_fc)),]  # Sort data based on log_fc values
  
  n = nrow(x)
  cutoff_n = round(n * 0.05)  # We set a 5% cutoff on the top end
  top_genes = tail(x, cutoff_n)

  # Print minimum log_fc of top genes and maximum log_fc of bottom genes
  print(paste("Minimum absolute log_fc of top genes: ", min(abs(top_genes$log_fc))))
  print(paste("Maximum absolute log_fc of top genes: ", max(abs(top_genes$log_fc))))
  c(top_genes$feature)  # Return the combined vector
})))

# build log2FC matrix 

deMat <- ldply(lapply(data.list,function(x) x[,c("feature","log_fc")]),data.frame)
deMat <- reshape(deMat, idvar = "feature", v.names = "log_fc",timevar = ".id", direction = "wide")
deMat <- deMat[!is.na(deMat$feature),]
row.names(deMat) <- deMat$feature
colnames(deMat) <- gsub("log_fc.","",colnames(deMat))
deMat <- deMat[,-1]
deMat[is.na(deMat)] <- 0
deMat <- deMat[, !colnames(deMat) %in% NT_COLS]
deMatsig <- deMat[intersect(deGenesAll,row.names(deMat)),]

write.table(deMat, paste(OUTPUT_DIR, "/deMat_topBot025Cutoff.txt", sep = ""))
write.table(deMatsig, paste(OUTPUT_DIR, "/deMatSig_topBot025Cutoff.txt", sep = ""))