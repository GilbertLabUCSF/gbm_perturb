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

INPUT_DIR = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/deseq/GL261_integrated_20230626_invitro_condNormalized_all"
OUTPUT_DIR = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/de_genes/GL261_integrated_20230626_invitro_condNormalized"
NT_COLS = c("CED_non-targeting_RT_non-targeting_noRT")
EXPR_CUTOFF = 0.5
ABS_CUTOFF = 0.01
TOP_BOT_CUTOFF = 0.01
SEED = 5220

##############################################################################
##############################################################################
# Code

# Import the DESeq outputs

file_names = list.files(INPUT_DIR)
data.list = list()
for (i in file_names) {
  name = gsub(".csv", "", i)
  df = read.table(paste(INPUT_DIR, "/", i, sep = ""))
  data.list[[name]] = df
}

# Filter for only noRT perturbations - this is for analysis development

df.names = names(data.list)
matched.names = df.names[grep("non-targeting_noRT", df.names)]
data.list = data.list[matched.names]

# Filter the dataframes. We apply a 1% cutoff to both the up- and down-regulated
# directions, remove all genes with no log fold changes measured, and admit
# only those genes in the top 50% of expression as measured by DESeq2's baseMean.

deGenesAll.topBot = unique(unlist(lapply(data.list, function(x) {
  x = x[!is.na(x$log_fc), ]
  
  n = nrow(x)
  cutoff_n = round(n * EXPR_CUTOFF)
  x = x[order(x$ave_expr),]
  x = tail(x, cutoff_n)
  print(paste("Minimum ave_expr of allowed genes: ", min(x$ave_expr)))
  
  n = nrow(x)
  cutoff_n = round(n * TOP_BOT_CUTOFF)
  x = x[order(x$log_fc),]
  top_genes = tail(x, cutoff_n)
  bottom_genes = head(x, cutoff_n)
  print(paste("Minimum log_fc of top genes: ", min(top_genes$log_fc)))
  print(paste("Maximum log_fc of bottom genes: ", max(bottom_genes$log_fc)))
  print(length(c(top_genes$feature, bottom_genes$feature)))
  c(top_genes$feature, bottom_genes$feature)  # Return the combined vector
})))

# Filter the dataframes based on only the absolute values.

deGenesAll.abs = unique(unlist(lapply(data.list, function(x) {
  x = x[!is.na(x$log_fc), ]
  
  n = nrow(x)
  cutoff_n = round(n * EXPR_CUTOFF)
  x = x[order(x$ave_expr),]
  x = tail(x, cutoff_n)
  print(paste("Minimum ave_expr of allowed genes: ", min(x$ave_expr)))
  
  n = nrow(x)
  cutoff_n = round(n * ABS_CUTOFF)
  print(cutoff_n)
  x = x[order(abs(x$log_fc)),]
  top_genes = tail(x, cutoff_n)
  print(paste("Minimum absolute log_fc of top genes: ", min(abs(top_genes$log_fc))))
  print(paste("Maximum absolute log_fc of top genes: ", max(abs(top_genes$log_fc))))
  print(quantile(top_genes$log_fc))
  print(length(top_genes$feature))
  c(top_genes$feature)  # Return the combined vector
  
  # x = x[abs(x$log_fc) > 1,]
  # print(dim(x))
  # c(x$feature)
})))

# Build log2FC matrices for both the topBot and Abs case

deMat <- ldply(lapply(data.list,function(x) x[,c("feature","log_fc")]),data.frame)
deMat <- reshape(deMat, idvar = "feature", v.names = "log_fc",timevar = ".id", direction = "wide")
deMat <- deMat[!is.na(deMat$feature),]
row.names(deMat) <- deMat$feature
colnames(deMat) <- gsub("log_fc.","",colnames(deMat))
deMat <- deMat[,-1]
deMat[is.na(deMat)] <- 0
deMat <- deMat[, !colnames(deMat) %in% NT_COLS]
deMatsig <- deMat[intersect(deGenesAll.topBot,row.names(deMat)),]

write.table(deMat, paste(OUTPUT_DIR, "/deMat_noRTOnly_topBot", TOP_BOT_CUTOFF, ".txt", sep = ""))
write.table(deMatsig, paste(OUTPUT_DIR, "/deMatSig_noRTOnly_topBot", TOP_BOT_CUTOFF, ".txt", sep = ""))

deMat <- ldply(lapply(data.list,function(x) x[,c("feature","log_fc")]),data.frame)
deMat <- reshape(deMat, idvar = "feature", v.names = "log_fc",timevar = ".id", direction = "wide")
deMat <- deMat[!is.na(deMat$feature),]
row.names(deMat) <- deMat$feature
colnames(deMat) <- gsub("log_fc.","",colnames(deMat))
deMat <- deMat[,-1]
deMat[is.na(deMat)] <- 0
deMat <- deMat[, !colnames(deMat) %in% NT_COLS]
deMatsig <- deMat[intersect(deGenesAll.abs,row.names(deMat)),]

write.table(deMat, paste(OUTPUT_DIR, "/deMat_noRTOnly_abs", ABS_CUTOFF, ".txt", sep = ""))
write.table(deMatsig, paste(OUTPUT_DIR, "/deMatSig_noRTOnly_abs", ABS_CUTOFF, ".txt", sep = ""))
