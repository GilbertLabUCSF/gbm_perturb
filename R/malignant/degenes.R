# Transform DElegate output into a matrix of log fold changes across all genes
# for differentially expressed genes. 
# 
# Author: Christopher Zou
# 
# Users should configure:
# - A list of input directories with DElegate CSV outputs
# - The output directory
# - Non-targeting cell columns (these should be shared with the DElegate
#   CSV filenames) to ignore
# - An expression level cutoff
# - A percentile cutoff for absolute value thresholding and top/bottom x%
#   thresholding
# - An adjusted p-value cutoff
# - A p-value cutoff
# - A log fold change cutoff.
# 
# Output: a design matrix with log fold changes across differentially expressed
# genes. The rows of this matrix are the DE genes and the columns are the 
# perturbation/context/condition cells.

##############################################################################
##############################################################################
# Dependencies:

library(Seurat)
library(DElegate)
library(plyr)
library(dplyr)

##############################################################################
##############################################################################
# Inputs:

INPUT_DIRS = c("")
OUTPUT_DIR = ""
NT_COLS = c("")
EXPR_CUTOFF = 0.5
ABS_CUTOFF = 0.01
TOP_BOT_CUTOFF = 0.01
ADJ_P_CUTOFF = 0.05
P_CUTOFF
LFC_CUTOFF = 0.1

##############################################################################
##############################################################################
# Code

# Import the DESeq outputs

data.list.all = list()
for (dir in INPUT_DIRS) {
  file_names = list.files(dir)
  for (i in file_names) {
    name = gsub(".csv", "", i)
    df = read.table(paste(dir, "/", i, sep = ""))
    data.list.all[[name]] = df
  }
}

# Filter the dataframe by a TOP_BOT_CUTOFF to the up- and down-regulated 
# directions of LFC. Remove all genes with no lfc, padj, or pvalue measured, 
# and admit only those genes in the top EXPR_CUTOFF of expression as measured 
# byDESeq2's baseMean.

deGenesAll.topBot = unique(unlist(lapply(data.list, function(x) {
  x = x[!is.na(x$log_fc), ]
  x = x[!is.na(x$pvalue), ]
  x = x[!is.na(x$padj), ]
  
  n = nrow(x)
  cutoff_n = round(n * EXPR_CUTOFF)
  x = x[order(x$ave_expr),]
  x = tail(x, cutoff_n)
  print(paste("Minimum ave_expr of allowed genes: ", min(x$ave_expr)))

  # Implement p/absolute-value threshold
  
  x$threshold_p = abs(x$log_fc) * -log10(x$pvalue)
  x$threshold_adjp = abs(x$log_fc) * -log10(x$padj)

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

# Filter the dataframe to include the top ABS_CUTOFF of genes by absolute value
# of log fold change. Admit only genes in the top EXPR_CUTOFF of expression
# as measured by DESeq2's baseMean.

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
  c(top_genes$feature)

})))

# Filter the dataframes based on ADJ_P_CUTOFF and LFC_CUTOFF. Include
# no baseMean cutoff.

deGenesAll.padj = unique(unlist(lapply(data.list, function(x) {
  x = x[!is.na(x$log_fc), ]
  x = x[!is.na(x$pvalue), ]
  x = x[!is.na(x$padj), ]
  x = x[x$padj < ADJ_P_CUTOFF, ]
  x = x[abs(x$log_fc) > LFC_CUTOFF, ]
  x$feature
})))

# Build log2FC matrices for the topBot, Abs, and LFC/adjp cutoffs.

deMat <- ldply(lapply(data.list,function(x) x[,c("feature","log_fc")]),data.frame)
deMat <- reshape(deMat, idvar = "feature", v.names = "log_fc",timevar = ".id", direction = "wide")
deMat <- deMat[!is.na(deMat$feature),]
row.names(deMat) <- deMat$feature
colnames(deMat) <- gsub("log_fc.","",colnames(deMat))
deMat <- deMat[,-1]
deMat[is.na(deMat)] <- 0
deMat <- deMat[, !colnames(deMat) %in% NT_COLS]
deMatsig <- deMat[intersect(deGenesAll.topBot,row.names(deMat)),]

write.table(deMat, paste(OUTPUT_DIR, "/deMat_topBot", TOP_BOT_CUTOFF, ".txt", sep = ""))
write.table(deMatsig, paste(OUTPUT_DIR, "/deMatSig_topBot", TOP_BOT_CUTOFF, ".txt", sep = ""))

deMat <- ldply(lapply(data.list,function(x) x[,c("feature","log_fc")]),data.frame)
deMat <- reshape(deMat, idvar = "feature", v.names = "log_fc",timevar = ".id", direction = "wide")
deMat <- deMat[!is.na(deMat$feature),]
row.names(deMat) <- deMat$feature
colnames(deMat) <- gsub("log_fc.","",colnames(deMat))
deMat <- deMat[,-1]
deMat[is.na(deMat)] <- 0
deMat <- deMat[, !colnames(deMat) %in% NT_COLS]
deMatsig <- deMat[intersect(deGenesAll.abs,row.names(deMat)),]

write.table(deMat, paste(OUTPUT_DIR, "/deMat_abs", ABS_CUTOFF, ".txt", sep = ""))
write.table(deMatsig, paste(OUTPUT_DIR, "/deMatSig_abs", ABS_CUTOFF, ".txt", sep = ""))

deMat <- ldply(lapply(data.list,function(x) x[,c("feature","log_fc")]),data.frame)
deMat <- reshape(deMat, idvar = "feature", v.names = "log_fc",timevar = ".id", direction = "wide")
deMat <- deMat[!is.na(deMat$feature),]
row.names(deMat) <- deMat$feature
colnames(deMat) <- gsub("log_fc.","",colnames(deMat))
deMat <- deMat[,-1]
deMat[is.na(deMat)] <- 0
deMat <- deMat[, !colnames(deMat) %in% NT_COLS]
deMatsig <- deMat[intersect(deGenesAll.padj,row.names(deMat)),]

write.table(deMat, paste(OUTPUT_DIR, "/deMat_adjp05_lfc1", ".txt", sep = ""))
write.table(deMatsig, paste(OUTPUT_DIR, "/deMatSig_adjp05_lfc1", ".txt", sep = ""))
