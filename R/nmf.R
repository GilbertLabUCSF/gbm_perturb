# Run non-negative matrix factorization given a list of differentially
# expressed genes.

##############################################################################
##############################################################################
# Dependencies:

library(Seurat)
library(DElegate)
library(plyr)
library(dplyr)
library(NMF)
set.seed(5220)

##############################################################################
##############################################################################
# Inputs:

PATH_TO_DE_GENES = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/de_genes/GL261_integrated_20230626_invitro_condNormalized/deMatSig_noRTOnly_abs0.01.txt"
PATH_TO_SEURAT_OBJECT = "/raleighlab/data1/liuj/gbm_perturb/analysis/GL261_integrated_20230619.Rds"
CONTEXT = "invitro"
RT = "noRT"
NT_COLS = c("CED_non-targeting_RT_non-targeting_noRT")
SEED = 5220

##############################################################################
##############################################################################
# Code

# Retrieve the expression data and subset it by context, RT and DE genes

data = readRDS(PATH_TO_SEURAT_OBJECT)
data = subset(data, source == CONTEXT)
data = subset(data, cond == RT)

deMatSig = read.table(PATH_TO_DE_GENES, header = TRUE)
de_genes = rownames(deMatSig)
data = subset(data, features = de_genes)

# Run NMF on the matrix. First, check to see if there are any negative values.

expr_data = data@assays$SCT@data
any_negative = any(expr_data < 0)
log_expr_matrix = log1p(expr_data)

# Determine an optimal rank. Note that NMF takes in a dense matrix. I also
# found an RccpL package that can run directly on sparse matrices; we will
# determine rank first though.

log_expr_matrix = as.matrix(log_expr_matrix)
res = nmfEstimateRank(log_expr_matrix, range = 2:30, nrun = 30)
saveRDS(res, file = "res.rds")
# plot(res)
# 
# # After choosing the optimal rank from the plot, run NMF
# 
# rank = 10
# nmf_result = nmf(log_expr_matrix, rank = rank, nrun = 20)
# 
# # Extract the basis and mixture matrices
# 
# w <- basis(nmf_result)
# h <- coef(nmf_result)