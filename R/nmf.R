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
log_expr_matrix = log1p(expr_matrix)

# Determine an optimal rank
res = nmfEstimateRank(log_expr_matrix, range = 2:30, nrun = 30)
plot(res)

# After choosing the optimal rank from the plot, run NMF
rank <- 10  # This should be replaced with the rank you have chosen
nmf_result <- nmf(log_expr_matrix, rank = rank, nrun = 20)

# Extract the basis and mixture matrices
w <- basis(nmf_result)  # W matrix
h <- coef(nmf_result)  # H matrix