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
library(Matrix)
library(RcppML)
library(parallel)
library(doParallel)
library(enrichR)
library(ComplexHeatmap)
library(colorRamp2)
library(ggplot2)
set.seed(NULL)

##############################################################################
##############################################################################
# Inputs:

PATH_TO_DE_GENES = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/de_genes/GL261_integrated_20230712_preinf_noRTNormalized_sorted/deMatSig_adjp05_lfc01.txt"
PATHS_TO_LFC = c("/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/deseq/GL261_integrated_20230705_preinf_noRTNormalized_sorted")
PATH_TO_SEURAT_OBJECT = "/raleighlab/data1/liuj/gbm_perturb/analysis/GL261_integrated_20230619.Rds"
PATH_TO_NORMALIZED_DATA = NULL
CONTEXTS = c("preinf")
RADIATION_CONDS = c("noRT", "RT")
NT_COLS = c("preinf_non-targeting_RT_non-targeting_noRT")
PERTURBS_TO_REMOVE = c("NA_RT", "NA_noRT", "non-targeting_B_RT", "non-targeting_B_noRT")
OUTPUT_DIR = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/nmf/RT_noRT/preinf_noRTNormalized"
OUTPUT_FILE_NAME = "preinf_res_list_log_noRTNormalized_ranks2_40.rds"
ZERO_INFLATION_PARAMETER = 0.01
MODULE_MEMBERSHIP_PERCENTILE_PARAMETER = 0.9

print(paste("PATH_TO_DE_GENES:", PATH_TO_DE_GENES))
print(paste("PATHS_TO_LFC:", PATHS_TO_LFC))
print(paste("PATH_TO_SEURAT_OBJECT:", PATH_TO_SEURAT_OBJECT))
print(paste("CONTEXTS:", CONTEXTS))
print(paste("RADIATION_CONDS:", RADIATION_CONDS))
print(paste("NT_COLS:", NT_COLS))
print(paste("PERTURBS_TO_REMOVE:", PERTURBS_TO_REMOVE))
print(paste("OUTPUT_DIR:", OUTPUT_DIR))
print(paste("OUTPUT_FILE_NAME:", OUTPUT_FILE_NAME))
print(paste("ZERO_INFLATION_PARAMETER:", ZERO_INFLATION_PARAMETER))
print(paste("MODULE_MEMBERSHIP_PERCENTILE_PARAMETER:", MODULE_MEMBERSHIP_PERCENTILE_PARAMETER))

##############################################################################
##############################################################################
# Code

# Retrieve the expression data and subset it by context, RT and DE genes. Note
# that because we eventually want to visualize expression with SCT, we only
# include genes that pass the SCT filters.

if (is.null(PATH_TO_NORMALIZED_DATA)) {
  data.main = readRDS(PATH_TO_SEURAT_OBJECT)
  DefaultAssay(data.main) = "RNA"
  data = NormalizeData(data.main, normalization.method = "LogNormalize", scale.factor = 1e6, verbose = TRUE)
} else {
  data = readRDS(PATH_TO_NORMALIZED_DATA)
}

data = subset(data, source %in% CONTEXTS)
data = subset(data, cond %in% RADIATION_CONDS)

sgRNAConds = unique(data$sgRNACond)
to_keep = setdiff(sgRNAConds, PERTURBS_TO_REMOVE)
data = subset(data, sgRNACond %in% to_keep)

deMatSig = read.table(PATH_TO_DE_GENES)
deMatSig = na.omit(deMatSig)

de_genes = rownames(deMatSig)
data = subset(data, features = de_genes)

# Create and then pseudobulk by sgRNAContextTreatment

data$sgRNAContextTreatment = paste(data$sgRNA, data$source, data$cond, sep = "_")

counts_data = as.matrix(GetAssayData(data))
metadata = data@meta.data
expr_df = as.data.frame(counts_data)
expr_df$gene = rownames(expr_df)
expr_long = tidyr::pivot_longer(expr_df, -gene, 
                                names_to = "cell_barcode",
                                values_to = "expression")
expr_long$sgRNAContextTreatment = metadata[expr_long$cell_barcode, "sgRNAContextTreatment"]
avg_data = expr_long %>%
  dplyr::group_by(gene, sgRNAContextTreatment) %>%
  dplyr::summarise(mean_expression = mean(expression, na.rm = TRUE))
avg_matrix = tidyr::pivot_wider(avg_data, names_from = sgRNAContextTreatment, values_from = mean_expression)
avg_matrix = as.data.frame(avg_matrix)
rownames(avg_matrix) = avg_matrix$gene
avg_matrix = as.matrix(avg_matrix[, -1])
avg_matrix = avg_matrix[apply(avg_matrix, 1, function(row) any(row != 0)), ]

# Obtain a rank estimate by getting NMF at many different ranks.

rank_range = 2:40
nrun = 10
ncores = 10

run_nmf = function(r) {
  model = NMF::nmf(avg_matrix, r, method = "brunet", nrun = 30)
}
res_list = mclapply(rank_range, run_nmf, mc.cores = ncores)
saveRDS(res_list, paste(OUTPUT_DIR, OUTPUT_FILE_NAME, sep = "/"))

# Examine and save the mean squared error and cophenetic correlation coefficient

mse_list = lapply(res_list, function(model) {
  reconstructed = fitted(model)
  mse = mean((avg_matrix - reconstructed)^2)
})
png(paste(OUTPUT_DIR, "mse_2-40.png", sep = "/"), height = 640, width = 640)
plot(1:length(mse_list), mse_list,
     xlab = "NMF Rank", ylab = "Mean Squared Error", pch = 19)
dev.off()

coph_corr_list = lapply(res_list, function(model) {
  consensus_matrix = consensus(model)
  coph_corr = cophcor(consensus_matrix)
})
png(paste(OUTPUT_DIR, "coph_2-40.png", sep = "/"), height = 640, width = 640)
plot(1:length(coph_corr_list), coph_corr_list,
     xlab = "NMF Rank", ylab = "Cophenetic Correlation", pch = 19)
dev.off()

# In general, the rank estimate for noRT invitro seems to hover around 20.
# Compute a stable NMF result for 15 and calculate gene modules.

rank = 20
res = res_list[[rank - 1]]
W = basis(res)

# Some genes are very zero inflated. We eliminate them before using the basis
# matrix.

row_means <- rowMeans(W)
W <- W[row_means >= ZERO_INFLATION_PARAMETER, ]

# Extract gene names into gene modules

gene_names = rownames(W)
W.df = as.data.frame(W)
gene_modules <- data.frame(gene = gene_names, modules = NA)
for (i in 1:nrow(W.df)) {
  gene_modules[i, "modules"] = which.max(W.df[i,])
  # threshold = MODULE_MEMBERSHIP_PERCENTILE_PARAMETER * max(W.df[i,])
  # belongs_to = which(W.df[i,] >= threshold)
  # gene_modules[i, "modules"] = paste(belongs_to, collapse = "_")
}

# Once we have modules, transform the modules into module groups. Throw out
# genes that aren't represented in SCT analysis or genes that are lncRNAs
# or Rik cDNA.

genes_to_remove = setdiff(rownames(data@assays$RNA@counts),
                          rownames(data@assays$SCT@counts))
module_groups = list()
for (i in 1:nrow(gene_modules)) {
  gene = gene_modules[i, "gene"]
  if (gene %in% genes_to_remove || grepl("Rik$|^Gm\\d+$", gene)) {
    print(gene)
    next
  }
  modules = strsplit(as.character(gene_modules[i, "modules"]), "_")[[1]]
  for (module in modules) {
    module_groups[[module]] = c(module_groups[[module]], gene)
  }
  # module = as.character(gene_modules[i, "modules"])
  # module_groups[[module]] = c(module_groups[[module]], gene)
}

lengths = lapply(module_groups, function(x) {length(x)})
flat_list = unlist(module_groups)
count_vec <- sapply(module_groups, function(x) sum(x %in% flat_list[duplicated(flat_list)]))

# Label the module groups using EnrichR

enrichr_groups = list()
ontology_groups = data.frame(msigdb_term = NA, go_term = NA, n_genes = NA,
                             msigdb_genes = NA, go_term_genes = NA, genes = NA)
for (i in 1:length(module_groups)) {
  setEnrichrSite("Enrichr")
  dbs <- c("MSigDB_Hallmark_2020","KEGG_2021_Human",
           "WikiPathway_2021_Human","GO_Biological_Process_2021")
  goi = module_groups[[i]]
  enrichROut = enrichr(goi, dbs)
  msigDB_terms = enrichROut$MSigDB_Hallmark_2020$Term
  GO_terms = enrichROut$GO_Biological_Process_2021$Term
  combined_term = paste(msigDB_terms[1], GO_terms[1], sep = "|")
  combined_term_and_genes = c(combined_term, goi)
  combined_term_and_genes = paste(combined_term_and_genes, collapse = "|")
  enrichr_groups[[combined_term_and_genes]] = goi

  ontology_groups[i, "msigdb_term"] = msigDB_terms[1]
  ontology_groups[i, "go_term"] = GO_terms[1]
  ontology_groups[i, "n_genes"] = length(goi)
  ontology_groups[i, "msigdb_genes"] = enrichROut$MSigDB_Hallmark_2020$Genes[1]
  ontology_groups[i, "go_term_genes"] = enrichROut$GO_Biological_Process_2021$Genes[1]
  ontology_groups[i, "genes"] = paste(goi, collapse = "|")
  Sys.sleep(3)
}
write.table(ontology_groups, paste(OUTPUT_DIR, "ontology_groups_rank20_MaxModules.txt", sep = "/"), sep = "\t")

# Generate the log fold change matrix. We keep only those rows
# that qualify and that have genes that are in the gene_names

data.list = list()
for (directory in PATHS_TO_LFC) {
  file_names = list.files(directory)
  for (i in file_names) {
    name = gsub(".csv", "", i)
    df = read.table(paste(directory, "/", i, sep = ""))
    data.list[[name]] = df
  }
}

nmfMat <- ldply(lapply(data.list,function(x) x[,c("feature","log_fc")]),data.frame)
nmfMat <- reshape(nmfMat, idvar = "feature", v.names = "log_fc",timevar = ".id", direction = "wide")
nmfMat <- nmfMat[!is.na(nmfMat$feature),]
row.names(nmfMat) <- nmfMat$feature
colnames(nmfMat) <- gsub("log_fc.","",colnames(nmfMat))
nmfMat <- nmfMat[,-1]
nmfMat[is.na(nmfMat)] <- 0
nmfMat <- nmfMat[, !colnames(nmfMat) %in% NT_COLS]
nmfMatsig <- nmfMat[intersect(gene_names,row.names(nmfMat)),]
nmfMatsig = as.matrix(nmfMatsig)

# Cluster the log fold changes into gene modules

module_means = lapply(enrichr_groups, function(g) {
  colMeans(nmfMatsig[g, , drop = FALSE])
})
expr_matrix = do.call(rbind, module_means)

# Plot a context separated heatmap with annotations.

phenoTbl <- read.table('/raleighlab/data1/liuj/gbm_perturb/analysis/GL261_1+2_results_table.txt',sep='\t',header=TRUE,row.names=1)
coverageTbl <- read.table('/raleighlab/data1/liuj/gbm_perturb/analysis/gbm_pdx_perturb_GL261_integrate_xfp/pdx_perturb_GL261_concordant_sgRNAs.txt',header=TRUE,sep='\t')
coverageTblCmp <- rbind(invitronoRT = apply(coverageTbl[c("GL261_48hit_noRT_1","GL261_48hit_noRT_2"),],2,sum),
                        invitroRT = apply(coverageTbl[c("GL261_48hit_RT_1","GL261_48hit_RT_2"),],2,sum),
                        preinfnoRT = apply(coverageTbl[c("GL261_noRT_preinf_MACSFACS_1","GL261_noRT_preinf_MACSFACS_2","GL261_noRT_preinf_MACSFACS_3"),],2,sum),
                        preinfRT = apply(coverageTbl[c("GL261_RT_preinf_MACSFACS_1","GL261_RT_preinf_MACSFACS_2","GL261_RT_preinf_MACSFACS_3"),],2,sum),
                        CEDnoRT = apply(coverageTbl[c("GL261_CED_pool","GL261_noRT_CED_MACSFACS"),],2,sum),
                        CEDRT = coverageTbl[c("GL261_RT_CED_MACSFACS"),])
colnames(coverageTblCmp) <- gsub("non.targeting","NTC",colnames(coverageTblCmp))
sgAnnot <- data.frame(row.names = colnames(coverageTblCmp), gamma = phenoTbl[colnames(coverageTblCmp),"gamma"], tau = phenoTbl[colnames(coverageTblCmp),"tau"], rho = phenoTbl[colnames(coverageTblCmp),"rho"])
sgAnnot["non.targeting",] <- c(0,0,0)
tblAnnot <- data.frame(row.names=colnames(nmfMatsig),
                       sgRNA = sapply(strsplit(colnames(nmfMatsig),"_"), `[`, 2),
                       source = sapply(strsplit(colnames(nmfMatsig),"_"), `[`, 1),
                       cond = sapply(strsplit(colnames(nmfMatsig),"_"), `[`, 3))
tblAnnot$sourceCond <- paste(tblAnnot$source,tblAnnot$cond,sep='')
tblAnnot$gamma <- sgAnnot[tblAnnot$sgRNA,"gamma"]
tblAnnot$tau <- sgAnnot[tblAnnot$sgRNA,"tau"]
tblAnnot$rho <- sgAnnot[tblAnnot$sgRNA,"rho"]
tblAnnot$coverage <- rep(0,nrow(tblAnnot))
for (i in 1:nrow(tblAnnot)){
  tryCatch({
    tblAnnot[i,"coverage"] <- coverageTblCmp[tblAnnot[i,"sourceCond"],tblAnnot[i,"sgRNA"]]
  }, error = function(e) {
    print(e)
  })
}
tblAnnot$sgRNACond = paste(tblAnnot$sgRNA,tblAnnot$cond,sep='_')
write.table(tblAnnot, paste(OUTPUT_DIR, "/", "annotation_table.txt", sep = ""), sep = "\t")

poi = colnames(expr_matrix)
ht = Heatmap(
  expr_matrix,
  name = "NMF Log Fold Change Matrix RT vs. noRT, CED",
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  top_annotation = HeatmapAnnotation(gamma = anno_barplot(as.numeric(tblAnnot[poi,"gamma"]), bar_width = 0.9),
                                     tau = anno_barplot(as.numeric(tblAnnot[poi,"tau"]), bar_width = 0.9),
                                     rho = anno_barplot(as.numeric(tblAnnot[poi,"rho"]), bar_width = 0.9),
                                     GO = tblAnnot[poi,"GO"],
                                     treatment = tblAnnot[poi,"cond"],
                                     context = tblAnnot[poi,"source"],
                                     coverage = as.numeric(tblAnnot[poi,"coverage"]),
                                     col = list(treatment = c("noRT" = "#2166AC", "RT" = "#B2182B", "control" = "green4"),
                                                context = c("invitro" = "#D53E4F", "preinf" = "#FC8D59","CED" = "#FEE08B"))),
  column_split = factor(tblAnnot[poi,"cond"], levels = c("noRT", "RT")),
  width = ncol(expr_matrix)*unit(4, "mm"),
  height = nrow(expr_matrix)*unit(4, "mm"),
)
pdf(paste(OUTPUT_DIR, "output_matrix_lognormalized_lfc_rank20_MaxModules_ced_noRTNormalized.pdf", sep = "/"), height = 30, width = 30)
draw(ht, heatmap_legend_side = "top", annotation_legend_side = "top")
dev.off()


ht <- Heatmap(
  expr_matrix,
  name = "Log Fold Changes Across Gene Modules and Perturbations",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  width = ncol(expr_matrix)*unit(7, "mm"),
  height = nrow(expr_matrix)*unit(7, "mm"),
)
pdf(paste(OUTPUT_DIR, "output_matrix_lognormalized_lfc_rank25.pdf", sep = "/"), height = 15, width = 60)
draw(ht, heatmap_legend_side = "top")
dev.off()
