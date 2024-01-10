## -------------------------------------------------------------------------------------------------
.libPaths(c("/raleighlab/data1/czou/gbm_perturb/gbm_perturb/renv/library/R-4.3/x86_64-pc-linux-gnu", .libPaths()))


## -------------------------------------------------------------------------------------------------
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
library(viridis)
library(knitr)
set.seed(NULL) # For NMF


## -------------------------------------------------------------------------------------------------
PATH.TO.SEURAT.OBJECT <- "/raleighlab/data1/liuj/gbm_perturb/analysis/SB28_micro_3_sgRNApos_20231203.Rds"
CELL.TYPE <- commandArgs(trailingOnly=TRUE)[1]
print(sprintf("Determined cell type is %s", CELL.TYPE))
PATH.TO.DESEQ.OUTPUT <- sprintf("/raleighlab/data1/czou/gbm_perturb/gbm_perturb_sb28_clean_outputs/deseq/SB28_integrated_20231207_ME/sgRNA_gene_sgNegCtrl3_knockouts/%s", CELL.TYPE)
OUTPUT.DIR <- sprintf("/raleighlab/data1/czou/gbm_perturb/gbm_perturb_sb28_clean_outputs/nmf/ME/prkdc_excluded_05pval/%s", CELL.TYPE)
RES.LIST.PATH <- sprintf("/raleighlab/data1/czou/gbm_perturb/gbm_perturb_sb28_clean_outputs/nmf/ME/prkdc_excluded_05pval/%s/%s_res_list_2-40_60.rds", CELL.TYPE, CELL.TYPE)


## -------------------------------------------------------------------------------------------------
# Import the DESeq data
deseq.list <- list()
file_names <- list.files(PATH.TO.DESEQ.OUTPUT)
for (i in file_names) {
  name <- gsub(".csv", "", i)
  df <- read.table(sprintf("%s/%s", 
                           PATH.TO.DESEQ.OUTPUT, i))
  deseq.list[[name]] <- df
}

# Filter for DE genes
deGenesAll.padj <- unique(unlist(lapply(deseq.list, function(x) {
  x = x[!is.na(x$log_fc), ]
  x = x[!is.na(x$pvalue), ]
  x = x[!is.na(x$padj), ]
  # x = x[x$padj < 0.05, ]
  x = x[x$pvalue < 0.05, ]
  x = x[abs(x$log_fc) > 0.1, ]
  x$feature
})))
print(sprintf("Found %s DE genes", length(deGenesAll.padj)))

# Build a Log 2FC matrix
deMat <- ldply(lapply(deseq.list,function(x) x[,c("feature","log_fc")]),data.frame)
deMat <- reshape(deMat, idvar = "feature", v.names = "log_fc",timevar = ".id", direction = "wide")
deMat <- deMat[!is.na(deMat$feature),]
row.names(deMat) <- deMat$feature
colnames(deMat) <- gsub("log_fc.","",colnames(deMat))
deMat <- deMat[,-1]
deMat[is.na(deMat)] <- 0
deMatsig <- deMat[intersect(deGenesAll.padj,row.names(deMat)),]


## -------------------------------------------------------------------------------------------------
# Import the data
data.main <- readRDS(PATH.TO.SEURAT.OBJECT)


## -------------------------------------------------------------------------------------------------
# Replace all of the oligodendrocyte progenitor mentions with OPCs
data.main$scMRMA_manual <- gsub("Oligodendrocyte progenitor cells", "OPCs", data.main$scMRMA_manual)

# Log-normalize the raw count data
DefaultAssay(data.main) <- "RNA"
data <- NormalizeData(data.main, normalization.method = "LogNormalize")

# Keep only the cells whose sgRNA_full guides fall into those we want to process
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

# Subset to the desired cell type
data <- subset(data, scMRMA_manual == CELL.TYPE)

# Subset to the union of DE genes across all perturbations for that cell type
de.genes <- rownames(deMatsig)
data <- subset(data, features = de.genes)

# Screen for minimum coverage and eliminate all perturbations that don't have more than 5 cells
sgRNA_gene.counts <- table(data$sgRNA_gene)
data = subset(data, sgRNA_gene %in% names(
  sgRNA_gene.counts[sgRNA_gene.counts > 5]
))


## -------------------------------------------------------------------------------------------------
# Pseudobulk by sgRNA_gene
counts_data = as.matrix(GetAssayData(data))
metadata = data@meta.data
expr_df = as.data.frame(counts_data)
expr_df$gene = rownames(expr_df)
expr_long = tidyr::pivot_longer(expr_df, -gene, 
                                names_to = "cell_barcode",
                                values_to = "expression")
expr_long$sgRNA_gene = metadata[expr_long$cell_barcode, "sgRNA_gene"]
avg_data = expr_long %>%
  dplyr::group_by(gene, sgRNA_gene) %>%
  dplyr::summarise(mean_expression = mean(expression, na.rm = TRUE))
avg_matrix = tidyr::pivot_wider(avg_data, names_from = sgRNA_gene, values_from = mean_expression)
avg_matrix = as.data.frame(avg_matrix)
rownames(avg_matrix) = avg_matrix$gene
avg_matrix = as.matrix(avg_matrix[, -1])


## -------------------------------------------------------------------------------------------------
# Filter for those perturbs we have DESeq output for. This means that we filter for only the perturbs that are included in colnames(deMatsig) AND the avg_matrix.
perturbs.in.deseq <- colnames(deMatsig)
perturbs.in.deseq <- gsub("_sgNegCtrl_3", "", 
                          perturbs.in.deseq)
avg_matrix <- avg_matrix[,intersect(perturbs.in.deseq, colnames(avg_matrix))]

# Remove Prkdc if needed
avg_matrix <- avg_matrix[, colnames(avg_matrix) != "Prkdc"]

# Ensure there are no zero or NULL rows
avg_matrix = avg_matrix[apply(avg_matrix, 1, function(row) any(row != 0)), ]

# Sync the deMatsig rows and columns with the avg_matrix columns
perturbations.in.deMatsig <- paste(colnames(avg_matrix), "sgNegCtrl_3", sep = "_")
deMatsig <- deMatsig[, perturbations.in.deMatsig]

print(sprintf("Evaluating %s genes", dim(avg_matrix)[1]))


## -------------------------------------------------------------------------------------------------
# Run NMF at many different ranks
print("Beginning NMF estimate")
rank.range <- 2:40
ncores <- 20
run.nmf = function(r) {
  model = NMF::nmf(avg_matrix, r, method = "brunet",
                   nrun = 60)
}
res.list = mclapply(rank.range, run.nmf, mc.cores = ncores)
saveRDS(res.list, sprintf("%s/%s_res_list_2-40_60.rds",
                          OUTPUT.DIR, CELL.TYPE))
print("Finished NMF estimate")


## -------------------------------------------------------------------------------------------------
# Load the res_list object
res_list = readRDS(RES.LIST.PATH)


## -------------------------------------------------------------------------------------------------
# Generate the mean squared error plot
mse_list = lapply(res_list, function(model) {
  reconstructed = fitted(model)
  mse = mean((avg_matrix - reconstructed)^2)
})
png(paste(OUTPUT.DIR, "mse_2-40_2.png", sep = "/"), height = 640, width = 640)
plot(1:length(mse_list), mse_list,
     xlab = "NMF Rank", ylab = "Mean Squared Error", pch = 19)
dev.off()

# Generate the cophenetic coefficient plot
coph_corr_list = lapply(res_list, function(model) {
  consensus_matrix = consensus(model)
  coph_corr = cophcor(consensus_matrix)
})
png(paste(OUTPUT.DIR, "coph_2-40_2.png", sep = "/"), height = 640, width = 640)
plot(1:length(coph_corr_list), coph_corr_list,
     xlab = "NMF Rank", ylab = "Cophenetic Correlation", pch = 19)
dev.off()


## -------------------------------------------------------------------------------------------------
RANK = 7
ZERO.INFLATION.PARAMETER = 0.01
res = res_list[[RANK - 1]]
W = basis(res)

# Eliminate very zero inflated genes
row.means <- rowMeans(W)
W <- W[row.means >= ZERO.INFLATION.PARAMETER, ]

# Assign each gene to a gene module
gene.names <- rownames(W)
W.df <- as.data.frame(W)
gene.modules <- data.frame(gene = gene.names, modules = NA)
for (i in 1:nrow(W.df)) {
  gene.modules[i, "modules"] <- which.max(W.df[i, ])
}

# Extract the assignments into module groups
genes.to.remove <- setdiff(rownames(data@assays$RNA@counts),
                          rownames(data@assays$SCT@counts))
module.groups <- list()
for (i in 1:nrow(gene.modules)) {
  gene = gene.modules[i, "gene"]
  if (gene %in% genes.to.remove || grepl("Rik$|^Gm\\d+$", gene)) {
    print(sprintf("Removing %s", gene))
    next
  }
  modules = strsplit(as.character(gene.modules[i, "modules"]), "_")[[1]]
  for (module in modules) {
    module.groups[[module]] = c(module.groups[[module]], gene)
  }
}


## -------------------------------------------------------------------------------------------------
# Label the groups using Enrichr
enrichr_groups = list()
ontology_groups = data.frame(msigdb_term = NA, go_term = NA, 
                             n_genes = NA, msigdb_genes = NA,
                             go_term_genes = NA, genes = NA)
for (i in 1:length(module.groups)) {
  setEnrichrSite("Enrichr")
  dbs <- c("MSigDB_Hallmark_2020","KEGG_2021_Human",
           "WikiPathway_2021_Human","GO_Biological_Process_2021")
  goi = module.groups[[i]]
  enrichROut = enrichr(goi, dbs)
  msigDB_terms = enrichROut$MSigDB_Hallmark_2020$Term
  GO_terms = enrichROut$GO_Biological_Process_2021$Term
  combined_term = paste(msigDB_terms[1], GO_terms[1], sep = "|")
  combined_term_and_genes = c(combined_term, goi)
  combined_term_and_genes = paste(combined_term_and_genes, collapse = "|")
  enrichr_groups[[combined_term]] = goi

  ontology_groups[i, "msigdb_term"] = msigDB_terms[1]
  ontology_groups[i, "go_term"] = GO_terms[1]
  ontology_groups[i, "n_genes"] = length(goi)
  ontology_groups[i, "msigdb_genes"] = enrichROut$MSigDB_Hallmark_2020$Genes[1]
  ontology_groups[i, "go_term_genes"] = enrichROut$GO_Biological_Process_2021$Genes[1]
  ontology_groups[i, "genes"] = paste(goi, collapse = "|")
  Sys.sleep(3)
}
write.table(ontology_groups, paste(OUTPUT.DIR, sprintf("ontology_groups_rank%s_MaxModules.txt", RANK), sep = "/"), sep = "\t")
ontology_groups = read.table(paste(OUTPUT.DIR, sprintf("ontology_groups_rank%s_MaxModules.txt", RANK), sep = "/"), sep = "\t")


## -------------------------------------------------------------------------------------------------
# Generate a list of all genes for later use
list_names = paste(ontology_groups$msigdb_term, ontology_groups$go_term, sep = "|")
list_items <- lapply(ontology_groups$genes, function(x) strsplit(x, split = "\\|")[[1]])
enrichr_groups = setNames(list_items, list_names)
all_genes = unlist(enrichr_groups)


## -------------------------------------------------------------------------------------------------
lfc.matrix <- deMatsig[intersect(all_genes,row.names(deMatsig)),]
module_means = lapply(enrichr_groups, function(g) {
  colMeans(lfc.matrix[g, , drop = FALSE])
})
expr_matrix = do.call(rbind, module_means)
colnames(expr_matrix) = sapply(strsplit(colnames(expr_matrix), "_"), function(x) paste(x[2], x[1], sep="_"))
write.table(expr_matrix, paste(OUTPUT.DIR, sprintf("output_matrix_lognormalized_lfc_rank%s_MaxModules.txt", RANK), sep = "/"), sep = "\t")


## -------------------------------------------------------------------------------------------------
# Generate a heatmap
poi = colnames(expr_matrix)
ht = Heatmap(expr_matrix,
        name = sprintf("Log Fold Change Matrix Clustered on Rank-%d NMF for %s", RANK, CELL.TYPE),
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        cluster_column_slices = FALSE,
        cluster_columns = TRUE,
        width = ncol(expr_matrix)*unit(4, "mm"),
        height = nrow(expr_matrix)*unit(4, "mm"),
)
pdf(sprintf("%s/output_matrix_lognormalized_lfc_rank%d_MaxModules_range.pdf", OUTPUT.DIR, RANK),
    width=35,
    height=19,
    useDingbats=FALSE)
draw(ht, heatmap_legend_side = "top", annotation_legend_side = "top")
dev.off()

print("Finished executing!")

