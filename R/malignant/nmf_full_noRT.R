# Run non-negative matrix factorization given a list of differentially
# expressed genes. This file provides a full example of our code run
# filtering for noRT only cells.
# 
# Author: Christopher Zou
# 
# Note that in addition to file names configured at the point of output,
# we configure:
# - The path to differentially expressed genes, the output of degenes.R
# - A list of paths to outputs of deseq.R, or DElegate output
# - The path to the base Seurat object, which should have source, cond, RT,
#   and sgRNACond metadata
# - The path to a downsampled Seurat object which will tell us which cells
#   to use
# - A list of contexts to consider
# - The condition, in this case radiation metadata desired
# - A list of non-targeting perturbations to omit in output
# - The desired output directory
# - The desired NMF output file name
# - The number of cores to use when running NMF
# - The zero inflation below which we should not consider genes in NMF output
# - The desired rank, once determined, of NMF
#
# Optionally, you may also include:
# - The path to an annotation table that includes coverage information
# - The path to a GO-term labeling table (this is used for color labeling
#   of go-terms)
# The script will run without these if you remove the annotations in the
# heatmap.
# 
# Once configured, this file can be run as a script
# 
# Output:
# - A list of NMF objects for ranks in RANK_RANGE, with NRUNS_PER_RANK runs
#   at each rank
# - Mean squared error and cophenetic coefficient plots for the rank of choice
# - The ontology-group labeled gene modules in a table for the rank of choice
# - The annotation table generated for all the perturbations considered
# - A heatmap of log fold changes, with genes condensed into NMF-informed
#   gene sets.
# - The underlying expression matrix for the heatmap

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
library(viridis)

##############################################################################
##############################################################################
# Inputs:

PATH_TO_DE_GENES = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/de_genes/GL261_integrated_20230819_3Context_noRTNormalized_noRT_downsampled/deMatSig_adjp05_lfc01.txt"
# PATHS_TO_LFC = c("/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/deseq/GL261_integrated_20230819_ced_noRTNormalized_downsampled_all",
#                  "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/deseq/GL261_integrated_20230819_invitro_noRTNormalized_downsampled_all",
#                  "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/deseq/GL261_integrated_20230819_preinf_noRTNormalized_downsampled_sorted")
PATHS_TO_LFC = c("/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/deseq/GL261_integrated_20230705_invitro_noRTNormalized_all",
                 "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/deseq/GL261_integrated_20230705_preinf_noRTNormalized_sorted",
                 "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/deseq/GL261_integrated_20230705_ced_noRTNormalized_all")
PATH_TO_SEURAT_OBJECT = "/raleighlab/data1/liuj/gbm_perturb/analysis/GL261_integrated_20230619.Rds"
PATH_TO_DOWNSAMPLED_SEURAT_OBJECT = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/downsampled_objects/GL261_integrated_downsampled_20230819.rds"
CONTEXTS = c("invitro", "preinf", "CED")
RT = "noRT"
NT_COLS = c("non-targeting_invitro_RT",
            "non-targeting_preinf_RT",
            "non-targeting_CED_RT")
PERTURBS_TO_REMOVE = c("NA_RT", "NA_noRT", "non-targeting_B_RT", "non-targeting_B_noRT")
OUTPUT_DIR = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/nmf/3Context_noRT_fullDEAndModules_downsampledLFC"
OUTPUT_FILE_NAME = "res_list_2-40.rds"
NCORES = 25
ZERO_INFLATION_PARAMETER = 0.01
MODULE_MEMBERSHIP_PERCENTILE_PARAMETER = 0.9
RANK = 25
RANK_RANGE = 2:40
NRUNS_PER_RANK = 60
SEED = 823
set.seed(SEED)
PHENOTBL = "/raleighlab/data1/liuj/gbm_perturb/analysis/GL261_1+2_results_table.txt"
GO_CATEGORIES = "/raleighlab/data1/liuj/gbm_perturb/analysis/features_48h_GL261_annot2.csv"
COVERAGE_FILTER = 5
DOWNSAMPLE = TRUE

cat("PATH_TO_DE_GENES: ", PATH_TO_DE_GENES, "\n")
cat("PATHS_TO_LFC: ", paste(PATHS_TO_LFC, collapse = ", "), "\n")
cat("PATH_TO_SEURAT_OBJECT: ", PATH_TO_SEURAT_OBJECT, "\n")
cat("CONTEXTS: ", paste(CONTEXTS, collapse = ", "), "\n")
cat("RT: ", RT, "\n")
cat("NT_COLS: ", paste(NT_COLS, collapse = ", "), "\n")
cat("PERTURBS_TO_REMOVE: ", paste(PERTURBS_TO_REMOVE, collapse = ", "), "\n")
cat("OUTPUT_DIR: ", OUTPUT_DIR, "\n")
cat("OUTPUT_FILE_NAME: ", OUTPUT_FILE_NAME, "\n")
cat("NCORES: ", NCORES, "\n")
cat("ZERO_INFLATION_PARAMETER: ", ZERO_INFLATION_PARAMETER, "\n")
cat("RANK: ", RANK, "\n")
cat("RANK RANGE: ", RANK_RANGE, "\n")
cat("NRUNS_PER_RANK: ", NRUNS_PER_RANK, "\n")
cat("PHENOTBL: ", PHENOTBL, "\n")
cat("GO CATEGORIES: ", GO_CATEGORIES, "\n")
cat("SEED: ", SEED, "\n")
cat("COVERAGE FILTER: ", COVERAGE_FILTER, "\n")

##############################################################################
##############################################################################
# Code

# Retrieve the expression data and subset it by context, RT and DE genes. Note
# that because we eventually want to visualize expression with SCT, we only
# include genes that pass the SCT filters.

data.main = readRDS(PATH_TO_SEURAT_OBJECT)
DefaultAssay(data.main) = "RNA"
data = NormalizeData(data.main, normalization.method = "LogNormalize", scale.factor = 1e6, verbose = TRUE)

data.downsampled = readRDS(PATH_TO_DOWNSAMPLED_SEURAT_OBJECT)
downsampled_cells = colnames(data.downsampled)
if (DOWNSAMPLE) {
  data = subset(data, cells = downsampled_cells)
}

data = subset(data, source %in% CONTEXTS)
data = subset(data, cond == RT)

sgRNAConds = unique(data$sgRNACond)
to_keep = setdiff(sgRNAConds, PERTURBS_TO_REMOVE)
data = subset(data, sgRNACond %in% to_keep)

deMatSig = read.table(PATH_TO_DE_GENES)
deMatSig = na.omit(deMatSig)
# rownames(deMatSig) = deMatSig$X
# deMatSig = deMatSig[,-1]

de_genes = rownames(deMatSig)
data = subset(data, features = de_genes)

# Filter for only those perturbations that have passed a coverage filter.

data$sgRNAContext = paste(data$sgRNA, data$source, sep = "_")
sgRNAContext_counts = table(data$sgRNAContext)
data = subset(data, sgRNAContext %in% 
                names(sgRNAContext_counts[sgRNAContext_counts > COVERAGE_FILTER]))

# Filter for only those perturbations that are represented in all three contexts.

all_sgRNAContexts = unique(data$sgRNAContext)
all_sgRNAs = unique(data$sgRNA)
qualified_sgRNAs = c()
for (sgRNA in all_sgRNAs) {
  necessary_sgRNAContexts = c(
    paste(sgRNA, "CED", sep = "_"),
    paste(sgRNA, "invitro", sep = "_"),
    paste(sgRNA, "preinf", sep = "_")
  )
  if (all(necessary_sgRNAContexts %in% all_sgRNAContexts)) {
    qualified_sgRNAs = unique(c(qualified_sgRNAs, sgRNA))
  }
}
data = subset(data, sgRNA %in% qualified_sgRNAs)

# Downsample each perturbation to its minimum among all three contexts

metadata_df = data@meta.data
all_sgRNAs = unique(data$sgRNA)
all_sgRNAContext_counts = table(data$sgRNAContext)
cells_to_include = c()
for (sgRNA in all_sgRNAs) {
  sgRNAContexts = c(
    paste(sgRNA, "CED", sep = "_"),
    paste(sgRNA, "invitro", sep = "_"),
    paste(sgRNA, "preinf", sep = "_")
  )
  min_count <- min(all_sgRNAContext_counts[sgRNAContexts])
  for (perturb_context in sgRNAContexts) {
    df = subset(metadata_df, sgRNAContext == perturb_context)
    cells = rownames(df)
    cells_to_add = sample(cells, size = min_count)
    cells_to_include = c(cells_to_include, cells_to_add)
  }
}
data = subset(data, cells = cells_to_include)

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

# Filter by those perturbs that we have deseq output for. This means they passed
# a coverage filter.

perturb_list = c()
for (directory in PATHS_TO_LFC) {
  file_names = list.files(directory)
  for (i in file_names) {
    split_string = strsplit(i, "_")[[1]]
    perturb = paste(split_string[2], split_string[1], split_string[3], sep = "_")
    perturb_list = c(perturb_list, perturb)
  }
}
noRT_indices <- grep("noRT", perturb_list, value = FALSE, ignore.case = TRUE)
noRT_list = perturb_list[noRT_indices]
perturb_list = intersect(noRT_list, colnames(avg_matrix))

perturb_list = c(perturb_list, paste("non-targeting", CONTEXTS[1], "noRT", sep = "_"))
perturb_list = c(perturb_list, paste("non-targeting", CONTEXTS[2], "noRT", sep = "_"))
perturb_list = c(perturb_list, paste("non-targeting", CONTEXTS[3], "noRT", sep = "_"))
avg_matrix = avg_matrix[,perturb_list]

# Obtain a rank estimate by getting NMF at many different ranks.

print("Beginning NMF")
rank_range = RANK_RANGE
ncores = NCORES

run_nmf = function(r) {
  model = NMF::nmf(avg_matrix, r, method = "brunet", nrun = NRUNS_PER_RANK, seed = SEED)
}
res_list = mclapply(rank_range, run_nmf, mc.cores = ncores)
saveRDS(res_list, paste(OUTPUT_DIR, OUTPUT_FILE_NAME, sep = "/"))

print("Finished executing")

# Load the res_list object

res_list = readRDS(RES_LIST_PATH)

if(sum(is.nan(avg_matrix)) > 0){
  print("Dataset has missing values.")
} else {
  print("Dataset doesn't have missing values.")
}

# Examine and save the mean squared error and cophenetic correlation coefficient

mse_list = lapply(res_list, function(model) {
  reconstructed = fitted(model)
  mse = mean((avg_matrix - reconstructed)^2)
})
png(paste(OUTPUT_DIR, "mse_2-40_2.png", sep = "/"), height = 640, width = 640)
plot(2:(length(mse_list) + 1), mse_list,
     xlab = "NMF Rank", ylab = "Mean Squared Error", pch = 19)
dev.off()

coph_corr_list = lapply(res_list, function(model) {
  consensus_matrix = consensus(model)
  coph_corr = cophcor(consensus_matrix)
})
png(paste(OUTPUT_DIR, "coph_2-40_2.png", sep = "/"), height = 640, width = 640)
plot(2:(length(mse_list) + 1), coph_corr_list,
     xlab = "NMF Rank", ylab = "Cophenetic Correlation", pch = 19)
dev.off()

# Compute a basis matrix for the rank chosen.

res = res_list[[RANK - 1]]
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
}
all_genes = unlist(module_groups)

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
  enrichr_groups[[combined_term]] = goi

  ontology_groups[i, "msigdb_term"] = msigDB_terms[1]
  ontology_groups[i, "go_term"] = GO_terms[1]
  ontology_groups[i, "n_genes"] = length(goi)
  ontology_groups[i, "msigdb_genes"] = enrichROut$MSigDB_Hallmark_2020$Genes[1]
  ontology_groups[i, "go_term_genes"] = enrichROut$GO_Biological_Process_2021$Genes[1]
  ontology_groups[i, "genes"] = paste(goi, collapse = "|")
  Sys.sleep(3)
}
write.table(ontology_groups, paste(OUTPUT_DIR, "ontology_groups_rank20_MaxModules_2.txt", sep = "/"), sep = "\t")
ontology_groups = read.table("ontology_groups_rank25_MaxModules.txt", sep = "\t")

list_names = paste(ontology_groups$msigdb_term, ontology_groups$go_term, sep = "|")
list_items <- lapply(ontology_groups$genes, function(x) strsplit(x, split = "\\|")[[1]])
enrichr_groups = setNames(list_items, list_names)

# Generate the log fold change matrix. We keep only those rows
# that qualify and that have genes that are in the gene_names

data.list = list()
for (directory in PATHS_TO_LFC) {
  file_names = list.files(directory)
  for (i in file_names) {
    split_string = strsplit(i, "_")[[1]]
    perturb = paste(split_string[2], split_string[1], split_string[3], sep = "_")
    df = read.table(paste(directory, "/", i, sep = ""))
    data.list[[perturb]] = df
  }
}
data.list.noRT = data.list[colnames(avg_matrix)]

nmfMat <- ldply(lapply(data.list.noRT,function(x) x[,c("feature","log_fc")]),data.frame)
nmfMat <- reshape(nmfMat, idvar = "feature", v.names = "log_fc",timevar = ".id", direction = "wide")
nmfMat <- nmfMat[!is.na(nmfMat$feature),]
row.names(nmfMat) <- nmfMat$feature
colnames(nmfMat) <- gsub("log_fc.","",colnames(nmfMat))
nmfMat <- nmfMat[,-1]
nmfMat[is.na(nmfMat)] <- 0
nmfMat <- nmfMat[, !colnames(nmfMat) %in% NT_COLS]
nmfMatsig <- nmfMat[intersect(all_genes,row.names(nmfMat)),]
nmfMatsig = as.matrix(nmfMatsig)

# Plot a context separated heatmap with annotations.

phenoTbl <- read.table(PHENOTBL,sep='\t',header=TRUE,row.names=1)
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
                       sgRNA = sapply(strsplit(colnames(nmfMatsig),"_"), `[`, 1),
                       source = sapply(strsplit(colnames(nmfMatsig),"_"), `[`, 2),
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

go_categories = read.table(GO_CATEGORIES, sep = ",")
colnames(go_categories) = c("Gene", "GO")
rownames(go_categories) = go_categories$Gene
tblAnnot$go_term = go_categories[tblAnnot$sgRNA, "GO"]
rownames(tblAnnot) = paste(tblAnnot$source, tblAnnot$sgRNA, sep = "_")

write.table(tblAnnot, paste(OUTPUT_DIR, "/", "noRT_annotation_table.txt", sep = ""), sep = "\t")

# Cluster the log fold changes into gene modules

module_means = lapply(enrichr_groups, function(g) {
  colMeans(nmfMatsig[g, , drop = FALSE])
})
expr_matrix = do.call(rbind, module_means)
colnames(expr_matrix) = sapply(strsplit(colnames(expr_matrix), "_"), function(x) paste(x[2], x[1], sep="_"))
write.table(expr_matrix, paste(OUTPUT_DIR, "output_matrix_lognormalized_lfc_rank25_MaxModules_full.txt", sep = "/"), sep = "\t")

# Plot the heatmap

poi = colnames(expr_matrix)
ht = Heatmap(expr_matrix,
        name = "Log Fold Change Matrix Clustered on Rank-25 NMF - FULL",
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        top_annotation = HeatmapAnnotation(gamma = anno_barplot(as.numeric(tblAnnot[poi,"gamma"]), bar_width = 0.9),
                                           tau = anno_barplot(as.numeric(tblAnnot[poi,"tau"]), bar_width = 0.9),
                                           rho = anno_barplot(as.numeric(tblAnnot[poi,"rho"]), bar_width = 0.9),
                                           GO = tblAnnot[poi,"go_term"],
                                           treatment = tblAnnot[poi,"cond"],
                                           context = tblAnnot[poi,"source"],
                                           coverage_log10 = log(as.numeric(tblAnnot[poi,"coverage"])),
                                           col = list(treatment = c("noRT" = "#2166AC", "RT" = "#B2182B","RTvsnoRT" = "green4"),
                                                      context = c("invitro" = "#D53E4F", "preinf" = "#FC8D59","CED" = "#FEE08B"),
                                                      GO = c("Acetyltransferase" = "#9E0142", "Ankyrin Repeat" = "#D53E4F",
                                                             "Cyclin-dependent kinase inhibitor" = "#F46D43", "Signaling" = "#FDAE61",
                                                             "DNA damage response" = "#FEE08B", "DNA Replication" = "#FFFFBF",
                                                             "LIM Domain" = "#E6F598", "Metabolism" = "#ABDDA4","Mitochondria" = "#66C2A5",
                                                             "Mitosis" = "#3288BD", "Proteasome" = "#5E4FA2", "Ras pathway" = "#1B9E77",
                                                             "Receptor" = "#D95F02", "RNA polymerase" = "#7570B3", "Telomere" = "#E7298A","Transcription" = "#66A61E","Translation" = "#E6AB02",
                                                             "Non-targeting" = "gray50"),
                                                      coverage_log10 = colorRamp2(quantile(log(as.numeric(tblAnnot[poi,"coverage"]),10), probs = seq(0, 1, by = 0.25)), viridis(5))
                                           )),
        column_split = factor(tblAnnot[poi,"source"], levels = c("invitro","preinf","CED")),
        cluster_column_slices = FALSE,
        cluster_columns = TRUE,
        width = ncol(expr_matrix)*unit(4, "mm"),
        height = nrow(expr_matrix)*unit(4, "mm"),
)
pdf(paste(OUTPUT_DIR, "output_matrix_lognormalized_lfc_rank25_MaxModules_range1_full_2.pdf", sep = "/"),
    width=35,
    height=19,
    useDingbats=FALSE)
draw(ht, heatmap_legend_side = "top", annotation_legend_side = "top")
dev.off()

print("Finished executing!")

# EOF
