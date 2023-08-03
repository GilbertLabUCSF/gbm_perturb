# Run DESeq2 on a Seurat object according to our filters and configurations
# 
# Users should configure:
# - The intput Seurat object
# - The identities of the sorted cells
# - The output directory
# - The experimental context to generate output for. We default to "iv"
# - A random seed
# - Whether we normalize to non-targeting noRT or normalize to the condition
# 
# Output: deseq output per perturbation. By default, we use 3 cells per
# pseudoreplicate.

##############################################################################
##############################################################################
# Dependencies:

library(Seurat)
library(DElegate)
library(dplyr)
library(NMF)
set.seed(5220)

##############################################################################
##############################################################################
# Inputs:

PATH_TO_SEURAT_OBJECT = "/raleighlab/data1/liuj/gbm_perturb/analysis/SB28_micro4_20230619.Rds"
NT_GUIDE = "sgNegCtrl"
DESEQ_OUTPUT_DIR = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_sb28_clean_outputs/deseq/sb28_integrated_micro4_macrophagesOnly"
DEGENES_OUTPUT_DIR = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_sb28_clean_outputs/de_genes/SB28_integrated_integrated_micro4_macrophagesOnly"
NMF_OUTPUT_DIR = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_sb28_clean_outputs/nmf/macrophageOnly"
NMF_OUTPUT_FILENAME = "res_list2-40.rds"
SEED = 5220
MINIMUM_COVERAGE = 5
DESEQ = TRUE
DEGENES = TRUE
RUN_NMF = TRUE

print(paste("PATH_TO_SEURAT_OBJECT =", PATH_TO_SEURAT_OBJECT))
print(paste("NT_GUIDES =", NT_GUIDE))
print(paste("DESEQ_OUTPUT_DIR =", DESEQ_OUTPUT_DIR))
print(paste("SEED =", SEED))
print(paste("MINIMUM_COVERAGE =", MINIMUM_COVERAGE))

##############################################################################
##############################################################################
# Code

# Seurat object setup

data = readRDS(PATH_TO_SEURAT_OBJECT)
DefaultAssay(data) = "RNA"
data = NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 1e6, verbose = TRUE)

# Subset the data by sgRNA guide positive cells

data.context = subset(data, sgRNA_binary == TRUE)

# Screen for macrophages only

data.context = subset(data.context, new.cluster.ids == "Macrophage")

# Screen for having a guide that's not NA

perturbs = unique(data.context$GeneA)
perturbs <- perturbs[!is.na(perturbs)]
data.context = subset(data.context, GeneA %in% perturbs)

# Screen for cells part of a group with coverage of > 5 cells

sgRNACond_counts = table(data.context$GeneA)
data.context = subset(data.context, GeneA %in% 
                        names(sgRNACond_counts[sgRNACond_counts >= MINIMUM_COVERAGE]))

# DESEQ

# Define functions that help us to find differential genes

if (DESEQ) {
  find_deseq_differential_genes = function(data.obj, seed, group1, group2, group_column=NULL) {
    # Takes in a normalized Seurat object with metadata
    # indicating group1 and group2 and a group_column argument
    # indicating the metadata column to be used. Sets a seed, then
    # finds differential expressed genes and outputs a data frame.
    # Note that log fold comparisons will be output as log2(group1/group2)
    set.seed(seed)
    
    # Set futures to greater than max capacity; you may need to tweak this
    options(future.globals.maxSize = 5000 * 1024^2)
    
    df = findDE(object = data.obj, group_column = group_column,
                compare = c(group1, group2), method = 'deseq')
    return(df)
  }
  
  build_filename = function(group1, group2) {
    filename = paste(group1, group2, sep = "_")
    directory = paste(DESEQ_OUTPUT_DIR, "/", sep = "")
    extension = ".csv"
    return(paste(directory, filename, extension, sep = ""))
  }
  
  # Find differential genes, normalizing to sgNegCtrl
  
  for (perturb in perturbs) {
    print(sprintf("Calculating diff genes for %s", perturb))
    found = FALSE
    tryCatch({
      df = find_deseq_differential_genes(data.context, SEED, perturb, NT_GUIDE, group_column = "GeneA")
      found = TRUE
    }, error = function(err) {
      print(paste("Failed to get model because of ", err, "in", perturb))
    })
    if (found) {
      filename = build_filename(perturb, NT_GUIDE)
      write.table(df, filename)
    }
  }
  
  print("Finished running DESeq!")
}

# Filter for DE genes

if (DEGENES) {
  data.list.all = list()
  for (dir in DESEQ_OUTPUT_DIR) {
    file_names = list.files(dir)
    for (i in file_names) {
      name = gsub(".csv", "", i)
      df = read.table(paste(dir, "/", i, sep = ""))
      data.list.all[[name]] = df
    }
  }
  
  deGenesAll.padj = unique(unlist(lapply(data.list.all, function(x) {
    x = x[!is.na(x$log_fc), ]
    x = x[!is.na(x$pvalue), ]
    x = x[!is.na(x$padj), ]
    x = x[x$pvalue < 0.05, ]
    x = x[abs(x$log_fc) > 0.1, ]
    x$feature
  })))
  
  deMat <- ldply(lapply(data.list.all,function(x) x[,c("feature","log_fc")]),data.frame)
  deMat <- reshape(deMat, idvar = "feature", v.names = "log_fc",timevar = ".id", direction = "wide")
  deMat <- deMat[!is.na(deMat$feature),]
  row.names(deMat) <- deMat$feature
  colnames(deMat) <- gsub("log_fc.","",colnames(deMat))
  deMat <- deMat[,-1]
  deMat[is.na(deMat)] <- 0
  deMat <- deMat[, !colnames(deMat) %in% NT_COLS]
  deMatsig <- deMat[intersect(deGenesAll.padj,row.names(deMat)),]
  
  write.table(deMat, paste(DEGENES_OUTPUT_DIR, "/deMat_p05_lfc01", ".txt", sep = ""))
  write.table(deMatsig, paste(DEGENES_OUTPUT_DIR, "/deMatSig_p05_lfc01", ".txt", sep = ""))
}


print("Finished generating DE genes!")

# Generate average matrix

deMatSig = read.table(paste(DEGENES_OUTPUT_DIR, "deMatSig_p05_lfc01.txt", sep = "/"))
deMatSig = na.omit(deMatSig)
de_genes = rownames(deMatSig)
data = subset(data, features = de_genes)

# Create and then pseudobulk by perturb

data$perturb = data$GeneA
counts_data = as.matrix(GetAssayData(data))
metadata = data@meta.data
expr_df = as.data.frame(counts_data)
expr_df$gene = rownames(expr_df)
expr_long = tidyr::pivot_longer(expr_df, -gene, 
                                names_to = "cell_barcode",
                                values_to = "expression")
expr_long$perturb = metadata[expr_long$cell_barcode, "perturb"]
avg_data = expr_long %>%
  dplyr::group_by(gene, perturb) %>%
  dplyr::summarise(mean_expression = mean(expression, na.rm = TRUE))
avg_matrix = tidyr::pivot_wider(avg_data, names_from = perturb, values_from = mean_expression)
avg_matrix = as.data.frame(avg_matrix)
rownames(avg_matrix) = avg_matrix$gene
avg_matrix = as.matrix(avg_matrix[, -1])
avg_matrix = avg_matrix[apply(avg_matrix, 1, function(row) any(row != 0)), ]

# Filter by those perturbs that we have deseq output for. This means they passed
# a coverage filter.

perturb_list = c()
for (directory in DESEQ_OUTPUT_DIR) {
  file_names = list.files(directory)
  for (i in file_names) {
    split_string = strsplit(i, "_")[[1]]
    perturb = split_string[1]
    perturb_list = c(perturb_list, perturb)
  }
}

avg_matrix = avg_matrix[,perturb_list]
avg_matrix = avg_matrix + 0.0001

# Run NMF

if (RUN_NMF) {
  print("Beginning NMF")
  rank_range = 2:40
  
  run_nmf = function(r) {
    model = NMF::nmf(avg_matrix, r, method = "brunet", nrun = 60)
  }
  res_list = mclapply(rank_range, run_nmf, mc.cores = 50)
  saveRDS(res_list, paste(NMF_OUTPUT_DIR, OUTPUT_FILE_NAME, sep = "/"))
}


mse_list = lapply(res_list, function(model) {
  reconstructed = fitted(model)
  mse = mean((avg_matrix - reconstructed)^2)
})
png(paste(NMF_OUTPUT_DIR, "mse_2-40.png", sep = "/"), height = 640, width = 640)
plot(1:length(mse_list), mse_list,
     xlab = "NMF Rank", ylab = "Mean Squared Error", pch = 19)
dev.off()

coph_corr_list = lapply(res_list, function(model) {
  consensus_matrix = consensus(model)
  coph_corr = cophcor(consensus_matrix)
})
png(paste(NMF_OUTPUT_DIR, "coph_2-40.png", sep = "/"), height = 640, width = 640)
plot(1:length(coph_corr_list), coph_corr_list,
     xlab = "NMF Rank", ylab = "Cophenetic Correlation", pch = 19)
dev.off()

# Compute a basis matrix for the rank(s) chosen.

RANKS = c(20, 25, 30)
for (RANK in RANKS) {
  print(sprintf("Calculating for Rank %s", RANK))
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
  write.table(ontology_groups, paste(NMF_OUTPUT_DIR, sprintf("ontology_groups_rank%s_MaxModules.txt", RANK), sep = "/"), sep = "\t")
  
  # Generate the log fold change matrix. We keep only those rows
  # that qualify and that have genes that are in the gene_names
  
  data.list = list()
  file_names = list.files(DESEQ_OUTPUT_DIR)
  for (i in file_names) {
    name = gsub(".csv", "", i)
    df = read.table(paste(DESEQ_OUTPUT_DIR, "/", i, sep = ""))
    data.list[[name]] = df
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
  
  # Plot the heatmap
  
  poi = colnames(expr_matrix)
  ht = Heatmap(expr_matrix, 
               name = sprintf("Log Fold Change Matrix Clustered on Rank-%s NMF", RANK),
               col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
               cluster_column_slices = FALSE,
               cluster_columns = TRUE,
               width = ncol(expr_matrix)*unit(4, "mm"),
               height = nrow(expr_matrix)*unit(4, "mm"),
  )
  pdf(paste(NMF_OUTPUT_DIR, sprintf("output_matrix_lognormalized_lfc_rank%s_MaxModules_range3.pdf", RANK), sep = "/"),
      width=35,
      height=19,
      useDingbats=FALSE)
  draw(ht, heatmap_legend_side = "top", annotation_legend_side = "top")
  dev.off()
}

