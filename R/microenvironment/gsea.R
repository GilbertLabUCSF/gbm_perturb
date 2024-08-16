# Microenvironment Gene Set Enrichment Analysis (GSEA)
# 
# This script contains code for generating and plotting GSEA results from
# microenvironment DESeq output stored in a directory called
# data/microenvironment/deseq.
# 
# In an attempt to provide the greatest possible degree of transparency and to
# fit with the example DESeq2 data provided, this code retains a significant 
# number of data wrangling steps.
# 
# Inputs:
# - DESeq output directory containing a CSV file per perturbation x cell type. 
#   We assume that the directory contains many different cell types, and show
#   one way to group the clustering according to cell type. This example code
#   assumes a format that is the same as the example CSV provided in
#   ~/example_data/microenvironment/deseq
# - An annotation table that contains manual annotations for which
#   perturbations fit in which GO categories.
# - An annotation table that contains screening phenotype data for each
#   perturbation.
# 
# Inputs:
# - DESeq output directory containing a CSV file per perturbation
# 
# Outputs:
# 
# For both all cell types grouped together and individual cell types:
# - Table of genes that don't convert to ENSEMBL
# - GSEA R object
# - Heatmap of normalized enrichment scores
# - Bubble plot of normalized enrichment scores

## -----------------------------------------------------------------------------------------------------------------------------------------------
# Import libraries:

library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(knitr)
library(msigdbr)
library(fgsea)
library(dplyr)
library(data.table)
library(ggplot2)
library(grid)
library(ComplexHeatmap)
library(colorRamp2)
library(tidyr)
library(patchwork)
library(ggtree)
library(gridExtra)
library(ape)
library(RColorBrewer)
library(reactome.db)
library(annotables)
library(biomaRt)
set.seed(5220)


## -----------------------------------------------------------------------------------------------------------------------------------------------
# Define constants:

INPUT.DIR.BASE <- "output/microenvironment/deseq"
INPUT.DIRS <- list.dirs(INPUT.DIR.BASE)[-1]
OUTPUT.DIR <- "output/gsea/Combined"


## -----------------------------------------------------------------------------------------------------------------------------------------------
# Import the DESeq output:

data.list.symbols <- list()
for (dir in INPUT.DIRS) {
  file_names <- list.files(dir)
  split.string <- strsplit(dir, "/")[[1]]
  cell.type <- tail(split.string, 1)
  for (i in file_names) {
    name <- gsub(".csv", "", i)
    name <- paste(cell.type, name, sep = "_")
    df <- read.table(paste(dir, "/", i, sep = ""))
    data.list.symbols[[name]] <- df
  }
}


## -----------------------------------------------------------------------------------------------------------------------------------------------
# Prepare the ENSEMBL gene ID list, get rid of low expression, and order:

# Prepare the ENSEMBL gene ID list
mapping.df <- annotables::grcm38
missing.genes <- list()
data.list.ensembl <- list()
for (name in names(data.list.symbols)) {
  deseq_table <- data.list.symbols[[name]]
  mapping.df <- mapping.df[!duplicated(mapping.df$symbol), ]
  mapping.table <- mapping.df[mapping.df$symbol %in% deseq_table$feature, ]
  mapping.table <- as.data.frame(mapping.table)
  rownames(mapping.table) <- mapping.table$symbol
  deseq_table$ensembl <- mapping.table[deseq_table$feature, "ensgene"]
  deseq_table.missing <- deseq_table[is.na(deseq_table$ensembl), ]
  missing.genes[[name]] <- deseq_table.missing$feature
  data.list.ensembl[[name]] <- deseq_table[!is.na(deseq_table$ensembl), ]
}
write.table(missing.genes[[1]], file = paste(OUTPUT.DIR, "missing_genes.txt", sep = "/"), row.names = FALSE, quote = FALSE)

# Get rid of low expression
processed.deseq.tables <- list()
for (name in names(data.list.ensembl)) {
  deseq_table <- data.list.ensembl[[name]]
  deseq_table <- deseq_table[!is.na(deseq_table$pvalue), ]
  deseq_table <- deseq_table[!is.na(deseq_table$padj), ]
  processed.deseq.tables[[name]] <- deseq_table
}
data.list <- processed.deseq.tables

# Create ordered lists
ranked_gene_set_list <- list()
for (name in names(processed.deseq.tables)) {
  deseq_table <- processed.deseq.tables[[name]]
  # Group by unique Ensembl ID and average
  deseq_table <- deseq_table %>%
    group_by(ensembl) %>%
    summarize(log_fc = mean(log_fc))
  deseq_table <- deseq_table %>% mutate(rank = rank(log_fc,  ties.method = "random"))
  deseq_table <- deseq_table[order(-deseq_table$rank),] # Rank deterministically
  
  # Generate a ranked list
  gene_list <- deseq_table$log_fc
  gene_names <- deseq_table$ensembl
  names(gene_list) <- gene_names
  ranked_gene_set_list[[name]] <- gene_list
}


## -----------------------------------------------------------------------------------------------------------------------------------------------
# Run GSEA

msigdbr_df <- msigdbr(species = "Mus musculus", category = "H")
msigdbr_list = split(x = msigdbr_df$ensembl_gene, f = msigdbr_df$gs_name)

fgsea.output.list <- list()
for (name in names(ranked_gene_set_list)) {
  gsea_results <- fgsea(pathways = msigdbr_list, stats = ranked_gene_set_list[[name]],
                        maxSize = 500, eps = 0.0, nproc = 1) # allow arbitrarily low p-values
  fgsea.output.list[[name]] <- gsea_results
}
saveRDS(fgsea.output.list, paste(OUTPUT.DIR, "fgsea_ensembl_ontology_list.rds", sep = "/"))


## -----------------------------------------------------------------------------------------------------------------------------------------------
# Generate a heatmap showing NES across perturbation x conditions

extracted_cols <- lapply(names(fgsea.output.list), function(name) {
  col <- fgsea.output.list[[name]]$NES
  names(col) <- fgsea.output.list[[name]]$pathway
  col
})
combined.df <- do.call(cbind, extracted_cols)
colnames(combined.df) <- names(fgsea.output.list)

ht = Heatmap(combined.df,
             name = "Enrichment Scores across Perturbations and Conditions",
             col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
             cluster_column_slices = FALSE,
             cluster_columns = FALSE,
             width = ncol(combined.df)*unit(4, "mm"),
             height = nrow(combined.df)*unit(4, "mm"),
)
pdf(paste(OUTPUT.DIR, "enrichment_combined_matrix_alphabetical.pdf", sep = "/"),
    width=35,
    height=19,
    useDingbats=FALSE)
draw(ht, heatmap_legend_side = "top", annotation_legend_side = "top")
dev.off()

## -----------------------------------------------------------------------------------------------------------------------------------------------
# Example bubble plot: cluster within cell types

combined.df <- bind_rows(fgsea.output.list, .id = "Perturbations")
combined.df$pathway <- gsub("HALLMARK_", "", combined.df$pathway)
combined.df$pathway <- gsub(" ", "", combined.df$pathway)
combined.df$log10pvalue <- -log10(combined.df$padj)
combined.df$Perturbations <- gsub("_sgNegCtrl_3", "", combined.df$Perturbations)
unfiltered.df <- combined.df

combined.df <- combined.df %>%
  group_by(pathway) %>%
  filter(any(padj < 0.05)) %>%
  ungroup()

es.df <- combined.df %>%
  dplyr::select(-pval, -padj, -log2err, -ES, -size, -leadingEdge, -log10pvalue)

es.df <- es.df %>%
  group_by(pathway, Perturbations) %>%
  pivot_wider(names_from = Perturbations, values_from = NES) %>%
  group_by(pathway)

es.df <- as.data.frame(es.df)
rownames(es.df) <- es.df$pathway
es.df <- es.df %>% dplyr::select(-pathway)

row_dist <- dist(es.df, method = "euclidean")
row_hclust <- hclust(row_dist, method = "complete")
row_order <- order.dendrogram(as.dendrogram(row_hclust))

# Group the cells within the cell types
ordering <- c()
col.hclusts <- list()
for (cell.type in c("Astrocytes", "Macrophages", "Microglia", "OPCs", "Oligodendrocytes")) {
  cols <- grep(sprintf("^%s_", cell.type), colnames(es.df), value = TRUE)
  col.df <- es.df[cols]
  col.dist <- dist(t(col.df), method = "euclidean")
  col.hclust <- hclust(col.dist, method = "complete")
  col.hclusts[[cell.type]] <- col.hclust
  col.order <- order.dendrogram(as.dendrogram(col.hclust))
  ordering <- c(ordering, cols[col.order])
}
es.df <- es.df[, ordering]

df.ordered <- combined.df %>%
  arrange(match(pathway, rownames(es.df)[row_order]), match(Perturbations, colnames(es.df)))

df.ordered$pathway <- factor(df.ordered$pathway, levels = rownames(es.df)[row_order])
df.ordered$Perturbations <- factor(df.ordered$Perturbations, levels = colnames(es.df))
unique_perturbations <- unique(df.ordered$Perturbations)

# Add in cell type fill bar

df.ordered$cell.type <- sapply(strsplit(as.character(df.ordered$Perturbations), "_"), `[`, 1)

rect.df <- data.frame(
  xmin = unique_perturbations,
  xmax = unique_perturbations,
  ymin = rep(max(as.numeric(df.ordered$pathway)) + 1, length(unique_perturbations)),
  ymax = rep(max(as.numeric(df.ordered$pathway)) + 2, length(unique_perturbations)),
  cell.type = df.ordered$cell.type[match(unique_perturbations, df.ordered$Perturbations)]
)
rect.df$numeric_xmin <- as.numeric(as.factor(rect.df$xmin)) - 0.5
rect.df$numeric_xmax <- as.numeric(as.factor(rect.df$xmax)) + 0.5

# Add in perturb fill bar

df.ordered$perturb <- sapply(strsplit(as.character(df.ordered$Perturbations), "_"), `[`, 2)

perturb.df <- data.frame(
  xmin = unique_perturbations,
  xmax = unique_perturbations,
  ymin = rep(max(as.numeric(df.ordered$pathway)) + 2, length(unique_perturbations)),
  ymax = rep(max(as.numeric(df.ordered$pathway)) + 3, length(unique_perturbations)),
  perturb = df.ordered$perturb[match(unique_perturbations, df.ordered$Perturbations)]
)
perturb.df$numeric_xmin <- as.numeric(as.factor(perturb.df$xmin)) - 0.5
perturb.df$numeric_xmax <- as.numeric(as.factor(perturb.df$xmax)) + 0.5

# Add significance layer
df.ordered$Significance <- ifelse(df.ordered$padj < 0.05, "adjp < 0.05", "NA")

# Clean up names
df.ordered$Pathways <- df.ordered$pathway
levels(df.ordered$Pathways) <- gsub("_", " ", levels(df.ordered$Pathways))

## -----------------------------------------------------------------------------------------------------------------------------------------------
# Make the plot:

color.map <- list(
  `OPCs` = "#A6CEE3",
  `Oligodendrocytes` = "#1F78B4",
  `Microglia` = "#B2DF8A",
  `Macrophages` = "#33A02C",
  `Astrocytes` = "#FB9A99",
  `Lyz2` = "#E31A1C",
  `Apoe` = "#FDBF6F",
  `Cd44` = "#FF7F00",
  `C1qa` = "#CAB2D6",
  `Ptprc` = "#6A3D9A",
  `Cd74` = "#FFFF99",
  `Prkdc` = "brown"
)
bubble.plot <- ggplot(df.ordered,
                      aes(x = Perturbations)) +
  
  # Add the bubble and outline significance layers
  geom_point(aes(y = Pathways, size = log10pvalue, color = NES), alpha = 1.0) +
  geom_point(data = subset(df.ordered, padj < 0.05),
             aes(y = Pathways, size = log10pvalue, shape = 'adjp < 0.05'),
             color = "black", fill = NA, alpha = 0.8) +
  geom_point(data = subset(df.ordered, is.na(pval) | is.na(padj) | is.na(log10pvalue)), 
             aes(y = Pathways), 
             color = "black", size = 1, alpha = 0.6) +
  scale_size(range = c(2, 7)) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "grey50"
  ) +
  scale_shape_manual(name = "Significance",
                     values = c(`adjp < 0.05` = 21),
                     guide = guide_legend(override.aes = list(color = "black", fill = NA))) +
  
  # Add the cell type bar layer
  geom_rect(data = rect.df, aes(xmin = numeric_xmin, xmax = numeric_xmax, ymin = ymin, ymax = ymax, fill = cell.type), inherit.aes = FALSE) +
  
  # Add the perturb bar layer
  geom_rect(data = perturb.df, aes(xmin = numeric_xmin, xmax = numeric_xmax, ymin = ymin,
                                   ymax = ymax, fill = perturb), inherit.aes = FALSE) +
  
  # Color things in:
  scale_fill_manual(values = color.map) +
  
  # Labels
  labs(
    color = "Effect size",
    size = "-Log10 padj"
  ) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.height = unit(0.3, 'cm'),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank())

perturb.dendrogram.OPCs <- ggtree(as.phylo(col.hclusts[["OPCs"]])) + theme_tree() + coord_flip() + scale_x_reverse()
perturb.dendrogram.Macrophages <- ggtree(as.phylo(col.hclusts[["Macrophages"]])) + theme_tree() + coord_flip() + scale_x_reverse()
perturb.dendrogram.Astrocytes <- ggtree(as.phylo(col.hclusts[["Astrocytes"]])) + theme_tree() + coord_flip() + scale_x_reverse()
perturb.dendrogram.Oligodendrocytes <- ggtree(as.phylo(col.hclusts[["Oligodendrocytes"]])) + theme_tree() + coord_flip() + scale_x_reverse()
perturb.dendrogram.Microglia <- ggtree(as.phylo(col.hclusts[["Microglia"]])) + theme_tree() + coord_flip() + scale_x_reverse()
module.dendrogram <- ggtree(as.phylo(row_hclust)) + theme_tree()
col.dendrograms <- (
  perturb.dendrogram.Astrocytes | perturb.dendrogram.Macrophages | perturb.dendrogram.Microglia | perturb.dendrogram.OPCs | perturb.dendrogram.Oligodendrocytes
)

vertical.plot <- col.dendrograms / bubble.plot + plot_layout(heights = c(0.5, 3))
dendrogram.vertical.plot <- plot_spacer() / module.dendrogram + plot_layout(heights = c(0.01, 1))
combined.plot <- dendrogram.vertical.plot | vertical.plot
combined.plot <- combined.plot + plot_layout(widths = c(0.2, 0.8)) + plot_annotation(
  title = "Enrichment of Hallmark pathways clustered within cell types",
  subtitle = "Prkdc excluded"
)
ggsave(paste(OUTPUT.DIR, "enrichment_bubble_plot_fgsea_clustered_within_cells.pdf", sep = "/"),
       combined.plot, height = 10, width = 9, device = pdf)

## -----------------------------------------------------------------------------------------------------------------------------------------------
# Overlay LFCs on the existing plot

processed_deseq_tables <- list()
for (name in names(data.list.ensembl)) {
  deseq_table <- data.list.ensembl[[name]]
  deseq_table <- deseq_table[!is.na(deseq_table$pvalue), ]
  deseq_table <- deseq_table[!is.na(deseq_table$padj), ]
  deseq_table <- deseq_table[deseq_table$padj < 0.05, ]
  deseq_table <- deseq_table[abs(deseq_table$log_fc) > 0.1, ]
  processed_deseq_tables[[name]] <- deseq_table
}
data_list <- processed_deseq_tables

# Get the union of all DE genes across perturbations x cell types
union_de_genes <- c()
for (name in names(data_list)) {
  deseq_table <- data_list[[name]]
  ensembl_ids <- deseq_table$ensembl
  union_de_genes <- unique(c(union_de_genes, ensembl_ids))
}

# Filter to DE genes
filter_to_de_genes <- function(deseq_table) {
  deseq_table <- subset(deseq_table, ensembl %in% union_de_genes)
}
de_gene_data_list <- lapply(data.list.ensembl, filter_to_de_genes)

msigdbr_df <- msigdbr(species = "Mus musculus", category = "H")
msigdbr_list = split(x = msigdbr_df$ensembl_gene, f = msigdbr_df$gs_name)

# For each perturbation x cell type x pathway, build a DESeq2 dataframe with
# the intersection of genes for that perturbation x cell type x pathway
celltype_perturb_pathway_dfs <- list()
for (celltype_perturb in names(de_gene_data_list)) {
  deseq_df <- de_gene_data_list[[celltype_perturb]]
  for (pathway in names(msigdbr_list)) {
    pathway_intersected_df <- subset(deseq_df, ensembl %in% msigdbr_list[[pathway]])
    celltype_perturb_pathway <- paste(celltype_perturb, pathway, sep = "_")
    celltype_perturb_pathway_dfs[[celltype_perturb_pathway]] <- pathway_intersected_df
  }
}

# Find means for each dataframe
calculate_average <- function(deseq_table) {
  deseq_table$log_fc[is.na(deseq_table$pvalue) & is.na(deseq_table$padj)] <- 0
  return(mean(deseq_table$log_fc))
}
celltype_perturb_pathway_means <- lapply(celltype_perturb_pathway_dfs, calculate_average)

# Map the means to the df.ordered table generated above
df.ordered$celltype_perturb_pathway <- 
  paste(df.ordered$Perturbations, df.ordered$Pathways, sep = "_")
df.ordered$celltype_perturb_pathway <- 
  gsub(" ", "_", df.ordered$celltype_perturb_pathway)

names(celltype_perturb_pathway_means) <- gsub("sgNegCtrl_3_HALLMARK_", "", 
                                              names(celltype_perturb_pathway_means))

df.ordered$celltype_perturb_pathway_means <- 
  unlist(celltype_perturb_pathway_means[df.ordered$celltype_perturb_pathway])

## -----------------------------------------------------------------------------------------------------------------------------------------------
# Make the plot on the cell-type grouped clusters

color.map <- list(
  `OPCs` = "#A6CEE3",
  `Oligodendrocytes` = "#1F78B4",
  `Microglia` = "#B2DF8A",
  `Macrophages` = "#33A02C",
  `Astrocytes` = "#FB9A99",
  `Lyz2` = "#E31A1C",
  `Apoe` = "#FDBF6F",
  `Cd44` = "#FF7F00",
  `C1qa` = "#CAB2D6",
  `Ptprc` = "#6A3D9A",
  `Cd74` = "#FFFF99",
  `Prkdc` = "brown"
)
bubble.plot <- ggplot(df.ordered,
                      aes(x = Perturbations)) +
  
  # Add the bubble and outline significance layers
  geom_point(aes(y = Pathways, size = log10pvalue, color = celltype_perturb_pathway_means), alpha = 1.0) +
  geom_point(data = subset(df.ordered, padj < 0.05),
             aes(y = Pathways, size = log10pvalue, shape = 'adjp < 0.05'),
             color = "black", fill = NA, alpha = 0.8) +
  geom_point(data = subset(df.ordered, is.na(pval) | is.na(padj) | is.na(log10pvalue)), 
             aes(y = Pathways), 
             color = "black", size = 1, alpha = 0.6) +
  scale_size(range = c(2, 7)) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "grey50",
    limits = c(-0.75, 0.75),
    oob = scales::squish
  ) +
  scale_shape_manual(name = "Significance",
                     values = c(`adjp < 0.05` = 21),
                     guide = guide_legend(override.aes = list(color = "black", fill = NA))) +
  
  # Add the cell type bar layer
  geom_rect(data = rect.df, aes(xmin = numeric_xmin, xmax = numeric_xmax, ymin = ymin, ymax = ymax, fill = cell.type), inherit.aes = FALSE) +
  
  # Add the perturb bar layer
  geom_rect(data = perturb.df, aes(xmin = numeric_xmin, xmax = numeric_xmax, ymin = ymin,
                                   ymax = ymax, fill = perturb), inherit.aes = FALSE) +
  
  # Color things in:
  scale_fill_manual(values = color.map) +
  
  # Labels
  labs(
    color = "Log fc",
    size = "-Log10 padj"
  ) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.height = unit(0.3, 'cm'),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank())

perturb.dendrogram.OPCs <- ggtree(as.phylo(col.hclusts[["OPCs"]])) + theme_tree() + coord_flip() + scale_x_reverse()
perturb.dendrogram.Macrophages <- ggtree(as.phylo(col.hclusts[["Macrophages"]])) + theme_tree() + coord_flip() + scale_x_reverse()
perturb.dendrogram.Astrocytes <- ggtree(as.phylo(col.hclusts[["Astrocytes"]])) + theme_tree() + coord_flip() + scale_x_reverse()
perturb.dendrogram.Oligodendrocytes <- ggtree(as.phylo(col.hclusts[["Oligodendrocytes"]])) + theme_tree() + coord_flip() + scale_x_reverse()
perturb.dendrogram.Microglia <- ggtree(as.phylo(col.hclusts[["Microglia"]])) + theme_tree() + coord_flip() + scale_x_reverse()
module.dendrogram <- ggtree(as.phylo(row_hclust)) + theme_tree()
col.dendrograms <- (
  perturb.dendrogram.Astrocytes | perturb.dendrogram.Macrophages | perturb.dendrogram.Microglia | perturb.dendrogram.OPCs | perturb.dendrogram.Oligodendrocytes
)

vertical.plot <- col.dendrograms / bubble.plot + plot_layout(heights = c(0.5, 3))
dendrogram.vertical.plot <- plot_spacer() / module.dendrogram + plot_layout(heights = c(0.01, 1))
combined.plot <- dendrogram.vertical.plot | vertical.plot
combined.plot <- combined.plot + plot_layout(widths = c(0.2, 0.8)) + plot_annotation(
  title = "Enrichment of Hallmark pathways clustered within cell types",
  subtitle = "Prkdc excluded"
)
ggsave(paste(OUTPUT.DIR, "enrichment_bubble_plot_fgsea_clustered_within_cells_logfc_overlay.pdf", sep = "/"),
       combined.plot, height = 10, width = 9, device = pdf)



