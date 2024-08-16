# GL261 Gene Set Enrichment Analysis (GSEA)
# 
# This script contains code for generating and plotting GSEA results on 
# GL261 DESeq output stored in a directory called data/malignant/deseq/invitro
# 
# In an attempt to provide the greatest possible degree of transparency and to
# fit with the example DESeq2 data provided, this code retains a significant 
# number of data wrangling steps.
# 
# Inputs:
# - DESeq output directory containing a CSV file per perturbation. We assume
#   that the directory contains both irradiated and non-irradiated cells, and
#   that we want to treat them as separate populations, but plot them together.
#   This example code assumes a "condNormalized" normalization scheme, and fits
#   with the example CSV provided in ~/example_data/malignant/deseq
# - An annotation table that contains manual annotations for which
#   perturbations fit in which GO categories.
# - An annotation table that contains screening phenotype data for each
#   perturbation.
# 
# Outputs: 
# - Table of genes that fail to convert to ENSEMBL
# - An R object that contains GSEA results for each perturbation x radiation
#   combination
# - Heatmap of normalized enrichment scores
# - Bubble plot of normalized enrichment scores
# - Counts of pathways in perturbations and perturbations in pathways

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

CONTEXT <- "invitro"
CONDITION <- "Combined"
INPUT_DIRS <- c(sprintf("output/malignant/deseq/%s", CONTEXT))
OUTPUT_DIR <- sprintf("output/malignant/gsea/%s/%s", CONTEXT, CONDITION)


## -----------------------------------------------------------------------------------------------------------------------------------------------
# Import the DESeq output:

data_list <- list()
for (dir in INPUT_DIRS) {
  file_names <- list.files(dir)
  for (i in file_names) {
    name <- gsub(".csv", "", i)
    df <- read.table(paste(dir, "/", i, sep = ""))
    data_list[[name]] <- df
  }
}


## -----------------------------------------------------------------------------------------------------------------------------------------------
# To plot separately, you could separate out the RT and noRT indices here.

if (CONDITION != "Combined") {
  noRT_indices <- grepl("_noRT_non-targeting_noRT", names(data_list))
  RT_indices <- grepl("_RT_non-targeting_RT", names(data_list))
  noRT_list <- data_list[noRT_indices]
  RT_list <- data_list[RT_indices]
  if (CONDITION == "RT") {
    data_list_symbols <- RT_list
  } else {
    data_list_symbols <- noRT_list
  }
} else {
  overall_list <- data_list
}

## -----------------------------------------------------------------------------------------------------------------------------------------------
# Prepare an ENSEMBL gene ID data list by converting using the `annotables` 
# package. Record the genes that end up going missing from ENSEMBL conversion.

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
write.table(missing.genes[[1]], file = paste(OUTPUT_DIR, "missing_genes.txt", sep = "/"), row.names = FALSE, quote = FALSE)


## -----------------------------------------------------------------------------------------------------------------------------------------------
# Remove any genes that don't pass DESeq's independent filtering, 
# indicated by having an NA p-value or adjusted p-value:

processed.deseq.tables <- list()
for (name in names(data.list.ensembl)) {
  deseq_table <- data.list.ensembl[[name]]
  deseq_table <- deseq_table[!is.na(deseq_table$pvalue), ]
  deseq_table <- deseq_table[!is.na(deseq_table$padj), ]
  processed.deseq.tables[[name]] <- deseq_table
}
data.list <- processed.deseq.tables


## -----------------------------------------------------------------------------------------------------------------------------------------------
# Generate an ordered list from the processed DESeq2 tables:

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
# Run GSEA using `msigdbr` and `fgsea`.

msigdbr_df <- msigdbr(species = "Mus musculus", category = "H")
msigdbr_list = split(x = msigdbr_df$ensembl_gene, f = msigdbr_df$gs_name)

set.seed(1) # required for deterministic performance
fgsea.output.list <- list()
for (name in names(ranked_gene_set_list)) {
  gsea_results <- fgsea(pathways = msigdbr_list, stats = ranked_gene_set_list[[name]],
                        maxSize = 500, eps = 0.0) # allow arbitrarily low p-values
  fgsea.output.list[[name]] <- gsea_results
}
saveRDS(fgsea.output.list, paste(OUTPUT_DIR, "fgsea_ensembl_ontology_list.rds", sep = "/"))


## -----------------------------------------------------------------------------------------------------------------------------------------------
# Plot the enrichment results on a heatmap:

extracted_cols <- lapply(names(fgsea.output.list), function(name) {
  col <- fgsea.output.list[[name]]$NES
  names(col) <- fgsea.output.list[[name]]$pathway
  col
})
combined.df <- do.call(cbind, extracted_cols)
colnames(combined.df) <- names(fgsea.output.list)

ht = Heatmap(combined.df,
        name = "Enrichment Scores across Perturbations and Conditions",
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        cluster_column_slices = FALSE,
        cluster_columns = TRUE,
        width = ncol(combined.df)*unit(4, "mm"),
        height = nrow(combined.df)*unit(4, "mm"),
)
pdf(paste(OUTPUT_DIR, "enrichment_combined_matrix.pdf", sep = "/"),
    width=35,
    height=19,
    useDingbats=FALSE)
draw(ht, heatmap_legend_side = "top", annotation_legend_side = "top")
dev.off()


## -----------------------------------------------------------------------------------------------------------------------------------------------
# Prepare the data for a bubble plot grouped within radiation conditions:

combined.df <- bind_rows(fgsea.output.list, .id = "Perturbations")
combined.df$pathway <- gsub("HALLMARK_", "", combined.df$pathway)
combined.df$pathway <- gsub(" ", "", combined.df$pathway)
combined.df$log10padj <- -log10(combined.df$padj)
combined.df$Perturbations <- gsub("_non-targeting_noRT", "", combined.df$Perturbations)
combined.df$Perturbations <- gsub("_non-targeting_RT", "", combined.df$Perturbations)
unfiltered.df <- combined.df

combined.df <- combined.df %>%
  group_by(pathway) %>%
  filter(any(padj < 0.05)) %>%
  ungroup()

nes.df <- combined.df %>%
  dplyr::select(-pval, -padj, -log2err, -ES, -size, -leadingEdge, -log10padj)

nes.df <- nes.df %>%
  group_by(pathway, Perturbations) %>%
  pivot_wider(names_from = Perturbations, values_from = NES) %>%
  group_by(pathway)

nes.df <- as.data.frame(nes.df)
rownames(nes.df) <- nes.df$pathway
nes.df <- nes.df %>% dplyr::select(-pathway)

row_dist <- dist(nes.df, method = "euclidean")
row_hclust <- hclust(row_dist, method = "complete")
row_order <- order.dendrogram(as.dendrogram(row_hclust))

ordering <- c()
col.hclusts <- list()
for (rt.condition in c("RT", "noRT")) {
  cols <- grep(sprintf("_%s$", rt.condition), colnames(nes.df), value = TRUE)
  col.df <- nes.df[cols]
  col.dist <- dist(t(col.df), method = "euclidean")
  col.hclust <- hclust(col.dist, method = "complete")
  col.hclusts[[rt.condition]] <- col.hclust
  col.order <- order.dendrogram(as.dendrogram(col.hclust))
  ordering <- c(ordering, cols[col.order])
}
nes.df <- nes.df[, ordering]

df.ordered <- combined.df %>%
  arrange(match(pathway, rownames(nes.df)[row_order]), match(Perturbations, colnames(nes.df)))

df.ordered$pathway <- factor(df.ordered$pathway, levels = rownames(nes.df)[row_order])
df.ordered$Perturbations <- factor(df.ordered$Perturbations, levels = colnames(nes.df))

# Add in RT fill bar
df.ordered$RT <- ifelse(grepl("_RT", df.ordered$Perturbations), "RT", "noRT")
unique_perturbations <- unique(df.ordered$Perturbations)

rt.df <- data.frame(
  xmin = unique_perturbations,
  xmax = unique_perturbations,
  ymin = rep(max(as.numeric(df.ordered$pathway)) + 1, length(unique_perturbations)),
  ymax = rep(max(as.numeric(df.ordered$pathway)) + 2, length(unique_perturbations)),
  RT = ifelse(df.ordered$RT[match(unique_perturbations, df.ordered$Perturbations)] == "RT", "RT", "noRT")
)
rt.df$numeric_xmin <- as.numeric(as.factor(rt.df$xmin)) - 0.5
rt.df$numeric_xmax <- as.numeric(as.factor(rt.df$xmax)) + 0.5

unique_perturbations <- unique(gsub("_noRT", "", unique_perturbations))
unique_perturbations <- unique(gsub("_RT", "", unique_perturbations))

# Add in GO fill bar
go_categories_RT = read.table("data/gbm_perturb/analysis/features_48h_GL261_annot2.csv", sep = ",")
colnames(go_categories_RT) = c("Gene", "GO")
rownames(go_categories_RT) = go_categories_RT$Gene
rownames(go_categories_RT) <- paste(CONTEXT, rownames(go_categories_RT), "RT", sep = "_")
rownames(go_categories_RT) <- gsub("NTC", "non-targeting", rownames(go_categories_RT))

go_categories_noRT = read.table("data/gbm_perturb/analysis/features_48h_GL261_annot2.csv", sep = ",")
colnames(go_categories_noRT) = c("Gene", "GO")
rownames(go_categories_noRT) = go_categories_noRT$Gene
rownames(go_categories_noRT) <- paste(CONTEXT, rownames(go_categories_noRT), "noRT", sep = "_")
rownames(go_categories_noRT) <- gsub("NTC", "non-targeting", rownames(go_categories_noRT))

go_categories <- rbind(go_categories_RT, go_categories_noRT)
df.ordered$GO <- go_categories[as.character(df.ordered$Perturbations), "GO"]
unique_perturbations <- unique(df.ordered$Perturbations)

rect.df <- data.frame(
  xmin = unique_perturbations,
  xmax = unique_perturbations,
  ymin = rep(max(as.numeric(df.ordered$pathway)) + 2, length(unique_perturbations)),
  ymax = rep(max(as.numeric(df.ordered$pathway)) + 3, length(unique_perturbations)),
  GO = go_categories[as.character(unique_perturbations), "GO"]
)
rect.df$numeric_xmin <- as.numeric(as.factor(rect.df$xmin)) - 0.5
rect.df$numeric_xmax <- as.numeric(as.factor(rect.df$xmax)) + 0.5

# Hard coded color mapping for consistency
term.color.map <- list(
  `Mitosis` = "#A6CEE3",
  `Translation` = "#1F78B4",
  `DNA damage response` = "#B2DF8A",
  `Metabolism` = "#33A02C",
  `LIM Domain` = "#FB9A99",
  `Mitochondria` = "#E31A1C",
  `Proteasome` = "#FDBF6F",
  `Signaling` = "#FF7F00",
  `Cyclin-dependent kinase inhibitor` = "#CAB2D6",
  `DNA Replication` = "#6A3D9A",
  `Telomere` = "#FFFF99",
  `Transcription` = "#B15928",
  `Non-targeting` = "gray",
  `RT` = "darkred",
  `noRT` = 'darkblue'
)

# Add in screening phenotype annotations
phenoTbl <- read.table("data/gbm_perturb/analysis/GL261_1+2_results_table.txt",sep='\t',header=TRUE,row.names=1)
coverageTbl <- read.table('data/gbm_perturb/analysis/gbm_pdx_perturb_GL261_integrate_xfp/pdx_perturb_GL261_concordant_sgRNAs.txt',header=TRUE,sep='\t')
coverageTblCmp <- rbind(invitronoRT = apply(coverageTbl[c("GL261_48hit_noRT_1","GL261_48hit_noRT_2"),],2,sum),
                        invitroRT = apply(coverageTbl[c("GL261_48hit_RT_1","GL261_48hit_RT_2"),],2,sum),
                        preinfnoRT = apply(coverageTbl[c("GL261_noRT_preinf_MACSFACS_1","GL261_noRT_preinf_MACSFACS_2","GL261_noRT_preinf_MACSFACS_3"),],2,sum),
                        preinfRT = apply(coverageTbl[c("GL261_RT_preinf_MACSFACS_1","GL261_RT_preinf_MACSFACS_2","GL261_RT_preinf_MACSFACS_3"),],2,sum),
                        CEDnoRT = apply(coverageTbl[c("GL261_CED_pool","GL261_noRT_CED_MACSFACS"),],2,sum),
                        CEDRT = coverageTbl[c("GL261_RT_CED_MACSFACS"),])

colnames(coverageTblCmp) <- gsub("non.targeting","NTC",colnames(coverageTblCmp))
sgAnnot <- data.frame(row.names = colnames(coverageTblCmp), gamma = phenoTbl[colnames(coverageTblCmp),"gamma"], tau = phenoTbl[colnames(coverageTblCmp),"tau"], rho = phenoTbl[colnames(coverageTblCmp),"rho"])
sgAnnot["non-targeting",] <- c(0,0,0)
sgAnnot <- sgAnnot[!(rownames(sgAnnot) %in% c("NTC", "NTC_B")),]

RT.sgAnnot <- sgAnnot
rownames(RT.sgAnnot) <- paste(rownames(RT.sgAnnot), "_RT", sep = "")

noRT.sgAnnot <- sgAnnot
rownames(noRT.sgAnnot) <- paste(rownames(noRT.sgAnnot), "_noRT", sep = "")

annotation.table <- rbind(RT.sgAnnot, noRT.sgAnnot)
annotation.table <- annotation.table[rownames(annotation.table) != "non-targeting_noRT", ]

rownames(annotation.table) <- paste(CONTEXT, rownames(annotation.table), sep = "_")
annotation.table <- subset(annotation.table, 
                           rownames(annotation.table) %in% unique(df.ordered$Perturbations))
annotation.table$Perturbations <- rownames(annotation.table)
annotation.table$Perturbations <- factor(annotation.table$Perturbations, 
                                         levels = levels(df.ordered$Perturbations))
annotation.table <- annotation.table %>%
  arrange(Perturbations)

# Add significance layer
df.ordered$Significance <- ifelse(df.ordered$padj < 0.05, "adjp < 0.05", "NA")

# Clean up names
df.ordered$Pathways <- df.ordered$pathway
levels(df.ordered$Pathways) <- levels(df.ordered$Pathways)


## -----------------------------------------------------------------------------------------------------------------------------------------------
# Plot the bubble plot:

bubble.plot <- ggplot(df.ordered,
                      aes(x = Perturbations)) +
  
  # Add the bubble and outline significance layers
  geom_point(aes(y = Pathways, size = log10padj, color = NES), alpha = 1.0) +
  geom_point(data = subset(df.ordered, padj < 0.05),
             aes(y = Pathways, size = log10padj, shape = 'adjp < 0.05'),
             color = "black", fill = NA, alpha = 0.8) +
  geom_point(data = subset(df.ordered, is.na(pval) | is.na(padj) | is.na(log10padj)), 
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
  
  # Add the GO layer
  geom_rect(data = rect.df, aes(xmin = numeric_xmin, xmax = numeric_xmax, ymin = ymin, ymax = ymax, fill = GO), inherit.aes = FALSE) +
  scale_fill_manual(values = term.color.map) +
  
  # Add the RT layer
  geom_rect(data = rt.df, aes(xmin = numeric_xmin, xmax = numeric_xmax, ymin = ymin, ymax = ymax, fill = RT), inherit.aes = FALSE) +
  
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

annotation.plot.theme <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_blank(),
  axis.title.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))
)
gamma.plot <- ggplot(annotation.table, aes(x = Perturbations, y = gamma)) +
  geom_col(width = 1.0, fill = "grey") + annotation.plot.theme
tau.plot <- ggplot(annotation.table, aes(x = Perturbations, y = tau)) +
  geom_col(width = 1.0, fill = "grey") + annotation.plot.theme
rho.plot <- ggplot(annotation.table, aes(x = Perturbations, y = rho)) +
  geom_col(width = 1.0, fill = "grey") + annotation.plot.theme

perturb.dendrogram.RT <- ggtree(as.phylo(col.hclusts[["RT"]])) + theme_tree() + coord_flip() + scale_x_reverse()
perturb.dendrogram.noRT <- ggtree(as.phylo(col.hclusts[["noRT"]])) + theme_tree() + coord_flip() + scale_x_reverse()
module.dendrogram <- ggtree(as.phylo(row_hclust)) + theme_tree()

vertical.plot <- gamma.plot / tau.plot / rho.plot / (perturb.dendrogram.RT | perturb.dendrogram.noRT) / bubble.plot + plot_layout(heights = c(0.5, 0.5, 0.5, 0.5, 3))
dendrogram.vertical.plot <- plot_spacer() / module.dendrogram + plot_layout(heights = c(0.85, 1))
combined.plot <- dendrogram.vertical.plot | vertical.plot
combined.plot <- combined.plot + plot_layout(widths = c(0.2, 0.8)) + plot_annotation(
  title = "Enrichment of Hallmark pathways across perturbations",
  subtitle = paste(CONTEXT, CONDITION, "condNormalized")
)
ggsave(paste(OUTPUT_DIR, "enrichment_bubble_plot_fgsea_nes_within_RT_all_clustered.pdf", sep = "/"),
       combined.plot, height = 13, 
       width = 6.7 + (length(unique(df.ordered$Perturbations)) / 27 * (11 - 6.7)) / 0.8, device = pdf)

## -----------------------------------------------------------------------------------------------------------------------------------------------
# Overlay mean LFCs on the existing plot

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

union_de_genes <- c()
for (name in names(data_list)) {
  deseq_table <- data_list[[name]]
  ensembl_ids <- deseq_table$ensembl
  union_de_genes <- unique(c(union_de_genes, ensembl_ids))
}

# Note that this list isn't exactly equal to the DE gene plots at the moment 
# because we also filter out all genes that don't have a valid ENSEMBL ID. 
# Next, we get the list of DESeq dataframes and filter them for these genes.
filter_to_de_genes <- function(deseq_table) {
  deseq_table <- subset(deseq_table, ensembl %in% union_de_genes)
}
de_gene_data_list <- lapply(data.list.ensembl, filter_to_de_genes)

# Generate an `msigdbr_df` that represents the genes representing each pathway.
msigdbr_df <- msigdbr(species = "Mus musculus", category = "H")
msigdbr_list = split(x = msigdbr_df$ensembl_gene, f = msigdbr_df$gs_name)

# For each perturbation x condition in the `de_gene_data_list`, 
# build a DESeq dataframe that includes only the intersection between the 
# pathway's genes and the current ones.
perturb_condition_pathway_dfs <- list()
for (perturb_condition in names(de_gene_data_list)) {
  deseq_df <- de_gene_data_list[[perturb_condition]]
  for (pathway in names(msigdbr_list)) {
    pathway_intersected_df <- subset(deseq_df, ensembl %in% msigdbr_list[[pathway]])
    perturb_condition_pathway <- paste(perturb_condition, pathway, sep = "_")
    perturb_condition_pathway_dfs[[perturb_condition_pathway]] <- pathway_intersected_df
  }
}

# Find the average LFC across dataframe in `perturb_condition_pathway_dfs`. 
# Note that for any gene that doesn't pass expression filters, we set its 
# contribution to the average to 0.
calculate_average <- function(deseq_table) {
  deseq_table$log_fc[is.na(deseq_table$pvalue) & is.na(deseq_table$padj)] <- 0
  return(mean(deseq_table$log_fc))
}
perturb_condition_pathway_means <- lapply(perturb_condition_pathway_dfs, 
                                          calculate_average)

# Map the averages to the `df.ordered` table generated above.
df.ordered$perturb_cond_pathway <- paste(df.ordered$Perturbations, 
                                         df.ordered$Pathways, sep = "_")
df.ordered$perturb_cond_pathway <- gsub(" ", "_", 
                                        df.ordered$perturb_cond_pathway)

names(perturb_condition_pathway_means) <- gsub(
  "non-targeting_RT_HALLMARK_|non-targeting_noRT_HALLMARK_|invitro_|CED_|preinf_", 
  "", names(perturb_condition_pathway_means))

df.ordered$perturb_cond_pathway_means <- 
  unlist(perturb_condition_pathway_means[df.ordered$perturb_cond_pathway])

## -----------------------------------------------------------------------------------------------------------------------------------------------
# Make the plot, but instead of using NES, use the mean:
bubble.plot <- ggplot(df.ordered,
                      aes(x = Perturbations)) +
  
  # Add the bubble and outline significance layers
  geom_point(aes(y = Pathways, size = log10padj, color = perturb_cond_pathway_means), alpha = 1.0) +
  geom_point(data = subset(df.ordered, padj < 0.05),
             aes(y = Pathways, size = log10padj, shape = 'adjp < 0.05'),
             color = "black", fill = NA, alpha = 0.8) +
  geom_point(data = subset(df.ordered, is.na(pval) | is.na(padj) | is.na(log10padj)), 
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
  
  # Add the GO layer
  geom_rect(data = rect.df, aes(xmin = numeric_xmin, xmax = numeric_xmax, ymin = ymin, ymax = ymax, fill = GO), inherit.aes = FALSE) +
  scale_fill_manual(values = term.color.map) +
  
  # Add the RT layer
  geom_rect(data = rt.df, aes(xmin = numeric_xmin, xmax = numeric_xmax, ymin = ymin, ymax = ymax, fill = RT), inherit.aes = FALSE) +
  
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

annotation.plot.theme <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text.x = element_blank(),
  axis.title.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))
)
gamma.plot <- ggplot(annotation.table, aes(x = Perturbations, y = gamma)) +
  geom_col(width = 1.0, fill = "grey") + annotation.plot.theme
tau.plot <- ggplot(annotation.table, aes(x = Perturbations, y = tau)) +
  geom_col(width = 1.0, fill = "grey") + annotation.plot.theme
rho.plot <- ggplot(annotation.table, aes(x = Perturbations, y = rho)) +
  geom_col(width = 1.0, fill = "grey") + annotation.plot.theme

perturb.dendrogram.RT <- ggtree(as.phylo(col.hclusts[["RT"]])) + theme_tree() + coord_flip() + scale_x_reverse()
perturb.dendrogram.noRT <- ggtree(as.phylo(col.hclusts[["noRT"]])) + theme_tree() + coord_flip() + scale_x_reverse()
module.dendrogram <- ggtree(as.phylo(row_hclust)) + theme_tree()

vertical.plot <- gamma.plot / tau.plot / rho.plot / (perturb.dendrogram.RT | perturb.dendrogram.noRT) / bubble.plot + plot_layout(heights = c(0.5, 0.5, 0.5, 0.5, 3))
dendrogram.vertical.plot <- plot_spacer() / module.dendrogram + plot_layout(heights = c(0.85, 1))
combined.plot <- dendrogram.vertical.plot | vertical.plot
combined.plot <- combined.plot + plot_layout(widths = c(0.2, 0.8)) + plot_annotation(
  title = "Enrichment of Hallmark pathways across perturbations",
  subtitle = paste(CONTEXT, CONDITION, "condNormalized")
)
ggsave(paste(OUTPUT_DIR, "enrichment_bubble_plot_fgsea_nes_within_RT_logfc_overlay.pdf", sep = "/"),
       combined.plot, height = 13, 
       width = 6.7 + (length(unique(df.ordered$Perturbations)) / 27 * (11 - 6.7)) / 0.8, device = pdf)

## -----------------------------------------------------------------------------------------------------------------------------------------------
# GSEA analysis output on perturbations and pathways:

analysis.df <- merge(unfiltered.df, annotation.table, on = "Perturbations")
analysis.df$leadingEdge <- sapply(analysis.df$leadingEdge, function(x) paste(x, collapse = ","))
write.table(analysis.df, file = paste(OUTPUT_DIR, "analysis_df.txt", sep = "/"), row.names = FALSE, quote = FALSE)

# For each perturbation, how many significant hallmark gene sets are there?
perturb.grouped <- analysis.df %>%
  group_by(Perturbations) %>%
  summarize(n_hallmark_significant = sum(padj < 0.05, na.rm = TRUE)) %>%
  complete(Perturbations, fill = list(n_hallmark_significant = 0))
write.table(perturb.grouped, file = paste(OUTPUT_DIR, "sig_counts_perturb_grouped.txt", sep = "/"), row.names = FALSE, quote = FALSE)

# For each Hallmark gene set, how many significant perturbations are there?
hallmark.grouped <- analysis.df %>%
  group_by(pathway) %>%
  summarize(n_perturb_significant = sum(padj < 0.05, na.rm = TRUE)) %>%
  complete(pathway, fill = list(n_perturb_significant = 0))
write.table(hallmark.grouped, file = paste(OUTPUT_DIR, "sig_counts_hallmark_grouped.txt", sep = "/"), row.names = FALSE, quote = FALSE)

