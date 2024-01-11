# GL261 Gene Set Enrichment Analysis (GSEA)
# 
# This script contains code for generating and plotting GSEA results on 
# GL261 DESeq output stored in a directory called data/malignant/deseq/invitro. 
# Note that this example runs on DESeq for cells labeled `noRT`, 
# meaning they didn't receive radiation.
# 
# Inputs:
# - DESeq output directory containing a CSV file per perturbation
# 
# Outputs: 
# - Table of genes that don't convert to ENSEMBL
# - GSEA R object
# - Heatmap of enrichment scores
# - Bubble plot of enrichment scores
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
set.seed(0)


## -----------------------------------------------------------------------------------------------------------------------------------------------
# Define constants:

INPUT_DIRS <- c("data/malignant/deseq/invitro")
CONTEXT <- "invitro"
CONDITION <- "RT"
OUTPUT_DIR <- sprintf("data/malignant/gsea/%s/%s", CONTEXT, CONDITION)


## -----------------------------------------------------------------------------------------------------------------------------------------------
# Import the DESeq output:

data.list <- list()
for (dir in INPUT_DIRS) {
  file_names <- list.files(dir)
  for (i in file_names) {
    name <- gsub(".csv", "", i)
    df <- read.table(paste(dir, "/", i, sep = ""))
    data.list[[name]] <- df
  }
}


## -----------------------------------------------------------------------------------------------------------------------------------------------
# Separate out RT and noRT use cases:

noRT.indices <- grepl("_noRT_non-targeting_noRT", names(data.list))
RT.indices <- grepl("_RT_non-targeting_noRT", names(data.list))
noRT.list <- data.list[noRT.indices]
RT.list <- data.list[RT.indices]
overall.list <- data.list

if (CONDITION == "RT") {
  data.list.symbols <- RT.list
} else {
  data.list.symbols <- noRT.list
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
  col <- fgsea.output.list[[name]]$ES
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
# Prepare the data for the bubble plot

combined.df <- bind_rows(fgsea.output.list, .id = "Perturbations")
combined.df$pathway <- gsub("HALLMARK_", "", combined.df$pathway)
combined.df$pathway <- gsub(" ", "", combined.df$pathway)
combined.df$log10pvalue <- -log10(combined.df$padj)
combined.df$Perturbations <- gsub("_non-targeting_noRT", "", combined.df$Perturbations)
unfiltered.df <- combined.df

combined.df <- combined.df %>%
  group_by(pathway) %>%
  filter(any(padj < 0.05)) %>%
  ungroup()

es.df <- combined.df %>%
  dplyr::select(-pval, -padj, -log2err, -NES, -size, -leadingEdge, -log10pvalue)

es.df <- es.df %>%
  group_by(pathway, Perturbations) %>%
  pivot_wider(names_from = Perturbations, values_from = ES) %>%
  group_by(pathway)

es.df <- as.data.frame(es.df)
rownames(es.df) <- es.df$pathway
es.df <- es.df %>% dplyr::select(-pathway)

row_dist <- dist(es.df, method = "euclidean")
row_hclust <- hclust(row_dist, method = "complete")
row_order <- order.dendrogram(as.dendrogram(row_hclust))

col_dist <- dist(t(es.df), method = "euclidean")
col_hclust <- hclust(col_dist, method = "complete")
col_order <- order.dendrogram(as.dendrogram(col_hclust))

df.ordered <- combined.df %>%
  arrange(match(pathway, rownames(es.df)[row_order]), match(Perturbations, colnames(es.df)[col_order]))

df.ordered$pathway <- factor(df.ordered$pathway, levels = rownames(es.df)[row_order])
df.ordered$Perturbations <- factor(df.ordered$Perturbations, levels = colnames(es.df)[col_order])

# Add in RT fill bar
df.ordered$RT <- ifelse(grepl("_RT", df.ordered$Perturbations), "RT", "noRT")
unique_perturbations <- unique(df.ordered$Perturbations)

rect.df <- data.frame(
  xmin = unique_perturbations,
  xmax = unique_perturbations,
  ymin = rep(max(as.numeric(df.ordered$pathway)) + 1, length(unique_perturbations)),
  ymax = rep(max(as.numeric(df.ordered$pathway)) + 2, length(unique_perturbations)),
  RT = ifelse(df.ordered$RT[match(unique_perturbations, df.ordered$Perturbations)] == "RT", "RT", "noRT")
)
rect.df$numeric_xmin <- as.numeric(as.factor(rect.df$xmin)) - 0.5
rect.df$numeric_xmax <- as.numeric(as.factor(rect.df$xmax)) + 0.5

unique.perturbations <- unique(gsub("_noRT", "", unique_perturbations))
unique.perturbations <- unique(gsub("_RT", "", unique.perturbations))

# Add in GO fill bar
go_categories = read.table("/raleighlab/data1/liuj/gbm_perturb/analysis/features_48h_GL261_annot2.csv", sep = ",")
colnames(go_categories) = c("Gene", "GO")
rownames(go_categories) = go_categories$Gene
rownames(go_categories) <- paste(CONTEXT, rownames(go_categories), CONDITION, sep = "_")
rownames(go_categories) <- gsub("NTC", "non-targeting", rownames(go_categories))

df.ordered$GO <- go_categories[as.character(df.ordered$Perturbations), "GO"]
unique_perturbations <- unique(df.ordered$Perturbations)

rect.df <- data.frame(
  xmin = unique_perturbations,
  xmax = unique_perturbations,
  ymin = rep(max(as.numeric(df.ordered$pathway)) + 1, length(unique_perturbations)),
  ymax = rep(max(as.numeric(df.ordered$pathway)) + 2, length(unique_perturbations)),
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
  `Non-targeting` = "gray"
)

# Add in screening phenotype annotations
phenoTbl <- read.table("/raleighlab/data1/liuj/gbm_perturb/analysis/GL261_1+2_results_table.txt",sep='\t',header=TRUE,row.names=1)
coverageTbl <- read.table('/raleighlab/data1/liuj/gbm_perturb/analysis/gbm_pdx_perturb_GL261_integrate_xfp/pdx_perturb_GL261_concordant_sgRNAs.txt',header=TRUE,sep='\t')
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
df.ordered$Pathways <- gsub("_", " ", df.ordered$Pathways)


## -----------------------------------------------------------------------------------------------------------------------------------------------
# Plot the bubble plot:

bubble.plot <- ggplot(df.ordered,
               aes(x = Perturbations)) +
  
  # Add the bubble and outline significance layers
  geom_point(aes(y = Pathways, size = log10pvalue, color = ES), alpha = 1.0) +
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
  
  # Add the bar layer
  geom_rect(data = rect.df, aes(xmin = numeric_xmin, xmax = numeric_xmax, ymin = ymin, ymax = ymax, fill = GO), inherit.aes = FALSE) +
  scale_fill_manual(values = term.color.map) +

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

perturb.dendrogram <- ggtree(as.phylo(col_hclust)) + theme_tree() + coord_flip() + scale_x_reverse()
module.dendrogram <- ggtree(as.phylo(row_hclust)) + theme_tree()

vertical.plot <- gamma.plot / tau.plot / rho.plot / perturb.dendrogram / bubble.plot + plot_layout(heights = c(0.5, 0.5, 0.5, 0.5, 3))
dendrogram.vertical.plot <- plot_spacer() / module.dendrogram + plot_layout(heights = c(0.85, 1))
combined.plot <- dendrogram.vertical.plot | vertical.plot
combined.plot <- combined.plot + plot_layout(widths = c(0.2, 0.8)) + plot_annotation(
  title = "Enrichment of Hallmark pathways across perturbations",
  subtitle = paste(CONTEXT, CONDITION, "noRTNormalized")
)
ggsave(paste(OUTPUT_DIR, "enrichment_bubble_plot_fgsea_3.pdf", sep = "/"),
       combined.plot, height = 13, width = 14, device = pdf)


## -----------------------------------------------------------------------------------------------------------------------------------------------
# GSEA analysis questions on perturbations and pathways:

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

