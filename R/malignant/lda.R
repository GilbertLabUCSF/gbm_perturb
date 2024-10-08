# GL261 Linear Discriminant Analysis (LDA)
# 
# This vignette contains code for visualizing perturbation space for GL261 
# cells stored in a Seurat object called `gl261_data.rds`. Specifically, we 
# show example anaylsis of cells whose `source == 'invitro'`, representing 
# the *in vitro* cells in our datasets.
# 
# In an attempt to provide the greatest possible degree of transparency and to
# fit with the preprocessing scripts provided, this code retains a significant 
# number of data wrangling steps.
# 
# Inputs:
# - A Seurat object that contains cells with source, cond and sgRNACond metadata
# 
# Outputs: 
# - PCA elbow plot
# - UMAP of Seurat object, separated by condition and based on perturbation 
#   signatures
# - Seurat object with LDA stored
# - Seurat objects with UMAP run on LDA results
# - UMAPs of RT and noRT perturbations using LDA as the dimensionality reduction.
# - UMAPs of RT and noRT perturbations individually highlighted using LDA as the
#   dimensionality reduction.

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
library(patchwork)
library(scales)
library(reshape2)
library(SCpubr)
library(enrichR)
library(RColorBrewer)
library(EnhancedVolcano)
library(stringr)


## -----------------------------------------------------------------------------------------------------------------------------------------------
# Define some constants:

CONTEXT <- "invitro"
OUTPUT.DIR <- paste("output/malignant/deseq/invitro/lda/%s", CONTEXT, sep = "/")
PATH.TO.SEURAT.OBJECT <- "gl261_data.rds"
IGNORE.GUIDES <- c("NA_RT", "NA_noRT", "non-targeting_B_RT", "non-targeting_B_noRT")
COVERAGE.FLOOR <- 5


## -----------------------------------------------------------------------------------------------------------------------------------------------
# Load up the R object. Then:
# 1.  Remove guides that don't pass quality checks or that don't have coverage more than the coverage floor
# 2.  Normalize the data using a centered log ratio
# 3.  Find variable features and scale, then run PCA

data.combined <- readRDS(PATH.TO.SEURAT.OBJECT)
data <- subset(data.combined, source == CONTEXT)

guides = unique(data$sgRNACond)
guides <- guides[!(guides %in% IGNORE.GUIDES)]
data = subset(data, sgRNACond %in% guides)

sgRNACond_counts = table(data$sgRNACond)
data = subset(data, sgRNACond %in% 
                        names(sgRNACond_counts[sgRNACond_counts > COVERAGE.FLOOR]))

data <- NormalizeData(object = data, assay = "RNA", normalization.method = "CLR", margin = 2)
DefaultAssay(object = data) <- "RNA"

data <- NormalizeData(object = data) %>% FindVariableFeatures() %>% ScaleData()
data <- RunPCA(object = data)

plot <- ElbowPlot(data, ndims = 50)
ggsave(paste(OUTPUT.DIR, sprintf("elbow_plot_%s.png", CONTEXT), sep = "/"), plot, width = 10,
         height = 10)
data <- RunUMAP(object = data, dims = 1:40)


## -----------------------------------------------------------------------------------------------------------------------------------------------
# Calculate perturbation signals:

data <- CalcPerturbSig(
  object = data,
  assay = "RNA",
  slot = "data",
  gd.class = "sgRNA",
  nt.cell.class = "non-targeting",
  reduction = "pca",
  ndims = 40,
  num.neighbors = 20,
  new.assay.name = "PRTB",
  split.by = "cond"
)
DefaultAssay(data) <- "PRTB"

VariableFeatures(data) <- VariableFeatures(data[["RNA"]])
data <- ScaleData(data, do.scale = FALSE, do.center = TRUE)

data <- RunPCA(object = data, reduction.key = "prtbpca", reduction.name = "prtbpca")

data <- RunUMAP(
  object = data, 
  dims = 1:40, 
  reduction = 'prtbpca', 
  reduction.key = 'prtbumap', 
  reduction.name = 'prtbumap')


## -----------------------------------------------------------------------------------------------------------------------------------------------
# Plot the UMAP based on perturbation signals:

plot <- SCpubr::do_DimPlot(
  data,
  group.by = "cond",
  split.by = "sgRNA",
  reduction = 'prtbumap',
  ncol = 5, pt.size = 3, label = TRUE, label.box = TRUE,
  label.fill = NULL, plot_cell_borders = FALSE, label.color = "white",
)
ggsave(paste(OUTPUT.DIR, "umap_prtb_RT_highlight.png", sep = "/"), plot, width = 40, height = 20)


## -----------------------------------------------------------------------------------------------------------------------------------------------
# Run LDA and plot the result on UMAP space. Note that the number of dimensions 
# we choose is equal to the number of linear discriminants found.

cond_objects = list()
for (condition in c("noRT", "RT")) {
  cond_object = subset(data, cond == condition)
  cond_object <- MixscapeLDA(
    object = cond_object,
    assay = "RNA",
    pc.assay = "PRTB",
    labels = "sgRNA",
    nt.label = "non-targeting",
    npcs = 10,
    logfc.threshold = 0.1,
    verbose = TRUE
  )
  cond_objects[[condition]] <- cond_object
  saveRDS(cond_object, file = sprintf("%s/%s_lda_data.rds", OUTPUT.DIR, condition))
}

for (condition in c("noRT", "RT")) {
  cond_object = cond_objects[[condition]]
  n.linear.discriminants <- min(length(cond_object@reductions$lda), length(unique(cond_object$sgRNA)) - 1)
  cond_object <- RunUMAP(
    object = cond_object,
    dims = 1:n.linear.discriminants,
    reduction = 'lda',
    min.dist = 0.05,
    reduction.key = 'ldaumap',
    reduction.name = 'ldaumap'
  )
  cond_objects[[condition]] <- cond_object
}

saveRDS(cond_objects, file = sprintf("%s/%s_lda_data_cond_objects_post_umap.rds", OUTPUT.DIR, condition))


## -----------------------------------------------------------------------------------------------------------------------------------------------
# Highlight each of the perturbations on the UMAP:

guides.of.interest.list <- list()
guides.of.interest.list[["RT"]] <- unique(cond_objects[["RT"]]$sgRNA)
guides.of.interest.list[["noRT"]] <- unique(cond_objects[["noRT"]]$sgRNA)

for (condition in c("noRT", "RT")) {
  cond.object <- cond_objects[[condition]]
  Idents(cond.object) <- "sgRNA"
  plot.array <- SCpubr::do_DimPlot(
    cond.object,
    reduction = "ldaumap",
    split.by = "sgRNA",
    idents.keep = guides.of.interest.list[[condition]],
    ncol = 6,
    na.value = "#F1F1F1",
    pt.size = 4,
    plot_cell_borders = FALSE
  )
  ggsave(sprintf("%s/mixscape_prtb_LDAumap_%s_perturb_array.pdf", OUTPUT.DIR, condition), 
         plot.array, width = 50, height = 40, limitsize = FALSE, device = "pdf")
  
  combined.plot <- SCpubr::do_DimPlot(cond.object, group.by = 'sgRNA', 
                                      pt.size = 2, reduction = "ldaumap",
               label = TRUE, label.box = TRUE, label.size = 3, 
               label.fill = NULL, label.color = "white", plot_cell_borders = FALSE) 
  ggsave(sprintf("%s/mixscape_prtb_LDAumap_%s.pdf", OUTPUT.DIR, condition), 
         combined.plot, width = 10, height = 10, limitsize = FALSE, device = "pdf")
}

