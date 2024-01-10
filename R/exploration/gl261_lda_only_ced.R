## -----------------------------------------------------------------------------------------------------------------------------------
.libPaths(c("/raleighlab/data1/czou/gbm_perturb/gbm_perturb/renv/library/R-4.3/x86_64-pc-linux-gnu", .libPaths()))


## -----------------------------------------------------------------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------------------------------------------------------------
CONTEXT <- "CED"
OUTPUT.DIR <- paste("/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/lda", CONTEXT, sep = "/")
PATH.TO.SEURAT.OBJECT <- "/raleighlab/data1/liuj/gbm_perturb/analysis/GL261_integrated_20230619.Rds"
IGNORE.GUIDES <- c("NA_RT", "NA_noRT", "non-targeting_B_RT", "non-targeting_B_noRT")
MINIMUM.COVERAGE <- 5


## -----------------------------------------------------------------------------------------------------------------------------------
data.combined <- readRDS(PATH.TO.SEURAT.OBJECT)
data <- subset(data.combined, source == CONTEXT)

guides = unique(data$sgRNACond)
guides <- guides[!(guides %in% IGNORE.GUIDES)]
data = subset(data, sgRNACond %in% guides)

sgRNACond_counts = table(data$sgRNACond)
data = subset(data, sgRNACond %in% 
                        names(sgRNACond_counts[sgRNACond_counts > MINIMUM.COVERAGE]))

data <- NormalizeData(object = data, assay = "RNA", normalization.method = "CLR", margin = 2)
DefaultAssay(object = data) <- "RNA"

data <- NormalizeData(object = data) %>% FindVariableFeatures() %>% ScaleData()
data <- RunPCA(object = data)

plot <- ElbowPlot(data, ndims = 50)
ggsave(paste(OUTPUT.DIR, sprintf("elbow_plot_%s.png", CONTEXT), sep = "/"), plot, width = 10,
         height = 10)
data <- RunUMAP(object = data, dims = 1:40)


## -----------------------------------------------------------------------------------------------------------------------------------
s.genes <- str_to_title(cc.genes$s.genes)
g2m.genes <- str_to_title(cc.genes$g2m.genes)
data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Plot the entire set of data, split by RT
plot <- SCpubr::do_DimPlot(data, 
                           group.by = 'Phase', 
                           split.by = "cond", 
                           reduction = "umap", 
                           ncol = 1, pt.size = 3, label = TRUE, label.box = TRUE, 
                           label.fill = NULL, plot_cell_borders = FALSE, label.color = "white")
plot1 <- SCpubr::do_DimPlot(data, 
                           group.by = 'sgRNA', 
                           split.by = "cond", 
                           reduction = "umap", 
                           ncol = 1, pt.size = 3, label = TRUE, label.box = TRUE, 
                           label.fill = NULL, plot_cell_borders = FALSE, label.color = "white")
ggsave(paste(OUTPUT.DIR, "umap_raw.png", sep = "/"), plot | plot1, width = 40, height = 20)

# Plot with splitting from RT and Phase
plot <- SCpubr::do_DimPlot(
  data,
  group.by = "cond",
  split.by = "sgRNA",
  reduction = 'umap',
  ncol = 5, pt.size = 3, label = TRUE, label.box = TRUE,
  label.fill = NULL, plot_cell_borders = FALSE, label.color = "white",
)
ggsave(paste(OUTPUT.DIR, "umap_raw_RT_highlight.png", sep = "/"), plot, width = 40, height = 20)
plot <- SCpubr::do_DimPlot(
  data,
  group.by = "Phase",
  split.by = "sgRNA",
  reduction = 'umap',
  ncol = 5, pt.size = 3, label = TRUE, label.box = TRUE,
  label.fill = NULL, plot_cell_borders = FALSE, label.color = "white",
)
ggsave(paste(OUTPUT.DIR, "umap_raw_phase_highlight.png", sep = "/"), plot, 
width = 40, height = 20)


## -----------------------------------------------------------------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------------------------------------------------------------
plot <- SCpubr::do_DimPlot(
  data,
  group.by = "cond",
  split.by = "sgRNA",
  reduction = 'prtbumap',
  ncol = 5, pt.size = 3, label = TRUE, label.box = TRUE,
  label.fill = NULL, plot_cell_borders = FALSE, label.color = "white",
)
ggsave(paste(OUTPUT.DIR, "umap_prtb_RT_highlight.png", sep = "/"), plot, width = 40, height = 20)
plot <- SCpubr::do_DimPlot(
  data,
  group.by = "Phase",
  split.by = "sgRNA",
  reduction = 'prtbumap',
  ncol = 5, pt.size = 3, label = TRUE, label.box = TRUE,
  label.fill = NULL, plot_cell_borders = FALSE, label.color = "white",
)
ggsave(paste(OUTPUT.DIR, "umap_prtb_phase_highlight.png", sep = "/"), plot, 
width = 40, height = 20)


## -----------------------------------------------------------------------------------------------------------------------------------
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



## -----------------------------------------------------------------------------------------------------------------------------------
# for (condition in c("noRT", "RT")) {
#   cond_object = cond_objects[[condition]]
#   plot1 <- SCpubr::do_DimPlot(data, 
#                          group.by = 'sgRNA', 
#                          split.by = "cond", 
#                          reduction = "ldaumap", 
#                          ncol = 1, pt.size = 3, label = TRUE, label.box = TRUE, 
#                          label.fill = NULL, plot_cell_borders = FALSE, label.color = "white")
#   ggsave(sprintf("%s/lda/mixscape_prtb_LDAumap_%s.png", OUTPUT.DIR, condition), 
#          plot1, width = 10, height = 10, limitsize = FALSE)
#   
#   iter = names(table(cond_object$sgRNA))
#   for (i in 1:length(table(cond_object$sgRNA))){
#     plot <- SCpubr::do_DimPlot(data, 
#                          group.by = 'sgRNA', 
#                          split.by = "cond", 
#                          reduction = "ldaumap", 
#                          ncol = 1, pt.size = 3, label = TRUE, label.box = TRUE, 
#                          label.fill = NULL, plot_cell_borders = FALSE, label.color = "white",
#                          cells.highlight = colnames(cond_object)[cond_object$sgRNA == iter[i]]) + ggtitle(iter[i])
#     ggsave(paste(OUTPUT.DIR, "/lda_perturb_umaps", "/mixscape_prtb_LDAumap_",
#                  iter[i], condition, ".png", sep=''), plot, width = 10, height = 8)
#   }
# 
#   # Make custom plot highlighting only select genes
#   genepass <- c('Nat2', 'Cdkn2a', 'Ifnar1', 'non-targeting')
#   cond_object$sgRNAPass <- cond_object$sgRNA
#   cond_object$sgRNAPass[!(cond_object$sgRNAPass %in% genepass)] = 'Other'
# 
#   # cell cycle vs select sgRNAs
#   # Set colors for each perturbation.
#   col = brewer.pal(n = 4, name = "Dark2")
#   col = c(col[1:2],'gray80','gray90',col[3:4])
# 
#   Idents(object = cond_object) <- "sgRNAPass"
#   plot <- DimPlot(object = cond_object, group.by = 'Phase', pt.size = 3, 
#                 reduction = "ldaumap", ncol = 1) 
#   plot1 <- DimPlot(object = cond_object, group.by = 'sgRNAPass', pt.size = 3, 
#                    reduction = "ldaumap", ncol = 1) + scale_color_manual(values=col, 
#                                                                          drop=FALSE) 
#   ggsave(sprintf("%s/lda/mixscape_prtb_LDAumap_passgenes_%s.png", OUTPUT.DIR, condition), 
#          plot | plot1, width = 30, height = 20, limitsize = FALSE)
# 
#   # cluster sgRNAs by perturbation 
#   cond_object.PRTB <- AverageExpression(cond_object, assay = "PRTB", 
#                                  slot = 'scale.data', group.by = "sgRNA")
#   corsPRTB <- cor(as.matrix(cond_object.PRTB$PRTB))
# 
#   png(sprintf("%s/lda/heatmap_PRTB_scaledata_averageExpr_%s.png",
#               OUTPUT.DIR, condition), width = 6000, height= 6000,res = 600)
#   heatmap(corsPRTB)
# }


## -----------------------------------------------------------------------------------------------------------------------------------
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
  ggsave(sprintf("%s/lda/mixscape_prtb_LDAumap_%s_perturb_array.pdf", OUTPUT.DIR, condition), 
         plot.array, width = 60, height = 70, limitsize = FALSE, device = "pdf")
  
  combined.plot <- SCpubr::do_DimPlot(cond.object, group.by = 'sgRNA', 
                                      pt.size = 2, reduction = "ldaumap",
               label = TRUE, label.box = TRUE, label.size = 3, 
               label.fill = NULL, label.color = "white", plot_cell_borders = FALSE) 
  ggsave(sprintf("%s/lda/mixscape_prtb_LDAumap_%s.pdf", OUTPUT.DIR, condition), 
         combined.plot, width = 10, height = 10, limitsize = FALSE, device = "pdf")
}