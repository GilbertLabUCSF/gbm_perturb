## ---------------------------------------------------------------------------------------------------------------------------------------------------
# Set libpaths
# .libPaths(c("/raleighlab/data1/czou/gbm_perturb/gbm_perturb/renv/library/R-4.3/x86_64-pc-linux-gnu", .libPaths()))


## ---------------------------------------------------------------------------------------------------------------------------------------------------
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
library(RColorBrewer)


## ---------------------------------------------------------------------------------------------------------------------------------------------------
OUTPUT_DIR <- "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_explore_outputs/heterogeneity"
PATH_TO_SEURAT_OBJECT <- "/raleighlab/data1/liuj/gbm_perturb/analysis/GBM43_1_malignant_only_annotated_20230817.Rds"


## ---------------------------------------------------------------------------------------------------------------------------------------------------
# data <- readRDS(PATH_TO_SEURAT_OBJECT)
# off_targets <- c("GRAMD4", "FBXL16", "OVCH1")
# include_genes <- setdiff(unique(data$sgRNA), off_targets)
# data <- subset(data, sgRNA %in% include_genes)
# 
# data <- NormalizeData(object = data, assay = "RNA", normalization.method = "CLR", margin = 2)
# DefaultAssay(object = data) <- "RNA"
# 
# data <- NormalizeData(object = data) %>% FindVariableFeatures() %>% ScaleData()
# data <- RunPCA(object = data)
# 
# ElbowPlot(data, ndims = 50)
# data <- RunUMAP(object = data, dims = 1:40)
# 
# 
# ## ---------------------------------------------------------------------------------------------------------------------------------------------------
# s.genes <- paste("GRCh38-", cc.genes$s.genes, sep = "")
# g2m.genes <- paste("GRCh38-", cc.genes$g2m.genes, sep = "")
# data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# 
# # Plot the entire set of data, split by RT
# plot <- DimPlot(object = data, group.by = 'Phase', split.by = "cond", reduction = "umap", ncol = 1, pt.size = 3, label = TRUE, label.box = TRUE, label.color = "white")
# plot1 <- DimPlot(object = data, group.by = 'sgRNA', split.by = "cond", reduction = "umap", ncol = 1, pt.size = 3, label = TRUE, label.box = TRUE, label.color = "white")
# ggsave("/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_explore_outputs/heterogeneity/gbm_pdx_perturb_mixscape/umap_cellcycle.png", plot | plot1, width = 40, height = 20)
# 
# # Plot data for individual perturbations and show them across cell cycle
# for (i in c("non-targeting", "PRKDC", "CENPT")) {
#   perturb_data = subset(data, sgRNA == i)
#   perturb_cells = colnames(perturb_data)
#   plot <- DimPlot(object = data, group.by = 'Phase', split.by = "cond", reduction = "umap",
#                   ncol = 1, pt.size = 3)
#   plot1 <- DimPlot(object = data, group.by = 'sgRNA', split.by = "cond", reduction = "umap",
#                    ncol = 1, cells.highlight = perturb_cells, pt.size = 3)
#   ggsave(sprintf("/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_explore_outputs/heterogeneity/gbm_pdx_perturb_mixscape/umap_cellcycle_highlight_%s.png", i), plot | plot1, width = 40,
#          height = 20)
# }
# 
# 
# ## ---------------------------------------------------------------------------------------------------------------------------------------------------
# data <- CalcPerturbSig(
#   object = data,
#   assay = "RNA",
#   slot = "data",
#   gd.class = "sgRNA",
#   nt.cell.class = "non-targeting",
#   reduction = "pca",
#   ndims = 50,
#   num.neighbors = 30,
#   new.assay.name = "PRTB",
#   split.by = "cond"
# )
# DefaultAssay(data) <- "PRTB"
# 
# VariableFeatures(data) <- VariableFeatures(data[["RNA"]])
# data <- ScaleData(data, do.scale = FALSE, do.center = TRUE)
# 
# data <- RunPCA(object = data, reduction.key = "prtbpca", reduction.name = "prtbpca")
# 
# data <- RunUMAP(
#   object = data, 
#   dims = 1:50, 
#   reduction = 'prtbpca', 
#   reduction.key = 'prtbumap', 
#   reduction.name = 'prtbumap')
# 
# 
# ## ---------------------------------------------------------------------------------------------------------------------------------------------------
# plot <- DimPlot(object = data, group.by = 'Phase', pt.size = 2, split.by = "cond", reduction = "prtbumap", ncol = 1) 
# plot1 <- DimPlot(object = data, group.by = 'sgRNA', pt.size = 2, split.by = "cond", reduction = "prtbumap", ncol = 1, label = TRUE, label.box = TRUE, label.color = "white")
# ggsave("/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_explore_outputs/heterogeneity/gbm_pdx_perturb_mixscape/umap_prtb_cellcycle.png", plot | plot1, width = 40, height = 20)
# 
# 
# ## ---------------------------------------------------------------------------------------------------------------------------------------------------
# data <- RunMixscape(
#   object = data,
#   assay = "PRTB",
#   slot = "scale.data",
#   labels = "sgRNA",
#   nt.class.name = "non-targeting",
#   min.de.genes = 5,
#   logfc.threshold = 0.1,
#   iter.num = 10,
#   de.assay = "RNA",
#   verbose = TRUE,
#   prtb.type = "KO",
#   split.by = "cond"
# )
# 
# df <- table(data$mixscape_class.global, data$mixscape_class)
# df2 <- reshape2::melt(df)
# df2$Var2 <- as.character(df2$Var2)
# test <- df2[which(df2$Var1 == "KO"),]
# test <- test[order(test$value, decreasing = TRUE),]
# new.levels <- test$Var2
# df2$Var2 <- factor(df2$Var2, levels = new.levels )
# df2$Var1 <- factor(df2$Var1, levels = c("non-targeting", "NP", "KO"))
# df2$gene <- sapply(as.character(df2$Var2), function(x) strsplit(x, split = " ")[[1]][1])
# df3 <- df2[-c(which(df2$gene == "non-targeting")),]
# 
# p1 <- ggplot(df3, aes(x = gene, y = value, fill= Var1)) +
#   geom_bar(stat= "identity") +
#   theme_classic()+
#   scale_fill_manual(values = c("grey49", "grey79","coral1")) + 
#   ylab("n cells") +
#   xlab("sgRNA")
# 
# plot <- p1 + theme(axis.text.x = element_blank(),
#            axis.text.y = element_text(size = 10), 
#            axis.title = element_text(size = 8), 
#            strip.text = element_text(size=8, face = "bold")) + 
#   facet_wrap(vars(gene),ncol = 5, scales = "free") +
#   labs(fill = "mixscape class") +theme(legend.title = element_text(size = 10),
#                                        legend.text = element_text(size = 8))
# 
# ggsave("/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_explore_outputs/heterogeneity/gbm_pdx_perturb_mixscape/mixscape_prtb_ncell_calls.png", plot, width = 12, height = 8)
# 
# saveRDS(data, "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_explore_outputs/heterogeneity/gbm_pdx_perturb_mixscape/mixscaped_data.rds")
# 
# 
# ## ---------------------------------------------------------------------------------------------------------------------------------------------------
# perturb_plot <- PlotPerturbScore(object = data, target.gene.ident = "PRKDC", target.gene.class = "sgRNA", mixscape.class = "mixscape_class", col = "coral2", split.by = "cond") + labs(fill = "mixscape class")
# ggsave("/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_explore_outputs/heterogeneity/gbm_pdx_perturb_mixscape/prkdc_perturb_plot.png", perturb_plot, width = 12, height = 8)
# 
# # For remaining plots, split by RT and noRT
# 
# for (condition in c("noRT", "RT")) {
#   cond_object = subset(data, cond == condition)
#   vln_plot <- VlnPlot(cond_object, "mixscape_class_p_ko", idents = c("non-targeting", "PRKDC KO", "PRKDC NP")) +
#   theme(axis.text.x = element_text(angle = 0, hjust = 0.5), 
#         axis.text = element_text(size = 16), 
#         plot.title = element_text(size = 20)) + 
#   NoLegend() +
#   ggtitle("mixscape posterior probabilities")
#   ggsave(sprintf("/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_explore_outputs/heterogeneity/gbm_pdx_perturb_mixscape/prkdc_vln_plot_%s.png", condition), vln_plot, width = 12,
#          height = 8)
#   
#   Idents(object = cond_object) <- "sgRNA"
#   genepass <- c('PRKDC','CENPT','RBX1')
# 
#   for (i in 1:length(genepass)) {
#     plot <- MixscapeHeatmap(object = cond_object, 
#                             ident.1 = "non-targeting", 
#                             ident.2 = genepass[i], 
#                             balanced = FALSE, 
#                             assay = "RNA", 
#                             max.genes = 30, angle = 0, 
#                             group.by = "mixscape_class", 
#                             max.cells.group = 300,
#                             size=6.5) + NoLegend() + 
#       theme(axis.text.y = element_text(size = 16))
#       ggsave(paste("/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_explore_outputs/heterogeneity/gbm_pdx_perturb_mixscape/de_genes/mixscape_prtb_",genepass[i],condition,"_DEgenes.png", sep=''),
#              plot, width = 15, height = 10)
#   }
# }
data <- readRDS("/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_explore_outputs/heterogeneity/gbm_pdx_perturb_mixscape/mixscaped_data.rds")

## ---------------------------------------------------------------------------------------------------------------------------------------------------
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
  
  cond_object <- RunUMAP(
    object = cond_object,
    dims = 1:30,
    reduction = 'lda',
    reduction.key = 'ldaumap',
    reduction.name = 'ldaumap'
  )
  
  plot <- DimPlot(object = cond_object, group.by = 'Phase', pt.size = 1, 
                reduction = "ldaumap", ncol = 1) 
  plot1 <- DimPlot(object = cond_object, group.by = 'sgRNA', pt.size = 1, 
               reduction = "ldaumap", ncol = 1,
               label = TRUE, label.box = TRUE, label.color = "white") 
  ggsave(sprintf("/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_explore_outputs/heterogeneity/gbm_pdx_perturb_mixscape/lda/mixscape_prtb_LDAumap_%s.png", condition), 
         plot | plot1, width = 30, height = 20, limitsize = FALSE)
  
  iter = names(table(cond_object$sgRNA))
  for (i in 1:length(table(cond_object$sgRNA))){
    plot <- DimPlot(object = cond_object, group.by = 'sgRNA', 
                    pt.size = 1, reduction = "ldaumap",
                    ncol = 1,
                    cells.highlight = colnames(cond_object)[cond_object$sgRNA == iter[i]]) +
      NoLegend() + ggtitle(iter[i])
    ggsave(paste("/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_explore_outputs/heterogeneity/gbm_pdx_perturb_mixscape/lda/highlight/mixscape_prtb_LDAumap_",
                 iter[i], condition, ".png", sep=''), plot, width = 10, height = 8)
  }

  # Make custom plot highlighting only select genes
  genepass <- c('PRKDC','CENPT','RBX1','non-targeting')
  cond_object$sgRNAPass <- cond_object$sgRNA
  cond_object$sgRNAPass[!(cond_object$sgRNAPass %in% genepass)] = 'Other'

  # cell cycle vs select sgRNAs
  # Set colors for each perturbation.
  col = brewer.pal(n = 4, name = "Dark2")
  col = c(col[1:2],'gray80','gray90',col[3:4])

  Idents(object = cond_object) <- "sgRNAPass"
  plot <- DimPlot(object = cond_object, group.by = 'Phase', pt.size = 3, 
                reduction = "ldaumap", ncol = 1) 
  plot1 <- DimPlot(object = cond_object, group.by = 'sgRNAPass', pt.size = 3, 
                   reduction = "ldaumap", ncol = 1) + scale_color_manual(values=col, 
                                                                         drop=FALSE) 
  ggsave(sprintf("/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_explore_outputs/heterogeneity/gbm_pdx_perturb_mixscape/lda/mixscape_prtb_LDAumap_passgenes_%s.png", condition), 
         plot | plot1, width = 30, height = 20, limitsize = FALSE)

  # cluster sgRNAs by perturbation 
  cond_object.PRTB <- AverageExpression(cond_object, assay = "PRTB", 
                                 slot = 'scale.data', group.by = "sgRNA")
  corsPRTB <- cor(cond_object.PRTB$PRTB)

  png(sprintf("/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_explore_outputs/heterogeneity/gbm_pdx_perturb_mixscape/lda/heatmap_PRTB_scaledata_averageExpr_%s.png",
              condition), width = 6000, height=6000,res = 600)
  heatmap(corsPRTB)
  saveRDS(cond_object, file = sprintf("/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_explore_outputs/heterogeneity/gbm_pdx_perturb_mixscape/%s_data.rds", condition))
}


## ---------------------------------------------------------------------------------------------------------------------------------------------------
# saveRDS(data, file = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_explore_outputs/heterogeneity/gbm_pdx_perturb_mixscape/overall_data.rds")
# knitr::purl("/raleighlab/data1/czou/gbm_perturb/gbm_perturb/R/exploration/heterogeneity.Rmd")

