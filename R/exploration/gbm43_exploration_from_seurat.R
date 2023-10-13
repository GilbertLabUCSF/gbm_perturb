library(Seurat)
library(SCpubr)

PATH_TO_SEURAT_OBJECT = "/raleighlab/data1/liuj/gbm_perturb/analysis/GBM43_1_malignant_only_annotated_20230817.Rds"
data = readRDS(PATH_TO_SEURAT_OBJECT)

# Screen for only human genes

DefaultAssay(data) = "RNA"
all_genes = rownames(data)
human_genes = all_genes[grep("GRCh38-", all_genes)]
data = subset(data, features = human_genes)

# Screen for sorted cells only (or other malignancy indicator)

if (SORTED_IDENTITES_ONLY) {
  data = subset(data, sorted %in% SORTED_IDENTITIES)
}

# Screen for cells part of a group with coverage of >= MINIMUM_COVERAGE cells

sgRNACond_counts = table(data$sgRNACond)
data = subset(data, sgRNACond %in% 
                names(sgRNACond_counts[sgRNACond_counts >= MINIMUM_COVERAGE]))

# Screen for just cells with a non-targeting guide

data = subset(data, sgRNA == "non-targeting")

# Rerun normalization and UMAP

data <- SCTransform(data, verbose = TRUE)
data <- RunPCA(data, verbose = TRUE)
mindist <- 0.7
data <- RunUMAP(data, reduction = "pca",dims = 1:30, min.dist = mindist)
DimPlot(data, reduction = "umap")
ggsave("umap_dimplot.png")
DimPlot(data, reduction = "umap", group.by = "cond")
ggsave("umap_dimplot_groupbyCond.png")

# Find clusters

data <- FindNeighbors(object = data, dims = 1:30) 
res <- 0.4
data <- FindClusters(object = data, reduction = "umap", resolution = res)

# Plot the clusters in UMAP space

plot <- DimPlot(data, reduction = "umap", group.by = paste("SCT_snn_res.", res, collapse = "", sep = ""),
                label = TRUE, repe = TRUE)
ggsave(paste("pheno_sct_sub_umap_mindist_",mindist, "_res", res, ".png", collapse = "", sep= "" ), 
       plot, width = 12, height = 9)


# Find markers

Idents(data) <- data$SCT_snn_res.0.4
markers.sub <- FindAllMarkers(object = data, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
markers.sub %>% group_by(cluster) %>% top_n(2, avg_log2FC)

# Build a marker heatmap

top10.sub <- markers.sub %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top25.sub <- markers.sub %>% group_by(cluster) %>% top_n(25, avg_log2FC)

write.table(top25.sub,paste("pheno_top25markers_res_sub",res,".txt",sep=''),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)
plot <- DoHeatmap(subset(data, downsample = 200), features = top10.sub$gene, size = 3) + NoLegend()
ggsave("pheno_sub_top10_heatmap.pdf", plot, width = 12, height = 24)

# Submit these to Enrichr for screening

top25.sub$gene = gsub("GRCh38-", "", top25.sub$gene)
cluster_genes <- list()
for (i in 0:6) {
  print(cluster)
  cluster.df <- subset(top25.sub, top25.sub$cluster == i)
  print(dim(cluster.df))
  genes = cluster.df$gene
  cluster_genes[[as.character(i)]] = genes
}

for(gene in cluster_genes[["6"]]) {
  cat(gene, "\n")
}
