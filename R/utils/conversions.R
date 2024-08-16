# Utility to convert Seurat objects into h5ad objects that can be loaded as AnnData objects

library(Seurat)
library(SeuratDisk)
options(Seurat.object.assay.version = 'v3')

data <- readRDS("~/code/genome-biology-methods/data/simplified_seurat.rds")
DefaultAssay(data) <- "RNA"
data[["RNA"]] <- as(object = data[["RNA"]], Class = "Assay")
SaveH5Seurat(data, filename = "~/code/genome-biology-methods/data/simplified_seurat.h5Seurat")
Convert("~/code/genome-biology-methods/data/simplified_seurat.h5Seurat", "~/code/genome-biology-methods/data/simplified_seurat.h5ad")