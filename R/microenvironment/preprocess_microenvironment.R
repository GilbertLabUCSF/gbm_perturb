## Perturb-seq of microenvironment cells - preprocessing analysis 

library(Seurat)
library(ggplot2)
library(CellChat)
library(circlize)
library(patchwork)
library(SCpubr)

options(stringsAsFactors = FALSE)

# load SB28 in vivo data
smpid <- dir(file.path("data/cellranger_out_XFP/gbm_pdx_micro_sb28_integrate/"))

dataGene.list <- list()
sgR.list <- list()

for (i in 1:length(smpid)){
  df <- Read10X(data.dir = paste("data/cellranger_out_XFP/gbm_pdx_micro_sb28_integrate/",smpid[i],"/outs/filtered_feature_bc_matrix/",sep=''))
  dataGene.list[[i]] <- CreateSeuratObject(counts=df$`Gene Expression`, project = smpid[i])
  sgR.list[[i]] <- CreateSeuratObject(counts = df$`CRISPR Guide Capture`, project = smpid[i])
  dataGene.list[[i]]$orig.ident <- smpid[i]
}

names(dataGene.list) <- smpid
names(sgR.list) <- smpid

dat.all <- merge(x = dataGene.list[[1]], y = dataGene.list[2:length(dataGene.list)])

# quality metrics for samples divided by ID 
mito.features <- grep(pattern = "^MT-", x = rownames(x = dat.all), value = TRUE, ignore.case = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = dat.all, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = dat.all, slot = 'counts'))
dat.all[['percent.mt']] <- percent.mito

# subset data and normalize
dat.sub <- subset(dat.all, subset = nFeature_RNA > 200)

# subset to high sgRNA coverage sample
Idents(dat.sub) <- dat.sub$orig.ident
dat.sub4 <- subset(dat.sub, idents = c("SB28_micro_4_unsorted"))

# subset data and normalize
dat.sub4 <- subset(dat.sub4, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 0.25)
dat.sub4 <- SCTransform(dat.sub4, verbose = FALSE)
dat.sub4 <- RunPCA(dat.sub4, verbose = FALSE)

# PCA dimensions
ElbowPlot(object = dat.sub4, ndims = 50) 

# Run UMAP embeddings
mindist=0.4
dat.sub4 <- RunUMAP(dat.sub4, reduction = "pca",dims = 1:30, min.dist = mindist)
SCpubr::do_DimPlot(dat.sub4, reduction = "umap", group.by = "orig.ident", label = TRUE, repel = TRUE, plot.axes = TRUE, legend.position = "none")

# find clusters, what is the heterogeneity due to?
res=0.1
dat.sub4 <- FindNeighbors(object = dat.sub4, dims = 1:30) 
dat.sub4 <- FindClusters(object = dat.sub4, reduction = "umap", resolution = res)

# plot clusters in UMAP space 
SCpubr::do_DimPlot(dat.sub4, reduction = "umap", pt.size=0.3, group.by = paste("SCT_snn_res.",res,collapse="",sep=""),label = TRUE, repel = TRUE, plot.axes = TRUE, legend.position = "none")

# Find cell markers for scMRMA
markers.sub <- FindAllMarkers(object = dat.sub4, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
markers.sub %>% group_by(cluster) %>% top_n(2, avg_log2FC)

# marker heatmap
top10.sub <- markers.sub %>% group_by(cluster) %>% top_n(10, avg_log2FC

# scMRMA annotation of cell types
result <- scMRMA(input = dat.sub4, species = "Mm")
dat.sub4[["scMRMA"]] <- result$multiR$annotationResult[colnames(dat.sub4),ncol(result$multiR$annotationResult)] # use the final assignment
SCpubr::do_DimPlot(dat.sub4, reduction = "umap", group.by = "scMRMA", label = TRUE, repel = TRUE,pt.size=0.3,plot.axes = TRUE, legend.position = "none")


##########################
# Generate sgRNA calls 
# Fill in sgRNA metadata 
sgR.all <- merge(x = sgR.list[[1]], y = sgR.list[2:length(sgR.list)])

sgR.sub <- subset(sgR.all, nFeature_RNA > 0)
sgR.sub <- GetAssayData(sgR.sub, slot="counts")
sgR.sub <- data.frame(sgR.sub)
sgR.sub <- sgR.sub[apply(sgR.sub,1,sum) > 0,]
colnames(sgR.sub) <- gsub("\\.","-",colnames(sgR.sub))

seqMetrics.list <- list()
sgRNATbl.list <- list()

for (i in 1:length(smpid)){
  seqMetrics.list[[i]] <- read.table(paste("data/cellranger_out_XFP/gbm_pdx_micro_sb28_integrate/",smpid[i],"/outs/metrics_summary.csv",sep=''),sep=',',header=TRUE)
  sgRNA_tags <- read.table(paste("data/cellranger_out_XFP/gbm_pdx_micro_sb28_integrate/",smpid[i],"/outs/crispr_analysis/protospacer_calls_per_cell.csv",sep=''),sep=',',header=TRUE)
  sgRNA_tags$cell_barcode <- gsub("-1","",sgRNA_tags$cell_barcode)
  row.names(sgRNA_tags) <- as.character(sgRNA_tags$cell_barcode)
  sgRNATbl.list[[i]] <- sgRNA_tags
}


# number of cells expressing sgRNAs 
sgRNATblConcordant <- sgRNATbl[sgRNATbl$num_features == 1,]
sgRNATblConcordant$GeneA <- sapply(strsplit(sgRNATblConcordant$feature_call,"_"), `[`, 1)

# Add total sgRNA UMIs and A and B counts 
sgRNATblConcordant$total_umis <- sgR.subTotal[row.names(sgRNATblConcordant)]

# Add cell type annotations to the sgRNA Tables
sgRNATblConcordant$cellType <- dat.sub$scMRMA_manual[row.names(sgRNATblConcordant)]


# Analysis of sgRNA UMIs 
sgR.subTotal <- apply(sgR.sub,2,sum)
dat.sub$sgRNA_UMI <- sgR.subTotal[colnames(dat.sub)]
dat.sub$sgRNA_logUMI <- log(sgR.subTotal[colnames(dat.sub)],10)
dat.sub$sgRNA_UMI[is.na(dat.sub$sgRNA_UMI)] <- 0
dat.sub$sgRNA_logUMI[is.na(dat.sub$sgRNA_logUMI)] <- 0
dat.sub$sgRNA_binary <- dat.sub$sgRNA_UMI > 0

# Get sgRNA target knockdown
# Normalize data to CPM
dat.sub4 <- NormalizeData(dat.sub4, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE)

# Normalize to sgNTCs but only keep cell types that have >2 sgNTCs 
coi <- unique(sgRNATblConcordant$cellType4)
coi <- coi[!is.na(coi)]

# build non targeting expression matrix 
cellTypes <- names(table(Idents(dat.sub4)))
targets <- names(table(sgRNATblConcordant$GeneA))
targets <- unique(targets)[grep("sgNeg",unique(targets),invert=TRUE)]
ntTargetTbl <- data.frame(row.names = targets, matrix(nrow = length(targets), ncol=length(cellTypes))) 
colnames(ntTargetTbl) <- cellTypes

# build bulk expression table for sgRNA on target cells 
kdTargetTbl <- data.frame(row.names = targets, matrix(nrow = length(targets), ncol=length(cellTypes))) 
colnames(kdTargetTbl) <- cellTypes

for (i in 1:ncol(kdTargetTbl)){
  ntcells = sgRNATblConcordant[sgRNATblConcordant$cellType4 == names(kdTargetTbl)[i] & sgRNATblConcordant$GeneA == "sgNegCtrl","full_barcode"]
  ntcells = ntcells[!is.na(ntcells)]
  df = GetAssayData(dat.all, slot = 'data')
  df = data.frame(df)
  dfNt = df[as.character(targets),gsub("-",".",ntcells)]
  ntTargetTbl[,i] <- apply(data.frame(dfNt),1,mean)
  
  for (j in 1:length(targets)){
    targetcells = sgRNATblConcordant[sgRNATblConcordant$cellType4 == names(kdTargetTbl)[i]  & sgRNATblConcordant$GeneA == targets[j],"full_barcode"]
    targetcells = targetcells[!is.na(targetcells)]
    dfTarget = df[targets[j],gsub("-",".",targetcells)]
    kdTargetTbl[j,i] <- mean(as.numeric(dfTarget),na.rm = TRUE)
  }
}

# get mean knockdown values
ntTargetTbl[ntTargetTbl==0] <- NaN
remainingRNATbl <- round((kdTargetTbl+1)/(ntTargetTbl+1),3)
remainingRNATbl <- remainingRNATbl[,intersect(coi,colnames(remainingRNATbl))]
remainingRNATbl <- remainingRNATbl[rowSums(is.na(remainingRNATbl)) != ncol(remainingRNATbl), ]
remainingRNATbl <- remainingRNATbl[,colSums(is.na(remainingRNATbl)) != nrow(remainingRNATbl)]

# generate ceiling of 1.0
remainingRNATbl[remainingRNATbl >= 1.0 ] = 1.0