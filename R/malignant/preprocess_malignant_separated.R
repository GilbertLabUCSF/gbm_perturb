## Analysis of in VIVO perturb-seq separated for analysis by context
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(SCpubr)
library(scMRMA)
library(ggpubr)

# load all GL261  data
smpid <- dir(file.path("/raleighlab/data1/liuj/gbm_perturb/cellranger_out_XFP/invivo_perturb_GL261_integrate"))
dataGene.list <- list()
sgR.list <- list()

for (i in 1:length(smpid)){
  df <- Read10X(data.dir = paste("/raleighlab/data1/liuj/gbm_perturb/cellranger_out_XFP/invivo_perturb_GL261_integrate/",smpid[i],"/outs/filtered_feature_bc_matrix/",sep=''))
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
dat.sub <- SCTransform(dat.sub, verbose = FALSE)
dat.sub <- RunPCA(dat.sub, verbose = FALSE)

##########################
# Analysis of sgRNA coverage 
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
  seqMetrics.list[[i]] <- read.table(paste("/raleighlab/data1/liuj/gbm_perturb/cellranger_out_XFP/invivo_perturb_GL261_integrate/",smpid[i],"/outs/metrics_summary.csv",sep=''),sep=',',header=TRUE)
  sgRNA_tags <- read.table(paste("/raleighlab/data1/liuj/gbm_perturb/cellranger_out_XFP/invivo_perturb_GL261_integrate/",smpid[i],"/outs/crispr_analysis/protospacer_calls_per_cell.csv",sep=''),sep=',',header=TRUE)
  sgRNA_tags$cell_barcode <- gsub("-1","",sgRNA_tags$cell_barcode)
  row.names(sgRNA_tags) <- as.character(sgRNA_tags$cell_barcode)
  sgRNATbl.list[[i]] <- sgRNA_tags
}

# histogram of sgRNA dectections
names(sgRNATbl.list) <- smpid
sgRNATbl <- ldply(sgRNATbl.list, data.frame)

# categorize sample types
sgRNATbl$source <- sgRNATbl$.id
sgRNATbl$sourceCond <- sgRNATbl$.id
smpTbl <- data.frame(orig.ident = dat.sub$orig.ident, sorted = dat.sub$sorted, source = dat.sub$source, cond = dat.sub$cond)
smpTbl <- unique(smpTbl)
row.names(smpTbl) <- as.character(smpTbl$orig.ident)
smpTbl$sourceCond <- paste(smpTbl$source,smpTbl$cond,sep="_")

sgRNATbl$source <- smpTbl[sgRNATbl$.id,"source"]
sgRNATbl$sourceCond <- smpTbl[sgRNATbl$.id,"sourceCond"]

# Tabulate sgRNA read data - what proportion of each sample has cells with sgRNAs?
sgR.subTotal <- apply(sgR.sub,2,sum)
dat.sub$sgRNA_UMI <- sgR.subTotal[colnames(dat.sub)]
dat.sub$sgRNA_logUMI <- log(sgR.subTotal[colnames(dat.sub)],10)
dat.sub$sgRNA_UMI[is.na(dat.sub$sgRNA_UMI)] <- 0
dat.sub$sgRNA_logUMI[is.na(dat.sub$sgRNA_logUMI)] <- 0
dat.sub$sgRNA_binary <- dat.sub$sgRNA_UMI > 0
dat.sub$source <- factor(dat.sub$source,levels=c("invitro","preinf","CED"))


SCpubr::do_FeaturePlot(dat.sub, features = "sgRNA_logUMI", plot.title = "GL261 sgRNA UMI", order = TRUE, split.by.idents ='TRUE')
SCpubr::do_FeaturePlot(dat.sub, features = "sgRNA_logUMI", plot.title = "GL261 sgRNA UMI", order = TRUE, split.by='source')
SCpubr::do_FeaturePlot(dat.sub, features = "sgRNA_logUMI", plot.title = "GL261 sgRNA UMI", order = TRUE, split.by='sorted')

# Which cell types are infected with sgRNA_UMI? #
cellTypeInfxTbl <- data.frame(table(data.frame(sgRNA = dat.sub$sgRNA_binary, cellType = dat.sub$scMRMA_manual, source =  dat.sub$source)))
cellTypeInfxTbl <- cellTypeInfxTbl[cellTypeInfxTbl$sgRNA == TRUE,]
cellTypeInfxTbl$source <- factor(cellTypeInfxTbl$source,levels=c("invitro","preinf","CED"))
cellTypeInfxTbl$cellType <- factor(cellTypeInfxTbl$cellType,levels=rev(as.character(unique(cellTypeInfxTbl$cellType))))


# number of cells with expression of concordant sgRNAs 
sgRNATblConcordant <- sgRNATbl[sgRNATbl$num_features == 2,]
sgRNATblConcordant <- sgRNATblConcordant %>% separate(feature_call, c("sgRNA_A", "sgRNA_B"), sep = "\\|")
sgRNATblConcordant$GeneA <- sapply(strsplit(sgRNATblConcordant$sgRNA_A,"_"), `[`, 1)
sgRNATblConcordant$GeneB <- sapply(strsplit(sgRNATblConcordant$sgRNA_B,"_"), `[`, 1)
sgRNATblConcordant <- sgRNATblConcordant[sgRNATblConcordant$GeneA == sgRNATblConcordant$GeneB,]

# identify negative controls
sgRNATblConcordant[sgRNATblConcordant$sgRNA_A == "non-targeting_00025 ","GeneA"] <- "non-targeting"

# Add cell names with barcode suffixes and use them as row names
smpidNumTbl <- data.frame(smpid = smpid, num = seq(1,length(smpid),1))
row.names(smpidNumTbl) <- as.character(smpidNumTbl$smpid)
sgRNATblConcordant$full_barcode <- paste(sgRNATblConcordant$cell_barcode,"-1_",smpidNumTbl[sgRNATblConcordant$.id,"num"],sep='')
row.names(sgRNATblConcordant) <- as.character(sgRNATblConcordant$full_barcode)

# Add total sgRNA UMIs and A and B counts
sgRNATblConcordant$total_umis <- sgR.subTotal[row.names(sgRNATblConcordant)]
sgRNATblConcordant$GeneA_umi <- as.numeric(sapply(strsplit(sgRNATblConcordant$num_umis,"\\|"), `[`, 1))
sgRNATblConcordant$GeneB_umi <- as.numeric(sapply(strsplit(sgRNATblConcordant$num_umis,"\\|"), `[`, 2))
sgRNATblConcordant$source <- factor(sgRNATblConcordant$source,levels=c("invitro","preinf","CED"))

#########################################################
# separate the samples by in vitro, preinfected and CED #
Idents(dat.sub) <- dat.sub$source
dat.preinf <- subset(dat.sub, ident="preinf")
dat.ced <- subset(dat.sub, ident="CED")

# Transform separately in UMAP space
dat.preinf <- SCTransform(dat.preinf, verbose = FALSE)
dat.preinf <- RunPCA(dat.preinf, verbose = FALSE)

dat.ced <- SCTransform(dat.ced, verbose = FALSE)
dat.ced <- RunPCA(dat.ced, verbose = FALSE)

# PCA dimensions
ElbowPlot(object = dat.preinf,ndims = 50) 
ElbowPlot(object = dat.ced,ndims = 50) 

# Run UMAP embeddings
colors <- c("invitro" = "#E64B35FF",
            "preinf" = "#4DBBD5FF",
            "CED" = "#E18727FF",
            "noRT" = "#3C5488FF",
            "unsorted" = "#5C88DAFF",
            "MACSFACS" = "#FFCD00FF",
            "RT" = "#DC0000FF")

mindist=0.7
dat.preinf <- RunUMAP(dat.preinf, reduction = "pca",dims = 1:30, min.dist = mindist)

SCpubr::do_DimPlot(dat.preinf, reduction = "umap", group.by = "orig.ident", label = FALSE, repel = TRUE, plot.axes = TRUE, legend.position = "right",pt.size = 0.5)
SCpubr::do_DimPlot(dat.preinf, reduction = "umap", group.by = "sorted", label = FALSE, repel = TRUE, plot.axes = TRUE, legend.position = "bottom",pt.size = 0.5,colors.use = colors)
SCpubr::do_DimPlot(dat.preinf, reduction = "umap", group.by = "source", label = FALSE, repel = TRUE, plot.axes = TRUE, legend.position = "bottom",pt.size = 0.5,colors.use = colors)
SCpubr::do_DimPlot(dat.preinf, reduction = "umap", group.by = "cond", label = FALSE, repel = TRUE, plot.axes = TRUE, legend.position = "bottom",pt.size = 0.5,colors.use = colors)

# CED umaps
dat.ced <- RunUMAP(dat.ced, reduction = "pca",dims = 1:30, min.dist = mindist)

SCpubr::do_DimPlot(dat.ced, reduction = "umap", group.by = "orig.ident", label = FALSE, repel = TRUE, plot.axes = TRUE, legend.position = "right",pt.size = 0.5)
SCpubr::do_DimPlot(dat.ced, reduction = "umap", group.by = "sorted", label = FALSE, repel = TRUE, plot.axes = TRUE, legend.position = "bottom",pt.size = 0.5,colors.use = colors)
SCpubr::do_DimPlot(dat.ced, reduction = "umap", group.by = "source", label = FALSE, repel = TRUE, plot.axes = TRUE, legend.position = "bottom",pt.size = 0.5,colors.use = colors)
gSCpubr::do_DimPlot(dat.ced, reduction = "umap", group.by = "cond", label = FALSE, repel = TRUE, plot.axes = TRUE, legend.position = "bottom",pt.size = 0.5,colors.use = colors)

# Combined UMAPS
SCpubr::do_DimPlot(dat.sub, reduction = "umap", group.by = "orig.ident", label = FALSE, repel = TRUE, plot.axes = TRUE, legend.position = "right",pt.size = 0.5)
SCpubr::do_DimPlot(dat.sub, reduction = "umap", group.by = "sorted", label = FALSE, repel = TRUE, plot.axes = TRUE, legend.position = "bottom",pt.size = 0.5,colors.use = colors)
SCpubr::do_DimPlot(dat.sub, reduction = "umap", group.by = "source", label = FALSE, repel = TRUE, plot.axes = TRUE, legend.position = "bottom",pt.size = 0.5,colors.use = colors)
SCpubr::do_DimPlot(dat.sub, reduction = "umap", group.by = "cond", label = FALSE, repel = TRUE, plot.axes = TRUE, legend.position = "bottom",pt.size = 0.5,colors.use = colors)

# Only show sgRNA UMI for concordant genes 
dat.preinf$sgRNA_logUMI_concordant <- dat.preinf$sgRNA_logUMI
dat.preinf$sgRNA_logUMI_concordant[is.na(dat.preinf$sgRNA)] <- 0 

dat.ced$sgRNA_logUMI_concordant <- dat.ced$sgRNA_logUMI
dat.ced$sgRNA_logUMI_concordant[is.na(dat.ced$sgRNA)] <- 0 

dat.sub$sgRNA_logUMI_concordant <- dat.sub$sgRNA_logUMI
dat.sub$sgRNA_logUMI_concordant[is.na(dat.sub$sgRNA)] <- 0 

# Find clusters
dat.preinf <- FindNeighbors(object = dat.preinf, dims = 1:30) 
res=0.4
dat.preinf <- FindClusters(object = dat.preinf, reduction = "umap", resolution = res)
Cpubr::do_DimPlot(dat.preinf, reduction = "umap", group.by = paste("SCT_snn_res.",res,collapse="",sep=""), label = TRUE, repel = TRUE,
                           plot.axes = TRUE, legend.position = "none")

dat.ced <- FindNeighbors(object = dat.ced, dims = 1:30) 
res=0.4
dat.ced <- FindClusters(object = dat.ced, reduction = "umap", resolution = res)
SCpubr::do_DimPlot(dat.ced, reduction = "umap", group.by = paste("SCT_snn_res.",res,collapse="",sep=""), label = TRUE, repel = TRUE,
                           plot.axes = TRUE, legend.position = "none")

# Featureplot of exogenous genes
goi <- c("GFP","BFP","KRAB","sgRNA_logUMI_concordant")
SCpubr::do_FeaturePlot(dat.preinf, features = goi, plot.title = "GL261 tumors exogenous genes",order = TRUE,ncol = 4)

goi <- c("GFP","BFP","KRAB","sgRNA_logUMI_concordant")
SCpubr::do_FeaturePlot(dat.ced, features = goi, plot.title = "GL261 tumors exogenous genes",order = TRUE,ncol = 4)

goi <- c("GFP","BFP","KRAB","sgRNA_logUMI_concordant")
SCpubr::do_FeaturePlot(dat.sub, features = goi, plot.title = "GL261 tumors exogenous genes",order = TRUE,ncol = 4)

# Dimplot of sgRNAs
SCpubr::do_DimPlot(dat.preinf[, dat.preinf$sgRNA %in% names(table(dat.preinf$sgRNA)) ], 
                   reduction = "umap", group.by = "sgRNA", label = FALSE, repel = TRUE, plot.axes = FALSE,shuffle=FALSE,na.value = NA)
SCpubr::do_DimPlot(dat.ced[, dat.ced$sgRNA %in% names(table(dat.ced$sgRNA)) ], 
                   reduction = "umap", group.by = "sgRNA", label = FALSE, repel = TRUE, plot.axes = FALSE,shuffle=FALSE,na.value = NA)
SCpubr::do_DimPlot(dat.sub[, dat.sub$sgRNA %in% names(table(dat.sub$sgRNA)) ], 
                   reduction = "umap", group.by = "sgRNA", label = FALSE, repel = TRUE, plot.axes = FALSE,shuffle=FALSE,na.value = NA)

# SCMRMA runs 
# scMRMA annotation of cell types
result.preinf <- scMRMA(input = dat.preinf, species = "Mm")
dat.preinf[["scMRMA"]] <- result.preinf$multiR$annotationResult[colnames(dat.preinf),ncol(result.preinf$multiR$annotationResult)] # use the final assignment

result.ced <- scMRMA(input = dat.ced, species = "Mm")
dat.ced[["scMRMA"]] <- result.ced$multiR$annotationResult[colnames(dat.ced),ncol(result.ced$multiR$annotationResult)] # use the final assignment

SCpubr::do_DimPlot(dat.preinf, reduction = "umap", group.by = "scMRMA", label = TRUE, repel = TRUE, plot.axes = TRUE, legend.position = "none")
SCpubr::do_DimPlot(dat.ced, reduction = "umap", group.by = "scMRMA", label = TRUE, repel = TRUE, plot.axes = TRUE, legend.position = "none")

# Manually adjust scMRMA annotations based on known marker expression and proximity to GL261 cells 
dat.preinf$scMRMA_manual <- dat.preinf$scMRMA
dat.preinf$scMRMA_manual <- gsub("Enteric glia cells", "GL261",dat.preinf$scMRMA_manual)
dat.preinf$scMRMA_manual <- gsub("Radial glial cell", "GL261",dat.preinf$scMRMA_manual)
dat.preinf$scMRMA_manual <- gsub("Epithelial cells", "GL261",dat.preinf$scMRMA_manual)
dat.preinf$scMRMA_manual <- gsub("Epithelial cell", "GL261",dat.preinf$scMRMA_manual)
dat.preinf$scMRMA_manual <- gsub("Nuocytes", "Microglia",dat.preinf$scMRMA_manual)
dat.preinf$scMRMA_manual <- gsub("Hepatic stellate cells", "Pericytes",dat.preinf$scMRMA_manual)
dat.preinf$scMRMA_manual <- gsub("Glutaminergic neurons", "Neurons",dat.preinf$scMRMA_manual)
dat.preinf$scMRMA_manual <- gsub("DC", "Dendritic cells",dat.preinf$scMRMA_manual)
dat.preinf$scMRMA_manual[dat.preinf$SCT_snn_res.0.4 == "14"] <- "Astrocyte" 

SCpubr::do_DimPlot(dat.preinf, reduction = "umap", group.by = "scMRMA_manual", label = TRUE, repel = TRUE,
                           plot.axes = FALSE, legend.position = "none")

dat.ced$scMRMA_manual <- dat.ced$scMRMA
dat.ced$scMRMA_manual <- gsub("Enteric glia cells", "GL261",dat.ced$scMRMA_manual)
dat.ced$scMRMA_manual <- gsub("Radial glial cell", "GL261",dat.ced$scMRMA_manual)
dat.ced$scMRMA_manual <- gsub("Epithelial cells", "GL261",dat.ced$scMRMA_manual)
dat.ced$scMRMA_manual <- gsub("Epithelial cell", "GL261",dat.ced$scMRMA_manual)
dat.ced$scMRMA_manual <- gsub("Bergmann glial cell", "Glial cell",dat.ced$scMRMA_manual)
dat.ced$scMRMA_manual <- gsub("Bergmann glia", "Glial cell",dat.ced$scMRMA_manual)
dat.ced$scMRMA_manual <- gsub("Nuocytes", "Microglia",dat.ced$scMRMA_manual)
dat.ced$scMRMA_manual <- gsub("Hepatic stellate cells", "Pericytes",dat.ced$scMRMA_manual)
dat.ced$scMRMA_manual <- gsub("Glutaminergic neurons", "Neurons",dat.ced$scMRMA_manual)
dat.ced$scMRMA_manual <- gsub("Plasmacytoid dendritic cells", "Dendritic cells",dat.ced$scMRMA_manual)
dat.ced$scMRMA_manual <- gsub("Astrocytes", "Astrocyte",dat.ced$scMRMA_manual)
dat.ced$scMRMA_manual[dat.ced$SCT_snn_res.0.4 == "0"] <- "GL261" 
dat.ced$scMRMA_manual[dat.ced$SCT_snn_res.0.4 == "1"] <- "GL261" 
dat.ced$scMRMA_manual[dat.ced$SCT_snn_res.0.4 == "5"] <- "GL261" 
dat.ced$scMRMA_manual[dat.ced$SCT_snn_res.0.4 == "6"] <- "GL261" 
dat.ced$scMRMA_manual[dat.ced$SCT_snn_res.0.4 == "16"] <- "GL261" 

SCpubr::do_DimPlot(dat.ced, reduction = "umap", group.by = "scMRMA_manual", label = TRUE, repel = TRUE,
                           plot.axes = FALSE, legend.position = "none")

dat.sub$scMRMA_manual <- gsub("in vivo", "",dat.sub$scMRMA_manual)
dat.sub$scMRMA_manual <- gsub("in vitro", "",dat.sub$scMRMA_manual)
dat.sub$scMRMA_manual <- gsub("Glutaminergic neurons", "Neurons",dat.sub$scMRMA_manual)
dat.sub$scMRMA_manual <- gsub("Plasmacytoid dendritic cells", "Dendritic cells",dat.sub$scMRMA_manual)
dat.sub$scMRMA_manual <- gsub("Chondrocytes", "Fibroblasts",dat.sub$scMRMA_manual)

SCpubr::do_DimPlot(dat.sub, reduction = "umap", group.by = "scMRMA_manual", label = TRUE, repel = TRUE,
                           plot.axes = FALSE, legend.position = "none")


#########################
# Analysis of knockdown #
# Get knockdown values
features <- read.table('../../meta/features_48h_NTs.csv',sep=',',header=TRUE,row.names=1)

# subset on only GL261 cells
dat.sub.GL261 <- subset(dat.sub, idents=c("GL261"))
dat.sub.GL261 <- subset(dat.sub.GL261, cells=row.names(sgRNATblConcordant))

# Normalize data to CPM
dat.sub.GL261Norm <- dat.sub.GL261
DefaultAssay(dat.sub.GL261Norm) <- "RNA"
dat.sub.GL261Norm <- NormalizeData(dat.sub.GL261Norm, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE)

# build non targeting expression matrix 
targets <- sapply(strsplit(as.character(features$name),"_"), `[`, 1)
targets <- unique(targets)[1:48]
ntTargetTbl <- data.frame(row.names = targets, matrix(nrow = length(targets), ncol=length(smpid))) 
colnames(ntTargetTbl) <- smpid

# build bulk expression table for sgRNA on target cells 
kdTargetTbl <- data.frame(row.names = targets, matrix(nrow = length(targets), ncol=length(smpid))) 
colnames(kdTargetTbl) <- smpid

df = GetAssayData(dat.sub.GL261Norm, slot = 'data')
colnames(df) = sapply(strsplit(colnames(df),"-"), `[`, 1)
df <- data.frame(df)

for (i in 1:length(smpid)){
  ntcells = sgRNATblConcordant[sgRNATblConcordant$.id == smpid[i] & sgRNATblConcordant$GeneA == "non-targeting","cell_barcode"]
  dfNt = df[targets,intersect(ntcells,colnames(df))]
  ntTargetTbl[,i] <- apply(data.frame(dfNt),1,mean)+0.1 # pseudocount
  
  for (j in 1:length(targets)){
    targetcells = sgRNATblConcordant[sgRNATblConcordant$.id == smpid[i] & sgRNATblConcordant$GeneA == targets[j],"cell_barcode"]
    dfTarget = df[targets[j],intersect(targetcells,colnames(df))]
    kdTargetTbl[j,i] <- mean(as.numeric(dfTarget),na.rm = TRUE)+0.1 # pseudocount
  }
}

# get mean knockdown values
ntTargetTbl[is.na(ntTargetTbl)] <- 0
remainingRNATbl <- round(kdTargetTbl/ntTargetTbl,3)

# generate ceiling of 1.0
remainingRNATbl[remainingRNATbl >= 1.0 ] = 1.0

