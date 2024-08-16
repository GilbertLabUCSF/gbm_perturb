## Perturb-seq of malignant cells - preprocessing analysis 

library(Seurat)
library(ggplot2)
library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(SCpubr)
library(scMRMA)
library(ggpubr)


# load all GL261  data
smpid <- dir(file.path("data/cellranger_out_XFP/invivo_perturb_GL261_integrate"))

dataGene.list <- list()
sgR.list <- list()

for (i in 1:length(smpid)){
  df <- Read10X(data.dir = paste("data/cellranger_out_XFP/invivo_perturb_GL261_integrate/",smpid[i],"/outs/filtered_feature_bc_matrix/",sep=''))
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

# PCA dimensions
ElbowPlot(object = dat.sub,ndims = 50) 

# assign annotations to samples
dat.sub$sorted <- dat.sub$orig.ident
Idents(dat.sub) <- dat.sub$sorted 
dat.sub <- RenameIdents(dat.sub,  'GL261_48hit_noRT_1' = 'invitro',
                                  'GL261_48hit_noRT_2' = 'invitro',
                                  'GL261_48hit_RT_1' = 'invitro',
                                  'GL261_48hit_RT_2' = 'invitro',
                                  'GL261_CED_pool' = 'unsorted',
                                  'GL261_CED_RT_pool' = 'unsorted',
                                  'GL261_noRT_CED' = 'unsorted',
                                  'GL261_noRT_CED_MACSFACS' = 'MACSFACS',
                                  'GL261_noRT_preinf' = 'unsorted',
                                  'GL261_noRT_preinf_MACSFACS_1' = 'MACSFACS',
                                  'GL261_noRT_preinf_MACSFACS_2' = 'MACSFACS',
                                  'GL261_noRT_preinf_MACSFACS_3' = 'MACSFACS',
                                  'GL261_RT_CED' = 'unsorted',
                                  'GL261_RT_CED_MACSFACS' = 'MACSFACS',
                                  'GL261_RT_preinf' = 'unsorted',
                                  'GL261_RT_preinf_MACSFACS_1' = 'MACSFACS',
                                  'GL261_RT_preinf_MACSFACS_2' = 'MACSFACS',
                                  'GL261_RT_preinf_MACSFACS_3' = 'MACSFACS')

dat.sub$sorted <- Idents(dat.sub)

dat.sub$source <- dat.sub$orig.ident
Idents(dat.sub) <- dat.sub$source 
dat.sub <- RenameIdents(dat.sub,  'GL261_48hit_noRT_1' = 'invitro',
                                  'GL261_48hit_noRT_2' = 'invitro',
                                  'GL261_48hit_RT_1' = 'invitro',
                                  'GL261_48hit_RT_2' = 'invitro',
                                  'GL261_CED_pool' = 'CED',
                                  'GL261_CED_RT_pool' = 'CED',
                                  'GL261_noRT_CED' = 'CED',
                                  'GL261_noRT_CED_MACSFACS' = 'CED',
                                  'GL261_noRT_preinf' = 'preinf',
                                  'GL261_noRT_preinf_MACSFACS_1' = 'preinf',
                                  'GL261_noRT_preinf_MACSFACS_2' = 'preinf',
                                  'GL261_noRT_preinf_MACSFACS_3' = 'preinf',
                                  'GL261_RT_CED' = 'CED',
                                  'GL261_RT_CED_MACSFACS' = 'CED',
                                  'GL261_RT_preinf' = 'preinf',
                                  'GL261_RT_preinf_MACSFACS_1' = 'preinf',
                                  'GL261_RT_preinf_MACSFACS_2' = 'preinf',
                                  'GL261_RT_preinf_MACSFACS_3' = 'preinf')

dat.sub$source <- Idents(dat.sub)

dat.sub$cond <- dat.sub$orig.ident
Idents(dat.sub) <- dat.sub$cond 
dat.sub <- RenameIdents(dat.sub, 'GL261_48hit_noRT_1' = 'noRT',
                        'GL261_48hit_noRT_2' = 'noRT',
                        'GL261_48hit_RT_1' = 'RT',
                        'GL261_48hit_RT_2' = 'RT',
                        'GL261_CED_pool' = 'noRT',
                        'GL261_CED_RT_pool' = 'noRT',
                        'GL261_noRT_CED' = 'noRT',
                        'GL261_noRT_CED_MACSFACS' = 'noRT',
                        'GL261_noRT_preinf' = 'noRT',
                        'GL261_noRT_preinf_MACSFACS_1' = 'noRT',
                        'GL261_noRT_preinf_MACSFACS_2' = 'noRT',
                        'GL261_noRT_preinf_MACSFACS_3' = 'noRT',
                        'GL261_RT_CED' = 'RT',
                        'GL261_RT_CED_MACSFACS' = 'RT',
                        'GL261_RT_preinf' = 'RT',
                        'GL261_RT_preinf_MACSFACS_1' = 'RT',
                        'GL261_RT_preinf_MACSFACS_2' = 'RT',
                        'GL261_RT_preinf_MACSFACS_3' = 'RT')
dat.sub$cond <- Idents(dat.sub)


# Run UMAP embeddings
colors <- c("invitro" = "#E64B35FF",
            "preinf" = "#4DBBD5FF",
            "CED" = "#E18727FF",
            "noRT" = "#3C5488FF",
            "unsorted" = "#5C88DAFF",
            "MACSFACS" = "#FFCD00FF",
            "RT" = "#DC0000FF")

mindist=0.7
dat.sub <- RunUMAP(dat.sub, reduction = "pca",dims = 1:30, min.dist = mindist)
SCpubr::do_DimPlot(dat.sub, reduction = "umap", group.by = "orig.ident", label = FALSE, repel = TRUE, plot.axes = TRUE, legend.position = "right",pt.size = 0.5)
SCpubr::do_DimPlot(dat.sub, reduction = "umap", group.by = "sorted", label = FALSE, repel = TRUE, plot.axes = TRUE, legend.position = "bottom",pt.size = 0.5,colors.use = colors)
SCpubr::do_DimPlot(dat.sub, reduction = "umap", group.by = "source", label = FALSE, repel = TRUE, plot.axes = TRUE, legend.position = "bottom",pt.size = 0.5,colors.use = colors)
SCpubr::do_DimPlot(dat.sub, reduction = "umap", group.by = "cond", label = FALSE, repel = TRUE, plot.axes = TRUE, legend.position = "bottom",pt.size = 0.5,colors.use = colors)


# scMRMA annotation of cell types
result <- scMRMA(input = dat.sub,
                 species = "Mm",
)

dat.sub[["scMRMA"]] <- result$multiR$annotationResult[colnames(dat.sub),ncol(result$multiR$annotationResult)] # use the final assignment

# Manually addend scMRMA annotations 
dat.sub$scMRMA_manual <- dat.sub$scMRMA
dat.sub$scMRMA_manual <- gsub("Radial glial cell", "GL261 in vivo",dat.sub$scMRMA_manual)
dat.sub$scMRMA_manual <- gsub("Enteric glia cells", "GL261 in vitro",dat.sub$scMRMA_manual)
dat.sub$scMRMA_manual <- gsub("Bergmann glia", "Astrocytes 2",dat.sub$scMRMA_manual)
dat.sub$scMRMA_manual <- gsub("Nuocytes", "Microglia 2",dat.sub$scMRMA_manual)

SCpubr::do_DimPlot(dat.sub, reduction = "umap", group.by = "scMRMA_manual", label = TRUE, repel = TRUE, plot.axes = TRUE, legend.position = "none",pt.size = 0.5)


# Get markers for scMRMA cell types 
Idents(dat.sub) <- dat.sub$scMRMA_manual
markers.sub <- FindAllMarkers(object = dat.sub, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
markers.sub %>% group_by(cluster) %>% top_n(2, avg_log2FC)

# marker heatmap
top10.sub <- markers.sub %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DoHeatmap(subset(dat.sub, downsample = 200), features = top10.sub$gene, size = 3) + NoLegend()


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
  seqMetrics.list[[i]] <- read.table(paste("data/cellranger_out_XFP/invivo_perturb_GL261_integrate/",smpid[i],"/outs/metrics_summary.csv",sep=''),sep=',',header=TRUE)
  sgRNA_tags <- read.table(paste("data/cellranger_out_XFP/invivo_perturb_GL261_integrate/",smpid[i],"/outs/crispr_analysis/protospacer_calls_per_cell.csv",sep=''),sep=',',header=TRUE)
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

###############################
# subset on only GL261 cells ##
Idents(dat.sub) <- dat.sub$scMRMA_manual
dat.sub.GL261 <- subset(dat.sub, idents=c("GL261 in vivo","GL261 in vitro"))
dat.sub.GL261 <- subset(dat.sub.GL261, cells=row.names(sgRNATblConcordant))
  
# rerun normalization and UMAP on subset and give to David Wu
dat.sub.GL261 <- SCTransform(dat.sub.GL261, verbose = FALSE)
dat.sub.GL261 <- RunPCA(dat.sub.GL261, verbose = FALSE)

# Run UMAP embeddings
mindist=0.7
dat.sub.GL261 <- RunUMAP(dat.sub.GL261, reduction = "pca",dims = 1:30, min.dist = mindist)
SCpubr::do_DimPlot(dat.sub.GL261, reduction = "umap", group.by = "orig.ident", label = FALSE, repel = TRUE, plot.axes = TRUE, legend.position = "right")

# summary of concordant hits - GL261 cells only
sgRNATblConcordantGL261 <- sgRNATblConcordant[colnames(dat.sub.GL261),] 
df <- as.matrix(table(sgRNATblConcordantGL261[,c(".id","GeneA")]))


#########################
# Analysis of knockdown #
# Get knockdown values
features <- read.table('../../meta/features_48h_NTs.csv',sep=',',header=TRUE,row.names=1)

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





