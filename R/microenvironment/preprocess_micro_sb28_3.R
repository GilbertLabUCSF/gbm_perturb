## Analysis of in VIVO perturb-seq - SB28 microenvironment H11
library(Seurat)
library(ggplot2)
library(CellChat)
library(circlize)
library(patchwork)
library(SCpubr)
options(stringsAsFactors = FALSE)


# load new SB28 in vivo data
smpid <- dir(file.path("/raleighlab/data1/liuj/gbm_perturb/cellranger_out_XFP/gbm_pdx_micro_sb28_3"))

dataGene.list <- list()
sgR.list <- list()

for (i in 1:length(smpid)){
  df <- Read10X(data.dir = paste("/raleighlab/data1/liuj/gbm_perturb/cellranger_out_XFP/gbm_pdx_micro_sb28_3/",smpid[i],"/outs/filtered_feature_bc_matrix/",sep=''))
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

# assign annotations to samples
dat.sub$sorted <- dat.sub$orig.ident
dat.sub$sorted <- gsub("SB28_micro_3_A_FACS","sort",dat.sub$sorted)
dat.sub$sorted <- gsub("SB28_micro_3_A_unsort","no_sort",dat.sub$sorted)
dat.sub$sorted <- gsub("SB28_micro_3_B_FACS","sort",dat.sub$sorted)
dat.sub$sorted <- gsub("SB28_micro_3_B_unsort","no_sort",dat.sub$sorted)
dat.sub$sorted <- gsub("SB28_micro_3_C_FACS","sort",dat.sub$sorted)
dat.sub$sorted <- gsub("SB28_micro_3_C_unsort","no_sort",dat.sub$sorted)
dat.sub$sorted <- gsub("SB28_micro_3_D_FACS","sort",dat.sub$sorted)
dat.sub$sorted <- gsub("SB28_micro_3_D_unsort","no_sort",dat.sub$sorted)        

# PCA dimensions
ElbowPlot(object = dat.sub,ndims = 50) 

# find clusters, what is the heterogeneity due to?
res=0.4
dat.sub <- FindNeighbors(object = dat.sub, dims = 1:30) 
dat.sub <- FindClusters(object = dat.sub, reduction = "umap", resolution = res)

# plot clusters in UMAP space 
DimPlot(object = dat.sub, reduction = "umap", pt.size=0.3, group.by = paste("SCT_snn_res.",res,collapse="",sep=""),label = TRUE,raster = FALSE)

# Run UMAP embeddings
mindist=0.7
dat.sub <- RunUMAP(dat.sub, reduction = "pca",dims = 1:30, min.dist = mindist)
DimPlot(object = dat.sub, group.by = 'orig.ident', reduction = "umap", pt.size=0.3,raster=FALSE)
DimPlot(object = dat.sub, group.by = 'sorted', reduction = "umap", pt.size=0.3,raster=FALSE)

colors <- c("no_sort" = "#5C88DAFF",
            "sort" = "#FFCD00FF")

SCpubr::do_DimPlot(dat.sub, reduction = "umap", group.by = "sorted", label = FALSE, repel = TRUE,
                           plot.axes = FALSE, legend.position = "bottom",colors.use = colors)

# Featureplot of exogenous genes
goi <- c("GFP","RFP","KRAB","Nras")
SCpubr::do_FeaturePlot(dat.sub, features = goi, plot.title = "SB28 microenvironment marker genes",order = TRUE,ncol = 4)

# scMRMA annotation of cell types
result <- scMRMA(input = dat.sub,species = "Mm")
dat.sub[["scMRMA"]] <- result$multiR$annotationResult[colnames(dat.sub),ncol(result$multiR$annotationResult)] # use the final assignment

SCpubr::do_DimPlot(dat.sub, reduction = "umap", group.by = "scMRMA", label = TRUE, repel = TRUE,
                           plot.axes = TRUE, legend.position = "none")
SCpubr::do_DimPlot(dat.sub, reduction = "umap", split.by = "scMRMA", label = TRUE, repel = TRUE, plot.axes = TRUE, legend.position = "none")


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
  seqMetrics.list[[i]] <- read.table(paste("/raleighlab/data1/liuj/gbm_perturb/cellranger_out_XFP/gbm_pdx_micro_sb28_3/",smpid[i],"/outs/metrics_summary.csv",sep=''),sep=',',header=TRUE)
  sgRNA_tags <- read.table(paste("/raleighlab/data1/liuj/gbm_perturb/cellranger_out_XFP/gbm_pdx_micro_sb28_3/",smpid[i],"/outs/crispr_analysis/protospacer_calls_per_cell.csv",sep=''),sep=',',header=TRUE)
  sgRNA_tags$cell_barcode <- gsub("-1","",sgRNA_tags$cell_barcode)
  row.names(sgRNA_tags) <- as.character(sgRNA_tags$cell_barcode)
  sgRNATbl.list[[i]] <- sgRNA_tags
}

# histogram of sgRNA dectections
names(sgRNATbl.list) <- smpid
sgRNATbl <- ldply (sgRNATbl.list, data.frame)

plot <- ggplot(sgRNATbl) + 
  geom_histogram(aes(x=num_features, fill=.id),color=NA,show.legend = FALSE,size=0.25,binwidth=1) +
  scale_x_continuous(limits = c(0,10),breaks = seq(0,10,1),expand=c(0.01,0.01)) +
  #scale_y_continuous(expand=c(0.01,0.01), limits = c(0,9000)) +
  labs(title = "Number of sgRNAs detected", x="sgRNAs", y="cells") +
  facet_wrap(~.id, scales="free", ncol = 4) +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0.5,vjust=1),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black',fill=NA, size=0.25),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave("pdx_perturb_micro_sb28_3_unique_sgRNAs_hist_sub.png", plot, width =6, height = 6)

# number of concordant sgRNAs - need to allow for the double sgRNAs to be included

sgRNATblConcordant <- sgRNATbl[sgRNATbl$num_features == 1,]
sgRNATblConcordant$GeneA <- sapply(strsplit(sgRNATblConcordant$feature_call,"_"), `[`, 1)

# Add cell names with barcode suffixes and use them as row names
smpidNumTbl <- data.frame(smpid = smpid, num = seq(1,length(smpid),1))
row.names(smpidNumTbl) <- as.character(smpidNumTbl$smpid)
sgRNATblConcordant$full_barcode <- paste(sgRNATblConcordant$cell_barcode,"-1_",smpidNumTbl[sgRNATblConcordant$.id,"num"],sep='')
row.names(sgRNATblConcordant) <- as.character(sgRNATblConcordant$full_barcode)

# Analysis of sgRNA UMIs 
sgR.subTotal <- apply(sgR.sub,2,sum)
dat.sub$sgRNA_UMI <- sgR.subTotal[colnames(dat.sub)]
dat.sub$sgRNA_logUMI <- log(sgR.subTotal[colnames(dat.sub)],10)
dat.sub$sgRNA_UMI[is.na(dat.sub$sgRNA_UMI)] <- 0
dat.sub$sgRNA_logUMI[is.na(dat.sub$sgRNA_logUMI)] <- 0
dat.sub$sgRNA_binary <- dat.sub$sgRNA_UMI > 0

SCpubr::do_FeaturePlot(dat.sub, features = "sgRNA_logUMI", plot.title = "SB28 micro sgRNA UMI", order = TRUE, split.by.idents ='TRUE')


# heatmap visualization of concordant sgRNA detections 
df <- as.matrix(table(sgRNATblConcordant[colnames(dat.subCore),c(".id","feature_call")]))
smpTbl <- data.frame(orig.ident = dat.sub$orig.ident, sorted = dat.sub$sorted)
smpTbl <- unique(smpTbl)
row.names(smpTbl) <- as.character(smpTbl$orig.ident)

Heatmap(log(df+1,10), name = "log10(Cells)",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.0f", df[i, j]), x, y, gp = gpar(fontsize = 5))},
        right_annotation = rowAnnotation(sort = smpTbl[row.names(df),"sorted"],
                                         col = list(sort = c("no_sort" = "#32BD9F","sort" = "#3288BD"))),
        #column_split = factor(sgRNAClassRT$cond, levels = c("NoRT","5x1.8")),
        cluster_columns = FALSE,
        width = ncol(df)*unit(4, "mm"),
        height = nrow(df)*unit(4, "mm"),
        col = colorRamp2(quantile(log(df+1,10), seq(0, 1, by = 0.25)), viridis(5)),
        
)

# Heatmap of which cell types got which sgRNAs
df <- as.matrix(table(sgRNATblConcordant[colnames(dat.subCore),c("feature_call","cellType")]))

Heatmap(log(df+1,10), name = "log10(Cells)",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.0f", df[i, j]), x, y, gp = gpar(fontsize = 5))},
        cluster_rows = FALSE,
        cluster_column_slice = FALSE,
        width = ncol(df)*unit(4, "mm"),
        height = nrow(df)*unit(4, "mm"),
        col = colorRamp2(quantile(log(df+1,10), seq(0, 1, by = 0.25)), viridis(5)),
        
)

df <- as.matrix(table(sgRNATblConcordant[colnames(dat.subCore),c("GeneA","cellType")]))
Heatmap(log(df+1,10), name = "log10(Cells)",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.0f", df[i, j]), x, y, gp = gpar(fontsize = 5))},
        cluster_rows = FALSE,
        cluster_column_slice = FALSE,
        width = ncol(df)*unit(4, "mm"),
        height = nrow(df)*unit(4, "mm"),
        col = colorRamp2(quantile(log(df+1,10), seq(0, 1, by = 0.25)), viridis(5)),
        
)

# Subset cells that are sorted and have sgRNA
dat.subCore <- subset(dat.sub,subset = sorted == "sort")
dat.subCore <- subset(dat.subCore,subset = sgRNA_binary == TRUE)
dat.subCore <- subset(dat.subCore,subset = sgRNA_full != "<NA>")


####################################################
# Analysis of knockdown at the per cell type level #
# Normalize data to CPM
DefaultAssay(dat.subCore) <- "RNA"
dat.subCore <- NormalizeData(dat.subCore, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE)

# Normalize to sgNTCs but only keep cell types that have >150 cells with sgRNAs 
coi <- names(table(sgRNATblConcordant$cellType)[table(sgRNATblConcordant$cellType) > 150])

# build non targeting expression matrix 
targets <- names(table(sgRNATblConcordant$feature_call))
targets <- unique(targets)[grep("sgNeg",unique(targets),invert=TRUE)]
ntTargetTbl <- data.frame(row.names = targets, matrix(nrow = length(targets), ncol=length(coi))) 
colnames(ntTargetTbl) <- coi

# build bulk expression table for sgRNA on target cells 
kdTargetTbl <- data.frame(row.names = targets, matrix(nrow = length(targets), ncol=length(ntTargetTbl))) 
colnames(kdTargetTbl) <- colnames(ntTargetTbl)

for (i in 1:ncol(kdTargetTbl)){
  ntcells = sgRNATblConcordant[sgRNATblConcordant$cellType == names(kdTargetTbl)[i] & sgRNATblConcordant$GeneA == "sgNegCtrl","full_barcode"]
  ntcells = ntcells[!is.na(ntcells)]
  df = GetAssayData(dat.subCore, slot = 'data')
  df = data.frame(df)
  dfNt = df[sapply(strsplit(targets,"_"), `[`, 1),intersect(names(df),gsub("-",".",ntcells))] # full sgRNA option
  ntTargetTbl[,i] <- apply(data.frame(dfNt),1,mean)
  
  for (j in 1:length(targets)){
    targetcells = sgRNATblConcordant[sgRNATblConcordant$cellType == names(kdTargetTbl)[i]  & sgRNATblConcordant$feature_call == targets[j],"full_barcode"]
    targetcells = targetcells[!is.na(targetcells)]
    dfTarget = df[sapply(strsplit(targets,"_"), `[`, 1)[j],intersect(names(df),gsub("-",".",targetcells))] # full sgRNA option
    
    kdTargetTbl[j,i] <- mean(as.numeric(dfTarget),na.rm = TRUE)
  }
}

# get mean knockdown values
remainingRNATbl <- round((kdTargetTbl+1)/(ntTargetTbl+1),3)
remainingRNATbl

# create heatmap of knockdown 
Heatmap(as.matrix(remainingRNATbl), name = "mRNA remaining",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", remainingRNATbl[i, j]), x, y, gp = gpar(fontsize = 5))},
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        cluster_column_slice = FALSE,
        width = ncol(remainingRNATbl)*unit(4, "mm"),
        height = nrow(remainingRNATbl)*unit(4, "mm"),
        col = colorRamp2(seq(0, 1, by = 0.25), rev(magma(5))),
        
)


