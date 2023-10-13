# Explore the GBM43 dataset in order to understand what's going on

library(Seurat)
library(ggplot2)
library(SCpubr)
library(scMRMA)

###############################################################################
# Files #

OUTPUT_DIR <- "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_explore_outputs"
setwd(OUTPUT_DIR)

CELLRANGER_INPUT_DIR <- "/raleighlab/data1/liuj/gbm_perturb/cellranger_out_combined/gbm_pdx_perturb_gbm43_1"
SCREEN_DATA <- "/raleighlab/data1/liuj/crispr_screens/hs_rt_integrated/GBM43_LN18_T98G_zim3_phenotable.txt"
SCREEN_HITS <- "/raleighlab/data1/liuj/gbm_perturb/scripts/features_GBM43_LN18_RT_hits.csv"

###############################################################################
# Parameters #

UMAP.MIN.DIST = 0.7
CLUSTER.RES = 0.4

###############################################################################
# Import data

screen_data.table <- read.table(SCREEN_DATA, sep='\t', header=TRUE, row.names=1)
screen_hits.table <- read.table(SCREEN_HITS, sep=',', header=TRUE, row.names=1)

# Import Cellranger data to create Seurat objects

dirs.list <- list.files(CELLRANGER_INPUT_DIR)

gene_expr.list <- list()
guide_capture.list <- list()

for (dir in dirs.list) {
  cat("Processing:", dir, "\n")
  cellranger_dir_appended <- paste(CELLRANGER_INPUT_DIR, dir, sep = "/")
  df <- Read10X(data.dir = paste(cellranger_dir_appended, 
                                 "/outs/filtered_feature_bc_matrix/",sep = ''))
  gene_expr.list[[dir]] <- CreateSeuratObject(counts = df$`Gene Expression`, 
                                              project = dir)
  guide_capture.list[[dir]] <- CreateSeuratObject(counts = df$`CRISPR Guide Capture`, 
                                                  project = dir)
}
seurat.data <- merge(x = gene_expr.list[[1]], 
                     y = gene_expr.list[2:length(gene_expr.list)])

# Remove from memory

rm(gene_expr.list)
gc()

###############################################################################
# Run analysis finding clusters

# QC

seurat.data$percent.mito = PercentageFeatureSet(seurat.data, pattern = "-MT-")
qc.plot = VlnPlot(seurat.data, features = c("nFeature_RNA", 
                                            "nCount_RNA", 
                                            "percent.mito"), ncol = 3)
ggsave(filename = "qc_vln_plot.png", qc.plot, height = 20, width = 20)

# Subset data, run SCTransform, and run PCA

data <- subset(seurat.data, subset = nFeature_RNA > 200 & percent.mito < 10)
data <- SCTransform(data, verbose = TRUE)
data <- RunPCA(data, verbose = TRUE)

# PCA dimensions

elbow.plot = ElbowPlot(object = data,ndims = 50)
ggsave(filename = "elbow_plot.png", elbow.plot, height = 20, width = 20)

# Assign annotations

data$sorted <- data$orig.ident
Idents(data) <- data$sorted
data <- RenameIdents(data,  
                        'GBM43_CED_noRT_FACS_A' = 'FACS',
                        'GBM43_CED_noRT_FACS_B' = 'FACS',
                        'GBM43_CED_noRT_FACS_C' = 'FACS',
                        'GBM43_CED_noRT_FACS_D' = 'FACS',
                        'GBM43_CED_noRT_unsort_A' = 'unsorted',
                        'GBM43_CED_noRT_unsort_B' = 'unsorted',
                        'GBM43_CED_noRT_unsort_C' = 'unsorted',
                        'GBM43_CED_noRT_unsort_D' = 'unsorted',
                        'GBM43_CED_RT_FACS_A1' = 'FACS',
                        'GBM43_CED_RT_FACS_A2' = 'FACS',
                        'GBM43_CED_RT_FACS_B' = 'FACS',
                        'GBM43_CED_RT_FACS_C1' = 'FACS',
                        'GBM43_CED_RT_FACS_C2' = 'FACS',
                        'GBM43_CED_RT_FACS_D' = 'FACS',
                        'GBM43_CED_RT_unsort_A' = 'unsorted',
                        'GBM43_CED_RT_unsort_C' = 'unsorted')
data$sorted <- Idents(data)
table(data$sorted)

data$cond <- data$orig.ident
Idents(data) <- data$cond 
data <- RenameIdents(data,  
                        'GBM43_CED_noRT_FACS_A' = 'noRT',
                        'GBM43_CED_noRT_FACS_B' = 'noRT',
                        'GBM43_CED_noRT_FACS_C' = 'noRT',
                        'GBM43_CED_noRT_FACS_D' = 'noRT',
                        'GBM43_CED_noRT_unsort_A' = 'noRT',
                        'GBM43_CED_noRT_unsort_B' = 'noRT',
                        'GBM43_CED_noRT_unsort_C' = 'noRT',
                        'GBM43_CED_noRT_unsort_D' = 'noRT',
                        'GBM43_CED_RT_FACS_A1' = 'RT',
                        'GBM43_CED_RT_FACS_A2' = 'RT',
                        'GBM43_CED_RT_FACS_B' = 'RT',
                        'GBM43_CED_RT_FACS_C1' = 'RT',
                        'GBM43_CED_RT_FACS_C2' = 'RT',
                        'GBM43_CED_RT_FACS_D' = 'RT',
                        'GBM43_CED_RT_unsort_A' = 'RT',
                        'GBM43_CED_RT_unsort_C' = 'RT')
data$cond <- Idents(data)

# Run and plot UMAP

data <- RunUMAP(data, reduction = "pca",dims = 1:30, min.dist = UMAP.MIN.DIST)
umap.plot <- SCpubr::do_DimPlot(data, reduction = "umap", group.by = "orig.ident", 
                   label = FALSE, repel = TRUE, plot.axes = TRUE, 
                   legend.position = "right",pt.size = 0.5)
ggsave("pheno_umap_sct_sub.png", umap.plot, width = 13, height = 9)

colors <- c("noRT" = "#3C5488FF",
            "unsorted" = "#5C88DAFF",
            "FACS" = "#FFCD00FF",
            "RT" = "#DC0000FF")

uamp.plot.sorted <- SCpubr::do_DimPlot(data, reduction = "umap", group.by = "sorted", 
                   label = FALSE, repel = TRUE, plot.axes = TRUE, 
                   legend.position = "bottom",pt.size = 0.5, colors.use = colors)
ggsave("pheno_umap_sct_sub_sorted.png", umap.plot.sorted, width = 10, height = 9)

umap.plot.cond <- SCpubr::do_DimPlot(data, reduction = "umap", group.by = "cond", 
                   label = FALSE, repel = TRUE, plot.axes = TRUE, 
                   legend.position = "bottom",pt.size = 0.5, colors.use = colors)
ggsave("pheno_umap_sct_sub_cond.png", umap.plot.cond, width = 10, height = 9)

# Find clusters

data <- FindNeighbors(data, dims = 1:30)
data <- FindClusters(data, resolution = CLUSTER.RES)

plot <- SCpubr::do_DimPlot(data, reduction = "umap", 
                           group.by = paste("SCT_snn_res.",res,collapse="",sep=""), 
                           label = TRUE, repel = TRUE,
                           plot.axes = TRUE, legend.position = "none")
ggsave(paste("pheno_sct_sub_umap_mindist_",mindist,"_res",res,".png",collapse="",sep=""), 
       plot, width =12, height = 9)

# Find markers

Idents(data) <- data$SCT_snn_res.0.4
markers.sub <- FindAllMarkers(object = data, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
markers.sub %>% group_by(cluster) %>% top_n(2, avg_log2FC)

# Generate a marker heatmap

top10.sub <- markers.sub %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top25.sub <- markers.sub %>% group_by(cluster) %>% top_n(25, avg_log2FC)
write.table(top25.sub,paste("pheno_top25markers_res_sub",res,".txt",sep=''),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)

plot <- DoHeatmap(subset(data, downsample = 200), features = top10.sub$gene, size = 3) + NoLegend()
ggsave("pheno_sub_top10_heatmap.pdf", plot, width = 12, height = 24)

###############################################################################
# scMRMA annotations

# Subset for only human genes

features.all <- rownames(data)
features.hs <- features.all[grep("GRCh38-", features.all)]
counts <- GetAssayData(data, assay = "RNA", slot = "counts")
counts <- counts[features.hs, ]
rownames(counts) <- gsub("GRCh38-", "", rownames(counts))
data.hs <- CreateSeuratObject(counts = counts, meta.data = data@meta.data)

# Rerun normalization followed by UMAP

data.hs <- SCTransform(data.hs, verbose = TRUE)
data.hs <- RunPCA(data.hs, verbose = TRUE)
data.hs <- RunUMAP(data.hs, reduction = "pca", 
                   dims = 1:30, min.dist = UMAP.MIN.DIST)

# Run scMRMA and plot the classifications on UMAP

result <- scMRMA(input = data,
                 species = "Hs",
)
data[["scMRMA"]] <- result$multiR$annotationResult[colnames(data), 
                                                   ncol(result$multiR$annotationResult)] # use the final assignment
SCpubr::do_DimPlot(data, reduction = "umap", group.by = "scMRMA_manual", 
                   label = TRUE, repel = TRUE, plot.axes = TRUE, 
                   legend.position = "none",pt.size = 0.5)

# Get markers for scMRMA cell types 

Idents(dat.sub) <- dat.sub$scMRMA_manual
markers.sub <- FindAllMarkers(object = dat.sub, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
markers.sub %>% group_by(cluster) %>% top_n(2, avg_log2FC)

top10.sub <- markers.sub %>% group_by(cluster) %>% top_n(10, avg_log2FC)
DoHeatmap(subset(dat.sub, downsample = 200), features = top10.sub$gene, size = 3) + NoLegend()

###############################################################################
# Run analysis to associate the cells with their appropriate sgRNAs

sgR.all <- merge(x = sgR.list[[1]], y = sgR.list[2:length(sgR.list)])

sgR.sub <- subset(sgR.all, nFeature_RNA > 0)
sgR.sub <- GetAssayData(sgR.sub, slot="counts")
sgR.sub <- data.frame(sgR.sub)
sgR.sub <- sgR.sub[apply(sgR.sub,1,sum) > 0,]
colnames(sgR.sub) <- gsub("\\.","-",colnames(sgR.sub))

seqMetrics.list <- list()
sgRNATbl.list <- list()

for (i in 1:length(smpid)){
  seqMetrics.list[[i]] <- read.table(paste("/raleighlab/data1/liuj/gbm_perturb/cellranger_out_combined/gbm_pdx_perturb_gbm43_1/",smpid[i],"/outs/metrics_summary.csv",sep=''),sep=',',header=TRUE)
  sgRNA_tags <- read.table(paste("/raleighlab/data1/liuj/gbm_perturb/cellranger_out_combined/gbm_pdx_perturb_gbm43_1/",smpid[i],"/outs/crispr_analysis/protospacer_calls_per_cell.csv",sep=''),sep=',',header=TRUE)
  sgRNA_tags$cell_barcode <- gsub("-1","",sgRNA_tags$cell_barcode)
  row.names(sgRNA_tags) <- as.character(sgRNA_tags$cell_barcode)
  sgRNATbl.list[[i]] <- sgRNA_tags
}

# histogram of sgRNA dectections
names(sgRNATbl.list) <- smpid
sgRNATbl <- ldply(sgRNATbl.list, data.frame)

# categorize sample types
sgRNATbl$sorted <- sgRNATbl$.id
sgRNATbl$sortedCond <- sgRNATbl$.id
smpTbl <- data.frame(orig.ident = dat.sub$orig.ident, sorted = dat.sub$sorted, cond = dat.sub$cond)
smpTbl <- unique(smpTbl)
row.names(smpTbl) <- as.character(smpTbl$orig.ident)
smpTbl$sortedCond <- paste(smpTbl$sorted,smpTbl$cond,sep="_")

sgRNATbl$sorted <- smpTbl[sgRNATbl$.id,"sorted"]
sgRNATbl$sortedCond <- smpTbl[sgRNATbl$.id,"sortedCond"]

plot <- ggplot(sgRNATbl) + 
  geom_histogram(aes(x=num_features, fill=.id),color=NA,show.legend = FALSE,size=0.25,binwidth=1) +
  scale_x_continuous(limits = c(0,10),breaks = seq(0,10,1),expand=c(0.01,0.01)) +
  #scale_y_continuous(expand=c(0.01,0.01), limits = c(0,9000)) +
  labs(title = "Number of sgRNAs detected", x="sgRNAs", y="cells") +
  facet_wrap(~.id, scales="free", ncol = 5) +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0.5,vjust=1),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black',fill=NA, size=0.25),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave("pdx_perturb_1_unique_sgRNAs_hist_sub.pdf", plot, width =8, height = 6)

plot <- ggplot(sgRNATbl) + 
  geom_histogram(aes(x=num_features, fill=sortedCond),color=NA,show.legend = FALSE,size=0.25,binwidth=1) +
  scale_x_continuous(limits = c(0,10),breaks = seq(0,10,1),expand=c(0.01,0.01)) +
  #scale_y_continuous(expand=c(0.01,0.01), limits = c(0,9000)) +
  labs(title = "Number of sgRNAs detected", x="sgRNAs", y="cells") +
  facet_wrap(~sortedCond, scales="free_y", ncol = 2) +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=0,hjust=0.5,vjust=1),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black',fill=NA, size=0.25),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave("pdx_perturb_1_unique_sgRNAs_hist_sortedCond.pdf", plot, width =4, height = 3.5)

# Tabulate sgRNA read data - what proportion of each sample has cells with sgRNAs?
sgR.subTotal <- apply(sgR.sub,2,sum)
dat.sub$sgRNA_UMI <- sgR.subTotal[colnames(dat.sub)]
dat.sub$sgRNA_logUMI <- log(sgR.subTotal[colnames(dat.sub)],10)
dat.sub$sgRNA_UMI[is.na(dat.sub$sgRNA_UMI)] <- 0
dat.sub$sgRNA_logUMI[is.na(dat.sub$sgRNA_logUMI)] <- 0
dat.sub$sgRNA_binary <- dat.sub$sgRNA_UMI > 0

plot <- SCpubr::do_FeaturePlot(dat.sub, features = "sgRNA_logUMI", plot.title = " sgRNA UMI", order = TRUE, split.by.idents ='TRUE')
ggsave("integrate_XFP_sgRNA_logUMI.png", plot = plot, width = 12, height = 9, units = 'in', dpi = 300)

plot <- SCpubr::do_FeaturePlot(dat.sub, features = "sgRNA_logUMI", plot.title = " sgRNA UMI", order = TRUE, split.by='sorted')
ggsave("integrate_XFP_sgRNA_logUMI_sorted.png", plot = plot, width = 24, height = 10, units = 'in', dpi = 300)

goi <- c("mm10---Cd44","mm10---Cd74","mm10---Apoe","mm10---C1qa","mm10---Gfap","mm10---Gzmb","mm10---Lyz2","GRCh38-NF1","GRCh38-EGFR","GRCh38-NRAS","GRCh38-RPS2","GRCh38-DDX5","GRCh38-COL1A1","GRCh38-RPLP1")
plot <- SCpubr::do_FeaturePlot(dat.sub, features = goi, plot.title = "GBM43 marker genes",order = TRUE,ncol = 5)
ggsave("integrate_marker_expression.png", plot = plot, width = 24, height = 21, units = 'in', dpi = 300)

# proportion of samples with sgRNAs detected 
df <- PrctCellExpringGene(dat.sub ,genes =c("sgRNA_UMI","GRCh38-EGFR","mm10---Rplp1"), group.by = "orig.ident")

# manually fill in sgRNA UMI proportion
Idents(dat.sub) <- dat.sub$orig.ident
factors = unique(df$Feature)

for(i in 1:length(factors)){
  sub <- subset(dat.sub,ident=factors[i])
  df[df$Feature == factors[i] & df$Markers == "sgRNA_UMI",2] <- sum(sub$sgRNA_UMI>0)/nrow(sub)
}

plot <- ggplot(df) +
  geom_bar(aes(x = Feature, y = Cell_proportion, fill = Markers),stat = 'identity') +
  labs(title = "", x="", y="Cell_proportion") +
  scale_fill_brewer(palette="Dark2") +
  facet_grid(cols = vars(df$Markers)) +
  theme_few() +
  coord_flip() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave("sgRNA_UMI_marker_proportion.pdf", plot, width = 6, height = 3, useDingbats = FALSE, limitsize=FALSE)

#####
# Which cell types are infected with sgRNA_UMI? #
cellTypeInfxTbl <- data.frame(table(data.frame(sgRNA = dat.sub$sgRNA_binary, cellType = dat.sub$SCT_snn_res.0.4, sorted =  dat.sub$sorted)))
cellTypeInfxTbl <- cellTypeInfxTbl[cellTypeInfxTbl$sgRNA == TRUE,]
cellTypeInfxTbl$sorted <- factor(cellTypeInfxTbl$sorted,levels=c("FACS","unsorted"))
cellTypeInfxTbl$cellType <- factor(cellTypeInfxTbl$cellType,levels=rev(as.character(unique(cellTypeInfxTbl$cellType))))

plot <- ggplot(cellTypeInfxTbl) +
  geom_bar(aes(x = cellType, y = Freq),fill="#3C5488FF",stat = 'identity') +
  labs(title = "", x="", y="Cells") +
  theme_few() +
  coord_flip() +
  facet_grid(cols = vars(cellTypeInfxTbl$sorted),scales="free") +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_rect(colour='black'))
ggsave("GBM43_pdx_infected_celltypes.pdf", plot, width = 2.2, height = 3, useDingbats = FALSE, limitsize=FALSE)

# number of concordant sgRNAs 
sgRNATblConcordant <- sgRNATbl[sgRNATbl$num_features == 2,]
sgRNATblConcordant <- sgRNATblConcordant %>% separate(feature_call, c("sgRNA_A", "sgRNA_B"), sep = "\\|")
sgRNATblConcordant$GeneA <- sapply(strsplit(sgRNATblConcordant$sgRNA_A,"_"), `[`, 1)
sgRNATblConcordant$GeneB <- sapply(strsplit(sgRNATblConcordant$sgRNA_B,"_"), `[`, 1)
sgRNATblConcordant <- sgRNATblConcordant[sgRNATblConcordant$GeneA == sgRNATblConcordant$GeneB,]

# Add cell names with barcode suffixes and use them as row names
smpidNumTbl <- data.frame(smpid = smpid, num = seq(1,length(smpid),1))
row.names(smpidNumTbl) <- as.character(smpidNumTbl$smpid)
sgRNATblConcordant$full_barcode <- paste(sgRNATblConcordant$cell_barcode,"-1_",smpidNumTbl[sgRNATblConcordant$.id,"num"],sep='')
row.names(sgRNATblConcordant) <- as.character(sgRNATblConcordant$full_barcode)

# Add total sgRNA UMIs and A and B counts - 5/4/23
sgRNATblConcordant$total_umis <- sgR.subTotal[row.names(sgRNATblConcordant)]
sgRNATblConcordant$GeneA_umi <- as.numeric(sapply(strsplit(sgRNATblConcordant$num_umis,"\\|"), `[`, 1))
sgRNATblConcordant$GeneB_umi <- as.numeric(sapply(strsplit(sgRNATblConcordant$num_umis,"\\|"), `[`, 2))
sgRNATblConcordant$sorted <- factor(sgRNATblConcordant$sorted,levels=c("FACS","unsorted"))

ggplot(sgRNATblConcordant)+
  geom_point(aes(x=log(GeneA_umi + GeneB_umi,10), y=log(total_umis,10), color = sorted),size=0.2,alpha=1,shape=16)+
  geom_abline(linetype = "dashed",linewidth=0.25)+
  #scale_color_brewer(palette="Dark2") + 
  scale_color_manual(values=c("#5C88DAFF","#FFCD00FF")) + 
  theme_few() +
  facet_grid(~sorted) + 
  stat_cor(aes(x=log(GeneA_umi + GeneB_umi,10), y=log(total_umis,10), color=sorted), label.x = 0, digits=5,size=1)+
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.ticks=element_line(colour='black',linewidth=0.25),panel.border=element_rect(colour='black')) +
  labs(title = "UMI Detection")
ggsave('GBM43_totalvAB_UMI_scatter.pdf',width=4,height=1.8,useDingbats=FALSE)

ggplot(sgRNATblConcordant)+
  geom_point(aes(x=log(GeneA_umi,10), y=log(GeneB_umi,10), color = sorted),size=0.2,alpha=0.5,shape=16)+
  geom_abline(linetype = "dashed",linewidth=0.25)+
  scale_color_brewer(palette="Dark2") + 
  theme_few() +
  facet_grid(~sorted) + 
  stat_cor(aes(x=log(GeneA_umi,10), y=log(GeneB_umi,10), color=sorted), label.x = 0, digits=5,size=1)+
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.ticks=element_line(colour='black',linewidth=0.25),panel.border=element_rect(colour='black')) +
  labs(title = "UMI Detection")
ggsave('GBM43_AvB_UMI_scatter.pdf',width=4,height=1.8,useDingbats=FALSE)



