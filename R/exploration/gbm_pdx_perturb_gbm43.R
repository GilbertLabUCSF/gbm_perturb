## Analysis of CED perturb-seq in GBM43 8/17/23
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(reshape2)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(SCpubr)
library(scMRMA)
library(ggpubr)
library(dplyr)
library(tidyr)
library(reshape2)

setwd('/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_explore_outputs/human_reference_only_analysis/')

# tabulate how many cells have x gene expression
PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}

# Load screen data
phenoTbl <- read.table('/raleighlab/data1/liuj/crispr_screens/hs_rt_integrated/GBM43_LN18_T98G_zim3_phenotable.txt',sep='\t',header=TRUE,row.names=1)

# load all GBM43  data
smpid <- dir(file.path("/raleighlab/data1/liuj/gbm_perturb/cellranger_out_combined/gbm_pdx_perturb_gbm43_human"))

dataGene.list <- list()
sgR.list <- list()

for (i in 1:length(smpid)){
  print(smpid[i])
  print(paste("/raleighlab/data1/liuj/gbm_perturb/cellranger_out_combined/gbm_pdx_perturb_gbm43_human/",smpid[i],"/outs/filtered_feature_bc_matrix/",sep=''))
  df <- Read10X(data.dir = paste("/raleighlab/data1/liuj/gbm_perturb/cellranger_out_combined/gbm_pdx_perturb_gbm43_human/",smpid[i],"/outs/filtered_feature_bc_matrix/",sep=''))
  # dataGene.list[[i]] <- CreateSeuratObject(counts=df$`Gene Expression`, project = smpid[i])
  # sgR.list[[i]] <- CreateSeuratObject(counts = df$`CRISPR Guide Capture`, project = smpid[i])
  # dataGene.list[[i]]$orig.ident <- smpid[i]
}

names(dataGene.list) <- smpid
names(sgR.list) <- smpid

dat.all <- merge(x = dataGene.list[[1]], y = dataGene.list[2:length(dataGene.list)])
rm(dataGene.list)
gc()

# quality metrics for samples divided by ID 
mito.features <- grep(pattern = "-MT-", x = rownames(x = dat.all), value = TRUE, ignore.case = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = dat.all, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = dat.all, slot = 'counts'))
dat.all[['percent.mt']] <- percent.mito

VlnPlot(dat.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
ggsave(filename = "qc_violin_sub.pdf",width = 20, height=10)

# subset data and normalize
dat.sub <- subset(dat.all, subset = nFeature_RNA > 200)
dat.sub <- SCTransform(dat.sub, verbose = FALSE)
dat.sub <- RunPCA(dat.sub, verbose = FALSE)

# PCA dimensions
ElbowPlot(object = dat.sub,ndims = 50) 

# assign annotations to samples
dat.sub$sorted <- dat.sub$orig.ident
Idents(dat.sub) <- dat.sub$sorted 
dat.sub <- RenameIdents(dat.sub,  
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
dat.sub$sorted <- Idents(dat.sub)
table(dat.sub$sorted)

dat.sub$cond <- dat.sub$orig.ident
Idents(dat.sub) <- dat.sub$cond 
dat.sub <- RenameIdents(dat.sub,  
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
dat.sub$cond <- Idents(dat.sub)

# Run UMAP embeddings
mindist=0.7
dat.sub <- RunUMAP(dat.sub, reduction = "pca",dims = 1:30, min.dist = mindist)
SCpubr::do_DimPlot(dat.sub, reduction = "umap", group.by = "orig.ident", label = FALSE, repel = TRUE, plot.axes = TRUE, legend.position = "right",pt.size = 0.5)
ggsave("pheno_umap_sct_sub.png", width =13, height = 9)

colors <- c("noRT" = "#3C5488FF",
            "unsorted" = "#5C88DAFF",
            "FACS" = "#FFCD00FF",
            "RT" = "#DC0000FF")

SCpubr::do_DimPlot(dat.sub, reduction = "umap", group.by = "sorted", label = FALSE, repel = TRUE, plot.axes = TRUE, legend.position = "bottom",pt.size = 0.5,colors.use = colors)
ggsave("pheno_umap_sct_sub_sorted.png", width =10, height = 9)

SCpubr::do_DimPlot(dat.sub, reduction = "umap", group.by = "cond", label = FALSE, repel = TRUE, plot.axes = TRUE, legend.position = "bottom",pt.size = 0.5,colors.use = colors)
ggsave("pheno_umap_sct_sub_cond.png", width =10, height = 9)

# find clusters, what is the heterogeneity due to?
dat.sub <- FindNeighbors(object = dat.sub, dims = 1:30) 

res=0.4
dat.sub <- FindClusters(object = dat.sub, reduction = "umap", resolution = res)

# plot clusters in UMAP space 
plot <- SCpubr::do_DimPlot(dat.sub, reduction = "umap", group.by = paste("SCT_snn_res.",res,collapse="",sep=""), label = TRUE, repel = TRUE,
                           plot.axes = TRUE, legend.position = "none")
ggsave(paste("pheno_sct_sub_umap_mindist_",mindist,"_res",res,".png",collapse="",sep=""), plot, width =12, height = 9)

# get markers
Idents(dat.sub) <- dat.sub$SCT_snn_res.0.4
markers.sub <- FindAllMarkers(object = dat.sub, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
markers.sub %>% group_by(cluster) %>% top_n(2, avg_log2FC)

# marker heatmap
top10.sub <- markers.sub %>% group_by(cluster) %>% top_n(10, avg_log2FC)
top25.sub <- markers.sub %>% group_by(cluster) %>% top_n(25, avg_log2FC)

write.table(top25.sub,paste("pheno_top25markers_res_sub",res,".txt",sep=''),row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)

# marker heatmap plot
plot <- DoHeatmap(subset(dat.sub, downsample = 200), features = top10.sub$gene, size = 3) + NoLegend()
ggsave("pheno_sub_top10_heatmap.pdf", plot, width = 12, height = 24)

# Analysis of sgRNA coverage 
##########################
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
  seqMetrics.list[[i]] <- read.table(paste("/raleighlab/data1/liuj/gbm_perturb/cellranger_out_combined/gbm_pdx_perturb_gbm43_human/",smpid[i],"/outs/metrics_summary.csv",sep=''),sep=',',header=TRUE)
  sgRNA_tags <- read.table(paste("/raleighlab/data1/liuj/gbm_perturb/cellranger_out_combined/gbm_pdx_perturb_gbm43_human/",smpid[i],"/outs/crispr_analysis/protospacer_calls_per_cell.csv",sep=''),sep=',',header=TRUE)
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


###############################
# subset on only GBM43 cells ##
Idents(dat.sub) <- dat.sub$sorted
dat.sub.GBM43 <- subset(dat.sub, idents=c("FACS"))
Idents(dat.sub.GBM43) <- dat.sub.GBM43$SCT_snn_res.0.4
dat.sub.GBM43 <- subset(dat.sub.GBM43, idents=c("1","2","6","9","12","17"))
dat.sub.GBM43 <- subset(dat.sub.GBM43, cells=row.names(sgRNATblConcordant))

# rerun normalization and UMAP
dat.sub.GBM43 <- SCTransform(dat.sub.GBM43, verbose = FALSE)
dat.sub.GBM43 <- RunPCA(dat.sub.GBM43, verbose = FALSE)

# Run UMAP embeddings
mindist=0.7
dat.sub.GBM43 <- RunUMAP(dat.sub.GBM43, reduction = "pca",dims = 1:30, min.dist = mindist)
SCpubr::do_DimPlot(dat.sub.GBM43, reduction = "umap", group.by = "orig.ident", label = FALSE, repel = TRUE, plot.axes = TRUE, legend.position = "right")
ggsave("GBM43_umap_sct_sub.png", width =13, height = 9)

SCpubr::do_DimPlot(dat.sub.GBM43, reduction = "umap", group.by = "sorted", label = FALSE, repel = TRUE, plot.axes = TRUE, legend.position = "bottom")
ggsave("GBM43_umap_sct_sub_sorted.png", width =10, height = 9)

SCpubr::do_DimPlot(dat.sub.GBM43, reduction = "umap", group.by = "cond", label = FALSE, repel = TRUE, plot.axes = TRUE, legend.position = "bottom")
ggsave("GBM43_umap_sct_sub_cond.png", width =10, height = 9)

# summary of concordant hits - GBM43 cells only
sgRNATblConcordantGBM43 <- sgRNATblConcordant[colnames(dat.sub.GBM43),]
write.table(table(sgRNATblConcordantGBM43[,c(".id","GeneA")]),"pdx_perturb_GBM43_concordant_sgRNAs.txt",sep='\t',row.names=TRUE,quote=FALSE)

# heatmap visualization of concordant sgRNA detections 
df <- as.matrix(table(sgRNATblConcordantGBM43[,c(".id","GeneA")]))

# Filter the phenotable
phenoTbl <- phenoTbl[grep(paste(colnames(df),collapse="$|^"),phenoTbl$GBM43.target),]
phenoTbl <- phenoTbl[order(phenoTbl$GBM43.gamma_pval),]
row.names(phenoTbl) <- make.unique(as.character(phenoTbl$GBM43.target))

sgAnnot <- data.frame(row.names = colnames(df), gamma = phenoTbl[colnames(df),"GBM43.gamma_avg"], tau = phenoTbl[colnames(df),"GBM43.tauRT_avg"], rho = phenoTbl[colnames(df),"GBM43.rhoRT_avg"])
sgAnnot["non-targeting",] <- c(0,0,0)

smpAnnot <- smpTbl[row.names(df),]

pdf("pdx_perturb_GBM43_concordant_sgRNAs_heatmap_sub_annot.pdf",width = 22, height=8)

Heatmap(log(df+1,10), name = "log10(Cells)",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.0f", df[i, j]), x, y, gp = gpar(fontsize = 5))},
        top_annotation = HeatmapAnnotation(gamma = anno_barplot(sgAnnot$gamma, bar_width = 0.9),
                                           tau = anno_barplot(sgAnnot$tau, bar_width = 0.9),
                                           rho = anno_barplot(sgAnnot$rho, bar_width = 0.9)),
        right_annotation = rowAnnotation(cond = smpAnnot$cond,
                                         source = smpAnnot$source,
                                         sort = smpAnnot$sorted,
                                         col = list(cond = c("noRT" = "#2166AC", "RT" = "#B2182B"),
                                                    sort = c("unsorted" = "#32BD9F","FACS" = "#FFCD00FF"))),
        #column_split = factor(sgRNAClassRT$cond, levels = c("NoRT","5x1.8")),
        row_split = smpAnnot$source,
        cluster_column_slice = FALSE,
        cluster_rows = FALSE,
        width = ncol(df)*unit(4, "mm"),
        height = nrow(df)*unit(4, "mm"),
        col = colorRamp2(quantile(log(df+1,10), seq(0, 1, by = 0.25)), viridis(5)),
        
)
dev.off()

# export seurat object for Exploration
dat.sub.GBM43$sgRNA <- sgRNATblConcordant[colnames(dat.sub.GBM43),"GeneA"]
dat.sub.GBM43$sgRNACond <- paste(dat.sub.GBM43$sgRNA,dat.sub.GBM43$cond,sep="_") 

# Export seurat object 
saveRDS(dat.sub.GBM43,file = "../GBM43_1_malignant_only_annotated_20230817.Rds")


#########################
# Analysis of knockdown #
# Get knockdown values
features <- read.table('/raleighlab/data1/liuj/gbm_perturb/scripts/features_GBM43_LN18_RT_hits.csv',sep=',',header=TRUE,row.names=1)

# Normalize data to CPM
dat.sub.GBM43Norm <- dat.sub.GBM43
DefaultAssay(dat.sub.GBM43Norm) <- "RNA"
dat.sub.GBM43Norm <- NormalizeData(dat.sub.GBM43Norm, normalization.method = "RC", scale.factor = 1e6, verbose = TRUE)

# build non targeting expression matrix 
targets <- sapply(strsplit(as.character(features$name),"_"), `[`, 1)
targets <- unique(targets)[grep("NTC",unique(targets),invert=TRUE)]
ntTargetTbl <- data.frame(row.names = targets, matrix(nrow = length(targets), ncol=length(smpid))) 
colnames(ntTargetTbl) <- smpid

# build bulk expression table for sgRNA on target cells 
kdTargetTbl <- data.frame(row.names = targets, matrix(nrow = length(targets), ncol=length(smpid))) 
colnames(kdTargetTbl) <- smpid

df = GetAssayData(dat.sub.GBM43Norm, slot = 'data')
colnames(df) = sapply(strsplit(colnames(df),"-"), `[`, 1)
df <- data.frame(df)
row.names(df) <- gsub("GRCh38-","",row.names(df))

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
#remainingRNATbl[remainingRNATbl >= 1.0 ] = 1.0

# Remove unsorted samples
remainingRNATbl <- remainingRNATbl[,grep("FACS",colnames(remainingRNATbl))]
remainingRNATbl <- remainingRNATbl[complete.cases(remainingRNATbl),]

# create heatmap of knockdown 
df <- melt(data.frame(gene=row.names(remainingRNATbl), remainingRNATbl))
df <- data.frame(df,smpTbl[as.character(df$variable),])

plot <- ggplot(df, aes(x = gene,y=variable)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = round(value,2)),size=1.5) +
  scale_fill_gradient(low = "white", high = "steelblue",limit=c(0,1),na.value="gray50") +
  labs(title = "Mean sgRNA knockdown - pseudobulk", x="", y="") +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_blank(),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave("pdx_perturb_meanKD_heatmap.pdf", plot, width =6, height = 2.5)

plot <- ggplot(df) +
  geom_histogram(aes(color = variable,x=value), fill="white",binwidth =0.05,show.legend = FALSE) +
  scale_x_continuous(expand=c(0.01,0.01), limits = c(-0.05,1.05), breaks = seq(0,1,0.25)) +
  labs(title = "Mean sgRNA knockdown - pseudobulk", x="mRNA remaining", y="Gene Targets") +
  facet_wrap(.~variable,ncol =6,) +
  theme_base() +
  theme(axis.text=element_text(size=6),text=element_text(size=6),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_blank(),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave("pdx_perturb_meanKD_histograms.pdf", plot, width =6, height = 4)

# Violin plots summary of KD values
plot <- ggplot(df, aes(x = cond,y=value,fill=cond)) +
  geom_violin(scale='width',position= position_dodge(1),linewidth=0.25)+
  geom_boxplot(color="white",width=0.1,outlier.colour = NA,position = position_dodge(1),linewidth=0.25)+
  scale_fill_manual(values=c("#3C5488FF", "#DC0000FF"))+
  labs(title = "Mean sgRNA knockdown - pseudobulk", x="", y="mRNA Remaining") +
  theme_base() +
  theme(axis.text=element_text(size=5),text=element_text(size=5),
        axis.ticks=element_line(colour='black',size=0.25),panel.border=element_blank(),panel.background=element_blank(),axis.line = element_line(size=0.25))
ggsave("pdx_perturb_meanKD_violin.pdf", plot, width =3.5, height = 2.5)

# quantify the knockdown
df %>% group_by(cond) %>% summarise(mean = mean(value, na.rm = TRUE), median=median(value, na.rm = TRUE))


