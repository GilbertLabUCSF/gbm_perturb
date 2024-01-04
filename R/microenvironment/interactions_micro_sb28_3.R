## Analysis of SB28 microenvironment using full set for cell-cell interactions 
library(patchwork)
library(Seurat)
library(EnhancedVolcano)
library(dplyr)
library(ggthemes)
library(ComplexHeatmap)
library(circlize)
options(stringsAsFactors = FALSE)


#############
# Cell chat #
cellchat <- createCellChat(object = dat.sub, group.by = "scMRMA_manual")
CellChatDB <- CellChatDB.mouse
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat,population.size = TRUE,type =  "truncatedMean", trim = 0.05)

cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)                           

groupSize <- as.numeric(table(cellchat@idents))

# Networks
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Look at the individual pathways 
cellchat@netP$pathways
netAnalysis_contribution(cellchat, signaling = cellchat@netP$pathways)
netVisual_bubble(cellchat, remove.isolate = TRUE)

#########################################
# Analysis of macrophages with sgPtprc #
cellTbl <- data.frame(barcode = colnames(dat.sub), celltype = dat.sub$scMRMA_manual, sgRNA_gene = dat.sub$sgRNA_gene, sgRNA_full = dat.sub$sgRNA_full)
cellTbl$celltype_sgRNA_full <- cellTbl$celltype
cellTbl$celltype_sgRNA_gene <- cellTbl$celltype
cellTbl <- cellTbl %>%
  mutate(celltype_sgRNA_full = case_when(
    celltype_sgRNA_full == "Macrophages" ~ paste(celltype_sgRNA_full, sgRNA_full, sep = "_"),
    TRUE ~ celltype_sgRNA_full
  ))

cellTbl <- cellTbl %>%
  mutate(celltype_sgRNA_gene = case_when(
    celltype_sgRNA_gene == "Macrophages" ~ paste(celltype_sgRNA_gene, sgRNA_gene, sep = "_"),
    TRUE ~ celltype_sgRNA_gene
  ))

table(cellTbl$celltype_sgRNA_gene)

keepCells <- unique(c(row.names(cellTbl[grep("Macrophages",cellTbl$celltype,invert=TRUE),]), row.names(cellTbl[grep("138175255.23-P1P2|138175267.23-P1P2|Macrophages_sgNegCtrl_3",cellTbl$celltype_sgRNA_full),])))
dat.sub.noPtprc <- subset(dat.sub, cells = keepCells)               
Idents(dat.sub.noPtprc) <- cellTbl[colnames(dat.sub.noPtprc),"celltype_sgRNA_gene"]
table(Idents(dat.sub.noPtprc))

cellchat2 <- createCellChat(object = dat.sub.noPtprc, group.by = "ident")
cellchat2@DB <- CellChatDB
cellchat2 <- subsetData(cellchat2)
cellchat2 <- identifyOverExpressedGenes(cellchat2)
cellchat2 <- identifyOverExpressedInteractions(cellchat2)
cellchat2 <- computeCommunProb(cellchat2,population.size = FALSE,type =  "truncatedMean", trim = 0.01)

cellchat2 <- filterCommunication(cellchat2, min.cells = 5)
df.net <- subsetCommunication(cellchat2)
cellchat2 <- computeCommunProbPathway(cellchat2)
cellchat2 <- aggregateNet(cellchat2)                           

groupSize <- as.numeric(table(cellchat2@idents))

########################################################################
# Calculate the differential interaction values in sender-receiver pairs
coi <- c("T cells","Neutrophils","Microglia","SB28","Astrocytes","Dendritic cells","Endothelial cells","Oligodendrocytes","Oligodendrocyte precursor cells")
dfSource <- netVisual_bubble(cellchat2, sources.use = c("Macrophages_sgNegCtrl","Macrophages_Ptprc"), targets.use = coi, remove.isolate = FALSE, sort.by.source = TRUE,grid.on = TRUE,thresh=1.0, return.data = TRUE)$communication
dfTarget <- netVisual_bubble(cellchat2, targets.use = c("Macrophages_sgNegCtrl","Macrophages_Ptprc"), sources.use = coi, remove.isolate = FALSE, sort.by.target = TRUE,grid.on = TRUE,thresh=1.0, return.data = TRUE)$communication

# Significant interactions only
interactionSource <- names(table(dfSource$interaction_name))[table(dfSource$interaction_name) >= 2]
interactionTarget <- names(table(dfTarget$interaction_name))[table(dfTarget$interaction_name) >= 2]

# Differential interaction table
diffTblSource <- data.frame(matrix(NA, nrow = length(interactionSource), ncol = length(coi)), row.names = interactionSource)
colnames(diffTblSource) <- coi
diffTblTarget <- data.frame(matrix(NA, nrow = length(interactionTarget), ncol = length(coi)), row.names = interactionTarget)
colnames(diffTblTarget) <- coi


# calculate differences between sgPtprc vs sgNegCtrl in interactions - macrophages as source
for (i in 1:length(interactionSource)){
  interaction = row.names(diffTblSource)[i]

  for (j in 1:length(coi)){
    probPerturb = dfSource[dfSource$source == "Macrophages_Ptprc" & dfSource$target == coi[j] & dfSource$interaction_name == interaction,"prob"][1]
    probCtrl = dfSource[dfSource$source == "Macrophages_sgNegCtrl" & dfSource$target == coi[j] & dfSource$interaction_name == interaction,"prob"][1]
    
    # only include significant interactions
    logical = dfSource[dfSource$source == "Macrophages_Ptprc" & dfSource$target == coi[j] & dfSource$interaction_name == interaction,"pval"][1] > 1 | 
      dfSource[dfSource$source == "Macrophages_sgNegCtrl" & dfSource$target == coi[j] & dfSource$interaction_name == interaction,"pval"][1] > 1
    
    if (is.na(logical)) {
      # Skip to the next iteration
      next
    }
    if(logical == TRUE){
      diffTblSource[i,j] <- (probPerturb - probCtrl)/probCtrl
    }
  }
} 

# calculate differences between sgPtprc vs sgNegCtrl in interactions - macrophages as target 
for (i in 1:length(interactionTarget)){
  interaction = row.names(diffTblTarget)[i]
  
  for (j in 1:length(coi)){
    probPerturb = dfTarget[dfTarget$target == "Macrophages_Ptprc" & dfTarget$source == coi[j] & dfTarget$interaction_name == interaction,"prob"][1]
    probCtrl = dfTarget[dfTarget$target == "Macrophages_sgNegCtrl" & dfTarget$source == coi[j] & dfTarget$interaction_name == interaction,"prob"][1]
    
    # only include significant interactions
    logical = dfTarget[dfTarget$target == "Macrophages_Ptprc" & dfTarget$source == coi[j] & dfTarget$interaction_name == interaction,"pval"][1] > 1 | 
      dfTarget[dfSource$target == "Macrophages_sgNegCtrl" & dfTarget$source == coi[j] & dfTarget$interaction_name == interaction,"pval"][1] > 1
    
    if (is.na(logical)) {
      # Skip to the next iteration
      next
    }
    if(logical == TRUE){
      diffTblTarget[i,j] <- (probPerturb - probCtrl)/probCtrl
    }
  }
} 

# Now filter the differntial interaction table by the DEGs from deseq
diffTblSourceDE <- diffTblSource[grep(paste(degs,collapse="|"),row.names(diffTblSource),ignore.case = TRUE),]
diffTblTargetDE <- diffTblTarget[grep(paste(degs,collapse="|"),row.names(diffTblTarget),ignore.case = TRUE),]

diffTblSourceDE <- diffTblSourceDE[,grep("precursor",colnames(diffTblSourceDE),invert=TRUE)]
diffTblTargetDE <- diffTblTargetDE[,grep("precursor",colnames(diffTblTargetDE),invert=TRUE)]

diffTblSourceDE[is.na(diffTblSourceDE)] <- 0
diffTblTargetDE[is.na(diffTblTargetDE)] <- 0

# Only keep high values
diffTblSourceDE <- diffTblSourceDE[apply(diffTblSourceDE,1,function(x) max(abs(x)) > 0.05),]
diffTblTargetDE <- diffTblTargetDE[apply(diffTblTargetDE,1,function(x) max(abs(x)) > 0.05),]

# Annotate
dfSourceConvert <- dfSource[,c("interaction_name","interaction_name_2")]
dfSourceConvert <- dfSourceConvert[complete.cases(dfSourceConvert),]
dfTargetConvert <- dfTarget[,c("interaction_name","interaction_name_2")]
dfTargetConvert <- dfTargetConvert[complete.cases(dfTargetConvert),]

row.names(dfSourceConvert) <- make.unique(as.character(dfSourceConvert$interaction_name))
row.names(dfTargetConvert) <- make.unique(as.character(dfTargetConvert$interaction_name))

row.names(diffTblSourceDE) <- dfSourceConvert[row.names(diffTblSourceDE),"interaction_name_2"]
row.names(diffTblTargetDE) <- dfTargetConvert[row.names(diffTblTargetDE),"interaction_name_2"]

Heatmap(diffTblSourceDE, name = "Differential Interaction",
        col = colorRamp2(c(-0.1, 0, 0.1), c("blue", "gray90", "red")),
        width = ncol(diffTblSourceDE)*unit(6, "mm"),
        height = nrow(diffTblSourceDE)*unit(3.5, "mm"),
)


Heatmap(diffTblTargetDE, name = "Differential Interaction",
        col = colorRamp2(c(-0.1, 0, 0.1), c("blue", "gray90", "red")),
        width = ncol(diffTblTargetDE)*unit(6, "mm"),
        height = nrow(diffTblTargetDE)*unit(3.5, "mm"),
)





