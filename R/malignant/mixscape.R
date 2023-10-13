## Dependencies

library(Seurat)
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(reshape2)
library(sctransform)

PATH_TO_SEURAT_OBJECT = "/raleighlab/data1/liuj/gbm_perturb/analysis/GBM43_1_malignant_only_annotated_20230817.Rds"
OUTPUT_DIR = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_clean_outputs/mixscape"
setwd(OUTPUT_DIR)

data.main = readRDS(PATH_TO_SEURAT_OBJECT)

plot = DimPlot(data.main, group.by = "cond")
ggsave("dimplot_radiation.png")

data = subset(data.main, cond == "noRT")
plot = DimPlot(data, group.by = "sgRNA")
ggsave("dimplot_noRT_sgRNA.png", width = 10, height = 10)

# Run Mixscape

data$NT
data = RunMixscape(data, assay = "SCT",
                   labels = "sgRNA",
                   nt.class.name = "non-targeting")

# Calculate percentage of KO cells for target gene classes

df = prop.table(table(data$mixscape_class.global, data$sgRNA), 2)

df2 <- reshape2::melt(df)
df2$Var2 <- as.character(df2$Var2)
test <- df2[which(df2$Var1 == "KO"),]
test <- test[order(test$value, decreasing = T),]
new.levels <- test$Var2
df2$Var2 <- factor(df2$Var2, levels = new.levels )
df2$Var1 <- factor(df2$Var1, levels = c("non-targeting", "NP", "KO"))
df2$gene <- sapply(as.character(df2$Var2), function(x) strsplit(x, split = "g")[[1]][1])
df3 <- df2[-c(which(df2$gene == "NT")),]

p1 <- ggplot(df3, aes(x = guide_number, y = value*100, fill= Var1)) +
  geom_bar(stat= "identity") +
  theme_classic()+
  scale_fill_manual(values = c("grey49", "grey79","coral1")) + 
  ylab("% of cells") +
  xlab("sgRNA")
ggsave("mixscape_plot.png", p1, height = 10, width = 10)

p1 + theme(axis.text.x = element_text(size = 18, hjust = 1), 
           axis.text.y = element_text(size = 18), 
           axis.title = element_text(size = 16), 
           strip.text = element_text(size=16, face = "bold")) + 
  facet_wrap(vars(gene),ncol = 5, scales = "free") +
  labs(fill = "mixscape class") +theme(legend.title = element_text(size = 14),
                                       legend.text = element_text(size = 12))
