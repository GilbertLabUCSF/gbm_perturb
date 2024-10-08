# GL261 DE gene counts
# 
# This script contains code for counting the number of differentially expressed 
# genes for GL261 DESeq output stored in a directory called 
# data/malignant/deseq/invitro.
#
# In an attempt to provide the greatest possible degree of transparency and to
# fit with the DESeq2 example provided, this code retains a significant 
# number of data wrangling steps.
# 
# Inputs:
# - DESeq output directory containing a CSV file per perturbation
# 
# Outputs: 
# - A table and plots of DE genes per perturbation separated by radiation 
#   condition

## -----------------------------------------------------------------------------------------------------------------------------------------------
# Import libraries:

library(Seurat)
library(clusterProfiler)
library(org.Mm.eg.db)
library(knitr)
library(msigdbr)
library(fgsea)
library(dplyr)
library(data.table)
library(ggplot2)
library(grid)
library(ComplexHeatmap)
library(colorRamp2)
library(tidyr)
library(patchwork)
library(RColorBrewer)

## -----------------------------------------------------------------------------------------------------------------------------------------------
# Import the differential expression data from a condition of choice:

CONTEXT <- "invitro"
OUTPUT_DIR <- sprintf("output/malignant/de_genes/%s", NORMALIZATION, CONTEXT)
INPUT_DIRS <- c(
  sprintf("output/malignant/deseq/%s", CONTEXT)
)
data.list <- list()
for (dir in INPUT_DIRS) {
  file_names <- list.files(dir)
  for (i in file_names) {
    name <- gsub(".csv", "", i)
    df <- read.table(paste(dir, "/", i, sep = ""))
    data.list[[name]] <- df
  }
}

## -----------------------------------------------------------------------------------------------------------------------------------------------
# Construct a dataframe that has RT and noRT numbers of DE genes in its columns. 
# Cutoffs for DE genes are as follows:
#   - padj < 0.05
#   - abs(LFC) > 0.1
# Along with expression filters (genes must have a valid log_fc, p-value, 
# and p-adjusted value):

RT.df <- data.frame()
noRT.df <- data.frame()
names(data.list) <- gsub("_non-targeting_noRT", "", names(data.list))
for (name in names(data.list)) {
  name.df <- data.list[[name]]
  name.df = name.df[!is.na(name.df$log_fc), ]
  name.df = name.df[!is.na(name.df$pvalue), ]
  name.df = name.df[!is.na(name.df$padj), ]
  name.df <- name.df[name.df$padj < 0.05, ]
  name.df <- name.df[abs(name.df$log_fc) > 0.1, ]
  n.degenes <- dim(name.df)[1]
  
  split.name <- strsplit(name, "_")
  context <- split.name[[1]][1]
  perturb <- split.name[[1]][2]
  condition <- split.name[[1]][3]
  if (condition == "noRT") {
    noRT.df <- rbind(noRT.df, c(perturb, n.degenes))
  } else {
    RT.df <- rbind(RT.df, c(perturb, n.degenes))
  }
}
colnames(RT.df) <- c("perturb", "n.degenes.RT")
colnames(noRT.df) <- c("perturb", "n.degenes.noRT")
all.df <- merge(RT.df, noRT.df, on = "perturb", all.x = TRUE, all.y = TRUE)
all.df[is.na(all.df)] <- 0
all.df$n.degenes.RT <- as.numeric(all.df$n.degenes.RT)
all.df$n.degenes.noRT <- as.numeric(all.df$n.degenes.noRT)
write.table(all.df, file = paste(OUTPUT_DIR, sprintf("%s_de_genes.txt", CONTEXT), sep = "/"), row.names = FALSE, quote = FALSE)


## -----------------------------------------------------------------------------------------------------------------------------------------------
# Create a plot that shows the differential expression across 
# perturbations across conditions:

degenes.plot <- ggplot(all.df, aes(x = perturb)) +
  # Add the RT layer
  geom_point(aes(y = n.degenes.RT, size = ifelse(n.degenes.RT == 0, NA, n.degenes.RT)), color = "red", alpha = 0.5) +
  # Add the noRT layer
  geom_point(aes(y = n.degenes.noRT, size = ifelse(n.degenes.noRT == 0, NA, n.degenes.noRT)), color = "blue", alpha = 0.5) +
  scale_size(range = c(0.5, 3)) +
  # Labels
  labs(
    title = "Number of DE Genes across treatment conditions",
    subtitle = paste(CONTEXT, NORMALIZATION, "LFC > 0.1, adjp < 0.05"),
    color = "Treatment condition",
    size = "Number of DE genes",
    y = "Number of DE Genes"
  ) + 
  geom_text(aes(label = "n = 100", x = Inf, y = 10), hjust = "inward", vjust = -1) +
  geom_hline(yintercept = 100, linetype = "dotted", color = "black") +
  scale_color_identity(name = "Treatment Condition", 
                     breaks = c("red", "blue"), 
                     labels = c("RT", "noRT"), 
                     guide = "legend") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.key = element_blank())
ggsave(paste(OUTPUT_DIR, "de_genes_bubble_plot.pdf", sep = "/"),
       degenes.plot, height = 3.5, width = 10, device = pdf)
