# Utility functions for file manipulation and other tasks in the gbm_perturb
# directory.
# 
# Author: Christopher Zou

##############################################################################
##############################################################################
# Dependencies:

library(readr)
library(writexl)
library(purrr)
library(dplyr)
library(xlsx)
library(EnhancedVolcano)

##############################################################################
##############################################################################
# Util functions:

# Write all DESeq .csv files from the directories in DIR_PATHS into tabs in an 
# xlsx file specified at OUTPUT_PATH. Pass in COLUMN_NAMES for the output.

write_deseq_directories_to_xlsx = function(dir_paths, output_path, column_names) {
  all_data_list = list()
  
  for (dir_path in dir_paths) {
    csv_files = list.files(path = dir_path, pattern = "*.csv")
    data_list = map(csv_files, ~{
      data = read_delim(file.path(dir_path, .), delim = " ", col_names = column_names)
      
      data = data[-1,]
      data$ave_expr = as.numeric(data$ave_expr)
      data$log_fc = as.numeric(data$log_fc)
      data$stat = as.numeric(data$stat)
      data$pvalue = as.numeric(data$pvalue)
      data$padj = as.numeric(data$padj)
      data$rate1 = NULL
      data$rate2 = NULL
      data$row_number = NULL
      data$groups1_2 = NULL
      
      data[] = lapply(data, function(x) if(is.numeric(x)) signif(x, 4) else x)
      data
    })
    names(data_list) = gsub("\\.csv$", "", csv_files)
    names(data_list) = gsub("non-targeting", "nt", names(data_list))
    all_data_list = c(all_data_list, data_list)
  }
  
  write_xlsx(all_data_list, path = output_path)
}

# Generate an enhanced volcano plot based off a DESeq .csv file specified by
# INPUT_PATH. Writes the volcano plot as a png to OUTPUT_PATH.

generate_enhanced_volcano_plot = function(input_path, output_path) {
  
  
}

input_path = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/deseq/GL261_integrated_20230812_crossContext/preinf_non-targeting_noRT_invitro_non-targeting_noRT.csv"
deseq_table <- read.table(input_path,header=TRUE)
row.names(deseq_table) <- as.character(deseq_table$feature)
df <- deseq_table
df = df[!is.na(df$log_fc), ]
df = df[!is.na(df$pvalue), ]
df = df[!is.na(df$padj), ]
pattern_rows <- grep("Rik$|^Gm\\d+$", df$feature)
df <- df[-pattern_rows, ]

# Label the union of the top 10 most significant and top 10 highest and lowest
# log fold changes.

lowest_adjp_features <- df[order(df$padj), "feature"][1:10]
highest_logfc_features <- df[order(df$log_fc, decreasing = TRUE), "feature"][1:10]
lowest_logfc_features <- df[order(df$log_fc), "feature"][1:10]
goi = unique(c(lowest_adjp_features, lowest_logfc_features, highest_logfc_features))

plot = EnhancedVolcano(df,
                lab = df$feature,
                x = 'log_fc',
                y = 'padj',
                subtitle = "adjp < 0.05, |log_fc| > 0.1",
                xlab = bquote(~Log[2]~ 'fold change'),
                title = "Preinf vs. invitro effects",
                selectLab = goi,
                pCutoff = 0.05,
                FCcutoff = 0.1,
                pointSize = 2.0,
                labSize = 4.0,
                colAlpha = 1,
                col=c('gray50', 'gray50', 'gray50', 'red3'),
                legendPosition = 'right',
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                cutoffLineType = 'blank',
                vline = 0,
                hline = 0.05,
                legendLabSize = 6,
                legendIconSize = 2.0,
                drawConnectors = TRUE,
                max.overlaps = 15,
                widthConnectors = 0.25,
                arrowheads = FALSE)
ggsave(plot = plot, paste('preinf_non-targeting_noRT_invitro_non-targeting_noRT','.pdf',sep=''),width=7,height=7,useDingbats=FALSE)

