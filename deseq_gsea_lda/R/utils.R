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

##############################################################################
##############################################################################
# Util functions:

# Write all .csv files from a directory into tabs in an xlsx file.

write_deseq_directory_to_xlsx = function(dir_path, output_path, column_names) {
  csv_files = list.files(path = dir_path, pattern = "*.csv")
  # Read full data
  data_list = map(csv_files, ~{
    # Read full data
    data = read_delim(file.path(INPUT_DIR, .), delim = " ", col_names = column_names)
    data = data[-1,]
    data
  })
  names(data_list) = gsub("\\.csv$", "", csv_files)
  write_xlsx(data_list, path = output_path)
}

nt_RT = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_clean_outputs/deseq/noRTNormalized/non-targeting_RT_non-targeting_noRT.csv"
nt_RT_df = read.table(nt_RT)
nt_RT_df <- nt_RT_df[!is.na(nt_RT_df$log_fc) & !is.na(nt_RT_df$pvalue) & !is.na(nt_RT_df$padj),] # & abs(nt_RT_df$log_fc) > 0.1 & nt_RT_df$padj < 0.05, ]

prkdc = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_clean_outputs/deseq/noRTNormalized/PRKDC_RT_non-targeting_noRT.csv"
prkdc_df = read.table(prkdc)
prkdc_df = prkdc_df[!is.na(prkdc_df$log_fc) & !is.na(prkdc_df$pvalue) & !is.na(prkdc_df$padj),]

