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
