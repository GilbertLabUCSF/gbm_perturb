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

##############################################################################
##############################################################################
# Util functions:

# Write all .csv files from a directory into tabs in an xlsx file.

write_directory_to_xlsx = function(dir_path, output_path) {
  csv_files <- list.files(path = dir_path, pattern = "*.csv")
  data_list <- map(csv_files, read_csv)
  names(data_list) <- gsub("\\.csv$", "", csv_files)
  write_xlsx(data_list, path = output_path)
}