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
  print(dir_path)
  csv_files <- list.files(path = dir_path, pattern = "*.csv")
  print(csv_files)
  data_list <- map(csv_files, ~read_delim(file.path("/path/to/your/directory", .), delim = " "))
  names(data_list) <- gsub("\\.csv$", "", csv_files)
  write_xlsx(data_list, path = output_path)
}

INPUT_DIR = "/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/deseq/GL261_integrated_20230705_preinf_noRTNormalized_sorted"
OUTPUT_PATH = "./raleighlab/data1/czou/gbm_perturb/gbm_perturb_gl261_clean_outputs/deseq_xlsx/preinf_noRTNormalized_sorted.xlsx"

write_directory_to_xlsx(INPUT_DIR, OUTPUT_PATH)
