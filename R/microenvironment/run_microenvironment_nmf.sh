#!/bin/bash

# Define an array of cell types
CELL_TYPES=("Astrocytes" "Macrophages" "Microglia" "OPCs" "Oligodendrocytes")

# Iterate over each cell type and run the R script
for CELL_TYPE in "${CELL_TYPES[@]}"
do
   echo "Running NMF for $CELL_TYPE"
   nohup Rscript "/raleighlab/data1/czou/gbm_perturb/gbm_perturb/R/microenvironment/nmf_microenvironment_prkdc_included.R" > "/raleighlab/data1/czou/gbm_perturb/gbm_perturb/R/microenvironment/nmf_microenvironment_prkdc_included_$CELL_TYPE.out" $CELL_TYPE &
done

echo "All jobs have been started."
