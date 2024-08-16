# Utility functions for DESeq operations

library(DElegate)

find_deseq_differential_genes <- function(data.obj, 
                                          group1, 
                                          group2, 
                                          group_column,
                                          seed = NULL,
                                          replicate_column = NULL) {
    # Takes in a normalized Seurat object with metadata
    # indicating group1 and group2 and a group_column argument
    # indicating the metadata column to be used. Sets a seed, then
    # finds differential expressed genes and outputs a data frame.
    # Note that log fold comparisons will be output as log2(group1/group2)
    if (!is.null(seed)) {
        set.seed(seed)
    }
    
    # Set up futures:
    future::plan(strategy = 'future::multicore', workers = 12)
    
    # Set futures to greater than max capacity; you may need to tweak this
    options(future.globals.maxSize = 8 * 10^9)
    
    df = findDE(object = data.obj, 
                group_column = group_column,
                compare = c(group1, group2),
                method = 'deseq',
                replicate_column = replicate_column)
    return(df)
}