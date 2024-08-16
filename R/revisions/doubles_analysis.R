# This analysis runs preprocessing on GL261 cells without filtering out
# the cells that had multiple guides. Note, however, that we do examine those
# cells' RNA counts to make sure that we aren't ending up with doublets.
# 
# The outputs from this analysis will be used for:
# 1. Gene interaction analysis
# 2. Synergy analysis
# 
# Here's how we'll do this:
# - Load in all of the data
# - Subset to the cells contained in the exported object and
#   in the doublet barcodes
# - Run DESeq2

library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(scDblFinder)

utils.deseq <- new.env()
source("R/utils/deseq.R", local = utils.deseq)

OUTPUT_DIR <- "output/doublet_analysis"

#=========================================================
# Load in data
#=========================================================
# Load in the data object
data_all <- readRDS("data/GL261_integrated_20230619.Rds")

# Load in the double cell info
double_cell_info <- read.table(
    "data/pdx_perturb_GL261_sgRNATblDoublePass.txt",
    header = TRUE
)

# Check to make sure that all the cell barcodes are in the data object
all_in <- all(double_cell_info$full_barcode %in% colnames(data_all))
double_barcodes <- double_cell_info$full_barcode

# Give these cells their double status assignments
data_all$double <- ifelse(colnames(data_all) %in% double_barcodes, "double", "single")

#=========================================================
# Do comparative QC on the "doubles" compared to the other cells. Here, we'll
# pool the cells from different conditions.
#=========================================================
QC_OUTPUT_DIR <- file.path(OUTPUT_DIR, "qc")

data_all[["percent_mt"]] <- PercentageFeatureSet(data_all, pattern = "^mt-")

for (condition in unique(data_all$source)) {
    source_subset <- subset(data_all, source == condition)
    # Visualize percent mt, nFeature RNA, nCount RNA
    vln_plot <- VlnPlot(
        source_subset,
        features = c("nFeature_RNA", "nCount_RNA", "percent_mt"),
        ncol = 3,
        group.by = "double"
    )
    ggsave(
        file.path(QC_OUTPUT_DIR, sprintf("%s_qc_vln_plot.pdf", condition)),
        plot = vln_plot,
        device = "pdf",
        height = 10,
        width = 12
    )
    # Visualize a feature scatter against nFeature/nCount
    Idents(source_subset) <- "double"
    count_feature_scatter <- FeatureScatter(
        source_subset,
        feature1 = "nCount_RNA",
        feature2 = "nFeature_RNA",
    )
    ggsave(
        file.path(QC_OUTPUT_DIR, sprintf("%s_count_feature_scatter.pdf", condition)),
        plot = count_feature_scatter,
        device = "pdf",
        height = 10,
        width = 12
    )
}

genes_to_compare <- unique(c(double_cell_info$GeneA, double_cell_info$GeneB, "non-targeting"))
relevant_cells <- subset(data_all, sgRNA %in% genes_to_compare | double == "double")

# Assign the double cells their genes as well
sgRNA_data <- double_cell_info %>% 
    dplyr::select(full_barcode, Gene_AB)
rownames(sgRNA_data) <- sgRNA_data$full_barcode
sgRNA_data <- sgRNA_data[-1]
colnames(sgRNA_data) <- c("sgRNA")
relevant_cells <- AddMetaData(relevant_cells, sgRNA_data, col.name = "sgRNA")

#=========================================================
# Export for later use
#=========================================================
# Download the count matrix and metadata
count_matrix <- GetAssayData(relevant_cells, layer = "counts")
metadata <- relevant_cells@meta.data

writeMM(count_matrix, file = file.path(OUTPUT_DIR, "data", "count_matrix.mtx"))
write.csv(metadata, file = file.path(OUTPUT_DIR, "data", "metadata.csv"), row.names = TRUE)

# Download a simplified Seurat object without dimensionality reductions
simplified_obj <- CreateSeuratObject(counts = GetAssayData(relevant_cells, assay = "RNA", slot = "counts"))
simplified_obj@meta.data <- relevant_cells@meta.data
saveRDS(simplified_obj, file.path(OUTPUT_DIR, "data", "simplified_seurat.rds"))

#=========================================================
# Run DESeq against the non-targeting cells for each condition
#=========================================================
DESEQ_OUTPUT_DIR <- file.path(OUTPUT_DIR, "deseq")

# Set up a DESeq loop for preinf_RT and preinf_noRT
preinf_cells <- subset(relevant_cells, source == "preinf")
for (condition in unique(preinf_cells$cond)) {
    subset_data <- subset(relevant_cells, cond == condition)
    perturbs <- unique(subset_data$sgRNA)
    perturbs <- perturbs[perturbs != "non-targeting"]
    for (sgRNA in perturbs) {
        deseq_df <- utils.deseq$find_deseq_differential_genes(
            subset_data,
            sgRNA,
            "non-targeting",
            "sgRNA",
            seed = 5220
        )
        file_name <- sprintf("preinf_%s_%s_non-targeting.csv", condition, sgRNA)
        write.table(deseq_df, file = file.path(
            DESEQ_OUTPUT_DIR, file_name
        ))
    }
}

