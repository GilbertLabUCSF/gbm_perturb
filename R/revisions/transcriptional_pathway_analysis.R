# Analyze transcriptional pathways for GL261 and SB28 DESeq based on the two
# sets: P53, ferroptosis, apoptosis, and mitochondrial genes AND RAS, P53,
# Stat3/Stat5, and inflammatory/cytokine pathways

# To do this, we will rerun GSEA and generate bubble plots with a custom
# blend of gene sets

library(dplyr)
library(fgsea)
library(msigdbr)
library(nichenetr)

OUTPUT_DIR <- "output/transcriptional_pathway_analysis"

utils.gsea <- new.env()
source("R/utils/gsea.R", local = utils.gsea)

#=========================================================
# Load DESeq data
#=========================================================
extract_deseq_list <- function(file_location) {
    # Given a file location, extract DESeq2 output from all .csv
    # files in that file location and put them in a list
    deseq_list <- list()
    for (file_name in list.files(file_location)) {
        list_name <- gsub(".csv", "", file_name)
        deseq_list[[list_name]] <- read.table(sprintf("%s/%s", file_location, file_name))
    }
    return(deseq_list)
}
raw_deseq_dfs <- extract_deseq_list("output/deseq/condNormalized_CED")

#=========================================================
# Load the gene sets of interest
#=========================================================
# MsigDB gene sets
msigdbr_df <- msigdbr(species = "Mus musculus", category = "H")
msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
gene_set_list <- msigdbr_list[
    c(
        "HALLMARK_APOPTOSIS",
        "HALLMARK_IL2_STAT5_SIGNALING",
        "HALLMARK_IL6_JAK_STAT3_SIGNALING",
        "HALLMARK_P53_PATHWAY",
        "HALLMARK_KRAS_SIGNALING_UP",
        "HALLMARK_KRAS_SIGNALING_DN"
    )
]

# Add in other gene sets, including human and mouse ones. For the human gene
# sets, we've converted them offline to mouse gene sets ahead of use.
human_gs_dir <- "data/transcriptional_pathway_analysis/gene_sets/human"
for (file_name in list.files(human_gs_dir)) {
    pathway <- fgsea::gmtPathways(file.path(human_gs_dir, file_name))
    pathway_name <- names(pathway)[1]
    human_genes <- pathway[[pathway_name]]
    mouse_genes <- nichenetr::convert_human_to_mouse_symbols(human_genes)
    mouse_genes <- mouse_genes[!is.na(mouse_genes)]
    mouse_genes <- mouse_genes[!duplicated(mouse_genes)]
    mouse_genes <- unname(mouse_genes)
    gene_set_list[[pathway_name]] <- mouse_genes
}

mouse_gs_dir <- "data/transcriptional_pathway_analysis/gene_sets/mouse"
for (file_name in list.files(mouse_gs_dir)) {
    pathway <- fgsea::gmtPathways(file.path(mouse_gs_dir, file_name))
    pathway_name <- names(pathway)[1]
    gene_set_list[[pathway_name]] <- pathway[[pathway_name]]
}

# Add in an extra mitochondrial gene set, which we will simply make using all
# the mitochondrial genes we have on file. This will be a proxy for the
# enrichment of mitochondrial genes that could be associated with cell death.
# all_genes <- raw_deseq_dfs[[names(raw_deseq_dfs)[1]]]$feature
# all_mt_genes <- all_genes[grepl("mt-", all_genes)]
# gene_set_list[["MITOCHONDRIAL"]] <- all_mt_genes

# The ferroptosis gene sets are mixed. We'll just use the driver one, but
# remove the overlaps with the suppressors
ferroptosis_dir <- "data/transcriptional_pathway_analysis/gene_sets/combined_ferroptosis"
drivers <- read.table(file.path(ferroptosis_dir, "ferroptosis_driver.txt"),
                     quote = "",
                     header = TRUE,
                     sep = "\t")
suppressors <- read.table(file.path(ferroptosis_dir, "ferroptosis_suppressor.txt"),
                          quote = "",
                          header = TRUE,
                          sep = "\t")
overlaps <- intersect(drivers$symbol, suppressors$symbol)
drivers <- drivers %>% filter(ensgstable != "_NA")
driver_genes <- nichenetr::convert_human_to_mouse_symbols(drivers$symbol)
driver_genes <- driver_genes[!is.na(driver_genes)]
driver_genes <- driver_genes[!duplicated(driver_genes)]
driver_genes <- unname(driver_genes)
gene_set_list[["FERR_DB_DRIVERS"]] <- driver_genes

for (file_name in list.files(ferroptosis_dir)) {
    pathway_name <- gsub(".txt", "", file_name)
    ferroptosis <- read.table(file.path(ferroptosis_dir, file_name), quote = "", header = TRUE, sep = "\t")
    ferroptosis <- ferroptosis %>% filter(ensgstable != "_NA_")
    genes <- ferroptosis$symbol
    mouse_genes <- nichenetr::convert_human_to_mouse_symbols(genes)
    mouse_genes <- mouse_genes[!is.na(mouse_genes)]
    mouse_genes <- mouse_genes[!duplicated(mouse_genes)]
    mouse_genes <- unname(mouse_genes)
    gene_set_list[[pathway_name]] <- mouse_genes
}

# Subset the gene set list based on which pathways we're interested in
set1 <- c(
    "HALLMARK_APOPTOSIS",
    "HALLMARK_P53_PATHWAY",
    "GOBP_FERROPTOSIS",
    "FERR_DB_DRIVERS",
    "KEGG_APOPTOSIS",
    "REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_CELL_DEATH_GENES",
    "WP_FERROPTOSIS",
    "REACTOME_APOPTOSIS",
    "REACTOME_PROGRAMMED_CELL_DEATH",
    "REACTOME_TRANSCRIPTIONAL_REGULATION_BY_TP53",
    "WP_APOPTOSIS",
    "GOBP_NECROPTOTIC_SIGNALING_PATHWAY"
)
set2 <- c(
    "HALLMARK_IL2_STAT5_SIGNALING",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_P53_PATHWAY",
    "HALLMARK_KRAS_SIGNALING_UP",
    "HALLMARK_KRAS_SIGNALING_DN",
    "PID_IL2_STAT5_PATHWAY",
    "REACTOME_SIGNALLING_TO_RAS",
    "REACTOME_STAT5_ACTIVATION",
    "REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_CELL_DEATH_GENES",
    "WP_IL6_SIGNALING_PATHWAY",
    "WP_RAS_SIGNALING",
    "REACTOME_TRANSCRIPTIONAL_REGULATION_BY_TP53"
)
gene_set_list <- gene_set_list[set1]

#=========================================================
# Run GSEA with these gene sets. Steps:
# - Convert the genes to ENSEMBL format, and remove any that don't convert
# - Get rid of anything that shows low expression
# - Generate ordered gene lists of the remaining genes
# - Run GSEA on symbols
#=========================================================
mapping_df <- read.table("data/transcriptional_pathway_analysis/deseq_table_ensembl.txt", header = TRUE)
rownames(mapping_df) <- mapping_df$feature
process_dfs <- function(df) {
    df <- df[!is.na(df$pvalue), ]
    df <- df[!is.na(df$padj), ]
    df <- df %>% filter(feature %in% rownames(mapping_df))
    df$ensembl <- mapping_df[df$feature, "ensembl"]
    return(df)
}
processed_deseq_dfs <- lapply(raw_deseq_dfs, process_dfs)

# Generate an ordered list from the processed DESeq2 tables.
set.seed(5220)
ranked_gene_set_list <- list()
for (name in names(processed_deseq_dfs)) {
    deseq_table <- processed_deseq_dfs[[name]]
    deseq_table <- deseq_table %>%
        group_by(ensembl) %>%
        mutate(log_fc = mean(log_fc)) %>%
        ungroup()
    deseq_table <- deseq_table %>% mutate(rank = rank(log_fc,  ties.method = "random"))
    deseq_table <- deseq_table[order(-deseq_table$rank),] # Rank deterministically
    gene_list <- deseq_table$log_fc
    gene_names <- deseq_table$feature
    names(gene_list) <- gene_names
    ranked_gene_set_list[[name]] <- gene_list
}

# Run GSEA on all of our gene sets
fgsea_output_list <- list()
for (name in names(ranked_gene_set_list)) {
    gsea_results <- fgsea(pathways = gene_set_list, stats = ranked_gene_set_list[[name]],
                          maxSize = 500, eps = 0.0) # allow arbitrarily low p-values
    fgsea_output_list[[name]] <- gsea_results
}
saveRDS(fgsea_output_list, paste(OUTPUT_DIR, "fgsea_ensembl_ontology_list_nes_cell_death.rds", sep = "/"))

#=========================================================
# Plot a bubble plot with all of the information above. We will split
# the bubble plot into two sets of gene sets. All of them will need to
# contain RT/noRT
#=========================================================
bubble_plot_df <- utils.gsea$generate_gsea_plottable_df(fgsea_output_list)$ordered_df
bubble_plot_df$leadingEdge <- sapply(bubble_plot_df$leadingEdge, function(x) paste(x, collapse = ","))
write.table(bubble_plot_df, file = file.path(OUTPUT_DIR, "enrichment_bubble_plot_df_ras_stat.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")
bubble_plot <- utils.gsea$make_gsea_bubble_plot_RT(
    fgsea_output_list,
    log_fc_overlay = FALSE
)
ggsave(file.path(OUTPUT_DIR, "enrichment_bubble_plot_fgsea_nes_within_RT_clustered_cell_death.pdf"),
       bubble_plot,
       height = 5,
       width = 6.7 + length(fgsea_output_list) / 27 * (11 - 6.7) / 0.8,
       device = "pdf")

#=========================================================
# Quantify the degree of overlap between the RAS and STAT pathways
#=========================================================
ras_stat_df <- read.table(file = file.path(OUTPUT_DIR, "enrichment_bubble_plot_df_ras_stat.txt"), sep = "\t",
                            header = TRUE)
# Calculate a Pearson R between the NES values across the perturbations
ras_signaling <- ras_stat_df %>% filter(pathway == "WP_RAS_SIGNALING")
stat3 <- ras_stat_df %>% filter(pathway == "HALLMARK_IL6_JAK_STAT3_SIGNALING")
stat5 <- ras_stat_df %>% filter(pathway == "HALLMARK_IL2_STAT5_SIGNALING")
cor(ras_signaling$NES, stat3$NES)
cor(ras_signaling$NES, stat5$NES)

#=========================================================
# Do a test across perturbations for mitochondrial gene content
#=========================================================
data_all <- readRDS("/raleighlab/data1/liuj/gbm_perturb/analysis/GL261_integrated_20230619.Rds")
data_all[["percent_mt"]] <- PercentageFeatureSet(data_all, pattern = "^mt-")
vln_plot <- VlnPlot(
    source_subset,
    features = "percent_mt",
    split.by = "sgRNACond",
    ncol = 3,
    group.by = "double"
)
ggsave(
    file.path(OUTPUT_DIR, sprintf("mt_vln_plot.pdf", condition)),
    plot = vln_plot,
    device = "pdf",
    height = 10,
    width = 12
)



