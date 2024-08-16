# Utils for GSEA processing and plotting

library(fgsea)
library(dplyr)
library(msigdbr)
library(tidyr)
library(ggplot2)
library(ggtree)
library(ape)
library(patchwork)

get_human_mapping_df <- function() {
  return(annotables::grch38)
}

prepare_deseq2_dfs_for_gsea <- function(data_list, 
                                        mapping_df, 
                                        missing_gene_destination,
                                        save_missing_genes = FALSE) {
  # Prepare a list of DESeq2 dataframes for GSEA by
  # generating an ordered list by log2FC for each, converting
  # to ENSEMBL using the annotables package, and removing genes
  # with low expression. Writes the genes dropped by ENSEMBL conversion 
  # to a destination of choice. Returns a list of ranked lists.
  # 
  # Note that this function doesn't set a seed! Set your seed beforehand
  # for consistent results.
  
  # Convert to ENSEMBL
  print("Converting to ENSEMBL...")
  mapping_df <- mapping_df[!duplicated(mapping_df$symbol), ]
  symbol_to_ensembl <- setNames(mapping_df$ensgene, mapping_df$symbol)
  
  missing_genes <- list()
  data_list_ensembl <- list()
  
  for (name in names(data_list)) {
    deseq_table <- data_list[[name]]
    ensembl_ids <- symbol_to_ensembl[deseq_table$feature]
    missing <- is.na(ensembl_ids)
    missing_genes[[name]] <- deseq_table$feature[missing]
    deseq_table$ensembl <- ensembl_ids
    data_list_ensembl[[name]] <- deseq_table[!missing, ]
  }
  if (save_missing_genes) {
    write.table(missing_genes[[1]], file = missing_gene_destination, 
                row.names = FALSE, quote = FALSE)
  }

  # Remove low expression
  print("Removing low expression...")
  processed_deseq_tables <- list()
  for (name in names(data_list_ensembl)) {
    deseq_table <- data_list_ensembl[[name]]
    deseq_table <- deseq_table[!is.na(deseq_table$pvalue), ]
    deseq_table <- deseq_table[!is.na(deseq_table$padj), ]
    processed_deseq_tables[[name]] <- deseq_table
  }
  
  # Generate ranked order lists
  print("Generating ranked order lists...")
  ranked_gene_set_list <- list()
  for (name in names(processed_deseq_tables)) {
    deseq_table <- processed_deseq_tables[[name]]
    # Group by unique Ensembl ID and average
    deseq_table <- deseq_table %>%
      group_by(ensembl) %>%
      summarize(log_fc = mean(log_fc))
    deseq_table <- deseq_table %>% mutate(rank = rank(log_fc,  ties.method = "random"))
    deseq_table <- deseq_table[order(-deseq_table$rank),]
    
    # Generate a ranked list
    gene_list <- deseq_table$log_fc
    gene_names <- deseq_table$ensembl
    names(gene_list) <- gene_names
    ranked_gene_set_list[[name]] <- gene_list
  }
  return(ranked_gene_set_list)
}

run_fgsea_on_ranked_lists <- function(ranked_lists,
                                      gene_sets,
                                      max_size_of_gene_set = 500,
                                      lowest_allowable_pvalue = 0.0
                                                  ) {
  # Run GSEA using fgsea on ranked ENSEMBL gene ID lists using the listed 
  # parameters and save the result to the result_destination if specified. 
  # Returns a list of GSEA results, one per ranked list.
  # 
  # Note that this function doesn't set a seed! Set your seed beforehand
  # for consistent results.
  fgsea_output_list <- list()
  for (name in names(ranked_lists)) {
    gsea_results <- fgsea(
      pathways = gene_sets,
      stats = ranked_lists[[name]],
      maxSize = max_size_of_gene_set,
      eps = lowest_allowable_pvalue
    )
    fgsea_output_list[[name]] <- gsea_results
  }
  return(fgsea_output_list)
}

make_gsea_bubble_plot_RT <- function(fgsea_output_list,
                                  log_fc_overlay = FALSE,
                                  data_list = NULL,
                                  mapping_df = NULL,
                                  normalization_scheme = "condNormalized") {
  # Given a list of GSEA results, generate a bubble plot for it. This function
  # currently supports plotting clusters within RT conditions for
  # condNormalized outputs.
  
  # Create an ordered dataframe and retrieve the necessary data
  combined_df <- bind_rows(fgsea_output_list, .id = "perturbations")
  # combined_df$pathway <- gsub("HALLMARK_", "", combined_df$pathway)
  # combined_df$pathway <- gsub(" ", "", combined_df$pathway)
  combined_df$log10padj <- -log10(combined_df$padj)

  # Set the normalization scheme
  if (normalization_scheme != "condNormalized" && 
      normalization_scheme != "noRTNormalized") {
    print("Normalization scheme not permitted!")
    return(NULL)
  }
  combined_df$perturbations <- gsub("_non-targeting_noRT", "", 
                                    combined_df$perturbations)
  if (normalization_scheme == "condNormalized") {
    combined_df$perturbations <- gsub("_non-targeting_RT", "",
                                      combined_df$perturbations)
  }
  
  combined_df <- combined_df %>%
    group_by(pathway) %>%
    filter(any(padj < 0.05)) %>%
    ungroup()
  
  # Ready clustering within RT conditions
  nes_df <- combined_df %>%
    dplyr::select(-pval, -padj, -log2err, -ES, -size, -leadingEdge, -log10padj)
  
  nes_df <- nes_df %>%
    group_by(pathway, perturbations) %>%
    pivot_wider(names_from = perturbations, values_from = NES) %>%
    group_by(pathway)
  
  nes_df <- as.data.frame(nes_df)
  rownames(nes_df) <- nes_df$pathway
  nes_df <- nes_df %>% dplyr::select(-pathway)
  
  row_dist <- dist(nes_df, method = "euclidean")
  row_heatmap_clustering <- hclust(row_dist, method = "complete")
  row_order <- order.dendrogram(as.dendrogram(row_heatmap_clustering))
  
  ordering <- c()
  col_heatmap_clusterings <- list()
  for (rt_condition in c("RT", "noRT")) {
    cols <- grep(sprintf("_%s$", rt_condition), colnames(nes_df), value = TRUE)
    col_df <- nes_df[cols]
    col_dist <- dist(t(col_df), method = "euclidean")
    col_hclust <- hclust(col_dist, method = "complete")
    col_heatmap_clusterings[[rt_condition]] <- col_hclust
    col_order <- order.dendrogram(as.dendrogram(col_hclust))
    ordering <- c(ordering, cols[col_order])
  }
  nes_df <- nes_df[, ordering]
  
  df_ordered <- combined_df %>%
    arrange(match(pathway, rownames(nes_df)[row_order]), 
            match(perturbations, colnames(nes_df)))
  df_ordered$pathway <- factor(df_ordered$pathway, 
                               levels = rownames(nes_df)[row_order])
  df_ordered$perturbations <- factor(df_ordered$perturbations,
                                     levels = colnames(nes_df))
  
  # Add the RT bar
  df_ordered$RT <- ifelse(grepl("_RT", df_ordered$perturbations), "RT", "noRT")
  unique_perturbations <- unique(df_ordered$perturbations)
  
  rt_df <- data.frame(
    xmin = unique_perturbations,
    xmax = unique_perturbations,
    ymin = rep(max(as.numeric(df_ordered$pathway)) + 1, length(unique_perturbations)),
    ymax = rep(max(as.numeric(df_ordered$pathway)) + 2, length(unique_perturbations)),
    RT = ifelse(df_ordered$RT[match(unique_perturbations, df_ordered$perturbations)] == "RT", "RT", "noRT")
  )
  rt_df$numeric_xmin <- as.numeric(as.factor(rt_df$xmin)) - 0.5
  rt_df$numeric_xmax <- as.numeric(as.factor(rt_df$xmax)) + 0.5
  
  # Hard coded color mapping for consistency
  term_color_map <- list(
    `RT` = "darkred",
    `noRT` = 'darkblue'
  )
  
  # Add significance layer
  df_ordered$significance <- ifelse(df_ordered$padj < 0.05, "adjp < 0.05", "NA")
  
  # If we want a log fold change overlay, fetch the overlay:
  if (log_fc_overlay && (is.null(mapping_df) || is.null(data_list))) {
    print("Cannot generate LFC overlay! Check mapping_df or data_list")
    return(NULL)
  } else {
    perturb_condition_pathway_means <- generate_log_fc_column(
      data_list,
      mapping_df
    )
    # Map the overlay to the dataframe
    df_ordered$perturb_cond_pathway <- paste(df_ordered$perturbations, df_ordered$pathway, sep = "_")

    # Generate the corresponding mapping for the means
    names(perturb_condition_pathway_means) <- gsub("non-targeting_RT_HALLMARK_|non-targeting_noRT_HALLMARK_|invitro_|CED_|preinf_", 
                                                   "", names(perturb_condition_pathway_means))
    df_ordered$perturb_cond_pathway_means <- unlist(perturb_condition_pathway_means[df_ordered$perturb_cond_pathway])
  }
  
  # Generate a bubble plot
  bubble_plot <- ggplot(df_ordered, aes(x = perturbations))
  
  # Overlay the log fold change if requested
  if (log_fc_overlay) {
    bubble_plot <- bubble_plot + 
      geom_point(aes(y = pathway,
                     size = log10padj, 
                     color = perturb_cond_pathway_means), alpha = 1.0) +
    scale_color_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      na.value = "grey50",
      limits = c(-0.75, 0.75),
      oob = scales::squish
    )
  } else {
    bubble_plot <- bubble_plot + geom_point(aes(y = pathway,
                                                size = log10padj,
                                                color = NES), alpha = 1.0) +
    scale_color_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      na.value = "grey50",
    )
  }
  bubble_plot <- bubble_plot +
    geom_point(data = subset(df_ordered, padj < 0.05),
               aes(y = pathway, size = log10padj, shape = 'adjp < 0.05'),
               color = "black", fill = NA, alpha = 0.8) +
    geom_point(data = subset(df_ordered, is.na(pval) | is.na(padj) | is.na(log10padj)), 
               aes(y = pathway), 
               color = "black", size = 1, alpha = 0.6) +
    scale_size(range = c(2, 7)) +
    scale_shape_manual(name = "Significance",
                       values = c(`adjp < 0.05` = 21),
                       guide = guide_legend(override.aes = list(color = "black", fill = NA))) +
    
    # Add the RT layer
    geom_rect(data = rt_df, aes(xmin = numeric_xmin, xmax = numeric_xmax, 
                                ymin = ymin, ymax = ymax, fill = factor(RT)), inherit.aes = FALSE) +
    scale_fill_manual(values = term_color_map) +

    # Labels
    labs(
      color = "Normalized enrichment score",
      size = "-Log10 padj"
    ) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key.height = unit(0.3, 'cm'),
          axis.line = element_line(colour = "black"),
          legend.key = element_blank())
  
  perturb_dendrogram_RT <- ggtree(as.phylo(col_heatmap_clusterings[["RT"]])) + 
    theme_tree() + coord_flip() + scale_x_reverse()
  perturb_dendrogram_noRT <- ggtree(as.phylo(col_heatmap_clusterings[["noRT"]])) + 
    theme_tree() + coord_flip() + scale_x_reverse()
  module_dendrogram <- ggtree(as.phylo(row_heatmap_clustering)) + theme_tree()
  
  vertical_plot <- (perturb_dendrogram_RT | perturb_dendrogram_noRT) / bubble_plot + plot_layout(heights = c(0.5, 3))
  dendrogram_vertical_plot <- plot_spacer() / module_dendrogram + plot_layout(heights = c(0.85, 1))
  combined_plot <- dendrogram_vertical_plot | vertical_plot
  combined_plot <- combined_plot + plot_layout(widths = c(0.2, 0.8)) + plot_annotation(
    title = "Enrichment of pathways across perturbations",
    subtitle = paste(normalization_scheme, "clustered within RT", sep = ", ")
  )
  return(combined_plot)
}

generate_gsea_plottable_df <- function(fgsea_output_list,
                                       log_fc_overlay = FALSE,
                                       log_fc_mapping_df = NULL,
                                       sample_cluster_labels = NULL) {
  # General function for generating a plottable dataframe.
  # Params:
  #   - fgsea_output_list, a list of outputs from fgsea
  #   - log_fc_overlay (opt.), whether to overlay log fold changes
  #   - log_fc_mapping_df (opt.), a list mapping sample names to avg. LFCs
  # Returns:
  #   - A GSEA bubble plot dataframe whose structure is ready for plotting. 
  #     Its names and specific numbers can be further modified if needed.
  
  combined_df <- bind_rows(fgsea_output_list, .id = "samples")
  combined_df$log10padj <- -log10(combined_df$padj)
  
  combined_df <- combined_df %>%
    group_by(pathway) %>%
    filter(any(padj < 0.05)) %>%
    ungroup()

  # Cluster the columns and rows
  nes_df <- combined_df %>%
    dplyr::select(-pval, -padj, -log2err, -ES, -size, -leadingEdge, -log10padj)
  
  nes_df <- nes_df %>%
    group_by(pathway, samples) %>%
    pivot_wider(names_from = samples, values_from = NES) %>%
    group_by(pathway)
  
  nes_df <- as.data.frame(nes_df)
  rownames(nes_df) <- nes_df$pathway
  nes_df <- nes_df %>% dplyr::select(-pathway)
  
  row_dist <- dist(nes_df, method = "euclidean")
  row_heatmap_clustering <- hclust(row_dist, method = "complete")
  row_order <- order.dendrogram(as.dendrogram(row_heatmap_clustering))
  
  if (is.null(sample_cluster_labels)) {
    col_dist <- dist(t(nes_df), method = "euclidean")
    col_heatmap_clustering <- hclust(col_dist, method = "complete")
    col_order <- order.dendrogram(as.dendrogram(col_heatmap_clustering))
    nes_df <- nes_df[, colnames(nes_df)[col_order]]
  } else {
    ordering <- c()
    col_heatmap_clustering <- list()
    for (sample_cluster_label in sample_cluster_labels) {
      cols <- grep(sprintf("_%s$", sample_cluster_label), 
                   colnames(nes_df), value = TRUE)
      col_df <- nes_df[cols]
      col_dist <- dist(t(col_df), method = "euclidean")
      col_hclust <- hclust(col_dist, method = "complete")
      col_heatmap_clustering[[rt_condition]] <- col_hclust
      col_order <- order.dendrogram(as.dendrogram(col_hclust))
      ordering <- c(ordering, cols[col_order])
    }
    nes_df <- nes_df[, ordering]
  }
  df_ordered <- combined_df %>%
    arrange(match(pathway, rownames(nes_df)[row_order]), 
            match(samples, colnames(nes_df)))
  df_ordered$samples <- factor(df_ordered$samples,
                               levels = colnames(nes_df))
  df_ordered$pathway <- factor(df_ordered$pathway,
                               levels = rownames(nes_df)[row_order])
  df_ordered$significance <- ifelse(df_ordered$padj < 0.05, "adjp < 0.05", "NA")
  
  # If we want a log fold change overlay, fetch the overlay:
  if (log_fc_overlay) {
    if (is.null(log_fc_mapping_df) || is.null(fgsea_output_list)) {
      print("Cannot generate LFC overlay! Check log_fc_mapping_df or fgsea_output_list")
      return(NULL)
    } else {
      sample_pathway_log_fc_means <- generate_log_fc_column(
        fgsea_output_list,
        log_fc_mapping_df
      )
      # Map the overlay to the dataframe
      df_ordered$sample_pathway <- paste(df_ordered$samples, df_ordered$pathway, sep = "_")
      
      # Generate the corresponding mapping for the means
      df_ordered$sample_pathway_log_fc_means <- 
        unlist(sample_pathway_log_fc_means[df_ordered$sample_pathway])
    }
  }
  return(list(
    "ordered_df" = df_ordered,
    "row_heatmap_clustering" = row_heatmap_clustering,
    "column_heatmap_clustering" = col_heatmap_clustering
  ))
}

plot_gsea_plottable_df <- function(df_ordered,
                                   row_heatmap_clustering,
                                   col_heatmap_clustering,
                                   title,
                                   subtitle,
                                   log_fc_overlay = FALSE,
                                   slice_columns = FALSE) {
  # General function for plotting GSEA ordered dataframes.
  
  # Generate a bubble plot
  bubble_plot <- ggplot(df_ordered, aes(x = samples))
  
  # Overlay the log fold change if requested
  if (log_fc_overlay) {
    bubble_plot <- bubble_plot + 
      geom_point(aes(y = pathway,
                     size = log10padj, 
                     color = perturb_cond_pathway_means), alpha = 1.0) +
      scale_color_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = 0,
        na.value = "grey50",
        limits = c(-0.75, 0.75),
        oob = scales::squish
      )
  } else {
    bubble_plot <- bubble_plot + geom_point(aes(y = pathway,
                                                size = log10padj,
                                                color = NES), alpha = 1.0) +
      scale_color_gradient2(
        low = "blue",
        mid = "white",
        high = "red",
        midpoint = 0,
        na.value = "grey50",
      )
  }
  bubble_plot <- bubble_plot +
    geom_point(data = subset(df_ordered, padj < 0.05),
               aes(y = pathway, size = log10padj, shape = 'adjp < 0.05'),
               color = "black", fill = NA, alpha = 0.8) +
    geom_point(data = subset(df_ordered, is.na(pval) | is.na(padj) | is.na(log10padj)), 
               aes(y = pathway), 
               color = "black", size = 1, alpha = 0.6) +
    scale_size(range = c(2, 7)) +
    scale_shape_manual(name = "Significance",
                       values = c(`adjp < 0.05` = 21),
                       guide = guide_legend(override.aes = list(color = "black", fill = NA))) +
    # Labels
    labs(
      color = "Effect size",
      size = "-Log10 padj"
    ) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.key.height = unit(0.3, 'cm'),
          axis.line = element_line(colour = "black"),
          legend.key = element_blank())
  
  
  if (slice_columns) {
    pathway_dendrogram <- NULL
    for (cluster_name in names(col_heatmap_clustering)) {
      ggtree_plot <- ggtree(as.phylo(col_heatmap_clustering[[cluster_name]])) +
        theme_tree() + coord_flip() + scale_x_reverse()
      if (is.null(pathway_dendrogram)) {
        pathway_dendrogram <- ggtree_plot
      } else {
        pathway_dendrogram <- pathway_dendrogram | ggtree_plot
      }
    }
  } else {
    pathway_dendrogram <- ggtree(as.phylo(col_heatmap_clustering)) + 
      theme_tree() + coord_flip() + scale_x_reverse()
  }
  sample_dendrogram <- ggtree(as.phylo(row_heatmap_clustering)) + theme_tree()
  
  vertical_plot <- pathway_dendrogram / bubble_plot + plot_layout(heights = c(0.5, 3))
  dendrogram_vertical_plot <- plot_spacer() / sample_dendrogram + plot_layout(heights = c(0.85, 1))
  combined_plot <- dendrogram_vertical_plot | vertical_plot
  combined_plot <- combined_plot + plot_layout(widths = c(0.2, 0.8)) + plot_annotation(
    title = title,
    subtitle = subtitle
  )
  return(combined_plot)
}

generate_log_fc_column <- function(data_list, mapping_df) {
  # Helper function to generate a log fold change column
  # with the average LFC for all genes in the intersection between
  # a particular pathway's HALLMARK gene set and the union of DE
  # genes across all perturbations and conditions, where a gene
  # is considered DE if it has |log_fc| > 0.1 and p-adj < 0.05.
  
  mapping_df <- mapping_df[!duplicated(mapping_df$symbol), ]
  symbol_to_ensembl <- setNames(mapping_df$ensgene, mapping_df$symbol)
  
  missing_genes <- list()
  data_list_ensembl <- list()
  
  for (name in names(data_list)) {
    deseq_table <- data_list[[name]]
    ensembl_ids <- symbol_to_ensembl[deseq_table$feature]
    missing <- is.na(ensembl_ids)
    missing_genes[[name]] <- deseq_table$feature[missing]
    deseq_table$ensembl <- ensembl_ids
    data_list_ensembl[[name]] <- deseq_table[!missing, ]
  }
  
  # Get the differentially expressed genes for each perturbation
  print("Getting DE genes")
  processed_deseq_tables <- list()
  for (name in names(data_list_ensembl)) {
    deseq_table <- data_list_ensembl[[name]]
    deseq_table <- deseq_table[!is.na(deseq_table$pvalue), ]
    deseq_table <- deseq_table[!is.na(deseq_table$padj), ]
    deseq_table <- deseq_table[deseq_table$padj < 0.05, ]
    deseq_table <- deseq_table[abs(deseq_table$log_fc) > 0.1, ]
    processed_deseq_tables[[name]] <- deseq_table
  }
  data_list <- processed_deseq_tables
  
  # Find the union of all those DE genes
  union_de_genes <- c()
  for (name in names(data_list)) {
    deseq_table <- data_list[[name]]
    ensembl_ids <- deseq_table$ensembl
    union_de_genes <- unique(c(union_de_genes, ensembl_ids))
  }
  
  # Filter each dataframe to include only DE genes in the union
  filter_to_de_genes <- function(deseq_table) {
    deseq_table <- subset(deseq_table, ensembl %in% union_de_genes)
  }
  de_gene_data_list <- lapply(data_list_ensembl, filter_to_de_genes)
  
  # Get the gene sets for each Hallmark pathway
  msigdbr_df <- msigdbr(species = "Homo sapiens", category = "H")
  msigdbr_list = split(x = msigdbr_df$ensembl_gene, f = msigdbr_df$gs_name)
  
  # For each perturbation x condition, build a DESeq dataframe that only
  # includes the intersection of genes between the gene set and DE genes
  print("Building DESeq mapping to intersection...")
  perturb_condition_pathway_dfs <- list()
  for (perturb_condition in names(de_gene_data_list)) {
    deseq_df <- de_gene_data_list[[perturb_condition]]
    for (pathway in names(msigdbr_list)) {
      pathway_intersected_df <- subset(deseq_df, ensembl %in% msigdbr_list[[pathway]])
      perturb_condition_pathway <- paste(perturb_condition, pathway, sep = "_")
      perturb_condition_pathway_dfs[[perturb_condition_pathway]] <- pathway_intersected_df
    }
  }
  
  # Find the average LFC across dataframe in `perturb_condition_pathway_dfs`. 
  # Note that for any gene that doesn't pass expression filters, 
  # we set its contribution to the average to 0.
  calculate_average <- function(deseq_table) {
    deseq_table$log_fc[is.na(deseq_table$pvalue) & is.na(deseq_table$padj)] <- 0
    return(mean(deseq_table$log_fc))
  }
  perturb_condition_pathway_means <- lapply(perturb_condition_pathway_dfs, calculate_average)
  return(perturb_condition_pathway_means)
}
