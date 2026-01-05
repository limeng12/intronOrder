#' Main pipeline for intron splicing order analysis
#'
#' @param bed_file Path to BED annotation file
#' @param iso_files Vector of paths to intron splicing order files
#' @param output_file Path to output file
#' @param gene_trans_map Optional file mapping transcripts to genes
#' @param read_count_threshold Read count threshold (default: 0)
#' @param trans_exp_file Optional file with expressed transcript IDs
#' @param read_cov_threshold Coverage threshold (default: 0.95)
#' @param trim_trans_id_by_dot Whether to trim transcript ID by dot (default: TRUE)
#' @param alpha Smoothing parameter for clustering (default: 0.1)
#' @return List containing analysis results
#' @export
run_iso_analysis <- function(bed_file, iso_files, output_file,
                            gene_trans_map = "",
                            read_count_threshold = 0,
                            trans_exp_file = "",
                            read_cov_threshold = 0.95,
                            trim_trans_id_by_dot = TRUE,
                            alpha = 0.1) {

  # 1. Extract introns from BED file
  cat("\n=== Step 1: Extracting introns from BED file ===\n")
  exon_intron_data <- extract_introns_from_bed(bed_file)
  intron_pos_mat_fr <- exon_intron_data$intron
  exon_pos_mat_fr <- exon_intron_data$exon

  # Filter by expressed transcripts if provided
  if(trans_exp_file != "") {
    express_trans_ids <- readLines(trans_exp_file)
    intron_pos_mat_fr <- intron_pos_mat_fr[intron_pos_mat_fr[, "trans_id"] %in% express_trans_ids, ]
    exon_pos_mat_fr <- exon_pos_mat_fr[exon_pos_mat_fr[, "trans_id"] %in% express_trans_ids, ]
  }

  # 2. Build isoform object
  cat("\n=== Step 2: Building isoform object ===\n")
  iso_final <- build_iso_object(iso_files, intron_pos_mat_fr,
                               trans_exp_file = trans_exp_file)

  # 3. Calculate isoform summary
  cat("\n=== Step 3: Calculating isoform summary ===\n")
  iso_summary <- get_iso_summary(iso_final, intron_pos_mat_fr)

  # 4. Build adjacency matrices
  cat("\n=== Step 4: Building adjacency matrices ===\n")
  igraph_list <- build_adjacency_matrices(iso_final, iso_summary,
                                         read_cov_threshold, read_count_threshold)

  # 5. Cluster introns
  cat("\n=== Step 5: Clustering introns ===\n")
  igraph_list <- cluster_introns(igraph_list, alpha)

  # 6. Add intron position index
  cat("\n=== Step 6: Adding intron position index ===\n")
  igraph_list <- add_intron_pos_index(igraph_list, intron_pos_mat_fr)

  # 7. Calculate most likely order
  cat("\n=== Step 7: Calculating most likely order ===\n")
  igraph_list <- calculate_most_likely_order(igraph_list, alpha)

  # 8. Calculate heterogeneity
  cat("\n=== Step 8: Calculating heterogeneity ===\n")
  igraph_list <- calculate_heterogeneity(igraph_list, alpha)

  # 9. Output results
  cat("\n=== Step 9: Outputting results ===\n")
  igraph_list <- output_mlo_results(igraph_list, gene_trans_map,
                                   intron_pos_mat_fr, output_file,
                                   trim_trans_id_by_dot)

  
  #generate_splicing_order_report(igraph_list, paste0(output_file,".splicing_order_report.html"));
  
  # draw_mlo_plot(igraph_list, paste0(output_file,".plot.pdf") );

  # draw_mol_table_plot(igraph_list, paste0(output_file,".table.pdf") );
  
  

  cat("\n=== Analysis Complete ===\n")
  print(paste("Results saved to:", output_file))
  print(paste("Total transcripts analyzed:", length(igraph_list)))

  return(list(
    intron_data = intron_pos_mat_fr,
    exon_data = exon_pos_mat_fr,
    iso_data = iso_final,
    iso_summary = iso_summary,
    key_re_list = igraph_list
  ))
}

#' Calculate most likely order for all transcripts
#'
#' @param t_igraph_list List of graph objects
#' @param t_alpha Smoothing parameter
#' @return Updated graph list with most likely orders
#' @export
calculate_most_likely_order <- function(t_igraph_list, t_alpha = 0.1) {

  cat("\n")
  print("Calculating most likely order: ")
  print(length(t_igraph_list))
  pb <- txtProgressBar(min = 0, max = length(t_igraph_list), initial = 0, width = 100, style = 3)

  for(i in 1:length(t_igraph_list)) {
    setTxtProgressBar(pb, i)

    adj_mat <- t_igraph_list[[i]]$adjacency_matrix

    best_order_ls <- find_path_global(adj_mat, t_alpha)

    t_igraph_list[[i]]$best_order <- best_order_ls$best_order
    t_igraph_list[[i]]$p_value <- best_order_ls$permut_p
    t_igraph_list[[i]]$chi_stat <- best_order_ls$chi_stat
    t_igraph_list[[i]]$permut_p <- best_order_ls$permut_p
    t_igraph_list[[i]]$entropy <- best_order_ls$entropy
    t_igraph_list[[i]]$number_of_maximum_order <- best_order_ls$number_of_maximum_order

    best_orders <- t_igraph_list[[i]]$best_order

    if(length(best_orders) == 1) {
      t_igraph_list[[i]]$spearman_p_value <- NA
      t_igraph_list[[i]]$spearman_rho <- NA
      t_igraph_list[[i]]$spearman_rho_abs <- NA
    } else {
      t_igraph_list[[i]]$spearman_p_value_greater <- cor.test(best_orders, 1:length(best_orders),
                                                              alternative = "greater",
                                                              method = "spearman")$p.value

      t_igraph_list[[i]]$spearman_p_value_less <- cor.test(best_orders, 1:length(best_orders),
                                                           alternative = "less",
                                                           method = "spearman")$p.value

      t_igraph_list[[i]]$spearman_rho <- cor(best_orders, 1:length(best_orders), method = "spearman")
      t_igraph_list[[i]]$spearman_rho_abs <- abs(cor(best_orders, 1:length(best_orders), method = "spearman"))
    }
  }

  close(pb)
  return(t_igraph_list)
}

#' Calculate heterogeneity of splicing patterns
#'
#' @param t_igraph_list List of graph objects
#' @param t_alpha Smoothing parameter
#' @return Updated graph list with heterogeneity measures
#' @export
calculate_heterogeneity <- function(t_igraph_list, t_alpha) {

  cat("\n")
  print("calculate intron splicing order matrix heterogeneity:")
  pb <- txtProgressBar(min = 1, max = length(t_igraph_list), initial = 0, width = 100, style = 3)

  for(i in 1:length(t_igraph_list)) {
    setTxtProgressBar(pb, i)

    t_adj_mat <- t_igraph_list[[i]]$adjacency_matrix
    t_best_orders <- t_igraph_list[[i]]$best_order
    t_adj_mat_order <- t_adj_mat[t_best_orders, t_best_orders]

    entropy_sum <- 0
    n <- 0
    t_adj_mat_li <- t_adj_mat_order

    for(ii in 1:nrow(t_adj_mat_li)) {
      for(jj in 1:ncol(t_adj_mat_li)) {
        if(ii < jj) {
          p <- (t_adj_mat_li[ii, jj] + t_alpha) / (t_adj_mat_li[ii, jj] + t_adj_mat_li[jj, ii] + 2 * t_alpha)
          entropy_sum <- entropy_sum + (-1 * p * log2(p))
          n <- n + 1
        }
      }
    }

    t_entropy_sum_normalized <- entropy_sum / n
    t_igraph_list[[i]]$entropy_sum_normalized <- t_entropy_sum_normalized
  }

  close(pb)
  return(t_igraph_list)
}
