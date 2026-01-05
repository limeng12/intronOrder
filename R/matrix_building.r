#' Build adjacency matrices for transcripts
#' 
#' @param t_iso_final Data frame of intron splicing order pairs
#' @param t_iso_slow_sumary Data frame with isoform summary statistics
#' @param intron_pair_cov_threshold Coverage threshold for intron pairs
#' @param t_read_count_threshold Read count threshold
#' @return List of adjacency matrices for each transcript
build_adjacency_matrices <- function(t_iso_final, t_iso_slow_sumary, 
                                    intron_pair_cov_threshold = 0.95, 
                                    t_read_count_threshold = 0) {
  
  isoform_num_produce <- sum(t_iso_slow_sumary[, "percent_intron_pair_coverage"] > intron_pair_cov_threshold)
  
  print(paste0("Number of multi introns transcripts that contain intron pairs >",
               sprintf("%1.2f%%", intron_pair_cov_threshold * 100), "= ", isoform_num_produce))
  
  igraph_list <- list()
  
  t_iso_slow_sumary <- t_iso_slow_sumary[t_iso_slow_sumary[, "percent_intron_pair_coverage"] > intron_pair_cov_threshold, ]
  t_iso_final <- t_iso_final[t_iso_final[, "id"] %in% t_iso_slow_sumary[, "trans_id"], ]
  
  cat("\n")
  print("Build adjacent matrix: ")
  pb <- txtProgressBar(min = 0, max = nrow(t_iso_slow_sumary), initial = 0, width = 100, style = 3)
  
  for(g in 1:nrow(t_iso_slow_sumary)) {
    setTxtProgressBar(pb, g)
    
    t_trans_id <- t_iso_slow_sumary[g, "trans_id"]
    iso_1 <- unique(t_iso_final[t_iso_final[, "id"] %in% t_iso_slow_sumary[g, "trans_id"], ])
    
    t_total_intron_count <- t_iso_slow_sumary[g, "intron_count"]
    
    adjacency_matrix <- matrix(0, nrow = t_total_intron_count, ncol = t_total_intron_count)
    rownames(adjacency_matrix) <- (1:t_total_intron_count)
    colnames(adjacency_matrix) <- (1:t_total_intron_count)
    
    jc_pair_matrix <- matrix(0, nrow = t_total_intron_count, ncol = t_total_intron_count)
    rownames(jc_pair_matrix) <- (1:t_total_intron_count)
    colnames(jc_pair_matrix) <- (1:t_total_intron_count)
    
    gencode_intron_o_first_number <- iso_1[, "intron_order_first"]
    gencode_intron_o_next_number <- iso_1[, "intron_order_next"]
    
    for(i in 1:nrow(iso_1)) {
      adjacency_matrix[gencode_intron_o_first_number[i], gencode_intron_o_next_number[i]] <- (iso_1[i, "read_count"])
      jc_pair_matrix[gencode_intron_o_first_number[i], gencode_intron_o_next_number[i]] <- (iso_1[i, "read_count_jc"])
    }
    
    jc_pair_matrix[is.na(jc_pair_matrix)] <- 0
    
    # for(ii in 1:nrow(adjacency_matrix)) {
    #   for(jj in 1:ncol(adjacency_matrix)) {
    #     if(adjacency_matrix[ii, jj] + adjacency_matrix[jj, ii] < t_read_count_threshold) {
    #       adjacency_matrix[ii, jj] <- adjacency_matrix[jj, ii] <- 0
    #     }
    #   }
    # }
    sum_mat <- adjacency_matrix + t(adjacency_matrix)
    
    # 找出所有满足 sum < threshold 的位置（返回逻辑矩阵）
    mask <- (sum_mat < t_read_count_threshold)
    
    # 将这些位置在原矩阵中设为 0（同时处理 [i,j] 和 [j,i]）
    adjacency_matrix[mask] <- 0
    
    
    percent_coverage_pair <- cal_intron_pair_cov(adjacency_matrix)
    
    if(percent_coverage_pair < intron_pair_cov_threshold) {
      next
    }
    
    node_size <- rep(10, t_total_intron_count)
    names(node_size) <- 1:t_total_intron_count
    
    node_color <- rep("green", t_total_intron_count)
    names(node_color) <- 1:t_total_intron_count
    
    igraph_list[[paste0(t_trans_id)]] <- list(
      gene_symbol = "",
      trans_id = t_trans_id,
      strand = iso_1[1, "strand"],
      adjacency_matrix = adjacency_matrix,
      jc_pair_matrix = jc_pair_matrix,
      percent_coverage_pair = percent_coverage_pair,
      intron_pair_count = t_iso_slow_sumary[g, "intron_pair_count"],
      node_size = node_size,
      node_color = node_color
    )
  }
  
  close(pb)
  
  igraph_list
}