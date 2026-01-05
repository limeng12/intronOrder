#' Add intron position index to graph list
#' 
#' @param t_igraph_list List of graph objects
#' @param t_intron_pos_mat_fr Data frame of intron positions
#' @return Updated graph list with intron position index
add_intron_pos_index <- function(t_igraph_list, t_intron_pos_mat_fr) {
  
  for(i in 1:length(t_igraph_list)) {
    
    intron_pos_mat_fr_one <- t_intron_pos_mat_fr[t_intron_pos_mat_fr$trans_id == t_igraph_list[[i]]$trans_id, ]
    intron_pos_mat_fr_one <- intron_pos_mat_fr_one[order(intron_pos_mat_fr_one$intron_order, decreasing = FALSE), ]
    
    intron_pos_index_fr_one <- data.frame(
      pos = paste(intron_pos_mat_fr_one[, "chr"], ":", 
                 intron_pos_mat_fr_one[, "start"], "-", 
                 intron_pos_mat_fr_one[, "end"]),
      index = intron_pos_mat_fr_one[, "intron_order"],
      stringsAsFactors = FALSE
    )
    
    t_igraph_list[[i]]$index_pos_fr <- intron_pos_index_fr_one
    
  }
  
  return(t_igraph_list)
}

#' Cluster introns based on splicing patterns
#' 
#' @param t_igraph_list List of graph objects
#' @param t_alpha Smoothing parameter for clustering
#' @return Updated graph list with cluster membership
cluster_introns <- function(t_igraph_list, t_alpha = 0.05) {
  
  for(i in 1:length(t_igraph_list)) {
    members <- cluster_introns_matrix(t_igraph_list[[i]]$adjacency_matrix, t_alpha)
    t_igraph_list[[i]]$members <- members
  }
  
  return(t_igraph_list)
}

#' Cluster introns from adjacency matrix
#' 
#' @param tt_adj_mat Adjacency matrix
#' @param t_alpha_v Smoothing parameter
#' @return Vector of cluster assignments
#' @importFrom dbscan optics extractDBSCAN
cluster_introns_matrix <- function(tt_adj_mat, t_alpha_v = 0.05) {
  
  if(nrow(tt_adj_mat) <= 2) {
    members <- rep(0, nrow(tt_adj_mat))
    names(members) <- colnames(tt_adj_mat)
    
    for(i in 1:length(members)) {
      if(members[i] == 0) {
        members[i] <- "not_in_unit(0)"
      } else {
        members[i] <- paste0("unit_", members[i])
      }
    }
    
    return(members)
  }
  
  t_adj_mat_p <- tt_adj_mat + t_alpha_v
  
  for(i in 1:nrow(t_adj_mat_p)) {
    for(j in 1:ncol(t_adj_mat_p)) {
      if(i == j) {
        t_adj_mat_p[i, j] <- 0.5
        next
      }
      if((t_adj_mat_p[i, j] == 0) && (t_adj_mat_p[j, i] == 0)) {
        next
      }
      
      if(i > j) {
        t_adj_mat_p[i, j] <- (t_adj_mat_p[i, j]) / (t_adj_mat_p[i, j] + t_adj_mat_p[j, i])
        t_adj_mat_p[j, i] <- 1 - t_adj_mat_p[i, j]
      }
    }
  }
  
  a <- optics(t_adj_mat_p, minPts = 2)
  eps_t <- 2 - 0.8
  
  if(!is.na(a$eps_cl)) {
    eps_t <- a$eps_cl
  }
  
  b <- extractDBSCAN(a, eps_cl = eps_t)
  members <- b$cluster
  names(members) <- colnames(t_adj_mat_p)
  
  for(i in 1:length(members)) {
    if(members[i] == 0) {
      members[i] <- "not_in_unit(0)"
    } else {
      members[i] <- paste0("unit_", members[i])
    }
  }
  
  members
}