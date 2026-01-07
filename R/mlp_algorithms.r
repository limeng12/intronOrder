#' Calculate probability of intron order
#' 
#' @param t_adj_mat Adjacency matrix
#' @param t_intron_order Vector of intron orders
#' @param t_alpha_v Smoothing parameter
#' @return Log probability of the order
calp2 <- function(t_adj_mat, t_intron_order, t_alpha_v) {
  
  t_adj_mat <- t_adj_mat + t_alpha_v
  
  p_all <- 0
  
  order_comb <- combn(t_intron_order, 2)
  
  for(i in 1:ncol(order_comb)) {
    
    order_first <- order_comb[1, i]
    order_next <- order_comb[2, i]
    
    if((t_adj_mat[order_first, order_next] == t_alpha_v) && 
       (t_adj_mat[order_next, order_first] == t_alpha_v)) {
      next
    }
    
    p_one <- log(t_adj_mat[order_first, order_next] /
                 (t_adj_mat[order_first, order_next] + t_adj_mat[order_next, order_first]))
    
    p_all <- p_all + p_one
  }
  
  p_all
}

#' Find most likely intron splicing order
#' 
#' @param t_adj_mat Adjacency matrix
#' @param t_alpha_v Smoothing parameter
#' @param is_verbose Whether to print verbose output
#' @return List containing best order and statistics
find_path_global <- function(t_adj_mat, t_alpha_v = 0.1, is_verbose = FALSE) {
  #print("finish yy")
  
  path_list <- NULL
  
  if(nrow(t_adj_mat) < 10) {
    path_list <- find_best_order_full2_c(t_adj_mat, t_alpha_v)
  } else {
    path_list <- lp_kenemy(t_adj_mat, t_alpha_v)
  }
  
  best_order <- path_list$best_order
  best_score <- calp2(t_adj_mat, best_order, t_alpha_v)
  permut_p <- path_list$permut_p
  entropy_one <- path_list$entropy
  
  if(!is.na(entropy_one)) {
    entropy_one <- entropy_one / log2(factorial(length(best_order)))
  }
  
  number_of_maximum_order <- path_list$number_of_maximum_order
  m_ini_intron_order <- colnames(t_adj_mat)[order(as.numeric(colnames(t_adj_mat)))]
  in_order_spliced_v <- calp2(t_adj_mat, m_ini_intron_order, t_alpha_v)
  
  relative_likelihood <- ((in_order_spliced_v - best_score))
  
  if(is_verbose && (!is.null(number_of_maximum_order))) {
    print(paste("Best path number: ", number_of_maximum_order))
  }
  
  if((!is.na(permut_p)) && (permut_p <= 0)) {
    permut_p <- 1 / 1000000
  }
  
  list(
    best_order = best_order,
    number_of_maximum_order = number_of_maximum_order,
    entropy = entropy_one,
    chi_stat = relative_likelihood,
    permut_p = log(permut_p)
  )
}

#' Linear programming approach for Kenemy problem
#' 
#' @param t_adj_mat Adjacency matrix
#' @param t_alpha_v Smoothing parameter
#' @param verbose_hill Whether to print verbose output
#' @return List containing best order
#' @importFrom igraph graph_from_adjacency_matrix topo_sort
#' @importFrom utils combn
lp_kenemy <- function(t_adj_mat, t_alpha_v, verbose_hill = FALSE) {
  #print("finish xx")
  
  t_adj_mat <- t_adj_mat + t_alpha_v
  t_adj_mat_li <- t_adj_mat
  
  for(i in 1:nrow(t_adj_mat_li)) {
    for(j in 1:ncol(t_adj_mat_li)) {
      if(i > j) {
        p <- t_adj_mat_li[i, j] / (t_adj_mat_li[i, j] + t_adj_mat_li[j, i])
        t_adj_mat_li[i, j] <- log(p)
        t_adj_mat_li[j, i] <- log(1 - p)
      }
      
      if(i == j) {
        t_adj_mat_li[i, j] <- 0
      }
    }
  }
  
  n <- nrow(t_adj_mat_li)
  f.obj <- as.vector(t(t_adj_mat_li))
  
  constraints <- matrix(0, nrow = (n * (n - 1) / 2 + choose(n, 3) * 2), ncol = n * n)
  
  t <- 1
  for(i in 1:nrow(t_adj_mat_li)) {
    for(j in 1:ncol(t_adj_mat_li)) {
      if(i > j) {
        constraints[t, (i - 1) * n + j] <- 1
        constraints[t, (j - 1) * n + i] <- 1
        t <- t + 1
      }
    }
  }
  #print("finish dd")
  
  all_combns <- combn(n, 3)
  
  for(i in 1:ncol(all_combns)) {
    one_combn <- all_combns[, i]
    
    constraints[t, (one_combn[2] - 1) * n + one_combn[1]] <- 1
    constraints[t, (one_combn[3] - 1) * n + one_combn[2]] <- 1
    constraints[t, (one_combn[1] - 1) * n + one_combn[3]] <- 1
    t <- t + 1
    
    one_combn <- rev(one_combn)
    constraints[t, (one_combn[2] - 1) * n + one_combn[1]] <- 1
    constraints[t, (one_combn[3] - 1) * n + one_combn[2]] <- 1
    constraints[t, (one_combn[1] - 1) * n + one_combn[3]] <- 1
    t <- t + 1
  }
  
  rm(all_combns)
  
  f.con <- constraints
  f.dir <- rep("==", n * (n - 1) / 2)
  f.dir <- c(f.dir, rep(">=", choose(n, 3) * 2))
  f.rhs <- rep(1, length(f.dir))
  
  #print("Got all coefficients")
  max <- TRUE
  
  #save(f.obj,f.con,f.dir,f.rhs,max,file="xx.Rd")
  #a <- lpsymphony_solve_LP(obj = f.obj, mat = f.con, dir = f.dir, 
  #                        rhs = f.rhs, max = max, types = "B")
  
  
  
  
  
  a <- Rsymphony_solve_LP2(
    obj = f.obj,
    mat = f.con,
    dir = f.dir,
    rhs = f.rhs,
    types = rep("B", length(f.obj)),
    max = TRUE
  )
  
  
  
  #print("finish LP")
  
  g <- graph_from_adjacency_matrix(matrix(a$solution, byrow = TRUE, nrow = n), 
                                   mode = "directed")
  
  best_order <- as.numeric(topo_sort(g))
  permut_p <- NA
  entropy <- NA
  
  list(best_order = best_order, permut_p = permut_p, 
       entropy = entropy, number_of_maximum_order = NA)
}


