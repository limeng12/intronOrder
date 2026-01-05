#' Output most likely order results
#'
#' @param t_igraph_list List of graph objects
#' @param t_gene_trans_id_map Optional gene-transcript mapping
#' @param t_intron_pos_mat_fr Intron position data frame
#' @param output_file Path to output file
#' @param trim_trans_id_by_dot Whether to trim transcript ID by dot
#' @return Updated graph list
#' @export
output_mlo_results <- function(t_igraph_list, t_gene_trans_id_map,
                               t_intron_pos_mat_fr, output_file,
                               trim_trans_id_by_dot) {

  if(!is.null(t_gene_trans_id_map) && t_gene_trans_id_map != "") {
    gene_map <- read.table(t_gene_trans_id_map, quote = "$",
                           header = FALSE, as.is = TRUE, sep = "\t")[, 1:3]
    colnames(gene_map) <- c("gene_id", "trans_id", "gene_symbol")
    t_gene_trans_id_map_v <- gene_map[, "gene_symbol"]
    names(t_gene_trans_id_map_v) <- gene_map[, "trans_id"]
  }

  cat("\n")
  print("Summarise and output to file:")
  pb <- txtProgressBar(min = 0, max = length(t_igraph_list), initial = 0, width = 100, style = 3)

  for(g in 1:length(t_igraph_list)) {
    setTxtProgressBar(pb, g)

    if(trim_trans_id_by_dot) {
      t_igraph_list[[g]]$trans_id <- sapply(strsplit(t_igraph_list[[g]]$trans_id, "\\."), "[", 1)
      names(t_igraph_list)[g] <- t_igraph_list[[g]]$trans_id
    }

    if(!is.null(t_gene_trans_id_map) && t_gene_trans_id_map != "") {
      t_igraph_list[[g]]$gene_symbol <- t_gene_trans_id_map_v[t_igraph_list[[g]]$trans_id]
    }
  }

  close(pb)

  unlink(output_file)
  cat(c("gene_symbol$transcript_id$spearman_p_value_one_side_less$best_order$number_of_orders_have_same_prob$percent_coverage_order_pair$normalized_entropy$spearman_rho$spearman_rho_abs$relative_likelihood\n"),
      file = paste0(output_file))

  for(i in 1:length(t_igraph_list)) {
    cat(
      paste0(t_igraph_list[[i]]$gene_symbol, "$",
             t_igraph_list[[i]]$trans_id, "$",
             (t_igraph_list[[i]]$spearman_p_value_less), "$",
             paste0(t_igraph_list[[i]]$best_order, collapse = ","), "$",
             (t_igraph_list[[i]]$number_of_maximum_order), "$",
             (t_igraph_list[[i]]$percent_coverage_pair), "$",
             (t_igraph_list[[i]]$entropy), "$",
             t_igraph_list[[i]]$spearman_rho, "$",
             t_igraph_list[[i]]$spearman_rho_abs, "$",
             t_igraph_list[[i]]$chi_stat,
             "\n"),
      append = TRUE, file = paste0(output_file))
  }

  t_igraph_list
}

#' Draw MLO plot for transcript order visualization
#'
#' @param tt_igraph_list List of graph objects
#' @param output_path Path to output PDF file
#' @param all_alt_se_introns Vector of SE intron positions
#' @param all_alt_a5_introns Vector of A5 intron positions
#' @param all_alt_a3_introns Vector of A3 intron positions
#' @param all_alt_ri_introns Vector of RI intron positions
#' @return NULL (creates PDF file)
#' @export
draw_mlo_plot <- function(tt_igraph_list, output_path,
                          all_alt_se_introns = NULL,
                          all_alt_a5_introns = NULL,
                          all_alt_a3_introns = NULL,
                          all_alt_ri_introns = NULL) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) stop("Please install 'igraph' package")
  if (!requireNamespace("ggraph", quietly = TRUE)) stop("Please install 'ggraph' package")
  if (!requireNamespace("tidygraph", quietly = TRUE)) stop("Please install 'tidygraph' package")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("Please install 'stringr' package")
  
  int_num <- sapply(tt_igraph_list, function(x) {
    length(x$best_order)
  })
  
  tt_igraph_list <- tt_igraph_list[order(int_num, decreasing = TRUE)]
  
  print("plot mlo to file:")
  pb <- txtProgressBar(min = 0, max = length(tt_igraph_list), initial = 0, width = 100, style = 3)
  
  
  pdf(output_path, width = 12, height = 8)
  
  for(k in 1:min(length(tt_igraph_list), .Machine$integer.max)) {
    setTxtProgressBar(pb, k)
    
    orders <- tt_igraph_list[[k]]$best_order
    index_pos_fr <- tt_igraph_list[[k]]$index_pos_fr
    
    #message(paste0(tt_igraph_list[[k]]$gene_symbol, ":", 
    #               tt_igraph_list[[k]]$trans_id, ":", k))
    
    ### sort intro position by order
    intron_pos <- integer(length(orders))
    
    for(i in seq_along(orders)) {
      if(orders[i] %in% index_pos_fr[, "index"]) {
        intron_pos[i] <- index_pos_fr[index_pos_fr[, "index"] == orders[i], "pos"]
      } else {
        intron_pos[i] <- NA_integer_
      }
    }
    
    # 创建节点列表
    NodeList <- data.frame(
      name = as.character(orders),
      label = as.character(orders),
      intron_pos = intron_pos,
      stringsAsFactors = FALSE
    )
    
    # 设置节点颜色和标签
    col <- character(nrow(NodeList))
    for(i in seq_len(nrow(NodeList))) {
      col[i] <- "green"
      
      if(!is.null(all_alt_se_introns) && NodeList[i, "intron_pos"] %in% all_alt_se_introns) {
        col[i] <- "yellow"
        NodeList[i, "label"] <- paste0(NodeList[i, "label"], ":SE")
      }
      if(!is.null(all_alt_a5_introns) && NodeList[i, "intron_pos"] %in% all_alt_a5_introns) {
        col[i] <- "cyan"
        NodeList[i, "label"] <- paste0(NodeList[i, "label"], ":A5")
      }
      if(!is.null(all_alt_a3_introns) && NodeList[i, "intron_pos"] %in% all_alt_a3_introns) {
        col[i] <- "wheat"
        NodeList[i, "label"] <- paste0(NodeList[i, "label"], ":A3")
      }
      if(!is.null(all_alt_ri_introns) && NodeList[i, "intron_pos"] %in% all_alt_ri_introns) {
        col[i] <- "white"
        NodeList[i, "label"] <- paste0(NodeList[i, "label"], ":RI")
      }
    }
    
    NodeList$color <- col
    
    # 创建边列表
    m_adj_matris <- tt_igraph_list[[k]]$adjacency_matrix
    from <- character(0)
    to <- character(0)
    weight <- numeric(0)
    
    for(i in seq_along(orders)) {
      for(j in seq_along(orders)) {
        if(i != j) {
          from <- c(from, as.character(orders[i]))
          to <- c(to, as.character(orders[j]))
          weight <- c(weight, m_adj_matris[orders[i], orders[j]])
        }
      }
    }
    
    EdgeList <- data.frame(from = from, to = to, weight = weight, 
                           stringsAsFactors = FALSE)
    
    # 创建图
    g3 <- igraph::graph_from_data_frame(d = EdgeList, vertices = NodeList, directed = TRUE)
    
    # 使用 suppressWarnings 来抑制特定警告
    p <- suppressWarnings({
      ggraph::ggraph(g3, layout = 'linear') +
        ggraph::geom_edge_arc(aes(label = weight), strength = 0.5,
                              arrow = arrow(length = unit(4, 'mm')),
                              end_cap = circle(5, 'mm'), alpha = 0.3) +
        ggraph::geom_node_point(aes(color = color), size = 8) +
        ggraph::geom_node_text(aes(label = label), size = 5) +
        ggraph::theme_graph(foreground = 'steelblue', fg_text_colour = 'white',
                            border = FALSE, base_size = 5,
                            title_size = 12, base_family = "Helvetica") +
        ggplot2::ggtitle(paste0(tt_igraph_list[[k]]$gene_symbol, "-",
                                tt_igraph_list[[k]]$trans_id)) +
        ggplot2::scale_color_identity()  # 使用 identity scale 来保持颜色
    })
    
    print(p)
  }
  
  dev.off()
  message("MLO plot saved to: ", output_path)
  message("Total plots generated: ", length(tt_igraph_list))
}

#' Draw probability table plot for order comparisons
#'
#' @param t_igraph_list List of graph objects
#' @param output_path Path to output PDF file
#' @return NULL (creates PDF file)
#' @export
draw_mol_table_plot <- function(t_igraph_list, output_path) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install 'ggplot2' package")
  if (!requireNamespace("reshape2", quietly = TRUE)) stop("Please install 'reshape2' package")

  int_num <- sapply(t_igraph_list, function(x) {
    length(x$best_order)
  })

  t_igraph_list <- t_igraph_list[order(int_num, decreasing = TRUE)]

  print("plot table to file:")
  pb <- txtProgressBar(min = 0, max = length(t_igraph_list), initial = 0, width = 100, style = 3)
  
  
  pdf(output_path, width = 12, height = 10)

  for(k in 1:min(length(t_igraph_list), .Machine$integer.max)) {
    setTxtProgressBar(pb, k)
    
    #print(paste0(t_igraph_list[[k]]$gene_symbol, ":",
    #             t_igraph_list[[k]]$trans_id, ":", k))

    adj <- t_igraph_list[[k]]$adjacency_matrix
    orders <- t_igraph_list[[k]]$best_order
    t_alpha_v <- 0.1

    t_adj_mat_li <- adj + t_alpha_v

    for(i in 1:nrow(t_adj_mat_li)) {
      for(j in 1:ncol(t_adj_mat_li)) {
        if(i > j) {
          p <- t_adj_mat_li[i, j] / (t_adj_mat_li[i, j] + t_adj_mat_li[j, i])
          t_adj_mat_li[i, j] <- p
          t_adj_mat_li[j, i] <- 1 - p
        }
        if(i == j) {
          t_adj_mat_li[i, j] <- 0
        }
      }
    }

    t_adj_mat_li_format_tmp <- t_adj_mat_li[(orders), orders]
    t_adj_mat_li <- t_adj_mat_li_format_tmp

    colnames(t_adj_mat_li) <- 1:length(orders)
    rownames(t_adj_mat_li) <- 1:length(orders)

    t_adj_mat_li_melt <- reshape2::melt(t_adj_mat_li)
    t_adj_mat_li_melt[, "fre"] <- as.factor(t_adj_mat_li_melt[, "value"] > 0.5)
    levels(t_adj_mat_li_melt[, "fre"]) <- c("<0.5", ">0.5")

    t_adj_mat_li_melt[, "value"] <- format(t_adj_mat_li_melt[, "value"], digits = 2)

    p <- ggplot2::ggplot(t_adj_mat_li_melt) +
      ggplot2::geom_tile(aes(x = Var2, y = Var1, fill = fre),
                         color = "black", size = 0.5, height = 1, width = 1) +
      ggplot2::geom_text(aes(x = Var2, y = Var1, label = value)) +
      ggplot2::scale_x_continuous(expand = c(0, 0),
                                  breaks = 1:length(orders),
                                  labels = orders,
                                  position = "top") +
      ggplot2::scale_y_continuous(expand = c(0, 0),
                                  trans = "reverse",
                                  breaks = 1:length(orders),
                                  labels = orders) +
      ggplot2::theme(
        plot.title = element_text(size = rel(1.2)),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.position = "right") +
      ggplot2::ggtitle(paste0(t_igraph_list[[k]]$gene_symbol, "-",
                              t_igraph_list[[k]]$trans_id)) +
      ggplot2::scale_fill_brewer(palette = "Set3")

    print(p)
  }

  dev.off()
  message("Table plot saved to: ", output_path)
}



