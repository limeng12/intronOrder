
#' Build isoform object from intron splicing order files
#' 
#' @param files_all Vector of file paths containing intron splicing order data
#' @param intron_anno Data frame of intron annotations
#' @return Data frame of intron splicing order pairs
#' @importFrom dplyr group_by summarise
#' @importFrom stringr str_c
build_iso_object <- function(files_all, intron_anno) {
  
  intron_anno[, "gencode_intron_region"] <- str_c(intron_anno[, "chr"], ":",
                                                intron_anno[, "start"], "-", intron_anno[, "end"])
  
  intron_anno[, "gencode_intron_o"] <- str_c(
    intron_anno[, "trans_id"],
    "_",
    "intron",
    "_",
    intron_anno[, "intron_order"]
  )
  
  intron_o_frame_first <- data.frame(
    gencode_intron_o_first = intron_anno[, "gencode_intron_o"],
    gencode_intron_region = intron_anno[, "gencode_intron_region"],
    trans_id = intron_anno[, "trans_id"],
    intron_order_first = intron_anno[, "intron_order"],
    stringsAsFactors = FALSE
  )
  
  intron_o_frame_next <- data.frame(
    gencode_intron_o_next = intron_anno[, "gencode_intron_o"],
    gencode_intron_region = intron_anno[, "gencode_intron_region"],
    trans_id = intron_anno[, "trans_id"],
    intron_order_next = intron_anno[, "intron_order"],
    stringsAsFactors = FALSE
  )
  
  
  
  print(paste0("load file: ", files_all[1]))
  iso_list <- list()
  
  for(i in seq_along(files_all)) {
    print(paste0("load file: ", files_all[i]))
    
    if(!file.exists(files_all[i]) || file.size(files_all[i]) == 0) {
      warning(paste("File not found or empty:", files_all[i]))
      next
    }
    
    iso_tmp <- read.table(files_all[i], header = FALSE, sep = "\t", 
                          as.is = TRUE, quote = "")
    
    # 统一列处理
    if(ncol(iso_tmp) == 7) {
      iso_tmp <- iso_tmp[, c(1:4, 6, 7)]
    } else {
      iso_tmp <- iso_tmp[, c(1:4, 6)]
      iso_tmp <- cbind(iso_tmp, rep(0, nrow(iso_tmp)))
    }
    
    colnames(iso_tmp) <- c("id", "nexti", "first", "strand", 
                           "read_count", "read_count_jc")
    
    iso_list[[i]] <- iso_tmp
  }
  
  # 合并所有数据
  iso_raw <- do.call(rbind, iso_list)
  
    
  # 假设你的数据叫 iso_raw
  iso <- as.data.frame(iso_raw %>% 
                         group_by(id, nexti, first, strand) %>% 
                         summarise(read_count = sum(read_count),
                                   read_count_jc = sum(read_count_jc) , .groups = 'drop'   ) )  ;
  


  iso <- iso[(iso[, "nexti"] %in% intron_anno[, "gencode_intron_region"]) &
             (iso[, "first"] %in% intron_anno[, "gencode_intron_region"]), ]
  
  iso_final <- inner_join(iso, intron_o_frame_next, 
                         by = c("nexti" = "gencode_intron_region", "id" = "trans_id")) %>%
               inner_join(intron_o_frame_first, 
                         by = c("first" = "gencode_intron_region", "id" = "trans_id"))
  
  iso_final <- iso_final[(!is.na(iso_final$nexti)) & (!is.na(iso_final$first)) & 
                         (!is.na(iso_final$id)) & (!is.na(iso_final$gencode_intron_o_first)) &
                         (!is.na(iso_final$gencode_intron_o_next)), ]
  
  iso_final <- iso_final[iso_final$gencode_intron_o_first != iso_final$gencode_intron_o_next, ]
  
  print(paste0("Total intron splicing order pairs: ", nrow(iso_final)))
  
  iso_final
}

#' Calculate isoform summary statistics
#' 
#' @param t_iso_final Data frame of intron splicing order pairs
#' @param t_anno_intron Data frame of intron annotations
#' @return Data frame with isoform summary statistics
#' @importFrom dplyr group_by summarise n_distinct inner_join
get_iso_summary <- function(t_iso_final, t_anno_intron) {
  
  iso_slow_sumary <- t_anno_intron %>% 
    group_by(trans_id) %>%
    summarise(intron_count = max(intron_order))
  
  iso_slow_sumary <- as.data.frame(iso_slow_sumary)
  
  small_intron <- apply(t_iso_final[c("nexti", "first")], 1, function(x) min(x))
  large_intron <- apply(t_iso_final[c("nexti", "first")], 1, function(x) max(x))
  
  t_iso_final[, "intron_pair"] <- paste(small_intron, large_intron)
  
  iso_edge_count <- as.data.frame(t_iso_final %>% 
                                   group_by(id) %>%
                                   summarise(intron_pair_count = n_distinct(intron_pair)))
  
  iso_slow_sumary <- inner_join(iso_slow_sumary, iso_edge_count, 
                               by = c("trans_id" = "id"))
  
  iso_slow_sumary[, "percent_intron_pair_coverage"] <- 
    (iso_slow_sumary$intron_pair_count) /
    (iso_slow_sumary[, "intron_count"] * (iso_slow_sumary[, "intron_count"] - 1) / 2)
  
  iso_slow_sumary <- iso_slow_sumary[order(iso_slow_sumary[, "percent_intron_pair_coverage"], 
                                          decreasing = TRUE), ]
  
  print(paste0("Number of multi introns transcripts detected =", nrow(iso_slow_sumary)))
  
  iso_slow_sumary
}