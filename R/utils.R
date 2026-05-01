
#' Calculate intron pair coverage
#' 
#' @param t_mat Adjacency matrix
#' @return Percentage coverage of intron pairs
cal_intron_pair_cov <- function(t_mat) {
  count_of_nonzero_ele <- 0
  
  for(i in 1:nrow(t_mat)) {
    for(j in 1:ncol(t_mat)) {
      if((t_mat[i, j] + t_mat[j, i] > 0) && (i > j)) {
        count_of_nonzero_ele <- count_of_nonzero_ele + 1
      }
    }
  }
  
  n <- nrow(t_mat)
  percent_coverage_pair <- count_of_nonzero_ele / (n * (n - 1) / 2 + 0.0)
  
  percent_coverage_pair
}


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

#' Extract exons and introns from BED file
#' 
#' @param t_bed_file Path to BED file
#' @param m_trim_trans_id_by_dot Whether to trim transcript ID by dot (default: TRUE)
#' @return List containing intron and exon data frames
#' @importFrom stringr str_split
extract_introns_from_bed <- function(t_bed_file, m_trim_trans_id_by_dot = TRUE) {
  
  bed_anno <- read.table(t_bed_file, header = FALSE, as.is = TRUE)
  colnames(bed_anno) <- c("chr", "start", "end", "trans_id", "score", "strand", 
                          "CDS_start", "CDS_end", "", "exon_count", "exon_len", "exon_start")
  
  bed_anno <- bed_anno[!duplicated(bed_anno$trans_id), ]
  
  exon_pos_mat <- matrix(nrow = sum(bed_anno[, "exon_count"]), ncol = 2)
  exon_pos_index <- 1
  
  intron_pos_mat <- matrix(nrow = sum(bed_anno[, "exon_count"] - 1), ncol = 9)
  intron_pos_index <- 1
  
  exon_list <- str_split(bed_anno[, "exon_len"], ",")
  exon_exon_list <- str_split(bed_anno[, "exon_start"], ",")
  
  cat("\n")
  print("Get exons and introns from BED file")
  pb <- txtProgressBar(min = 1, max = nrow(bed_anno), initial = 0, width = 100, style = 3)
  
  for(row_num in 1:nrow(bed_anno)) {
    setTxtProgressBar(pb, row_num)
    
    exon <- as.numeric(exon_list[[row_num]])
    exon_pos <- as.numeric(exon_exon_list[[row_num]])
    
    exon <- exon[-1 * length(exon)]
    exon_pos <- exon_pos[-1 * length(exon_pos)]
    
    exon_count <- length(exon)
    
    if(exon_count <= 2) {
      next
    }
    
    transcript_id <- bed_anno[row_num, "trans_id"]
    transcript_start_site <- bed_anno[row_num, "start"]
    strand <- bed_anno[row_num, "strand"]
    chr <- bed_anno[row_num, "chr"]
    
    # Process introns
    for(i in 1:(length(exon) - 1)) {
      start <- transcript_start_site + (exon_pos[i]) + (exon[i]) + 1
      end <- (exon_pos[i + 1]) + transcript_start_site + 1 - 1
      
      if(strand == "+") {
        id <- paste0(transcript_id, "_", i)
        intron_index <- i
      } else {
        id <- paste0(transcript_id, "_", (length(exon) - i))
        intron_index <- (length(exon) - i)
      }
      
      intron_pos_mat[intron_pos_index, ] <- c(chr, start, end, id, -1, strand, 
                                              transcript_id, intron_index, 
                                              paste0(chr, ":", start, "-", end))
      intron_pos_index <- intron_pos_index + 1
    }
    
    # Process exons
    for(i in 1:length(exon)) {
      end <- transcript_start_site + (exon_pos[i]) + (exon[i]) + 1 - 1
      start <- (exon_pos[i]) + transcript_start_site + 1 - 1 + 1
      
      exon_pos_mat[exon_pos_index, ] <- c(transcript_id, paste0(chr, ":", start, "-", end))
      exon_pos_index <- exon_pos_index + 1
    }
  }
  
  close(pb)
  
  exon_pos_mat_fr <- as.data.frame(exon_pos_mat, stringsAsFactors = FALSE)
  colnames(exon_pos_mat_fr) <- c("trans_id", "exon_pos")
  
  cat("\n")
  print(paste0("Total number of introns in the annotation: ", intron_pos_index))
  
  intron_pos_mat <- intron_pos_mat[1:(intron_pos_index - 1), ]
  intron_pos_mat_fr <- as.data.frame(intron_pos_mat, stringsAsFactors = FALSE)
  
  colnames(intron_pos_mat_fr) <- c("chr", "start", "end", "id", "score", "strand", 
                                   "trans_id", "intron_order", "intron_pos")
  
  intron_pos_mat_fr$intron_order <- as.numeric(intron_pos_mat_fr$intron_order)
  intron_pos_mat_fr$start <- as.numeric(intron_pos_mat_fr$start)
  intron_pos_mat_fr$end <- as.numeric(intron_pos_mat_fr$end)
  
  list(intron = intron_pos_mat_fr, exons = exon_pos_mat_fr)
}
