

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

#' Build isoform object from intron splicing order files
#' 
#' @param files_all Vector of file paths containing intron splicing order data
#' @param intron_anno Data frame of intron annotations
#' @param trans_exp_file Optional file with expressed transcript IDs
#' @return Data frame of intron splicing order pairs
#' @importFrom dplyr group_by summarise
#' @importFrom stringr str_c
build_iso_object <- function(files_all, intron_anno, trans_exp_file = "") {
  
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
  iso_raw <- read.table(files_all[1], header = FALSE, sep = "\t", as.is = TRUE)
  
  if(ncol(iso_raw) == 7) {
    iso_raw <- iso_raw[, c(1:4, 6, 7)]
  } else {
    iso_raw <- iso_raw[, c(1:4, 6)]
    iso_raw <- cbind(iso_raw, rep(0, nrow(iso_raw)))
  }
  
  colnames(iso_raw) <- c("id", "nexti", "first", "strand", "read_count", "read_count_jc")

  iso_tmp<-"dd"
    
  for(i in files_all[-1]) {
    print(paste0("load file: ", i))
    if(file.size(i) == 0 || (!file.exists(i))) {
      next
    }
    
    iso_tmp <- read.table(i, header = FALSE, sep = "\t", as.is = TRUE)
    if(ncol(iso_tmp) == 7) {
      iso_tmp <- iso_tmp[, c(1:4, 6, 7)]
    } else {
      iso_tmp <- iso_tmp[, c(1:4, 6)]
      iso_tmp <- cbind(iso_tmp, rep(0, nrow(iso_tmp)))
    }
    
    colnames(iso_tmp) <- c("id", "nexti", "first", "strand", "read_count", "read_count_jc")
    iso_raw <- rbind(iso_raw, iso_tmp)
  }
  
    
  # 假设你的数据叫 iso_raw
  iso <- as.data.frame(iso_raw %>% 
                         group_by(id, nexti, first, strand) %>% 
                         summarise(read_count = sum(read_count),
                                   read_count_jc = sum(read_count_jc) , .groups = 'drop'   ) )  ;
  
  #  save(iso,file="xx.Rd")
  
  if(trans_exp_file != "") {
    exp_trans <- readLines(trans_exp_file)
    iso[, "trim_id"] <- sapply(strsplit(iso[, "id"], "\\."), "[", 1)
    iso <- iso[iso$id %in% exp_trans, ]
    print(paste0("number of expressed transcripts: ", length(unique(iso$id))))
  }

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