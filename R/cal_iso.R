# iso_analysis_package.R
# 完整的R包代码，用于内含子剪接顺序分析


# ====================== 转录本数据处理 ======================

#' 从BED文件准备转录本数据
#'
#' @param bed_file BED文件路径
#' @param unique_intron 是否只使用唯一内含子
#' @return 转录本列表
prepare_transcript_data <- function(bed_file, unique_intron = FALSE) {
  cat("正在解析BED文件...\n")
  
  bed_lines <- readLines(bed_file)
  transcripts <- list()
  processed_introns <- new.env(hash = TRUE)
  valid_transcripts <- 0
  
  for (line in bed_lines) {
    # 跳过空行和注释
    if (line == "" || substr(line, 1, 1) == "#") next
    
    fields <- strsplit(line, "\t")[[1]]
    if (length(fields) < 12) next
    
    # 解析BED行
    chr <- fields[1]
    tx_start <- as.integer(fields[2]) + 1  # 转换为1-based
    tx_end <- as.integer(fields[3])
    tx_name <- fields[4]
    
    # 链信息
    strand <- "+"
    if (length(fields) > 5 && fields[6] != "") {
      strand <- fields[6]
      if (!strand %in% c("+", "-")) strand <- "+"
    }
    
    # 解析外显子块
    block_count <- as.integer(fields[10])
    if (block_count < 3) next  # 至少需要3个外显子
    
    block_sizes <- as.integer(strsplit(fields[11], ",")[[1]])
    block_starts <- as.integer(strsplit(fields[12], ",")[[1]])
    
    # 清理NA值
    block_sizes <- block_sizes[!is.na(block_sizes)]
    block_starts <- block_starts[!is.na(block_starts)]
    
    if (length(block_sizes) != block_count || length(block_starts) != block_count) {
      next
    }
    
    # 计算外显子坐标
    exons <- list()
    for (i in 1:block_count) {
      exon_start <- tx_start + block_starts[i]
      exon_end <- exon_start + block_sizes[i] - 1
      exons[[i]] <- list(start = exon_start, end = exon_end)
    }
    
    # 按起始位置排序
    exon_starts <- sapply(exons, function(x) x$start)
    exons <- exons[order(exon_starts)]
    
    # 生成内含子信息
    introns <- list()
    if (length(exons) >= 2) {
      for (i in 1:(length(exons) - 1)) {
        left_exon_end <- exons[[i]]$end
        right_exon_start <- exons[[i + 1]]$start
        
        intron_start <- left_exon_end + 1
        intron_end <- right_exon_start - 1
        
        if (intron_start < intron_end) {
          intron_key <- paste(chr, intron_start, intron_end, sep = ":")
          
          if (unique_intron) {
            if (exists(intron_key, envir = processed_introns)) next
            assign(intron_key, TRUE, envir = processed_introns)
          }
          
          introns[[length(introns) + 1]] <- list(
            start = intron_start,
            end = intron_end,
            left_exon_end = left_exon_end,
            right_exon_start = right_exon_start
          )
        }
      }
    }
    
    # 只保留至少有2个内含子的转录本
    if (length(introns) >= 2) {
      valid_transcripts <- valid_transcripts + 1
      
      transcripts[[tx_name]] <- list(
        name = tx_name,
        chr = chr,
        strand = strand,
        start = tx_start,
        end = tx_end,
        exons = exons,
        introns = introns,
        n_introns = length(introns)
      )
    }
  }
  
  cat(sprintf("解析完成，找到 %d 个有效转录本\n", valid_transcripts))
  return(transcripts)
}

#' 获取转录本的BAM数据
#'
#' @param bam_file BAM文件路径
#' @param transcript 转录本信息
#' @param padding 区域扩展大小
#' @param min_mapq 最小MAPQ值
#' @return BAM数据列表
get_bam_data_for_transcript <- function(bam_file, transcript, 
                                        padding = 1000, min_mapq = 1) {
  chr <- transcript$chr
  tx_start <- transcript$start
  tx_end <- transcript$end
  
  # 扩展查询区域
  query_start <- max(1, tx_start - padding)
  query_end <- tx_end + padding
  
  # 读取BAM数据
  param <- Rsamtools::ScanBamParam(
    what = c("qname", "rname", "pos", "qwidth", "cigar", "mapq"),
    tag = c("NH", "jI"),
    which = GRanges(chr, IRanges(query_start, query_end)),
    mapqFilter = min_mapq
  )
  
  bam <- Rsamtools::scanBam(bam_file, param = param)[[1]]
  
  if (length(bam$qname) == 0) {
    return(NULL)
  }
  
  # 计算read结束位置
  read_ends <- bam$pos + cigarWidthAlongReferenceSpace(bam$cigar) - 1
  
  # 处理标签
  nh_tags <- bam$tag$NH
  if (is.null(nh_tags)) nh_tags <- vector("list", length(bam$qname))
  
  ji_tags <- bam$tag$jI
  if (is.null(ji_tags)) ji_tags <- vector("list", length(bam$qname))
  
  # 返回格式化的数据
  return(list(
    qname = as.character(bam$qname),
    pos = as.integer(bam$pos),
    end = as.integer(read_ends),
    cigar = as.character(bam$cigar),
    mapq = as.integer(bam$mapq),
    nh = nh_tags,
    ji = ji_tags,
    n_reads = length(bam$qname)
  ))
}

# ====================== 主分析函数 ======================

#' 从BAM文件生成iso_file的主函数
#'
#' @param bed_file BED文件路径
#' @param bam_file BAM文件路径
#' @param output_file 输出文件路径
#' @param intron_flank_threshold 内含子侧翼阈值
#' @param consider_exon_in_intron 是否考虑外显子在内含子中
#' @param unique_intron 是否只使用唯一内含子
#' @param min_mapq 最小MAPQ值
#' @param n_threads 线程数
#' @param padding 区域扩展大小
#' @return 分析结果列表
#' @export
getIsoFromBam <- function(bed_file, bam_file, output_file,
                          intron_flank_threshold = 90,
                          consider_exon_in_intron = TRUE,
                          unique_intron = FALSE,
                          min_mapq = 1,
                          n_threads = 1,
                          padding = 1000) {
  
  cat("\n========================================\n")
  cat("从BAM文件分析内含子剪接顺序\n")
  cat("========================================\n")
  
  # 检查文件
  if (!file.exists(bed_file)) stop("BED文件不存在: ", bed_file)
  if (!file.exists(bam_file)) stop("BAM文件不存在: ", bam_file)
  
  # 检查依赖
  #check_dependencies()
  
  # 编译C++代码
  #if (file.exists("iso.cpp")) {
  #  Rcpp::sourceCpp("iso.cpp")
  #  cat("C++代码编译成功\n")
  #} else {
  #  stop("找不到iso.cpp文件")
  #}
  
  # 准备转录本数据
  cat("\n1. 准备转录本数据...\n")
  transcripts <- prepare_transcript_data(bed_file, unique_intron)
  
  if (length(transcripts) == 0) {
    stop("没有找到有效的转录本")
  }
  
  cat(sprintf("   找到 %d 个转录本\n", length(transcripts)))
  
  # 分析转录本
  cat("\n2. 分析转录本...\n")
  results <- list()
  result_cover_left_names<-c();
  
  tx_names <- names(transcripts)
  n_tx <- length(tx_names)
  
  pb <- txtProgressBar(min = 0, max = n_tx, style = 3)
  
  for (i in 1:n_tx) {
    setTxtProgressBar(pb, i)
    
    transcript <- transcripts[[i]]
    bam_data <- get_bam_data_for_transcript(bam_file, transcript, padding, min_mapq)
    
    if (is.null(bam_data) || bam_data$n_reads == 0) next
    
    tryCatch({
      result_list <- analyze_transcript_java_exact(
        transcript_info = transcript,
        bam_data = bam_data,
        intron_flank_threshold = intron_flank_threshold,
        consider_exon_in_intron = consider_exon_in_intron
      )
      
      result<-result_list[[1]];
      result_cover_left_names<-c(result_cover_left_names, result_list[[2]]);
      
      
      if (length(result$transcript) > 0) {
        results[[length(results) + 1]] <- as.data.frame(result)
      }
    }, error = function(e) {
      # 静默处理错误
    })
    
    if (i %% 50 == 0) gc()  # 定期清理内存
  }
  
  close(pb)
  
  # 合并结果
  cat("\n3. 合并结果...\n")
  if (length(results) == 0) {
    cat("   没有找到任何内含子对\n")
    final_df <- data.frame(
      transcript = character(),
      left_intron = character(),
      right_junction = character(),
      strand = character(),
      
      cover_count = integer(),
      junction_count = integer()
    )
  } else {
    
    
    final_df <- do.call(rbind, results);
    
    cat(sprintf("   总共找到 %d 个内含子对\n", nrow(final_df)))
    
    # 排序并保存
    final_df <- final_df[order(final_df$transcript, final_df$left_intron), ]
    
    final_df$xx <- rep("false",nrow(final_df)     )      # 先加到末尾
    col_order <- c(1:4, ncol(final_df), 5:(ncol(final_df)-1))  # 构造新顺序
    final_df <- final_df[col_order]
    
    cat(result_cover_left_names,file="qname.txt",sep = "\n");
    
    data.table::fwrite(final_df, output_file, sep = "\t", quote = FALSE,col.names = FALSE)
    cat(sprintf("   结果已保存到: %s\n", output_file))
  }
  
  cat("\n========================================\n")
  cat("分析完成\n")
  cat("========================================\n")
  

  
  return(list(
    iso_file = output_file,
    results = final_df,
    n_transcripts = n_tx,
    n_intron_pairs = nrow(final_df)
  ))
}

# ====================== 调试函数 ======================

#' 调试单个转录本
#'
#' 从BED文件读取第一个转录本进行详细调试
#'
#' @param bed_file BED文件路径
#' @param bam_file BAM文件路径
#' @param threshold 阈值（默认：90）
#' @param transcript_index 转录本索引（默认：1，即第一行）
#' @return 详细调试结果
#' @export
debug_single_transcript <- function(bed_file, bam_file, 
                                    threshold = 90,
                                    transcript_index = 1) {
  
  cat("\n========================================\n")
  cat("调试单个转录本\n")
  cat("========================================\n")
  
  # 检查文件
  if (!file.exists(bed_file)) stop("BED文件不存在: ", bed_file)
  if (!file.exists(bam_file)) stop("BAM文件不存在: ", bam_file)
  
  # 检查依赖
  check_dependencies()
  
  # 编译C++代码
  if (file.exists("iso.cpp")) {
    Rcpp::sourceCpp("iso.cpp")
  } else {
    stop("找不到iso.cpp文件")
  }
  
  # 读取转录本
  cat("\n1. 读取转录本数据...\n")
  transcripts <- prepare_transcript_data(bed_file, FALSE)
  
  if (length(transcripts) == 0) {
    stop("BED文件中没有找到有效的转录本")
  }
  
  # 检查索引
  if (transcript_index < 1 || transcript_index > length(transcripts)) {
    cat(sprintf("警告: 转录本索引 %d 超出范围，使用第一个转录本\n", transcript_index))
    transcript_index <- 1
  }
  
  # 获取转录本
  tx_name <- names(transcripts)[transcript_index]
  transcript <- transcripts[[tx_name]]
  
  cat(sprintf("   选择转录本: %s (索引: %d/%d)\n", 
              tx_name, transcript_index, length(transcripts)))
  cat(sprintf("   染色体: %s, 链: %s\n", transcript$chr, transcript$strand))
  cat(sprintf("   位置: %d - %d\n", transcript$start, transcript$end))
  cat(sprintf("   外显子数: %d, 内含子数: %d\n", 
              length(transcript$exons), length(transcript$introns)))
  
  # 获取BAM数据
  cat("\n2. 获取BAM数据...\n")
  bam_data <- get_bam_data_for_transcript(bam_file, transcript)
  
  if (is.null(bam_data) || bam_data$n_reads == 0) {
    cat("   没有reads数据，无法进行分析\n")
    return(NULL)
  }
  
  cat(sprintf("   读取到 %d 条reads\n", bam_data$n_reads))
  
  # 显示BAM统计信息
  cat("\n3. BAM数据统计:\n")
  cat(sprintf("   MAPQ >= 30: %d (%.1f%%)\n",
              sum(bam_data$mapq >= 30),
              mean(bam_data$mapq >= 30) * 100))
  cat(sprintf("   MAPQ >= 20: %d (%.1f%%)\n",
              sum(bam_data$mapq >= 20),
              mean(bam_data$mapq >= 20) * 100))
  
  has_splice <- sapply(bam_data$cigar, function(x) grepl("N", x))
  cat(sprintf("   包含剪接的reads: %d (%.1f%%)\n",
              sum(has_splice),
              mean(has_splice) * 100))
  
  # 显示内含子信息
  cat("\n4. 内含子信息:\n")
  for (i in 1:min(3, length(transcript$introns))) {
    intron <- transcript$introns[[i]]
    cat(sprintf("   %d: %s:%d-%d (长度: %d bp)\n",
                i, transcript$chr, intron$start, intron$end,
                intron$end - intron$start + 1))
  }
  
  # 运行分析
  cat("\n5. 运行分析...\n")
  result_list <- analyze_transcript_java_exact_debug(
    transcript_info = transcript,
    bam_data = bam_data,
    intron_flank_threshold = threshold,
    consider_exon_in_intron = TRUE
  )
  result<-result_list[[1]];
  
  # 显示结果
  cat("\n6. 分析结果:\n")
  if (length(result$transcript) > 0) {
    cat(sprintf("   找到 %d 个内含子对:\n", length(result$transcript)))
    
    for (i in 1:min(5, length(result$transcript))) {
      cat(sprintf("   %d: %s | %s -> %s | 覆盖=%d, 剪接=%d\n",
                  i, result$transcript[i], result$left_intron[i],
                  result$right_junction[i], result$cover_count[i],
                  result$junction_count[i]))
    }
    
    if (length(result$transcript) > 5) {
      cat(sprintf("   ... 还有 %d 个结果未显示\n", length(result$transcript) - 5))
    }
    
    cat(sprintf("\n   统计信息:\n"))
    cat(sprintf("     总内含子对: %d\n", length(result$transcript)))
    cat(sprintf("     平均覆盖数: %.1f\n", mean(result$cover_count)))
    cat(sprintf("     最大覆盖数: %d\n", max(result$cover_count)))
    
  } else {
    cat("   没有找到内含子对\n")
    
    # 提供一些调试建议
    cat("\n   调试建议:\n")
    cat("     1. 检查BAM文件是否包含该转录本区域的reads\n")
    cat("     2. 检查reads是否覆盖内含子区域\n")
    cat("     3. 尝试降低阈值或调整参数\n")
  }
  
  cat("\n========================================\n")
  cat("调试完成\n")
  cat("========================================\n")
  
  # 返回结果
  return(list(
    transcript_info = transcript,
    bam_summary = list(
      total_reads = bam_data$n_reads,
      high_quality_reads = sum(bam_data$mapq >= 20),
      splice_reads = sum(has_splice)
    ),
    analysis_result = if (length(result$transcript) > 0) as.data.frame(result) else NULL
  ))
}

# ====================== 使用示例 ======================

#' 检查文件格式
#'
#' @param file_path 文件路径
#' @param file_type 文件类型（"bed"或"bam"）
check_file_format <- function(file_path, file_type = "bed") {
  if (!file.exists(file_path)) {
    return(FALSE)
  }
  
  if (file_type == "bed") {
    # 简单的BED格式检查
    lines <- readLines(file_path, n = 5)
    if (length(lines) == 0) return(FALSE)
    
    # 检查是否有至少一行包含足够字段
    for (line in lines) {
      if (line == "" || substr(line, 1, 1) == "#") next
      fields <- strsplit(line, "\t")[[1]]
      if (length(fields) >= 12) return(TRUE)
    }
    return(FALSE)
    
  } else if (file_type == "bam") {
    # 检查BAM文件是否存在并可能创建索引
    if (!file.exists(paste0(file_path, ".bai"))) {
      cat("BAM索引文件不存在，可能需要创建...\n")
      cat("使用命令: indexBam(\"", file_path, "\")\n", sep = "")
    }
    return(TRUE)
  }
  
  return(FALSE)
}

# ====================== 包初始化 ======================
