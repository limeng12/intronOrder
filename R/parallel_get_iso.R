#' 并行获取isoform信息（追加模式）
#' @param bed_file BED文件路径
#' @param bam_file BAM文件路径
#' @param output_file 输出文件路径（如：SRR10538401_iso.tsv）
#' @param n_cores 并行核心数，默认为总核心数-1
#' @param chunks 拆分份数，默认为10
#' @export

parallel_get_iso <- function(bed_file, bam_file, output_file,
                             n_cores =8,
                             chunks = 40) {
  
  cat("开始参数检查...\n")
  
  # 1. 检查BED文件
  if (!file.exists(bed_file)) {
    stop(sprintf("BED文件不存在: %s", bed_file))
  }
  
  if (!grepl("\\.bed$", bed_file, ignore.case = TRUE)) {
    warning(sprintf("BED文件可能不是BED格式（缺少.bed扩展名）: %s", bed_file))
  }
  
  # 尝试读取前几行检查格式
  bed_test <- tryCatch({
    read.table(bed_file, sep = "\t", header = FALSE, nrows = 5, stringsAsFactors = FALSE)
  }, error = function(e) {
    stop(sprintf("无法读取BED文件（可能不是制表符分隔）: %s", e$message))
  })
  
  if (ncol(bed_test) < 1) {
    stop(sprintf("BED文件格式错误：至少需要3列，当前有%d列", ncol(bed_test)))
  }
  
  # 2. 检查BAM文件
  if (!file.exists(bam_file)) {
    stop(sprintf("BAM文件不存在: %s", bam_file))
  }
  
  if (!grepl("\\.bam$", bam_file, ignore.case = TRUE)) {
    warning(sprintf("BAM文件可能不是BAM格式（缺少.bam扩展名）: %s", bam_file))
  }
  
  # 检查BAM索引文件
  bai_files <- c(
    paste0(bam_file, ".bai"),
    paste0(bam_file, ".Bai"),
    paste0(sub("\\.bam$", "", bam_file, ignore.case = TRUE), ".bai"),
    paste0(sub("\\.bam$", "", bam_file, ignore.case = TRUE), ".Bai")
  )
  
  has_index <- any(file.exists(bai_files))
  if (!has_index) {
    warning(sprintf("BAM文件没有找到对应的.bai索引文件，处理可能失败或速度较慢: %s", bam_file))
  }
  
  # 3. 检查输出文件路径
  if (!is.character(output_file) || length(output_file) != 1) {
    stop(sprintf("output_file必须是单个字符串，当前为: %s", 
                 paste(class(output_file), collapse = ", ")))
  }
  
  
  
  
  
  
  
  cat(sprintf("并行处理：%s -> %s\n", basename(bam_file), output_file))
  cat(sprintf("使用 %d 个核心，拆分 %d 份\n", n_cores, chunks))
  
  # 1. 读取并拆分BED
  cat("读取BED文件...\n")
  bed_data <- read.table(bed_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  cat(sprintf("BED文件包含 %d 行记录\n", nrow(bed_data)))
  
  # 均匀拆分
  indices <- rep(1:chunks, length.out = nrow(bed_data))
  bed_splits <- split(bed_data, indices)
  cat(sprintf("拆分为 %d 份，每份约 %d-%d 行\n", 
              length(bed_splits), 
              min(sapply(bed_splits, nrow)),
              max(sapply(bed_splits, nrow))))
  
  # 2. 创建临时输出目录
  output_dir <- tempfile("parallel_output_")
  dir.create(output_dir, recursive = TRUE)
  cat(sprintf("临时文件目录：%s\n", output_dir))
  
  # 3. 清空最终输出文件
  write.table(NULL, output_file, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # 4. 并行处理
  cat("启动并行处理...\n")
  cl <- parallel::makeCluster(n_cores)
  
  # 处理单个chunk的函数
  process_chunk <- function(chunk_idx, bed_split, bam_file, output_dir) {
    # 加载所需包
    
    # 创建临时BED文件
    temp_bed <- tempfile(pattern = sprintf("chunk%d_", chunk_idx), fileext = ".bed")
    write.table(bed_split, temp_bed, sep = "\t", col.names = FALSE, 
                row.names = FALSE, quote = FALSE)
    
    # 临时输出文件（用于getIsoFromBam）
    temp_output <- tempfile(pattern = sprintf("temp_result%d_", chunk_idx), fileext = ".tsv")
    
    # 调用原始函数
    getIsoFromBam(temp_bed, bam_file, temp_output)
    
    # 结果文件路径（用于最终合并）
    chunk_output <- NULL
    
    # 如果有结果，复制到临时目录
    if (file.exists(temp_output) && file.size(temp_output) > 0) {
      chunk_output <- file.path(output_dir, sprintf("result_chunk_%04d.tsv", chunk_idx))
      file.copy(temp_output, chunk_output, overwrite = TRUE)
      result_size <- file.size(chunk_output)
    } else {
      result_size <- 0
    }
    
    # 清理临时文件
    unlink(c(temp_bed, temp_output))
    
    # 返回结果信息
    return(list(
      chunk_idx = chunk_idx,
      output_file = chunk_output,
      result_size = result_size,
      n_rows = ifelse(file.exists(chunk_output) && result_size > 0, 
                      length(readLines(chunk_output)), 0)
    ))
  };
  
  
  # 导出函数到集群
  parallel::clusterExport(cl, "getIsoFromBam", envir = globalenv())
  
  # 进度跟踪函数
  #pb <- txtProgressBar(max = min(chunks, length(bed_splits)), style = 3)
  #progress <- 0
  
  # 并行运行
  results <- parallel::parLapply(cl, 1:min(chunks, length(bed_splits)), function(i) {
    library(intronOrder)
    
    result <- process_chunk(i, bed_splits[[i]], bam_file, output_dir)
    
    # 更新进度条（需要在主进程中执行）
    # 这里我们只返回结果，进度条在主进程中更新
    return(result)
  })
  
  # 在主进程中更新进度条
  # for (i in 1:length(results)) {
  #   #progress <- progress + 1
  #   #setTxtProgressBar(pb, progress)
  #   
  #   # 打印每个chunk的信息
  #   if (!is.null(results[[i]]$output_file) && results[[i]]$result_size > 0) {
  #     cat(sprintf("\nChunk %d: 输出 %d 行", 
  #                 results[[i]]$chunk_idx, 
  #                 results[[i]]$n_rows))
  #   }
  # }
  #close(pb)
  
  # 停止集群
  parallel::stopCluster(cl)
  
  # 5. 合并所有临时文件
  cat("\n合并结果文件...\n")
  
  # 提取所有输出文件路径
  temp_files <- sapply(results, function(x) x$output_file)
  temp_files <- temp_files[!sapply(temp_files, is.null)]
  
  if (length(temp_files) == 0) {
    cat("警告：没有生成任何结果文件\n")
    unlink(output_dir, recursive = TRUE)
    return(NULL)
  }
  
  # 按chunk顺序排序
  temp_files <- temp_files[order(as.numeric(
    sapply(temp_files, function(f) {
      gsub(".*chunk_(\\d+).*", "\\1", basename(f))
    })
  ))]
  
  cat(sprintf("需要合并 %d 个文件\n", length(temp_files)))
  
  # 合并文件
  header_written <- FALSE
  total_rows <- 0
  
  for (i in seq_along(temp_files)) {
    temp_file <- temp_files[[i]]
    
    if (file.exists(temp_file) && file.size(temp_file) > 0) {
      # 读取临时文件内容
      lines <- readLines(temp_file)
      
      if (length(lines) > 0) {
        if (!header_written) {
          # 第一个文件：写入完整内容（包含表头）
          writeLines(lines, output_file)
          header_written <- TRUE
          total_rows <- length(lines)
          cat(sprintf("  文件 %d: 写入 %d 行（含表头）\n", i, length(lines)))
        } else {
          # 后续文件：只追加数据行（跳过表头）
          if (length(lines) > 1) {
            write(lines, output_file, append = TRUE)
            total_rows <- total_rows + length(lines) 
            cat(sprintf("  文件 %d: 追加 %d 行数据\n", i, length(lines) ))
          }
        }
      }
    }
  }
  
  # 6. 清理临时目录
  unlink(output_dir, recursive = TRUE)
  cat(sprintf("清理临时目录：%s\n", output_dir))
  
  # 7. 读取并返回最终结果
  if (file.exists(output_file) && file.size(output_file) > 0) {
    final_result <- tryCatch({
      df <- read.table(output_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
      cat(sprintf("完成！输出 %d 行数据\n", nrow(df)))
      df
    }, error = function(e) {
      cat(sprintf("读取输出文件时出错：%s\n", e$message))
      NULL
    })
    
    return(final_result)
  } else {
    cat("警告：输出文件为空\n")
    return(NULL)
  }
}
# 使用示例
if (FALSE) {
  results <- parallel_get_iso(
    "./Arabidopsis_thaliana.TAIR10.47.chr_no_thick.bed",
    "./SRR10538401_1.bam",
    "SRR10538401_iso.tsv"
  )
  
  
  bedfile <-  system.file("extdata", "Schizosaccharomyces_pombe.ASM294v2.43.chr_nothick.bed", package = "intronOrder");
  bamfile <-  system.file("extdata", "SRR6144325_junction_only.bam", package = "intronOrder");
  results <- parallel_get_iso(
    bedfile,
    bamfile,
    "SRR10538401_iso.tsv"
  )
  
  
}
