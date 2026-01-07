#' Generate Static HTML Report for Intron Splicing Order Analysis
#' 
#' @description
#' This function generates a comprehensive static HTML report from intron splicing
#' order analysis results. The report includes interactive visualizations of
#' Most Likely Order (MLO) networks, adjacency matrices, frequency matrices,
#' and detailed statistical summaries.
#' 
#' @details
#' The function processes pre-computed intron splicing order data and creates
#' a self-contained HTML file with the following features:
#' \itemize{
#'   \item Interactive MLO network visualization (ggraph-style with directional arcs)
#'   \item Read count and frequency matrices with sorting options
#'   \item Detailed statistical analysis of splicing order
#'   \item Gene and transcript browsing with search functionality
#'   \item All visualizations are interactive and require no backend server
#' }
#' 
#' The MLO network visualization mimics the appearance of ggplot2's ggraph plots,
#' with linear layouts, directional arcs, and interactive controls for adjusting
#' visualization parameters.
#'
#' @param t_igraph_list A named list containing intron splicing order analysis
#'   results. Each element should correspond to a transcript and contain:
#'   \itemize{
#'     \item \code{adjacency_matrix}: Adjacency matrix of intron splicing counts
#'     \item \code{best_order}: Vector of most likely splicing order
#'     \item \code{gene_symbol}: Gene symbol (optional)
#'     \item \code{trans_id}: Transcript ID
#'     \item \code{strand}: Strand information (optional)
#'     \item \code{statistics}: List of statistical measures including:
#'       \itemize{
#'         \item \code{spearman_rho}, \code{spearman_p_value_less}, \code{spearman_p_value_greater}
#'         \item \code{p_value}, \code{chi_stat}, \code{entropy}
#'         \item \code{percent_coverage_pair}, \code{number_of_maximum_order}
#'       }
#'   }
#' 
#' @param output_file Character string specifying the output file path.
#'   Defaults to "splicing_order_report.html" in the current working directory.
#' 
#' @param ... Additional arguments (currently not used, reserved for future extensions).
#'
#' @return Invisibly returns the processed data used to generate the report.
#'   The function primarily creates an HTML file as a side effect.
#'
#' @examples
#' \dontrun{
#' # Load your intron splicing order data
#' load("path/to/your/t_igraph_list_parent.Rd")
#' 
#' # Generate HTML report
#' generate_splicing_order_report(t_igraph_list_parent$key_re_list)
#' 
#' # Or specify custom output file
#' generate_splicing_order_report(
#'   t_igraph_list = t_igraph_list_parent$key_re_list,
#'   output_file = "my_splicing_report.html"
#' )
#' }
#'
#' @references
#' For more information about intron splicing order analysis:
#' \itemize{
#'   \item Methodological details: [Reference paper or documentation]
#'   \item Visualization techniques: Wickham H. (2016) ggplot2: Elegant Graphics for Data Analysis
#' }
#'
#' @seealso
#' Related functions in the intronOrder package:
#' \itemize{
#'   \item \code{\link{calculate_splicing_order}}: Calculate most likely splicing order
#'   \item \code{\link{plot_mlo_network}}: Interactive MLO network plotting
#'   \item \code{\link{extract_splicing_statistics}}: Extract splicing statistics
#' }
#'
#' @export
generate_splicing_order_report <- function(t_igraph_list, output_file) {
  cat("Generating splicing order HTML report...\n")
  
  # 步骤1: 提取数据
  cat("Step 1: Extracting data...\n")
  all_data <- extract_splicing_data(t_igraph_list)
  
  if (length(all_data) == 0) {
    cat("❌ No valid data found. Check your input.\n")
    return(NULL)
  }
  
  # 步骤2: 生成HTML文件
  cat("Step 2: Generating HTML file...\n")
  generate_html_page(all_data, output_file);
  
  cat("\n✅ Report generation complete!\n")
  cat("📁 Output file: splicing_order_report.html\n")
  
  # 显示统计信息
  trans_count <- length(all_data)
  genes <- unique(sapply(all_data, function(x) x$gene_symbol))
  
  cat("📈 Total transcripts:", trans_count, "\n")
  cat("🧬 Unique genes:", length(genes), "\n")
  
  # 显示数据分布
  intron_counts <- sapply(all_data, function(x) x$n_introns)
  cat("\n📊 Intron distribution:\n")
  cat("  Min:", min(intron_counts), "\n")
  cat("  Max:", max(intron_counts), "\n")
  cat("  Mean:", round(mean(intron_counts), 2), "\n")
  
  # 显示示例数据
  if(trans_count > 0) {
    first_trans <- names(all_data)[1]
    first_data <- all_data[[first_trans]]
    cat("\n📋 Example transcript:\n")
    cat("  Gene:", first_data$gene_symbol, "\n")
    cat("  Transcript:", first_trans, "\n")
    cat("  Best order:", paste(first_data$best_order, collapse = " → "), "\n")
    cat("  Introns:", first_data$n_introns, "\n")
    cat("  Spearman ρ:", format(first_data$statistics$spearman_cor, digits = 3), "\n")
  }
  
  cat("\n🎉 Report is ready! Open the HTML file in any web browser.\n")
  cat("Enhancements:\n")
  cat("  • Fixed R syntax error (replaced JavaScript ternary operator with R ifelse)\n")
  cat("  • MLO Net lines are thicker with Link Width control\n")
  cat("  • Added directional arrows: Blue for left-to-right, Orange for right-to-left\n")
  cat("  • Left panel transcripts sorted by gene name, then transcript ID\n")
  cat("  • Each transcript displayed with gene name, transcript ID, and statistics\n")
  cat("  • Interactive controls for Node Radius, Arc Height, Link Width, Font Size\n")
  cat("  • Color-coded edges with directional arrows\n")
  cat("  • Improved hover effects with thicker highlighted lines\n")
  cat("  • Legend showing arrow color meanings\n")
  
  return(all_data)
}




# 主函数：提取并处理数据
extract_splicing_data <- function(t_igraph_list) {
  all_data <- list()
  
  cat("Processing", length(t_igraph_list), "transcripts...\n")
  
  # 为每个转录本预计算可视化数据
  for(trans_id in names(t_igraph_list)) {
    data <- t_igraph_list[[trans_id]]
    
    if(is.null(data$adjacency_matrix)) {
      cat("Skipping", trans_id, "(no adjacency matrix)\n")
      next
    }
    
    # 获取最佳顺序
    orders <- data$best_order
    n <- length(orders)
    
    if(n == 0) {
      cat("Skipping", trans_id, "(no best order)\n")
      next
    }
    
    # 计算统计信息
    spearman_cor <- data$spearman_rho
    spearman_p_less <- data$spearman_p_value_less
    spearman_p_greater <- data$spearman_p_value_greater
    
    # 获取原始矩阵
    orig_adj_matrix <- data$adjacency_matrix
    orig_n <- nrow(orig_adj_matrix)
    
    # 按最佳顺序重新排序矩阵（类似app.R中的order by MLO）
    best_order_indices <- match(orders, rownames(orig_adj_matrix))
    adj_matrix_best <- orig_adj_matrix[best_order_indices, best_order_indices]
    
    adj_matrix_freq <- orig_adj_matrix
    
    # 准备频率矩阵数据（最佳顺序）
    freq_matrix_best <- adj_matrix_best
    if(n > 1) {
      for(i in 1:n) {
        for(j in 1:n) {
          if(i > j) {
            total <- adj_matrix_best[i, j] + adj_matrix_best[j, i]
            if(total > 0) {
              freq_matrix_best[i, j] <- adj_matrix_best[i, j] / total
              freq_matrix_best[j, i] <- adj_matrix_best[j, i] / total
            }
          }
        }
      }
      diag(freq_matrix_best) <- 0
    }
    
    # 准备频率矩阵数据（频率顺序）
    freq_matrix_freq <- adj_matrix_freq
    if(n > 1) {
      for(i in 1:n) {
        for(j in 1:n) {
          if(i > j) {
            total <- adj_matrix_freq[i, j] + adj_matrix_freq[j, i]
            if(total > 0) {
              freq_matrix_freq[i, j] <- adj_matrix_freq[i, j] / total
              freq_matrix_freq[j, i] <- adj_matrix_freq[j, i] / total
            }
          }
        }
      }
      diag(freq_matrix_freq) <- 0
    }
    
    # 为ggraph风格的MLO net准备数据（类似app.R中的线性图）
    # 创建节点数据
    mlo_nodes <- list()
    for(i in 1:n) {
      # 线性布局：节点在x轴上均匀分布
      x_pos <- (i - 1) * 100  # 类似app.R中的线性布局
      y_pos <- 0
      
      mlo_nodes[[i]] <- list(
        id = i - 1,
        label = as.character(orders[i]),
        x = x_pos,
        y = y_pos,
        radius = 15,
        color = "#4CAF50",
        order_idx = i - 1,
        original_label = as.character(orders[i])
      )
    }
    
    # 创建边数据（用于ggraph风格的弧线）
    mlo_edges <- list()
    edge_idx <- 1
    
    # 生成所有可能的边（包括双向）
    for(i in 1:n) {
      for(j in 1:n) {
        if(i != j && adj_matrix_best[i, j] > 0) {
          # 计算方向：从左到右还是从右到左
          is_forward <- i < j  # i在j左边就是从左到右
          direction <- ifelse(is_forward, "forward", "backward")  # 修复这里：使用ifelse而不是三元运算符
          
          mlo_edges[[edge_idx]] <- list(
            source = i - 1,
            target = j - 1,
            value = adj_matrix_best[i, j],
            weight = format(adj_matrix_best[i, j], digits = 2),
            source_label = as.character(orders[i]),
            target_label = as.character(orders[j]),
            direction = direction,
            is_forward = is_forward
          )
          edge_idx <- edge_idx + 1
        }
      }
    }
    
    # 存储所有数据
    all_data[[trans_id]] <- list(
      gene_symbol = ifelse(is.null(data$gene_symbol), trans_id, 
                           ifelse(is.numeric(data$gene_symbol), 
                                  as.character(data$gene_symbol), 
                                  data$gene_symbol)),
      trans_id = trans_id,
      strand = ifelse(is.null(data$strand), "unknown", data$strand),
      best_order = as.character(orders),
      n_introns = n,
      statistics = list(
        p_value = ifelse(is.null(data$p_value), NA, data$p_value),
        chi_stat = ifelse(is.null(data$chi_stat), NA, data$chi_stat),
        entropy = ifelse(is.null(data$entropy), NA, data$entropy),
        spearman_cor = spearman_cor,
        spearman_p_less = spearman_p_less,
        spearman_p_greater = spearman_p_greater,
        percent_coverage_pair = ifelse(is.null(data$percent_coverage_pair), NA, data$percent_coverage_pair),
        number_of_maximum_order = ifelse(is.null(data$number_of_maximum_order), NA, data$number_of_maximum_order),
        spearman_rho = ifelse(is.null(data$spearman_rho), NA, data$spearman_rho),
        spearman_rho_abs = ifelse(is.null(data$spearman_rho_abs), NA, data$spearman_rho_abs),
        entropy_sum_normalized = ifelse(is.null(data$entropy_sum_normalized), NA, data$entropy_sum_normalized)
      ),
      # 最佳顺序的矩阵
      adjacency_matrix_best = as.matrix(adj_matrix_best),
      frequency_matrix_best = as.matrix(freq_matrix_best),
      
      adjacency_matrix_freq = as.matrix(adj_matrix_freq),
      frequency_matrix_freq = as.matrix(freq_matrix_freq),
      # ggraph风格的MLO net数据
      mlo_nodes = mlo_nodes,
      mlo_edges = mlo_edges,
      # 原始顺序标签
      original_labels = rownames(orig_adj_matrix)
    )
  }
  
  cat("Successfully processed", length(all_data), "transcripts\n")
  return(all_data)
}


# 生成HTML页面 - 使用文件分段写入


# 生成HTML页面 - 使用文件分段写入
generate_html_page <- function(all_data, output_file = "splicing_order_report.html") {
  
  cat("Generating HTML file...\n")
  
  con <- file(output_file, "w", encoding = "UTF-8")
  
  writeLines('<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Intron Splicing Order Analysis</title>
    <link rel="icon" type="image/x-icon" href="https://img.icons8.com/color/48/000000/dna-helix.png">
<link rel="shortcut icon" type="image/x-icon" href="https://img.icons8.com/color/48/000000/dna-helix.png">
    <script src="https://d3js.org/d3.v7.min.js"></script>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
    
    <style>
        body { padding: 20px; font-family: Arial, sans-serif; background-color: #f5f5f5; }
        .header {
            background: linear-gradient(135deg, #4CAF50 0%, #2E7D32 100%);
            color: white; padding: 30px; border-radius: 10px; margin-bottom: 20px;
        }
        .mlo-svg {
            width: 100%; height: 700px; background: white; border-radius: 5px;
            overflow: visible !important;
        }
        .mlo-node {
            fill: #4CAF50; stroke: white; stroke-width: 2px;
            transition: all 0.3s ease;
        }
        .mlo-node:hover {
            stroke: #FF5722; stroke-width: 4px; transform: scale(1.1);
        }
        .mlo-node-text {
            font-size: 12px; font-weight: bold; fill: white;
            text-anchor: middle; pointer-events: none;
            text-shadow: 1px 1px 2px rgba(0,0,0,0.5);
        }
        .mlo-link { fill: none; stroke-opacity: 0.7; transition: all 0.3s ease; }
        .forward-link, .backward-link { stroke: #2196F3; }
        .zero-link {
            stroke: #999 !important; stroke-dasharray: 4,4 !important;
            stroke-opacity: 0.7 !important; stroke-width: 1.5 !important;
        }
        .mlo-link-text {
            font-size: 10px; font-weight: bold; fill: #333;
            text-anchor: middle; pointer-events: none;
            text-shadow: 1px 1px 2px rgba(255,255,255,0.8);
        }
        .zero-link-text { fill: #666 !important; font-size: 9px !important; }
        .transcript-item {
            padding: 5px 10px; border-bottom: 1px solid #eee; cursor: pointer;
            transition: background-color 0.2s;
        }
        .transcript-item:hover { background-color: #f5f5f5; }
        .transcript-item.active {
            background-color: #e3f2fd; border-left: 3px solid #2196F3;
        }
        .text-bg {
            fill: white; fill-opacity: 0.9; stroke-width: 0.5; stroke-opacity: 0.7;
        }
        .zero-text-bg { fill: #f5f5f5 !important; fill-opacity: 0.8 !important; }
        .stat-box {
            background: white; padding: 15px; border-radius: 8px;
            margin: 10px 0; border-left: 4px solid #4CAF50;
        }
        .order-sequence {
            font-family: "Courier New", monospace; background: #f1f8e9;
            padding: 10px; border-radius: 5px; margin: 10px 0;
        }
        .matrix-table {
            font-size: 11px; width: 100%;
        }
        .matrix-table th, .matrix-table td {
            text-align: center; padding: 4px; border: 1px solid #dee2e6;
        }
        .matrix-table th {
            background-color: #f8f9fa; font-weight: 600;
        }
        .tab-content {
            background: white; padding: 20px; border-radius: 0 0 8px 8px;
            border: 1px solid #dee2e6; border-top: none;
        }
        .nav-tabs {
            border-bottom: 1px solid #dee2e6;
        }
        .nav-tabs .nav-link {
            border: none; color: #6c757d; font-weight: 500;
        }
        .nav-tabs .nav-link.active {
            color: #4CAF50; border-bottom: 3px solid #4CAF50;
            background: transparent;
        }
        .mlo-net-container {
            background: white; border-radius: 8px; padding: 15px;
            border: 1px solid #dee2e6; margin-bottom: 20px;
        }
        .matrix-controls { margin-bottom: 15px; }
        .mlo-controls {
            background: #f8f9fa; padding: 10px; border-radius: 5px;
            margin-bottom: 15px;
        }
        .slider-container { margin: 5px 0; }
        .slider-label {
            font-size: 12px; margin-bottom: 2px; font-weight: 500;
        }
        .gene-name { font-weight: bold; color: #333; }
        .transcript-id {
            font-size: 0.85em; color: #666; font-family: monospace;
        }
        .transcript-stats { font-size: 0.8em; color: #888; }
    </style>
</head>
<body>
    <div class="container-fluid">
        <div class="header">
            <h3>Intron Splicing Order Analysis</h3>', con)
  
  writeLines(sprintf('<p class="lead mb-0">Transcripts: %d</p>
        </div>', length(all_data)), con)
  
  writeLines('        
        <div class="row">
            <div class="col-md-3">
                <div class="card">
                    <div class="card-header">
                        <h5 class="mb-0">Transcripts</h5>
                        <div class="form-group mt-2">
                            <input type="text" class="form-control" id="searchInput" placeholder="Search genes or transcripts...">
                        </div>
                    </div>
                    <div class="card-body" style="max-height: 600px; overflow-y: auto;">
                        <div id="transcriptList"></div>
                    </div>
                </div>
            </div>
            
            <div class="col-md-9">
                <div id="detailPanel">
                    <div class="card">
                        <div class="card-body text-center py-5">
                            <div class="display-1 text-muted mb-3">📊</div>
                            <h4 class="text-muted mb-3">Select a transcript to view details</h4>
                            <p class="text-muted">Click on a transcript from the left panel</p>
                        </div>
                    </div>
                </div>
            </div>
        </div>
        
        <div class="row mt-4">
            <div class="col-12 text-center text-muted">
                <p>Report generated on ', con)
  
  writeLines(as.character(Sys.Date()), con)
  
  writeLines(' | Single HTML file - No backend required</p>
            </div>
        </div>
    </div>', con)
  
  writeLines('    
    <script>
        const allData = ', con)
  
  json_str <- toJSON(all_data, auto_unbox = TRUE, force = TRUE)
  writeLines(json_str, con)
  
  writeLines(';
        let currentData = null;
        let currentMatrixOrder = "best";
        let svg = null;
        
        let mloParams = {
            nodeRadius: 10,
            nodeSpacing: 60,
            arcHeight: 200,
            fontSize: 12,
            linkOpacity: 0.7,
            linkWidthMultiplier: 0.3,
            width: 1000,
            height: 600
        };
        
        // ========== 辅助函数 ==========
        function formatNumber(num) {
            if (num === null || num === undefined || isNaN(num)) return "N/A";
            if (Math.abs(num) < 0.001 && num !== 0) return num.toExponential(2);
            return Number(num).toFixed(3);
        }
        
        function formatPValue(p) {
            if (p === null || p === undefined || isNaN(p)) return "N/A";
            if (p < 0.001) return p.toExponential(2);
            return Number(p).toFixed(3);
        }
        
        function sortTranscripts(transcripts) {
            return transcripts.sort((a, b) => {
                const dataA = allData[a];
                const dataB = allData[b];
                const geneCompare = dataA.gene_symbol.localeCompare(dataB.gene_symbol);
                if (geneCompare !== 0) return geneCompare;
                return a.localeCompare(b);
            });
        }
        
        function getCurrentMatrix(type) {
            if (!currentData) return [];
            if (currentMatrixOrder === "frequency") {
                if (type === "read") {
                    return currentData.adjacency_matrix_freq || [];
                } else if (type === "frequency") {
                    return currentData.frequency_matrix_freq || [];
                }
            } else {
                if (type === "read") {
                    return currentData.adjacency_matrix_best || [];
                } else if (type === "frequency") {
                    return currentData.frequency_matrix_best || [];
                }
            }
            return [];
        }
        
        function getCurrentOrder() {
            if (!currentData) return [];
            if (currentMatrixOrder === "frequency") {
                const matrix = currentData.adjacency_matrix_freq;
                const n = matrix ? matrix.length : 0;
                return Array.from({length: n}, (_, i) => `Intron ${i + 1}`);
            } else {
                return currentData.best_order || [];
            }
        }
        
        function renderMatrixTable(matrix, order, type) {
            if (!matrix || matrix.length === 0 || !order || order.length === 0) {
                return `<div class="text-center text-muted py-3">No matrix data available</div>`;
            }
            
            let html = `<table class="matrix-table"><thead><tr><th>From \\\\ To</th>`;
            
            order.forEach(label => {
                html += `<th>${label}</th>`;
            });
            html += `</tr></thead><tbody>`;
            
            for (let i = 0; i < order.length; i++) {
                html += `<tr><th>${order[i]}</th>`;
                for (let j = 0; j < order.length; j++) {
                    const value = matrix[i][j];
                    let displayValue = "0";
                    let cellStyle = "";
                    
                    if (type === "frequency" && i === j) {
                        displayValue = "-";
                        cellStyle = "background-color: #f8f9fa;";
                    } else if (value > 0) {
                        displayValue = formatNumber(value);
                        if (type === "read") {
                            const maxVal = Math.max(...matrix.flat());
                            const intensity = maxVal > 0 ? Math.min(value / maxVal, 0.7) : 0;
                            cellStyle = `background-color: rgba(76, 175, 80, ${intensity});`;
                        } else {
                            const intensity = Math.min(value, 0.7);
                            cellStyle = `background-color: rgba(76, 175, 80, ${intensity});`;
                        }
                    }
                    
                    html += `<td style="${cellStyle}">${displayValue}</td>`;
                }
                html += `</tr>`;
            }
            
            html += `</tbody></table>`;
            return html;
        }
        
        function renderStatsTab(data) {
            const stats = data.statistics;
            let statsHtml = "";
            
            if (stats.spearman_cor !== null && !isNaN(stats.spearman_cor)) {
                statsHtml += `<tr><td>Spearman ρ:</td><td class="text-end">${formatNumber(stats.spearman_cor)}</td></tr>`;
            }
            if (stats.spearman_rho_abs !== null && !isNaN(stats.spearman_rho_abs)) {
                statsHtml += `<tr><td>Spearman |ρ|:</td><td class="text-end">${formatNumber(stats.spearman_rho_abs)}</td></tr>`;
            }
            if (stats.spearman_p_less !== null && !isNaN(stats.spearman_p_less)) {
                statsHtml += `<tr><td>P-value (less):</td><td class="text-end">${formatPValue(stats.spearman_p_less)}</td></tr>`;
            }
            if (stats.spearman_p_greater !== null && !isNaN(stats.spearman_p_greater)) {
                statsHtml += `<tr><td>P-value (greater):</td><td class="text-end">${formatPValue(stats.spearman_p_greater)}</td></tr>`;
            }
            if (stats.chi_stat !== null && !isNaN(stats.chi_stat)) {
                statsHtml += `<tr><td>Log Likelihood:</td><td class="text-end">${formatNumber(stats.chi_stat)}</td></tr>`;
            }
            if (stats.p_value !== null && !isNaN(stats.p_value)) {
                statsHtml += `<tr><td>P-value (Chi):</td><td class="text-end">${formatPValue(stats.p_value)}</td></tr>`;
            }
            
            let entropyHtml = "";
            if (stats.entropy !== null && !isNaN(stats.entropy)) {
                entropyHtml += `<tr><td>Entropy:</td><td class="text-end">${formatNumber(stats.entropy)}</td></tr>`;
            }
            if (stats.entropy_sum_normalized !== null && !isNaN(stats.entropy_sum_normalized)) {
                entropyHtml += `<tr><td>Normalized Entropy:</td><td class="text-end">${formatNumber(stats.entropy_sum_normalized)}</td></tr>`;
            }
            if (stats.percent_coverage_pair !== null && !isNaN(stats.percent_coverage_pair)) {
                entropyHtml += `<tr><td>Pair Coverage:</td><td class="text-end">${formatNumber(stats.percent_coverage_pair * 100)}%</td></tr>`;
            }
            if (stats.number_of_maximum_order !== null) {
                entropyHtml += `<tr><td>Max Orders:</td><td class="text-end">${stats.number_of_maximum_order || "N/A"}</td></tr>`;
            }
            if (data.n_introns > 1) {
                entropyHtml += `<tr><td>Intron Pairs:</td><td class="text-end">${data.n_introns * (data.n_introns - 1) / 2}</td></tr>`;
            }
            
            return `
                <div class="row">
                    <div class="col-md-6">
                        <div class="stat-box">
                            <h5>Order Statistics</h5>
                            <table class="table table-sm">${statsHtml}</table>
                        </div>
                    </div>
                    <div class="col-md-6">
                        <div class="stat-box">
                            <h5>Entropy & Coverage</h5>
                            <table class="table table-sm">${entropyHtml}</table>
                        </div>
                    </div>
                </div>
            `;
        }
        
        // ========== 核心函数：标签切换 ==========
        function showTab(tabName, event) {
            if (event) event.preventDefault();
            
            // 更新标签按钮状态
            document.querySelectorAll("#detailTabs .nav-link").forEach(link => {
                link.classList.remove("active");
            });
            if (event && event.target) {
                event.target.classList.add("active");
            }
            
            // 隐藏所有标签内容
            ["overviewTab", "mlo_netTab", "matrixTab", "frequencyTab", "statsTab"].forEach(tabId => {
                const tab = document.getElementById(tabId);
                if (tab) tab.style.display = "none";
            });
            
            // 显示目标标签内容
            const targetTab = document.getElementById(tabName + "Tab");
            if (targetTab) {
                targetTab.style.display = "block";
                
                // 如果是MLO Net标签，渲染图形
                if (tabName === "mlo_net" && currentData) {
                    setTimeout(() => renderMLONet(), 50);
                }
                // 如果是Matrix标签，渲染矩阵
                else if (tabName === "matrix" && currentData) {
                    const order = getCurrentOrder();
                    const matrix = getCurrentMatrix("read");
                    const matrixContent = document.getElementById("readMatrixContent");
                    if (matrixContent) {
                        matrixContent.innerHTML = renderMatrixTable(matrix, order, "read");
                    }
                }
                // 如果是Frequency标签，渲染频率矩阵
                else if (tabName === "frequency" && currentData) {
                    const order = getCurrentOrder();
                    const matrix = getCurrentMatrix("frequency");
                    const matrixContent = document.getElementById("frequencyMatrixContent");
                    if (matrixContent) {
                        matrixContent.innerHTML = renderMatrixTable(matrix, order, "frequency");
                    }
                }
                // 如果是Stats标签，渲染统计
                else if (tabName === "stats" && currentData) {
                    const statsContent = document.getElementById("statsTab");
                    if (statsContent) {
                        statsContent.innerHTML = renderStatsTab(currentData);
                    }
                }
            }
        }
        
        // ========== 核心函数：矩阵排序切换 ==========
        function setMatrixOrder(orderType) {
            currentMatrixOrder = orderType;
            
            // 更新按钮状态
            document.querySelectorAll(".matrix-controls .btn").forEach(btn => {
                const btnText = btn.textContent.trim();
                if (btnText.includes("Best Order")) {
                    btn.className = currentMatrixOrder === "best" ? 
                                   "btn btn-sm btn-primary" : 
                                   "btn btn-sm btn-outline-primary";
                } else if (btnText.includes("Order by default")) {
                    btn.className = currentMatrixOrder === "frequency" ? 
                                   "btn btn-sm btn-primary" : 
                                   "btn btn-sm btn-outline-primary";
                }
            });
            

            
            // 重新渲染当前激活的标签页
            const activeTab = document.querySelector(".nav-link.active");
            if (activeTab) {
                const tabText = activeTab.textContent.toLowerCase();
                const order = getCurrentOrder();
                
                if (tabText.includes("matrix") || activeTab.textContent.includes("Read Matrix")) {
                    const matrix = getCurrentMatrix("read");
                    const matrixContent = document.getElementById("readMatrixContent");
                    if (matrixContent) {
                        matrixContent.innerHTML = renderMatrixTable(matrix, order, "read");
                    }
                } else if (tabText.includes("frequency")) {
                    const matrix = getCurrentMatrix("frequency");
                    const matrixContent = document.getElementById("frequencyMatrixContent");
                    if (matrixContent) {
                        matrixContent.innerHTML = renderMatrixTable(matrix, order, "frequency");
                    }
                }
            }
        }
        
            function calculateDynamicNodeSpacing(nIntrons) {
              // 基础公式：节点越多，间距越小
              // 最小间距 40，最大间距 120
              const minSpacing = 20;
              const maxSpacing = 400;
              const baseSpacing = 80;
              
              if (nIntrons <= 3) return maxSpacing;
              
              if (nIntrons <= 5) return 200;


              if (nIntrons <= 10) return 100;

              
              if (nIntrons <= 15) return 80;
              
              
              if (nIntrons <= 30) return 60;
              
              return minSpacing;
              // 线性递减：每增加一个节点，减少 4 个像素间距
              return Math.max(minSpacing, maxSpacing - (nIntrons - 3) * 4);
            }        
        
        
        // ========== MLO网络图渲染函数 ==========
        function renderMLONet() {
            if (!currentData) return;
            
            const container = document.getElementById("mloNetGraph");
            if (!container) return;
            
            container.innerHTML = "";
            
            const rawNodes = currentData.mlo_nodes || [];
            const edges = currentData.mlo_edges || [];
            const bestOrder = currentData.best_order || [];
            
            if (rawNodes.length === 0) {
                container.innerHTML = "<div class=\'text-center text-muted py-5\'>No network data available</div>";
                return;
            }
            
            const nIntrons = bestOrder.length;
                
                
            if (typeof mloParams.nodeSpacingUserSet === undefined || mloParams.nodeSpacingUserSet === false) {
                const nIntrons = currentData.best_order.length;
                if (nIntrons > 0) {
                    const dynamicSpacing = calculateDynamicNodeSpacing(nIntrons);
                    mloParams.nodeSpacing = dynamicSpacing;
                    
                    // 更新界面控件显示
                    const nodeSpacingSlider = document.getElementById("nodeSpacingSlider");
                    const nodeSpacingValueSpan = document.getElementById("nodeSpacingValue");
                    if (nodeSpacingSlider) nodeSpacingSlider.value = dynamicSpacing;
                    if (nodeSpacingValueSpan) nodeSpacingValueSpan.textContent = dynamicSpacing;
                }
            }
    
            
            
            
            // 1. 按best_order顺序排序节点
            const sortedNodes = [];
            const nodeMap = {};
            
            rawNodes.forEach(node => {
                const label = node.label || node.original_label || "";
                nodeMap[label] = node;
            });
            
            bestOrder.forEach((label, index) => {
                if (nodeMap[label]) {
                    const node = {...nodeMap[label]};
                    node.x = 0;
                    node.y = 0;
                    node.radius = mloParams.nodeRadius;
                    node.bestOrderIndex = index;
                    node.displayLabel = label;
                    sortedNodes.push(node);
                }
            });
            
            // 处理best_order中不存在的节点
            rawNodes.forEach(node => {
                const label = node.label || node.original_label || "";
                if (!bestOrder.includes(label) && nodeMap[label]) {
                    const nodeCopy = {...nodeMap[label]};
                    nodeCopy.x = 0;
                    nodeCopy.y = 0;
                    nodeCopy.radius = mloParams.nodeRadius;
                    nodeCopy.bestOrderIndex = 999;
                    nodeCopy.displayLabel = label;
                    sortedNodes.push(nodeCopy);
                }
            });
            
            // 2. 计算布局
            const totalWidth = sortedNodes.length * mloParams.nodeSpacing;
            const startX = 20;
            const centerY = mloParams.height / 2;
            
            sortedNodes.forEach((node, i) => {
                node.x = startX + i * mloParams.nodeSpacing;
                node.y = centerY;
                node.id = i;
            });
            
            // 3. 创建SVG
            svg = d3.select("#mloNetGraph")
                .append("svg")
                .attr("width", mloParams.width)
                .attr("height", mloParams.height)
                .attr("viewBox", `0 0 ${mloParams.width} ${mloParams.height}`)
                .attr("preserveAspectRatio", "xMidYMid meet")
                .attr("overflow", "visible");
            
            // 4. 定义箭头
            const defs = svg.append("defs");
            
            defs.append("marker").attr("id", "forward-arrow")
                .attr("viewBox", "0 -5 10 10").attr("refX", 8).attr("refY", 0)
                .attr("markerWidth", 4).attr("markerHeight", 4).attr("orient", "auto")
                .append("path").attr("d", "M0,-5L10,0L0,5")
                .attr("class", "forward-arrow").attr("fill", "#2196F3");
            
            defs.append("marker").attr("id", "backward-arrow")
                .attr("viewBox", "0 -5 10 10").attr("refX", 2).attr("refY", 0)
                .attr("markerWidth", 4).attr("markerHeight", 4).attr("orient", "auto")
                .append("path").attr("d", "M10,-5L0,0L10,5")
                .attr("class", "backward-arrow").attr("fill", "#2196F3");
            
            defs.append("marker").attr("id", "zero-arrow")
                .attr("viewBox", "0 -5 10 10").attr("refX", 8).attr("refY", 0)
                .attr("markerWidth", 3.5).attr("markerHeight", 3.5).attr("orient", "auto")
                .append("path").attr("d", "M0,-5L10,0L0,5")
                .attr("class", "zero-arrow").attr("fill", "#999");
            
            // 5. 创建弧线路径
            function createArcPath(source, target, isForward) {
                const x1 = source.x, y1 = source.y;
                const x2 = target.x, y2 = target.y;
                const distance = Math.abs(x2 - x1);
                
                let baseArcHeight = mloParams.arcHeight;
                const distanceRatio = distance / mloParams.nodeSpacing;
                let dynamicFactor = 1;
                
                if (distanceRatio > 1) dynamicFactor = 1 + Math.log(distanceRatio) * 0.8;
                const dynamicArcHeight = Math.min(baseArcHeight * dynamicFactor, baseArcHeight * 3);
                const midX = (x1 + x2) / 2;
                const arcDirection = isForward ? 1 : -1;
                const arcY = centerY - dynamicArcHeight * arcDirection;
                
                return `M ${x1} ${y1} Q ${midX} ${arcY} ${x2} ${y2}`;
            }
            
            // 6. 处理边数据
            const allEdges = [];
            const n = sortedNodes.length;
            const labelToNewIndex = {};
            
            sortedNodes.forEach((node, idx) => {
                labelToNewIndex[node.displayLabel] = idx;
            });
            
            edges.forEach(edge => {
                const sourceLabel = edge.source_label;
                const targetLabel = edge.target_label;
                
                if (labelToNewIndex[sourceLabel] !== undefined && 
                    labelToNewIndex[targetLabel] !== undefined) {
                    
                    const newSourceIdx = labelToNewIndex[sourceLabel];
                    const newTargetIdx = labelToNewIndex[targetLabel];
                    const isForward = newSourceIdx < newTargetIdx;
                    
                    allEdges.push({
                        source: newSourceIdx,
                        target: newTargetIdx,
                        value: edge.value,
                        weight: edge.weight,
                        source_label: sourceLabel,
                        target_label: targetLabel,
                        direction: isForward ? "forward" : "backward",
                        is_forward: isForward
                    });
                }
            });
            
            // 添加0值的边
            for (let i = 0; i < n; i++) {
                for (let j = 0; j < n; j++) {
                    if (i !== j) {
                        const existingEdge = allEdges.find(e => e.source === i && e.target === j);
                        const isForward = i < j;
                        
                        if (!existingEdge) {
                            allEdges.push({
                                source: i,
                                target: j,
                                value: 0,
                                weight: "0",
                                source_label: sortedNodes[i].displayLabel,
                                target_label: sortedNodes[j].displayLabel,
                                direction: isForward ? "forward" : "backward",
                                is_forward: isForward
                            });
                        }
                    }
                }
            }
            
            // 7. 创建链接
            const linkGroup = svg.append("g").attr("class", "mlo-links");
            
            const forwardLinks = allEdges.filter(d => d.is_forward);
            const backwardLinks = allEdges.filter(d => !d.is_forward);
            
            linkGroup.selectAll("path.forward-link")
                .data(forwardLinks.filter(d => d.value > 0))
                .enter().append("path")
                .attr("class", "mlo-link forward-link")
                .attr("d", d => createArcPath(sortedNodes[d.source], sortedNodes[d.target], true))
                .attr("stroke-width", d => Math.max(1, Math.sqrt(d.value) * mloParams.linkWidthMultiplier + 1))
                .attr("stroke-opacity", mloParams.linkOpacity)
                .attr("marker-end", "url(#forward-arrow)");
            
            linkGroup.selectAll("path.forward-zero-link")
                .data(forwardLinks.filter(d => d.value === 0))
                .enter().append("path")
                .attr("class", "mlo-link forward-link zero-link")
                .attr("d", d => createArcPath(sortedNodes[d.source], sortedNodes[d.target], true))
                .attr("stroke-width", 1.5).attr("stroke-opacity", 0.7)
                .attr("marker-end", "url(#zero-arrow)");
            
            linkGroup.selectAll("path.backward-link")
                .data(backwardLinks.filter(d => d.value > 0))
                .enter().append("path")
                .attr("class", "mlo-link backward-link")
                .attr("d", d => createArcPath(sortedNodes[d.source], sortedNodes[d.target], false))
                .attr("stroke-width", d => Math.max(1, Math.sqrt(d.value) * mloParams.linkWidthMultiplier + 1))
                .attr("stroke-opacity", mloParams.linkOpacity)
                .attr("marker-end", "url(#backward-arrow)");
            
            linkGroup.selectAll("path.backward-zero-link")
                .data(backwardLinks.filter(d => d.value === 0))
                .enter().append("path")
                .attr("class", "mlo-link backward-link zero-link")
                .attr("d", d => createArcPath(sortedNodes[d.source], sortedNodes[d.target], false))
                .attr("stroke-width", 1.5).attr("stroke-opacity", 0.7)
                .attr("marker-end", "url(#zero-arrow)");
            
            // 8. 创建链接文本
            const linkTextGroup = linkGroup.append("g").attr("class", "mlo-link-texts");
            
            linkTextGroup.selectAll("text")
                .data(allEdges)
                .enter().append("text")
                .attr("class", d => `mlo-link-text ${d.value === 0 ? "zero-link-text" : ""}`)
                .text(d => d.weight || d.value)
                .attr("font-size", d => d.value === 0 ? "9px" : "11px")
                .attr("font-weight", "bold")
                .attr("fill", d => d.value === 0 ? "#666" : "#2196F3");
            
            // 9. 创建节点
            const nodeGroup = svg.append("g").attr("class", "mlo-nodes");
            
            nodeGroup.selectAll("circle")
                .data(sortedNodes)
                .enter().append("circle")
                .attr("class", "mlo-node")
                .attr("cx", d => d.x).attr("cy", d => d.y)
                .attr("font-size", mloParams.fontSize + "px")  // 关键：使用参数而不是固定值
                .attr("r", d => d.radius)
                .attr("fill", "#4CAF50").attr("stroke", "white")
                .attr("stroke-width", 3);
            
            nodeGroup.selectAll("text.mlo-node-text")
                .data(sortedNodes)
                .enter().append("text")
                .attr("class", "mlo-node-text")
                .text(d => d.displayLabel)
                .attr("x", d => d.x).attr("y", d => d.y)
                .attr("font-size", mloParams.fontSize + "px")
                .attr("dy", "0.35em").attr("text-anchor", "middle")
                .style("pointer-events", "none");
            
            // 10. 更新链接文本位置
            setTimeout(() => {
                updateLinkTextPositions(sortedNodes, allEdges, svg, mloParams, centerY);
            }, 100);
            
            // 11. 添加图例和标题
            const legendX = 20, legendY = 40;
            const legend = svg.append("g").attr("transform", `translate(${legendX}, ${legendY})`);
            
            legend.append("text").attr("x", 0).attr("y", -15)
                .attr("font-size", "12px").attr("font-weight", "bold").attr("fill", "#333")
                .text("Connection Legend:");
            
            legend.append("line").attr("x1", 0).attr("y1", 0).attr("x2", 20).attr("y2", 0)
                .attr("stroke", "#2196F3").attr("stroke-width", 2)
                .attr("marker-end", "url(#forward-arrow)");
            legend.append("text").attr("x", 25).attr("y", 0).attr("dy", "0.35em")
                .attr("font-size", "11px").attr("fill", "#333")
                .text("Splicing direction (non-zero)");
            
            legend.append("line").attr("x1", 0).attr("y1", 20).attr("x2", 20).attr("y2", 20)
                .attr("stroke", "#999").attr("stroke-width", 1.5)
                .attr("stroke-dasharray", "4,4").attr("stroke-opacity", 0.7)
                .attr("marker-end", "url(#zero-arrow)");
            legend.append("text").attr("x", 25).attr("y", 20).attr("dy", "0.35em")
                .attr("font-size", "11px").attr("fill", "#666")
                .text("Zero read count");
            
            svg.append("text").attr("x", mloParams.width / 2).attr("y", 25)
                .attr("text-anchor", "middle").attr("font-size", "16px")
                .attr("font-weight", "bold").attr("fill", "#333")
                .text("");

                //.text("Most Likely Order Network");
            
            svg.append("text").attr("x", mloParams.width / 2).attr("y", mloParams.height - 15)
                .attr("text-anchor", "middle").attr("font-size", "11px")
                .attr("fill", "#666")
                .text("");
                //.text("Nodes arranged in best order sequence • Hover to highlight connections");
                
        }
        
        function updateLinkTextPositions(nodes, edges, svg, mloParams, centerY) {
            const linkTexts = svg.selectAll(".mlo-link-text");
            
            linkTexts.each(function(d) {
                const source = nodes[d.source];
                const target = nodes[d.target];
                if (!source || !target) return;
                
                const midX = (source.x + target.x) / 2;
                const distance = Math.abs(target.x - source.x);
                
                let baseArcHeight = mloParams.arcHeight;
                const distanceRatio = distance / mloParams.nodeSpacing;
                let dynamicFactor = 1;
                
                if (distanceRatio > 1) dynamicFactor = 1 + Math.log(distanceRatio) * 0.8;
                const dynamicArcHeight = Math.min(baseArcHeight * dynamicFactor, baseArcHeight * 3);
                const arcDirection = d.is_forward ? 1 : -1;
                const textY = centerY - dynamicArcHeight * arcDirection / 2;
                
                const fontSize = Math.min(Math.max(9, Math.sqrt(d.value) * 1.5), 13);
                const textElement = d3.select(this);
                
                textElement.attr("x", midX).attr("y", textY)
                    .attr("font-size", fontSize + "px").attr("dy", "0.35em")
                    .attr("text-anchor", "middle").attr("alignment-baseline", "middle");
                
                const bbox = this.getBBox();
                textElement.selectAll("rect").remove();
                
                const bgClass = d.value === 0 ? "zero-text-bg" : "text-bg";
                const bgFill = d.value === 0 ? "#f5f5f5" : "white";
                
                textElement.insert("rect", ":first-child")
                    .attr("x", bbox.x - 3).attr("y", bbox.y - 1)
                    .attr("width", bbox.width + 6).attr("height", bbox.height + 2)
                    .attr("rx", 3).attr("ry", 3).attr("class", bgClass)
                    .attr("fill", bgFill).attr("fill-opacity", d.value === 0 ? 0.8 : 0.9)
                    .attr("stroke", "#2196F3").attr("stroke-width", 0.5);
            });
        }
        
        // ========== 页面初始化函数 ==========
        //# 在 renderTranscriptList 函数中修改转录本条目的HTML
      function renderTranscriptList() {
          const listDiv = document.getElementById("transcriptList");
          const transcripts = Object.keys(allData);
          const sortedTranscripts = sortTranscripts(transcripts);
          
          listDiv.innerHTML = "";
          
          sortedTranscripts.forEach(transId => {
              const data = allData[transId];
              const item = document.createElement("div");
              item.className = "transcript-item";
              item.innerHTML = `
                  <div>
                      <div class="gene-name">${data.gene_symbol}</div>
                      <div class="transcript-id">${transId}</div>
                      <div class="transcript-stats">
                          <span class="badge bg-light text-dark">${data.n_introns} introns</span>
                          <span class="badge bg-light text-dark ms-1">ρ=${formatNumber(data.statistics.spearman_cor)}</span>
                      </div>
                      <!-- 新增：显示best orders -->
                      <div class="best-order mt-1 small text-muted" style="font-family: monospace; word-break: break-all;">
                          ${data.best_order.join(" → ")}
                      </div>
                  </div>
              `;
              item.addEventListener("click", function() {
                  loadTranscriptDetail(transId);
                  document.querySelectorAll(".transcript-item").forEach(i => i.classList.remove("active"));
                  item.classList.add("active");
              });
              listDiv.appendChild(item);
          });
          
          if (sortedTranscripts.length > 0) {
              loadTranscriptDetail(sortedTranscripts[0]);
              const firstItem = listDiv.querySelector(".transcript-item");
              if (firstItem) firstItem.classList.add("active");
          }
      }
        
        function loadTranscriptDetail(transId) {
            const data = allData[transId];
            currentData = data;
            currentMatrixOrder = "best";
            
            const detailHTML = `
                <div class="card">
                    <div class="card-header">
                        <h4 class="mb-0">${data.gene_symbol} - ${transId}
                            <span class="badge bg-secondary float-end">${data.strand}</span>
                        </h4>
                    </div>
                    <div class="card-body">
                        <ul class="nav nav-tabs mb-3" id="detailTabs">
                            <li class="nav-item"><button class="nav-link active" onclick="showTab(\'overview\', event)">Overview</button></li>
                            <li class="nav-item"><button class="nav-link" onclick="showTab(\'mlo_net\', event)">MLO Net</button></li>
                            <li class="nav-item"><button class="nav-link" onclick="showTab(\'matrix\', event)">Read Matrix</button></li>
                            <li class="nav-item"><button class="nav-link" onclick="showTab(\'frequency\', event)">Frequency</button></li>
                            <li class="nav-item"><button class="nav-link" onclick="showTab(\'stats\', event)">Detailed Stats</button></li>
                        </ul>
                        <div id="tabContent">
                            <div id="overviewTab">
                                <div class="row">
                                    <div class="col-12">
                                        <div class="stat-box">
                                            <h5>Most Likely Splicing Order</h5>
                                            <div class="order-sequence">${data.best_order.join(" → ")}</div>
                                            <small class="text-muted">${data.n_introns} introns in optimal order</small>
                                        </div>
                                    </div>
                                </div>
                                <div class="row mt-3">
                                    <div class="col-md-6">
                                        <div class="stat-box">
                                            <h5>Key Statistics</h5>
                                            <table class="table table-sm">
                                                <tr><td>Spearman ρ:</td><td class="text-end">${formatNumber(data.statistics.spearman_cor)}</td></tr>
                                                <tr><td>P-value (less):</td><td class="text-end">${formatPValue(data.statistics.spearman_p_less)}</td></tr>
                                                <tr><td>Entropy:</td><td class="text-end">${formatNumber(data.statistics.entropy)}</td></tr>
                                            </table>
                                        </div>
                                    </div>
                                    <div class="col-md-6">
                                        <div class="stat-box">
                                            <h5>Additional Information</h5>
                                            <table class="table table-sm">
                                                <tr><td>Transcript ID:</td><td class="text-end">${data.trans_id}</td></tr>
                                                <tr><td>Gene Symbol:</td><td class="text-end">${data.gene_symbol}</td></tr>
                                                <tr><td>Strand:</td><td class="text-end">${data.strand}</td></tr>
                                                <tr><td>Number of Introns:</td><td class="text-end">${data.n_introns}</td></tr>
                                            </table>
                                        </div>
                                    </div>
                                </div>
                            </div>
                            <div id="mlo_netTab" style="display: none;">
                                <div class="mlo-net-container">
                                    <h5>MLO Network (ggraph style)</h5>
                                    <p class="text-muted">Linear layout with arc connections - Blue arrows: splicing direction, Gray dashed: zero read count</p>
                                    <div class="mlo-controls">
                                        <div class="row">
                                            <div class="col-md-2">
                                                <div class="slider-label">Node Radius: <span id="nodeRadiusValue">${mloParams.nodeRadius}</span></div>
                                                <input type="range" class="form-range" min="5" max="30" step="1" value="${mloParams.nodeRadius}" oninput="updateNodeRadius(this.value)">
                                            </div>
                                            <div class="col-md-2">
                                                <div class="slider-label">Arc Height: <span id="arcHeightValue">${mloParams.arcHeight}</span></div>
                                                <input type="range" class="form-range" min="10" max="500" step="10" value="${mloParams.arcHeight}" oninput="updateArcHeight(this.value)">
                                            </div>
                                            <div class="col-md-2">
                                                <div class="slider-label">Link Width: <span id="linkWidthValue">${mloParams.linkWidthMultiplier}</span></div>
                                                <input type="range" class="form-range" min="0.1" max="3" step="0.1" value="${mloParams.linkWidthMultiplier}" oninput="updateLinkWidth(this.value)">
                                            </div>
                                            <div class="col-md-2">
                                                <div class="slider-label">Node Spacing: <span id="nodeSpacingValue">${mloParams.nodeSpacing}</span></div>
                                                <input type="range" class="form-range" min="50" max="500" step="10" value="${mloParams.nodeSpacing}" oninput="updateNodeSpacing(this.value)">
                                            </div>

                                            <div class="col-md-3">
                                                <div class="slider-label">Font Size: <span id="fontSizeValue">${mloParams.fontSize}</span></div>
                                                <input type="range" class="form-range" min="8" max="20" step="1" value="${mloParams.fontSize}" oninput="updateFontSize(this.value)">
                                            </div>
                                        </div>
                                    </div>
                                    <div id="mloNetGraph" class="mlo-svg"></div>
                                </div>
                            </div>
                            <div id="matrixTab" style="display: none;">
                                <div class="matrix-controls">
                                    <div class="btn-group" role="group">
                                        <button type="button" class="btn btn-sm btn-primary" onclick="setMatrixOrder(\'best\')">Best Order</button>
                                        <button type="button" class="btn btn-sm btn-outline-primary" onclick="setMatrixOrder(\'frequency\')">Order by default</button>
                                    </div>
                                </div>
                                <div id="readMatrixContent">Loading matrix...</div>
                            </div>
                            <div id="frequencyTab" style="display: none;">
                                <div class="matrix-controls">
                                    <div class="btn-group" role="group">
                                        <button type="button" class="btn btn-sm btn-primary" onclick="setMatrixOrder(\'best\')">Best Order</button>
                                        <button type="button" class="btn btn-sm btn-outline-primary" onclick="setMatrixOrder(\'frequency\')">Order by default</button>
                                    </div>
                                </div>
                                <div id="frequencyMatrixContent">Loading frequency matrix...</div>
                            </div>
                            <div id="statsTab" style="display: none;">Loading detailed statistics...</div>
                        </div>
                    </div>
                </div>
            `;
            
            document.getElementById("detailPanel").innerHTML = detailHTML;
            
            // 初始化Overview标签为激活状态
            setTimeout(() => {
                showTab("overview");
            }, 100);
        }
        
        // ========== MLO控制函数 ==========
        function updateNodeRadius(value) {
            mloParams.nodeRadius = parseInt(value);
            document.getElementById("nodeRadiusValue").textContent = value;
            renderMLONet();
        }
        
        function updateArcHeight(value) {
            mloParams.arcHeight = parseInt(value);
            document.getElementById("arcHeightValue").textContent = value;
            renderMLONet();
        }
        
        function updateLinkWidth(value) {
            mloParams.linkWidthMultiplier = parseFloat(value);
            document.getElementById("linkWidthValue").textContent = value;
            renderMLONet();
        }
        
        function updateFontSize(value) {
            mloParams.fontSize = parseInt(value);
            document.getElementById("fontSizeValue").textContent = value;
            
            // 直接更新已存在的文本元素，避免完全重新渲染
            //if (svg) {
            //    svg.selectAll(".mlo-node-text")
            //        .attr("font-size", mloParams.fontSize + "px");
            //} else {
                renderMLONet();  // 如果没有SVG，则重新渲染
            //}
        }
        
        function updateNodeSpacing(value) {
            mloParams.nodeSpacing = parseInt(value);
            mloParams.nodeSpacingUserSet = true;  // 标记为用户手动设置

            document.getElementById("nodeSpacingValue").textContent = value;
            renderMLONet();
        }
        
        // ========== 页面初始化 ==========
        document.addEventListener("DOMContentLoaded", function() {
            renderTranscriptList();
            
            document.getElementById("searchInput").addEventListener("input", function(e) {
                const searchTerm = e.target.value.toLowerCase();
                document.querySelectorAll(".transcript-item").forEach(item => {
                    const text = item.textContent.toLowerCase();
                    item.style.display = text.includes(searchTerm) ? "block" : "none";
                });
            });
        });
    </script>
</body>
</html>', con)
  
  close(con)
  cat("HTML report generated:", output_file, "\n")
}





