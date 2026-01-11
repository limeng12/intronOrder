// iso_compatible.cpp - 兼容现有R函数的版本

#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <fstream>

using namespace Rcpp;
using namespace std;

// ====================== 调试系统 ======================

// 调试级别控制
// 0: 无调试信息
// 1: 基本信息
// 2: 详细信息
// 3: 完整调试信息
#define DEBUG_LEVEL 2

// 全局调试日志文件流
static std::ofstream g_debug_log_file;
static bool g_debug_log_open = false;

// 初始化调试日志文件
void init_debug_log(const std::string& filename) {
  if (!g_debug_log_open) {
    g_debug_log_file.open(filename, std::ios::out | std::ios::trunc);
    if (g_debug_log_file.is_open()) {
      g_debug_log_open = true;
    } else {
      Rcpp::Rcerr << "Warning: Failed to open debug log file: " << filename << std::endl;
    }
  }
}

// 关闭日志文件
void close_debug_log() {
  if (g_debug_log_open) {
    g_debug_log_file.close();
    g_debug_log_open = false;
  }
}

// 调试打印函数
inline void debug_print(int level, const std::string& message) {
  if (level <= DEBUG_LEVEL && g_debug_log_open) {
    g_debug_log_file << "[DEBUG " << level << "] " << message << std::endl;
    g_debug_log_file.flush();
  }
}

// ====================== 辅助函数 ======================

inline vector<int> get_splice_junction_pos(const string& cigar, int read_start) {
  debug_print(3, "进入 get_splice_junction_pos");
  debug_print(3, "cigar: " + cigar);
  debug_print(3, "read_start: " + to_string(read_start));
  
  vector<int> junctions;
  int current_pos = read_start;
  string num_str;
  
  for (size_t i = 0; i < cigar.size(); i++) {
    char c = cigar[i];
    
    if (isdigit(c)) {
      num_str += c;
    } else {
      if (!num_str.empty()) {
        int length = stoi(num_str);
        
        if (c == 'N' && length > 10) {
          int intron_start = current_pos;
          int intron_end = current_pos + length - 1;
          
          junctions.push_back(intron_start);
          junctions.push_back(intron_end);
          
          debug_print(3, "找到内含子: " + to_string(intron_start) + "-" + to_string(intron_end));
        }
        
        if (c == 'M' || c == 'D' || c == 'N' || c == '=' || c == 'X') {
          current_pos += length;
        }
        
        num_str.clear();
      }
    }
  }
  
  debug_print(3, "离开 get_splice_junction_pos, 找到 " + to_string(junctions.size()) + " 个剪接点");
  return junctions;
}

inline bool all_junctions_in_intron(const std::vector<int>& junctions,
                             const std::unordered_set<std::string>& transcript_introns,
                             string chr,
                             const std::unordered_set<std::string>& all_introns) {
  
  // 如果 junctions 为空、长度不足2、或为奇数，说明没有有效剪接点对
  // Java 的逻辑是：这种情况下不进行任何过滤，视为"兼容"
  if (junctions.size() < 2 || junctions.size() % 2 != 0) {
    debug_print(3, "剪接点列表为空、不足或格式不完整，跳过检查（视为兼容）");
    return true; // 👈 关键：返回 true，表示"无冲突"
  }
  
  // 遍历所有有效的剪接点对 (donor, acceptor)
  for (size_t i = 0; i + 1 < junctions.size(); i += 2) {
    int donor = junctions[i];
    int acceptor = junctions[i + 1];
    
    // 构造内含子键（注意：Java 中可能是 "donor-acceptor"）
    std::string junction_key =chr+":"+ std::to_string(donor) + "-" + std::to_string(acceptor);
    
    // 检查该剪接点是否不在当前转录本中，但却存在于全局注释中
    if (transcript_introns.find(junction_key) == transcript_introns.end() &&
        all_introns.find(junction_key) != all_introns.end()) {
      // 发现一个"外来但已知"的剪接点 → 该 read 与当前转录本不兼容
      debug_print(3, "发现外来剪接点: " + junction_key);
      return false;
    }
  }
  return true;
}


inline bool check_read_coverage(int read_start, int read_end,
                                int left_intron_start, int left_intron_end,
                                int left_exon_end, int right_exon_start,
                                int intron_flank_threshold,
                                bool consider_exon_in_intron,
                                const vector<int>& junctions,
                                bool& has_junction_at_this_intron) {
  
  debug_print(3, "进入 check_read_coverage");
  debug_print(3, "read_start: " + to_string(read_start) + ", read_end: " + to_string(read_end));
  debug_print(3, "内含子: " + to_string(left_intron_start) + "-" + to_string(left_intron_end));
  
  has_junction_at_this_intron = false;
  
  // 检查是否有junction正好是这个内含子，如果是的话，再简单判断一下就可以返回false了。
  for (size_t j = 0; j < junctions.size(); j += 2) {
    if (j + 1 >= junctions.size()) break;
    
    if(junctions[j]>left_intron_end) break;
    
    if (junctions[j] == left_intron_start && junctions[j+1] == left_intron_end) {
      has_junction_at_this_intron = true;
      debug_print(3, "找到正好是左侧内含子的junction");
      
      bool covers_enough = (read_start < left_intron_start - intron_flank_threshold &&
                            read_end > left_intron_end + intron_flank_threshold);
      debug_print(3, "read覆盖足够区域: " + string(covers_enough ? "true" : "false"));
      
      
      has_junction_at_this_intron=has_junction_at_this_intron && covers_enough;
      return false;
    }
  }
  
  // 如果有junction的splice site在这个内含子区域内，返回false
  for (size_t j = 0; j < junctions.size(); j++) {
    if (junctions[j] >= left_intron_start && junctions[j] <= left_intron_end) {
      debug_print(3, "有junction在左侧内含子区域内，返回false");
      return false;
    }
    if(junctions[j]>left_intron_end) break;
    
  }

  
  
  // 检查是否有junction跨越这个内含子
  for (size_t j = 0; j < junctions.size(); j += 2) {
    if (j + 1 >= junctions.size()) break;
    
    if (junctions[j] <= left_intron_start && junctions[j+1] >= left_intron_end) {
      debug_print(3, "有junction跨越左侧内含子，返回false");
      return false;
    }
    if(junctions[j]>left_intron_end) break;
    
  }
  
  // 计算anchor_region_len
  int anchor_region_len = 100;
  if (intron_flank_threshold > anchor_region_len) {
    anchor_region_len = intron_flank_threshold + 1;
  }
  
  int intron_length = right_exon_start - left_exon_end - 1;
  if (intron_length <= 0) {
    debug_print(3, "内含子长度<=0，返回false");
    return false;
  }
  
  if (anchor_region_len > intron_length) {
    anchor_region_len = intron_length;
  }
  // if (anchor_region_len > intron_length) {
  //   anchor_region_len = intron_length;
  // }
  
  int current_threshold = intron_flank_threshold;
  if (current_threshold > intron_length) {
    current_threshold = intron_length - 1;
    if (current_threshold < 1) {
      current_threshold = 1;
    }
  }
  
  debug_print(3, "计算参数: anchor_region_len=" + to_string(anchor_region_len) + 
    ", intron_length=" + to_string(intron_length) + 
    ", current_threshold=" + to_string(current_threshold));
  
  bool result = false;
  
  if (consider_exon_in_intron) {//更严格
    bool condition1 = (read_start <= (right_exon_start - current_threshold)) &&
      (read_end >= (right_exon_start - anchor_region_len + current_threshold));//覆盖右侧外显子边界
    
    bool condition2 = (read_start <= (left_exon_end + anchor_region_len - current_threshold)) &&
      (read_end >= (left_exon_end + current_threshold));
    
    debug_print(3, "condition1: " + string(condition1 ? "true" : "false"));
    debug_print(3, "condition2: " + string(condition2 ? "true" : "false"));//覆盖左侧外显子边界
    
    result = (condition1 || condition2);
  } else {
    bool condition1_simple = (read_start <= (right_exon_start - current_threshold));
    bool condition2_simple = (read_end >= (left_exon_end + current_threshold));
    
    debug_print(3, "简化condition1: " + string(condition1_simple ? "true" : "false"));
    debug_print(3, "简化condition2: " + string(condition2_simple ? "true" : "false"));
    
    result = condition1_simple && condition2_simple;
  }
  
  debug_print(3, "覆盖检查结果: " + string(result ? "true" : "false"));
  debug_print(3, "离开 check_read_coverage");
  
  return result;
}

// ====================== 主分析函数 ======================

// 兼容 func2.R 中的 analyze_transcript_java_exact 函数
// [[Rcpp::export(rng=false)]]
List analyze_transcript_java_exact(List transcript_info,
                                   List bam_data,
                                   int intron_flank_threshold = 90,
                                   bool consider_exon_in_intron = true) {
  
  // 初始化调试日志
  init_debug_log("iso_debug.log");
  debug_print(1, "========================================");
  debug_print(1, "开始分析转录本");
  debug_print(1, "========================================");
  
  try {
    // 解析转录本信息（与 func2.R 中的格式匹配）
    string transcript_name = as<string>(transcript_info["name"]);
    string chr = as<string>(transcript_info["chr"]);
    string strand = as<string>(transcript_info["strand"]);
    int tx_start = as<int>(transcript_info["start"]);
    int tx_end = as<int>(transcript_info["end"]);
    List transcript_introns_list = transcript_info["introns"];
    
    debug_print(1, "转录本: " + transcript_name);
    debug_print(1, "染色体: " + chr + ", 链: " + strand);
    debug_print(1, "起始: " + to_string(tx_start) + ", 结束: " + to_string(tx_end));
    debug_print(1, "阈值: " + to_string(intron_flank_threshold));
    
    int n_introns = transcript_introns_list.size();
    debug_print(1, "内含子数量: " + to_string(n_introns));
    
    if (n_introns < 2) {
      debug_print(1, "内含子数量不足，跳过");
      close_debug_log();
      return List::create();
    }
    
    // 构建转录本内含子集合
    unordered_set<string> transcript_introns;
    vector<vector<int>> intron_info;  // [start, end, left_exon_end, right_exon_start]
    
    for (int i = 0; i < n_introns; i++) {
      List intron = transcript_introns_list[i];
      int intron_start = as<int>(intron["start"]);
      int intron_end = as<int>(intron["end"]);
      int left_exon_end = as<int>(intron["left_exon_end"]);
      int right_exon_start = as<int>(intron["right_exon_start"]);
      
      string intron_key = chr + ":" + to_string(intron_start) + "-" + to_string(intron_end);
      transcript_introns.insert(intron_key);
      
      vector<int> info = {intron_start, intron_end, left_exon_end, right_exon_start};
      intron_info.push_back(info);
      
      debug_print(3, "内含子 " + to_string(i) + ": " + intron_key);
    }
    
    // 解析BAM数据 - 使用C++标准数据结构以提高效率
    CharacterVector read_names_r = bam_data["qname"];
    IntegerVector read_starts_r = bam_data["pos"];
    IntegerVector read_ends_r = bam_data["end"];
    IntegerVector mapqs_r = bam_data["mapq"];
    CharacterVector cigars_r = bam_data["cigar"];
    List nh_tags = bam_data["nh"];
    List ji_tags = bam_data["ji"];
    
    // 转换为C++标准容器
    vector<string> read_names;
    vector<int> read_starts;
    vector<int> read_ends;
    vector<int> mapqs;
    vector<string> cigars;
    
    // 转换数据
    int n_reads = read_names_r.size();
    read_names.reserve(n_reads);
    read_starts.reserve(n_reads);
    read_ends.reserve(n_reads);
    mapqs.reserve(n_reads);
    cigars.reserve(n_reads);
    
    for (int i = 0; i < n_reads; i++) {
      read_names.push_back(as<string>(read_names_r[i]));
      read_starts.push_back(read_starts_r[i]);
      read_ends.push_back(read_ends_r[i]);
      mapqs.push_back(mapqs_r[i]);
      cigars.push_back(as<string>(cigars_r[i]));
    }
    
    debug_print(1, "读取到 " + to_string(n_reads) + " 条reads");
    
    if (n_reads == 0) {
      debug_print(1, "没有reads数据");
      close_debug_log();
      return List::create();
    }
    
    // 使用转录本内含子作为所有内含子（简化）
    unordered_set<string> all_introns = transcript_introns;
    
    // 存储结果
    vector<string> result_transcripts;
    vector<string> result_left_introns;
    vector<string> result_right_junctions;
    vector<string> result_strands;
    vector<int> result_cover_counts;
    vector<int> result_junction_counts;
    
    vector<string> result_cover_left_names;
    
    
    debug_print(2, "开始逐内含子分析...");
    
    std::unordered_map<std::string, vector<int> > junction_pos;
    std::unordered_map<std::string, bool> junction_in_intron;
    
    
    // 逐内含子处理
    for (int left_idx = 0; left_idx < n_introns; left_idx++) {
      debug_print(2, "处理左侧内含子 " + to_string(left_idx));
      
      int left_intron_start = intron_info[left_idx][0];
      int left_intron_end = intron_info[left_idx][1];
      int left_exon_end = intron_info[left_idx][2];
      int right_exon_start = intron_info[left_idx][3];
      
      string left_intron_key = chr + ":" + 
        to_string(left_intron_start) + "-" + 
        to_string(left_intron_end);
      
      debug_print(3, "左侧内含子: " + left_intron_key);
      debug_print(3, "外显子边界: " + to_string(left_exon_end) + " - " + to_string(right_exon_start));
      
      // 查找覆盖该内含子的reads
      unordered_set<string> cover_left_names;
      unordered_set<string> left_jc_names;
      
      
      debug_print(3, "查找覆盖左侧内含子的reads");
      int checked_reads = 0;
      int passed_reads = 0;
      
      for (int r = 0; r < n_reads; r++) {
        checked_reads++;
        
        // 基本过滤
        if (mapqs[r] < 1) {
          continue;
        }
        
        // 检查NH标签
        // if (false && !Rf_isNull(nh_tags[r])) {
        //   SEXP nh_val = nh_tags[r];
        //   if (TYPEOF(nh_val) == INTSXP) {
        //     int nh = INTEGER(nh_val)[0];
        //     if (nh > 1) {
        //       continue;
        //     }
        //   }
        // }
        
        int read_start = read_starts[r];
        int read_end = read_ends[r];
        
        // 快速位置检查
        if (read_start > left_intron_end || read_end < left_intron_start) {
          continue;
        }
        
        string read_name = read_names[r];
        string cigar = cigars[r];
        
        // 获取剪接点信息
        vector<int> junctions;
        // if (false && !Rf_isNull(ji_tags[r])) {
        //   SEXP ji_val = ji_tags[r];
        //   if (TYPEOF(ji_val) == INTSXP) {
        //     int* ji_data = INTEGER(ji_val);
        //     int ji_len = LENGTH(ji_val);
        //     for (int j = 0; j < ji_len; j++) {
        //       junctions.push_back(ji_data[j]);
        //     }
        //   }
        // } else 
          if (cigar.find('N') != string::npos) {
          
          if(junction_pos.find(to_string(read_start)+cigar)== junction_pos.end() ){
            junctions = get_splice_junction_pos(cigar, read_start);
            junction_pos.insert({to_string(read_start)+cigar, junctions});
            
          }else{
            junctions=junction_pos.at(to_string(read_start)+cigar);
          }
          
        }
        
        // 检查剪接点是否在已知内含子内
        if (!junctions.empty()) {
          
          bool j_i_in=false;
          if(junction_in_intron.find(to_string(read_start)+cigar)==junction_in_intron.end()){
            j_i_in=all_junctions_in_intron(junctions, transcript_introns, chr, all_introns);
            
            junction_in_intron.insert({to_string(read_start)+cigar, j_i_in});
            
          }else{
            j_i_in=junction_in_intron.at(to_string(read_start)+cigar);
 
          }
          if(j_i_in==false)
            continue;
          // 
          // if (!all_junctions_in_intron(junctions, transcript_introns, chr, all_introns)) {
          //   continue;
          // }
        }
        
        
        
        // 检查覆盖
        bool has_junction_at_this_intron = false;
        bool covers_intron = check_read_coverage(
          read_start, read_end,
          left_intron_start, left_intron_end,
          left_exon_end, right_exon_start,
          intron_flank_threshold,
          consider_exon_in_intron,
          junctions,
          has_junction_at_this_intron
        );//has_junction_at_this_intron按照地址传入
        
        if (has_junction_at_this_intron) {
          left_jc_names.insert(read_name);
          passed_reads++;
        } else if (covers_intron) {
          cover_left_names.insert(read_name);
          passed_reads++;
        }
        
        
      }
      
      debug_print(2, "检查 " + to_string(checked_reads) + " 条reads，通过 " + to_string(passed_reads) + " 条");
      debug_print(2, "覆盖左侧内含子的reads: " + to_string(cover_left_names.size()));
      debug_print(2, "有junction在左侧内含子的reads: " + to_string(left_jc_names.size()));
      
      if (cover_left_names.empty() && left_jc_names.empty()) {
        debug_print(2, "没有reads覆盖此内含子，并且没有reads在这个内含子有junction，跳过");
        continue;
      }
      
      
      
      
      // 在整个转录本区域查找配对的reads
      unordered_map<string, int> interest_iso;
      unordered_map<string, int> read_support_jc_pair;
      unordered_set<string> processed_pairs_cover;
      unordered_set<string> processed_pairs_jc;
      
      int paired_reads_checked = 0;
      int paired_reads_passed = 0;
      
      debug_print(3, "查找配对的reads");
      
      for (int r = 0; r < n_reads; r++) {
        paired_reads_checked++;
        
        if (mapqs[r] < 1) continue;
        
        // if (false&& !Rf_isNull(nh_tags[r])) {
        //   SEXP nh_val = nh_tags[r];
        //   if (TYPEOF(nh_val) == INTSXP) {
        //     int nh = INTEGER(nh_val)[0];
        //     if (nh > 1) continue;
        //   }
        // }
        
        string read_name = read_names[r];
        
        bool is_cover = cover_left_names.count(read_name) > 0;
        bool is_left_jc = left_jc_names.count(read_name) > 0;
        
        if (!is_cover && !is_left_jc) continue;
        
        int read_start = read_starts[r];
        int read_end = read_ends[r];
        
        if (read_start > tx_end || read_end < tx_start) continue;//不在转录本区域内
        string cigar = cigars[r];
        
        // 获取剪接点
        vector<int> junctions;
        // if (false && !Rf_isNull(ji_tags[r])) {
        //   SEXP ji_val = ji_tags[r];
        //   if (TYPEOF(ji_val) == INTSXP) {
        //     int* ji_data = INTEGER(ji_val);
        //     int ji_len = LENGTH(ji_val);
        //     for (int j = 0; j < ji_len; j++) {
        //       junctions.push_back(ji_data[j]);
        //     }
        //   }
        // } else {
          if (cigar.find('N') != string::npos) {
            //junctions = get_splice_junction_pos(cigar, read_start);
            if(junction_pos.find(to_string(read_start)+cigar)== junction_pos.end() ){
              junctions = get_splice_junction_pos(cigar, read_start);
              junction_pos.insert({to_string(read_start)+cigar, junctions});
              
            }else{
              junctions=junction_pos.at(to_string(read_start)+cigar);
            
            }
          }
        //}
        
        if (junctions.size() < 2) continue;//splice site少于2
        
        // if (!all_junctions_in_intron(junctions, transcript_introns, chr, all_introns)) {
        //   continue;
        // }
        
        if (!junctions.empty()) {
          
          bool j_i_in=false;
          if(junction_in_intron.find(to_string(read_start)+cigar)==junction_in_intron.end()){
            j_i_in=all_junctions_in_intron(junctions, transcript_introns, chr, all_introns);
            
            junction_in_intron.insert({to_string(read_start)+cigar, j_i_in});
            
          }else{
            j_i_in=junction_in_intron.at(to_string(read_start)+cigar);
            
          }
          if(j_i_in==false)
            continue;
        }
        
        
        
        for (size_t j = 0; j < junctions.size(); j += 2) {
          if (j + 1 >= junctions.size()) break;
          
          int junc_start = junctions[j];
          int junc_end = junctions[j+1];
          
          string right_jc = chr + ":" + 
            to_string(junc_start) + "-" + 
            to_string(junc_end);
          
          if (right_jc == left_intron_key) continue;
          
          string region_key = transcript_name + "\t" + 
            left_intron_key + "\t" + 
            right_jc + "\t" + 
            strand + "\t" + 
            "false";
          
          debug_print(3, "找到region: " + region_key);
          
          if (is_cover) {
            string pair_key = read_name + "-" + right_jc;
            if (processed_pairs_cover.count(pair_key) == 0) {//每个read或者read pair只能检测到一对intron pair
              processed_pairs_cover.insert(pair_key);
              interest_iso[region_key]++;
              paired_reads_passed++;
            }
          }
          
          if (is_left_jc) {
            string pair_key = read_name + "-" + right_jc;
            if (processed_pairs_jc.count(pair_key) == 0) {
              processed_pairs_jc.insert(pair_key);
              read_support_jc_pair[region_key]++;
              paired_reads_passed++;
            }
          }
        }
      }
      
      debug_print(2, "检查 " + to_string(paired_reads_checked) + " 条配对reads，找到 " + to_string(paired_reads_passed) + " 个有效对");
      debug_print(2, "找到 " + to_string(interest_iso.size()) + " 个不同的region");
      
      
      result_cover_left_names.insert(result_cover_left_names.end(), cover_left_names.begin(), cover_left_names.end());
      result_cover_left_names.insert(result_cover_left_names.end(), left_jc_names.begin(), left_jc_names.end());
      
      // 保存结果
      for (const auto& entry : interest_iso) {
        vector<string> parts;
        string part;
        istringstream tokenStream(entry.first);
        while (getline(tokenStream, part, '\t')) {
          parts.push_back(part);
        }
        
        if (parts.size() < 5) continue;
        
        result_transcripts.push_back(parts[0]);
        result_left_introns.push_back(parts[1]);
        result_right_junctions.push_back(parts[2]);
        result_strands.push_back(parts[3]);
        result_cover_counts.push_back(entry.second);
        
        int jc_count = 0;
        auto it = read_support_jc_pair.find(entry.first);
        if (it != read_support_jc_pair.end()) {
          jc_count = it->second;
        }
        result_junction_counts.push_back(jc_count);
        
        debug_print(3, "保存结果: " + entry.first + 
          " (覆盖=" + to_string(entry.second) + 
          ", 剪接=" + to_string(jc_count) + ")");
      }
    }
    
    
    debug_print(1, "分析完成，找到 " + to_string(result_transcripts.size()) + " 个内含子对");
    debug_print(1, "========================================");
    
    // 返回结果（与R函数期望的格式匹配）
    List result_list_1 = List::create(
      Named("transcript") = result_transcripts,
      Named("left_intron") = result_left_introns,
      Named("right_junction") = result_right_junctions,
      Named("strand") = result_strands,
      Named("cover_count") = result_cover_counts,
      Named("junction_count") = result_junction_counts
    
    );
    
    
    List result_list = List::create(      Named("result_list_1") = result_list_1,
                                          Named("result_cover_left_names")=result_cover_left_names
                                            
    );
    
    close_debug_log();
    
    return result_list;
    
  } catch (const std::exception& e) {
    debug_print(0, "错误: " + string(e.what()));
    close_debug_log();
    Rcpp::stop("分析过程中发生错误: " + string(e.what()));
  } catch (...) {
    debug_print(0, "未知错误");
    close_debug_log();
    Rcpp::stop("分析过程中发生未知错误");
  }
}

// 调试版本，兼容 func.r 中的 analyze_transcript_java_exact_debug 函数
//[[Rcpp::export]]
List analyze_transcript_java_exact_debug(List transcript_info,
                                         List bam_data,
                                         int intron_flank_threshold = 90,
                                         bool consider_exon_in_intron = true) {
  // 使用与正式版本相同的实现，但可能有额外的调试输出
  init_debug_log("iso_debug_detailed.log");
  debug_print(1, "运行调试版本的分析");
  
  List result = analyze_transcript_java_exact(
    transcript_info, bam_data, intron_flank_threshold, consider_exon_in_intron
  );
  
  debug_print(1, "调试版本分析完成");
  close_debug_log();
  
  return result;
}

// 测试覆盖逻辑函数，兼容 func.r
//[[Rcpp::export]]
List test_coverage_logic_debug(int read_start, int read_end,
                               int left_exon_end, int right_exon_start,
                               int intron_flank_threshold,
                               bool consider_exon_in_intron) {
  
  init_debug_log("coverage_test.log");
  debug_print(1, "测试覆盖逻辑");
  
  // 计算anchor_region_len
  int anchor_region_len = 100;
  if (intron_flank_threshold > anchor_region_len) {
    anchor_region_len = intron_flank_threshold + 1;
  }
  
  int intron_length = right_exon_start - left_exon_end - 1;
  if (intron_length <= 0) {
    intron_length = 1;
  }
  
  // if (anchor_region_len > intron_length) {
  //   anchor_region_len = intron_length;
  // }
  
  int current_threshold = intron_flank_threshold;
  if (current_threshold > intron_length) {
    current_threshold = intron_length - 1;
    if (current_threshold < 1) {
      current_threshold = 1;
    }
  }
  
  bool condition1 = false;
  bool condition2 = false;
  bool result = false;
  
  if (consider_exon_in_intron) {
    condition1 = (read_start <= (right_exon_start - current_threshold)) &&
      (read_end >= (right_exon_start - anchor_region_len + current_threshold));
    
    condition2 = (read_start <= (left_exon_end + anchor_region_len - current_threshold)) &&
      (read_end >= (left_exon_end + current_threshold));
    
    result = (condition1 || condition2);
  } else {
    bool condition1_simple = (read_start <= (right_exon_start - current_threshold));
    bool condition2_simple = (read_end >= (left_exon_end + current_threshold));
    
    result = condition1_simple && condition2_simple;
  }
  
  debug_print(1, "测试参数:");
  debug_print(1, "  read_start: " + to_string(read_start));
  debug_print(1, "  read_end: " + to_string(read_end));
  debug_print(1, "  left_exon_end: " + to_string(left_exon_end));
  debug_print(1, "  right_exon_start: " + to_string(right_exon_start));
  debug_print(1, "  intron_flank_threshold: " + to_string(intron_flank_threshold));
  debug_print(1, "  consider_exon_in_intron: " + string(consider_exon_in_intron ? "true" : "false"));
  debug_print(1, "  result: " + string(result ? "true" : "false"));
  
  close_debug_log();
  
  return List::create(
    Named("read_start") = read_start,
    Named("read_end") = read_end,
    Named("left_exon_end") = left_exon_end,
    Named("right_exon_start") = right_exon_start,
    Named("intron_flank_threshold") = intron_flank_threshold,
    Named("consider_exon_in_intron") = consider_exon_in_intron,
    Named("condition1") = condition1,
    Named("condition2") = condition2,
    Named("result") = result,
    Named("anchor_region_len") = anchor_region_len,
    Named("current_threshold") = current_threshold,
    Named("intron_length") = intron_length
  );
}

// 测试函数，兼容 func.r
//[[Rcpp::export]]
List test_debug_output() {
  init_debug_log("test_output.log");
  
  debug_print(0, "调试级别0信息");
  debug_print(1, "调试级别1信息");
  debug_print(2, "调试级别2信息");
  debug_print(3, "调试级别3信息");
  
  vector<int> test_vec = {1, 2, 3, 4, 5};
  
  unordered_set<string> test_set;
  test_set.insert("item1");
  test_set.insert("item2");
  test_set.insert("item3");
  
  debug_print(1, "测试完成");
  close_debug_log();
  
  return List::create(
    Named("message") = "调试输出测试完成",
    Named("status") = "success"
  );
}

// 批量处理函数（可选，用于提高效率）
//[[Rcpp::export]]
List analyze_multiple_transcripts(List transcript_list,
                                  List bam_data_list,
                                  int intron_flank_threshold = 90,
                                  bool consider_exon_in_intron = true) {
  
  init_debug_log("batch_analysis.log");
  debug_print(1, "开始批量分析多个转录本");
  
  int n_transcripts = transcript_list.size();
  debug_print(1, "总共 " + to_string(n_transcripts) + " 个转录本");
  
  // 存储所有结果
  vector<string> all_transcripts;
  vector<string> all_left_introns;
  vector<string> all_right_junctions;
  vector<string> all_strands;
  vector<int> all_cover_counts;
  vector<int> all_junction_counts;
  
  for (int i = 0; i < n_transcripts; i++) {
    debug_print(1, "分析转录本 " + to_string(i+1) + "/" + to_string(n_transcripts));
    
    List transcript = transcript_list[i];
    List bam_data = bam_data_list[i];
    
    List result = analyze_transcript_java_exact(
      transcript, bam_data, intron_flank_threshold, consider_exon_in_intron
    );
    
    CharacterVector tx_names = result["transcript"];
    CharacterVector left_introns = result["left_intron"];
    CharacterVector right_junctions = result["right_junction"];
    CharacterVector strands = result["strand"];
    IntegerVector cover_counts = result["cover_count"];
    IntegerVector junction_counts = result["junction_count"];
    
    for (int j = 0; j < tx_names.size(); j++) {
      all_transcripts.push_back(as<string>(tx_names[j]));
      all_left_introns.push_back(as<string>(left_introns[j]));
      all_right_junctions.push_back(as<string>(right_junctions[j]));
      all_strands.push_back(as<string>(strands[j]));
      all_cover_counts.push_back(cover_counts[j]);
      all_junction_counts.push_back(junction_counts[j]);
    }
    
    if ((i + 1) % 10 == 0) {
      debug_print(1, "进度: " + to_string(i+1) + "/" + to_string(n_transcripts));
    }
  }
  
  debug_print(1, "批量分析完成，找到 " + to_string(all_transcripts.size()) + " 个内含子对");
  close_debug_log();
  
  return List::create(
    Named("transcript") = all_transcripts,
    Named("left_intron") = all_left_introns,
    Named("right_junction") = all_right_junctions,
    Named("strand") = all_strands,
    Named("cover_count") = all_cover_counts,
    Named("junction_count") = all_junction_counts
  );
}