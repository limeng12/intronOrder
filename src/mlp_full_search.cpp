// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <vector>
#include <set>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <climits>
#include <limits>

using namespace Rcpp;

// Calculate probability for order
inline float calp2_c(float **t_read_count_mat_li, std::vector<int>& order_arr, int t_dim) {
  float p_sum = 0;
  
  for(int i = 0; i < t_dim - 1; i++) {
    for(int j = i + 1; j < t_dim; j++) {
      p_sum = p_sum + t_read_count_mat_li[order_arr[i]][order_arr[j]];
    }
  }
  
  return p_sum;
}

// [[Rcpp::export]]
List find_best_order_full2_c(NumericMatrix t_read_count_mat, double t_alpha_v = 0.1) {
  int dim = t_read_count_mat.nrow();
  
  // Allocate and initialize matrix
  float **read_count_mat_li = new float *[dim];
  for(int i = 0; i < dim; i++)
    read_count_mat_li[i] = new float[dim];
  
  for(int i = 0; i < dim; i++) {
    for(int j = 0; j < dim; j++) {
      if(i == j) {
        read_count_mat_li[i][j] = read_count_mat_li[j][i] = 0;
        continue;
      }
      
      if((t_read_count_mat(i, j) + t_read_count_mat(j, i)) == 0) {
        read_count_mat_li[i][j] = read_count_mat_li[j][i] = log(0.5);
        continue;
      }
      
      if(i > j) {
        float p = (t_read_count_mat(i, j) + t_alpha_v) / 
                  (t_read_count_mat(i, j) + t_read_count_mat(j, i) + 2 * t_alpha_v);
        read_count_mat_li[i][j] = log(p);
        read_count_mat_li[j][i] = log(1 - p);
      }
    }
  }
  
  std::vector<float> li;
  li.reserve(pow(2, dim));
  
  std::vector<int> init_order(dim);
  std::vector<int> best_order(dim);
  float max_li = -FLT_MAX;
  
  for(int i = 0; i < dim; i++) {
    init_order[i] = i;
  }
  
  float in_order_li = calp2_c(read_count_mat_li, init_order, dim);
  std::sort(init_order.begin(), init_order.end());
  
  // Exhaustive search through all permutations
  int number_of_maximum_order = 0;
  double sum_of_less_than_in_order = 0;
  double entropy = 0;
  double prob_sum = 0;
  
  do {
    float tmp_li = calp2_c(read_count_mat_li, init_order, dim);
    li.push_back(tmp_li);
    
    if(tmp_li > max_li) {
      max_li = tmp_li;
      std::copy(init_order.begin(), init_order.end(), best_order.begin());
    }
  } while(std::next_permutation(init_order.begin(), init_order.end()));
  
  // Calculate statistics
  for(auto li_d = li.begin(); li_d != li.end(); li_d++) {
    prob_sum += std::exp(*li_d);
  }
  
  for(auto li_d = li.begin(); li_d != li.end(); li_d++) {
    if(std::abs(*li_d - max_li) <= std::numeric_limits<double>::epsilon()) {
      number_of_maximum_order++;
    }
    
    double prob_one = std::exp(*li_d) / prob_sum;
    double entropy_one = prob_one * log2(prob_one);
    entropy += entropy_one;
    
    if(*li_d <= in_order_li) {
      sum_of_less_than_in_order += prob_one;
    }
  }
  
  double permut_p = sum_of_less_than_in_order;
  
  for(int i = 0; i < dim; i++) {
    best_order[i] = best_order[i] + 1;
  }
  
  // Clean up memory
  for(int i = 0; i < dim; i++)
    delete[] read_count_mat_li[i];
  delete[] read_count_mat_li;
  
  List L = List::create(
    Named("best_order") = best_order,
    _["entropy"] = -1 * entropy,
    _["best_score"] = max_li,
    _["number_of_maximum_order"] = number_of_maximum_order,
    _["permut_p"] = permut_p
  );
  
  return L;
}