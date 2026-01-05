// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <vector>
#include <set>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <climits>

using namespace Rcpp;

// Calculate probability for order
double calp2_c(double **t_read_count_mat_li, std::vector<int>& order_arr) {
  double p_sum = 0;
  int dim = order_arr.size();
  
  for(int i = 0; i < dim - 1; i++) {
    for(int j = i + 1; j < dim; j++) {
      p_sum = p_sum + t_read_count_mat_li[order_arr[i]][order_arr[j]];
    }
  }
  
  return p_sum;
}

// Hill climbing iteration
std::vector<int> hill_iter_c(double **t_read_count_mat_li, std::vector<int> init_order, std::set<int> full_order) {
  double best_score = -DBL_MAX;
  std::vector<int> best_order;
  
  std::set<int> candidates = full_order;
  for(auto it_v = init_order.begin(); it_v != init_order.end(); it_v++) {
    candidates.erase(*it_v);
  }
  
  for(auto it = candidates.begin(); it != candidates.end(); it++) {
    for(int j = 0; j < (init_order.size() + 1); j++) {
      std::vector<int> tmp_order = init_order;
      tmp_order.insert(tmp_order.begin() + j, *it);
      
      double tmp_score = calp2_c(t_read_count_mat_li, tmp_order);
      
      if(tmp_score > best_score) {
        best_score = tmp_score;
        best_order = tmp_order;
      }
    }
  }
  
  return best_order;
}

// [[Rcpp::export]]
List hill_c(NumericMatrix t_read_count_mat, double t_alpha_v = 0.1) {
  int dim = t_read_count_mat.nrow();
  
  // Allocate and initialize matrix
  double **read_count_mat_li = new double *[dim];
  for(int i = 0; i < dim; i++)
    read_count_mat_li[i] = new double[dim];
  
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
        double p = (t_read_count_mat(i, j) + t_alpha_v) / 
                   (t_read_count_mat(i, j) + t_read_count_mat(j, i) + 2 * t_alpha_v);
        read_count_mat_li[i][j] = log(p);
        read_count_mat_li[j][i] = log(1 - p);
      }
    }
  }
  
  std::vector<int> init_order(2);
  std::set<int> full_order;
  double max_li = -DBL_MAX;
  
  for(int i = 0; i < dim; i++)
    full_order.insert(i);
  
  // Find initial best pair
  for(int i = 0; i < dim; i++) {
    for(int j = 0; j < dim; j++) {
      if((read_count_mat_li[i][j] > max_li) && (i != j)) {
        max_li = read_count_mat_li[i][j];
        init_order[0] = i;
        init_order[1] = j;
      }
    }
  }
  
  // Hill climbing iterations
  while(init_order.size() < dim) {
    init_order = hill_iter_c(read_count_mat_li, init_order, full_order);
  }
  
  max_li = calp2_c(read_count_mat_li, init_order);
  
  for(int i = 0; i < dim; i++)
    init_order[i] = init_order[i] + 1;
  
  // Clean up memory
  for(int i = 0; i < dim; i++)
    delete[] read_count_mat_li[i];
  delete[] read_count_mat_li;
  
  List L = List::create(
    Named("best_order") = init_order,
    _["entropy"] = NA_REAL,
    _["best_score"] = max_li,
    _["number_of_maximum_order"] = NA_INTEGER,
    _["permut_p"] = NA_REAL
  );
  
  return L;
}