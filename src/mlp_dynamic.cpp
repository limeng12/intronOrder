// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cstdint>
#include <map>
#include <vector>
#include <set>
#include <bitset>
#include <math.h>
#include <unordered_map>
#include <algorithm>

using namespace Rcpp;

// Bit manipulation macros
#define BIT_SET(a,b) ((a) |= (1ULL<<(b)))
#define BIT_CLEAR(a,b) ((a) &= ~(1ULL<<(b)))
#define BIT_FLIP(a,b) ((a) ^= (1ULL<<(b)))
#define BIT_CHECK(a,b) (!!((a) & (1ULL<<(b))))

// Count number of set bits
inline int numberOfSetBits(uint32_t i) {
  i = i - ((i >> 1) & 0x55555555);
  i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
  return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

// Compare function for sets
bool compare_two_set(const uint32_t lhs, const uint32_t rhs) {
  if(numberOfSetBits(lhs) != numberOfSetBits(rhs)) {
    return numberOfSetBits(lhs) < numberOfSetBits(rhs);
  }
  return lhs < rhs;
}

// Calculate probability for a given order
double calp2_c(float **t_read_count_mat_li, std::vector<int>& order_arr) {
  double p_sum = 0;
  int dim = order_arr.size();
  
  for(int i = 0; i < dim - 1; i++) {
    for(int j = i + 1; j < dim; j++) {
      p_sum = p_sum + t_read_count_mat_li[order_arr[i]][order_arr[j]];
    }
  }
  
  return p_sum;
}

// Calculate permutation p-value
double calp2_p_v(float **t_read_count_mat_li, int t_dim, int sim_times = 1000000) {
  std::vector<int> init_order_p_v(t_dim);
  for(int i = 0; i < t_dim; i++)
    init_order_p_v[i] = i;
  
  std::sort(init_order_p_v.begin(), init_order_p_v.end());
  float in_order_li = calp2_c(t_read_count_mat_li, init_order_p_v);
  
  double sum_of_less_than_in_order = 0;
  double sim_sum = 0;
  
  for(int i = 0; i < sim_times; i++) {
    std::random_shuffle(init_order_p_v.begin(), init_order_p_v.end());
    double tmp_li = calp2_c(t_read_count_mat_li, init_order_p_v);
    sim_sum += std::exp(tmp_li);
    
    if(tmp_li <= in_order_li) {
      sum_of_less_than_in_order += std::exp(tmp_li);
    }
  }
  
  double permut_p = sum_of_less_than_in_order / sim_sum;
  return permut_p;
}

// Get log sum for a sink node
inline float get_log_sum_c_bit(int t_sink, uint32_t t_parents, float **mat_li, int t_dim) {
  float t_log_sum = 0;
  
  for(std::size_t i = 0; i < t_dim; ++i) {
    if(BIT_CHECK(t_parents, i) && (mat_li[i][t_sink] != 0)) {
      t_log_sum = t_log_sum + (mat_li[i][t_sink]);
    }
  }
  
  return t_log_sum;
}

// Generate all combinations
void combinations(uint32_t set, int at, int r, int n, std::vector<uint32_t>& subsets) {
  int elementsLeftToPick = n - at;
  if(elementsLeftToPick < r) {
    return;
  }
  
  if(r == 0) {
    subsets.push_back(set);
  } else {
    for(int i = at; i < n; i++) {
      set ^= (1 << i);
      combinations(set, i + 1, r - 1, n, subsets);
      set ^= (1 << i);
    }
  }
}

// Get all combinations
std::vector<uint32_t> get_all_combinations(int n) {
  std::vector<uint32_t> subsets;
  subsets.reserve(pow(2, n));
  
  for(int r = 1; r <= n; r++) {
    combinations(0, 0, r, n, subsets);
  }
  return subsets;
}

// Main dynamic programming algorithm
std::vector<double> find_opti_dynam_c_bit(float **read_count_mat_li, double t_alpha_v, int dim) {
  std::vector<uint32_t> v_set_power_vec = get_all_combinations(dim);
  std::sort(v_set_power_vec.begin(), v_set_power_vec.end(), compare_two_set);
  v_set_power_vec.erase(std::unique(v_set_power_vec.begin(), v_set_power_vec.end()), v_set_power_vec.end());
  
  std::vector<float> scores(pow(2, dim), -FLT_MAX);
  std::vector<unsigned char> sinks(pow(2, dim), -1);
  
  for(auto it = v_set_power_vec.begin(); it != v_set_power_vec.end(); ++it) {
    uint32_t one_sub = *it;
    uint32_t number_bit = numberOfSetBits(one_sub);
    
    if(number_bit == 0) {
      continue;
    }
    
    if(number_bit > 2) {
      for(std::size_t i = 0; i < dim; ++i) {
        if(!BIT_CHECK(one_sub, i)) {
          continue;
        }
        
        BIT_CLEAR(one_sub, i);
        float skore = scores[one_sub];
        skore = skore + get_log_sum_c_bit(i, one_sub, read_count_mat_li, dim);
        BIT_SET(one_sub, i);
        
        if((skore > scores[one_sub]) || (sinks[one_sub] == -1)) {
          scores[one_sub] = skore;
          sinks[one_sub] = i;
        }
      }
    } else if(number_bit == 2) {
      int index_one = 0;
      int index_two = 0;
      bool a = true;
      
      for(std::size_t i = 0; i < dim; ++i) {
        if(a && BIT_CHECK(one_sub, i)) {
          index_one = i;
          a = false;
        }
        
        if((!a) && BIT_CHECK(one_sub, i)) {
          index_two = i;
        }
      }
      
      if(read_count_mat_li[index_one][index_two] > read_count_mat_li[index_two][index_one]) {
        scores[one_sub] = read_count_mat_li[index_one][index_two];
        sinks[one_sub] = index_two;
      } else {
        scores[one_sub] = read_count_mat_li[index_two][index_one];
        sinks[one_sub] = index_one;
      }
    } else {
      scores[one_sub] = 0;
      for(std::size_t i = 0; i < dim; ++i) {
        if(BIT_CHECK(one_sub, i)) {
          sinks[one_sub] = i;
        }
      }
    }
  }
  
  std::vector<int> ord(dim + 1);
  uint32_t left = (pow(2, dim) - 1);
  uint32_t full_set = (pow(2, dim) - 1);
  float best_score = scores[full_set];
  
  for(int i = dim - 1; i >= 0; i--) {
    ord[i] = sinks[left];
    BIT_CLEAR(left, ord[i]);
  }
  
  ord[dim] = best_score;
  
  std::vector<double> ord_double;
  for(int i = 0; i < ord.size(); i++) {
    ord_double.push_back(ord[i]);
  }
  
  return ord_double;
}

// Rcpp export function
// [[Rcpp::export]]
List find_opti_dynam_r_cpp_bit(NumericMatrix t_read_count_mat, double t_alpha_v = 0.1) {
  int dim = t_read_count_mat.nrow();
  
  // Allocate memory
  double **m_read_count_mat = new double *[dim];
  for(int i = 0; i < dim; i++)
    m_read_count_mat[i] = new double[dim];
  
  for(int i = 0; i < dim; i++) {
    for(int j = 0; j < dim; j++) {
      m_read_count_mat[i][j] = t_read_count_mat(i, j);
    }
  }
  
  float **m_read_count_mat_li = new float *[dim];
  for(int i = 0; i < dim; i++)
    m_read_count_mat_li[i] = new float[dim];
  
  for(int i = 0; i < dim; i++) {
    for(int j = 0; j < dim; j++) {
      if(i == j) {
        m_read_count_mat_li[i][j] = m_read_count_mat_li[j][i] = 0;
        continue;
      }
      
      if((m_read_count_mat[i][j] + m_read_count_mat[j][i]) == 0) {
        m_read_count_mat_li[i][j] = m_read_count_mat_li[j][i] = log(0.5);
        continue;
      }
      
      if(i > j) {
        float p = (m_read_count_mat[i][j] + t_alpha_v) / 
                  (m_read_count_mat[i][j] + m_read_count_mat[j][i] + 2 * t_alpha_v);
        m_read_count_mat_li[i][j] = log(p);
        m_read_count_mat_li[j][i] = log(1 - p);
      }
    }
  }
  
  std::vector<double> ord = find_opti_dynam_c_bit(m_read_count_mat_li, t_alpha_v, dim);
  NumericVector out(dim);
  
  for(int i = 0; i < dim; i++) {
    out[i] = ord[i] + 1;
  }
  
  double permut_p = calp2_p_v(m_read_count_mat_li, dim);
  
  // Clean up memory
  for(int i = 0; i < dim; i++) {
    delete[] m_read_count_mat[i];
    delete[] m_read_count_mat_li[i];
  }
  delete[] m_read_count_mat;
  delete[] m_read_count_mat_li;
  
  List L = List::create(
    Named("best_order") = out,
    _["entropy"] = NA_REAL,
    _["best_score"] = ord[dim],
    _["number_of_maximum_order"] = NA_INTEGER,
    _["permut_p"] = permut_p
  );
  
  return L;
}