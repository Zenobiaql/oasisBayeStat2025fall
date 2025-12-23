#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;


struct Pool {
  double v; 
  double w; 
  int t;    
  int l;    
};

//' (Internal) C++ Implementation of OASIS AR1 Deconvolution
//' 
//' This is the backend function that performs the heavy lifting for the 
//' deconvolution. It implements the Pool Adjacent Violators Algorithm (PAVA) 
//' for the AR(1) model logic.
//' 
//' @param y Numeric vector. The observed fluorescence trace.
//' @param g Double. The autoregressive decay factor (gamma). 
//'   Range: [0, 1). Closer to 1 implies slower decay.
//' @param lam Double. Sparsity penalty parameter (lambda). 
//'   Larger values result in fewer detected spikes (sparser solution).
//' @param s_min Double. Minimum spike size constraint.
//' 
//' @return A List containing two vectors:
//' \itemize{
//'   \item \code{c}: The inferred denoised calcium traces (same length as y).
//'   \item \code{s}: The inferred spike train (same length as y).
//' }
//' 
//' @section Note:
//' This function is not intended to be called directly by users. 
//' Please use the R wrapper function \code{oasis_ar1()}.
//' 
//' @noRd
// [[Rcpp::export]]
List oasisAR1_cpp(NumericVector y, double g, double lam, double s_min){
  int T = y.size();
  double lg = std::log(g);
   
  std::vector<Pool> P;
  P.reserve(T);
   
  // Initialize first pool
  Pool newpool;
  newpool.v = y[0] - lam * (1 - g);
  newpool.w = 1.0;
  newpool.t = 0;
  newpool.l = 1;
  P.push_back(newpool);
   
  int i = 0; 
   
  // Main loop
  for(int t = 1; t < T; ++t){
    double penalty_factor = (t == T - 1) ? 1.0 : (1 - g);
    newpool.v = y[t] - lam * penalty_factor;
    newpool.w = 1.0;
    newpool.t = t;
    newpool.l = 1;
    P.push_back(newpool);
     
    i++; 
     
    // Backtracking
    while(i > 0){
      // P[i-1] is the previous pool, P[i] is the current (last) pool
      double prev_mean = P[i-1].v / P[i-1].w;
      double curr_mean = P[i].v / P[i].w;
       
      // Check constraint violation
      if(prev_mean * std::exp(lg * P[i-1].l) + s_min > curr_mean){
        i--;
         
        // Merge logic
        double shift_v = std::exp(lg * P[i].l); // P[i] here is the 'prev' one after decrementing i
        double shift_w = std::exp(lg * 2 * P[i].l);
           
        P[i].v += P[i+1].v * shift_v;
        P[i].w += P[i+1].w * shift_w;
        P[i].l += P[i+1].l;
           
        P.pop_back();
      }
      else
          break;
    }
  }
   
  // Construct c
  NumericVector c(T);
  for(size_t j = 0; j < P.size(); ++j){
    double tmp = P[j].v / P[j].w;
     
    if (tmp < 0) {
      tmp = 0;
    }
     
  for(int k = 0; k < P[j].l; ++k){
      c[k + P[j].t] = tmp;
      tmp *= g;
    }
  }
   
  // Construct s
  NumericVector s(T);
  s[0] = 0;
  for(int k = 1; k < T; ++k){
    s[k] = c[k] - g * c[k-1];
  }
   
  return List::create(_["c"] = c, _["s"] = s);
}


//' (Internal) Calculate Median of a Vector
//'
//' A helper function that computes the median of a vector using 
//' std::nth_element for efficient O(N) performance.
//'
//' @param v A std::vector<double> containing the data.
//' @return The median value as a double.
//' @section Warning:
//' This function modifies the input vector \code{v} in-place (partial sorting).
//' Do not use if the order of \code{v} must be preserved.
//' 
//' @noRd
double get_median_cpp(std::vector<double>& v){
  size_t n = v.size();
  if (n == 0) return 0.0;
  size_t target = n / 2;
  std::nth_element(v.begin(), v.begin() + target, v.end());
   
  return v[target];
}


//' (Internal) C++ Implementation of Parameter Estimation
//' 
//' This is the backend function that estimates the noise standard deviation 
//' (sn) and the AR(1) decay constant (g) from the raw fluorescence data.
//' 
//' @param y Numeric vector. The observed fluorescence trace.
//' 
//' @details 
//' The noise standard deviation (\code{sn}) is estimated using the Median 
//' Absolute Deviation (MAD) of the first differences of the signal, 
//' which is robust to large spikes. This assumes the underlying noise is 
//' Gaussian and robustly scales the MAD factor (1.4826) accordingly.
//' 
//' The decay constant (\code{g}) is estimated using the Lag-1 Autocorrelation 
//' of the signal.
//' 
//' @return A List containing two named elements:
//' \itemize{
//'   \item \code{g}: The estimated autoregressive decay factor.
//'   \item \code{sn}: The estimated noise standard deviation.
//' }
//' 
//' @section Note:
//' This function is not intended to be called directly by users. 
//' Please use the R wrapper function \code{estimate_parameters()}.
//' 
//' @noRd
// [[Rcpp::export]]
List estimate_parameters_cpp(NumericVector y){
  int T = y.size();
   
  if (T < 2){
    stop("Time series too short to estimate parameters.");
  }
   
  std::vector<double> abs_diffs;
  abs_diffs.reserve(T - 1);
   
  for(int t = 0; t < T - 1; ++t){
    abs_diffs.push_back(std::abs(y[t+1] - y[t]));
  }
   
  double med_val = get_median_cpp(abs_diffs);
  double estimated_sn = med_val / 0.6745 / std::sqrt(2.0);
   
  double sum_y = 0.0;
  for(double val : y) sum_y += val;
  double mean_y = sum_y / T;
   
  double num = 0.0; 
  double den = 0.0; 
   
  for (int t = 0; t < T - 1; ++t) {
    double centered_curr = y[t] - mean_y;
    double centered_next = y[t+1] - mean_y;
     
    num += centered_curr * centered_next;
    den += centered_curr * centered_curr;
  }
  den += std::pow(y[T-1] - mean_y, 2);
   
  double estimated_g = num / den;
  
  if(estimated_g < 0)
    estimated_g = 0;
  
  if(estimated_g >= 1)
    estimated_g = 0.99;
   
  return List::create(
    Named("g") = estimated_g,
    Named("sn") = estimated_sn
  );
}


//' (Internal) Estimate noise using PSD (C++ version)
//' 
//' @param y Input vector
//' @param range_ff Frequency range [low, high], e.g. c(0.25, 0.5)
//' @param method "mean" or "median"
//' 
//' @noRd
// [[Rcpp::export]]
double GetSn_cpp(arma::vec y, NumericVector range_ff, std::string method) {
  int n = y.n_elem;
  
  arma::cx_vec fy = arma::fft(y);
  
  int n_half = n / 2;
  
  arma::vec Pxx = arma::abs(fy.head(n_half)); 

  Pxx = arma::square(Pxx) / n; 
   
  int idx_start = std::ceil(range_ff[0] * n);
  int idx_end   = std::floor(range_ff[1] * n);
   
  if (idx_start < 0) 
    idx_start = 0;
  if (idx_end > n_half) 
    idx_end = n_half;
   
  std::vector<double> Pxx_ind;
   
  if(idx_end > idx_start)
    Pxx_ind.reserve(idx_end - idx_start);
   
  for(int i = idx_start; i < idx_end; ++i)
    Pxx_ind.push_back(Pxx[i]);
   
  if(Pxx_ind.empty()) 
    return 0.0;
   
  double result = 0.0;
   
  if(method == "mean"){
    double sum = 0.0;
    for (double val : Pxx_ind) sum += val;
    result = std::sqrt( (sum / Pxx_ind.size()) / 2.0 );
  }
  
  else if(method == "median"){
    size_t sz = Pxx_ind.size();
    size_t target = sz / 2;
    std::nth_element(Pxx_ind.begin(), Pxx_ind.begin() + target, Pxx_ind.end());
    double med = Pxx_ind[target];
     
    if(sz % 2 == 0){
      std::sort(Pxx_ind.begin(), Pxx_ind.end());
      med = (Pxx_ind[sz/2 - 1] + Pxx_ind[sz/2]) / 2.0;
    }
     
    result = std::sqrt(med / 2.0);
  }
   
  return result;
}


//' (Internal) Estimate AR(1) Decay Factor with Noise Correction
//' 
//' @param y Numeric vector.
//' @param sn Double.
//' 
//' @noRd
// [[Rcpp::export]]
double estimate_g_corrected_cpp(NumericVector y, double sn){
  int T = y.size();
  if (T < 2) return 0.0;
   
  double sum_y = 0.0;
  for(double val : y) sum_y += val;
  double mean_y = sum_y / T;
   
  double num = 0.0; 
  double den = 0.0; 
   
  for (int t = 0; t < T - 1; ++t) {
    double centered_curr = y[t] - mean_y;
    double centered_next = y[t+1] - mean_y;
     
    num += centered_curr * centered_next;
    den += centered_curr * centered_curr;
  }
  den += std::pow(y[T-1] - mean_y, 2);
   
  double noise_variance_term = T * sn * sn;
  double corrected_den = den - noise_variance_term;
   
  if (corrected_den <= 0.0001) {
    corrected_den = den; 
  }
   
  double g = num / corrected_den;
   
  if (g < 0) g = 0;
  if (g >= 1) g = 0.99; 
   
  return g;
}