#include <Rcpp.h>
using namespace Rcpp;

//' (Internal) C++ Implementation of Calcium Data Simulation
//'
//' This is the backend function that performs the heavy lifting for the 
//' simulation. It implements the AR process and Poisson spiking logic.
//'
//' @param T Integer. Total number of time frames.
//' @param N Integer. Number of neurons.
//' @param g Numeric vector. AR coefficients (length 1 for AR1, length 2 for AR2).
//' @param sn Double. Standard deviation of Gaussian noise.
//' @param b Double. Baseline fluorescence.
//' @param firerate Double. Target firing rate in Hz.
//' @param framerate Integer. Sampling rate in Hz.
//' 
//' @return A List containing three matrices:
//' \itemize{
//'   \item \code{Y}: Noisy fluorescence data (N x T).
//'   \item \code{c}: True calcium traces (N x T).
//'   \item \code{s}: Binary spike trains (N x T).
//' }
//' 
//' @section Note:
//' This function is not intended to be called directly by users. 
//' Please use the R wrapper function \code{generate_data()}.
//'
//' @noRd
// [[Rcpp::export]]
List generate_data_impl(int T, int N, 
                   NumericVector g, 
                   double sn, double b, 
                   double firerate, 
                   int framerate){
  
  // 1. initialize output: Y, fluorescence data; 
  //                c, ground truth calcium trace; 
  //                s, ground truth spike train.
  
  NumericMatrix Y(N, T);
  NumericMatrix c(N, T);
  NumericMatrix s(N, T);
  
  double threshold = firerate / (double)framerate;
  
  // 2. generate ground truth spike
  
  for(int t = 0; t < T; t++){
    for(int n = 0; n < N; n++){
      if (R::runif(0, 1) < threshold)
        s(n ,t) = 1;
      else
        s(n, t) = 0;
    }
  }
  
  // 3. generate ground truth calcium trace
  
  int p = g.size();
  int start_t = (p == 2) ? 2 : 1;
  
  
  for(int n = 0; n < N; n++){
      if(p == 2){
        c(n, 0) = s(n, 0);
        c(n, 1) = s(n, 1) + g[0] * c(n, 0);
      }
      else if(p == 1)
        c(n, 0) = s(n, 0);
    }
    
  for(int t = start_t; t < T; t++){
    for(int n = 0; n < N; n++){
      double past_val = 0.0;
      
      if(p == 2){
        past_val = g[0] * c(n, t-1)  + g[1] * c(n, t-2);
      }
      else if(p == 1){
        past_val = g[0] * c(n, t-1);
      }
      
      c(n, t) = past_val + s(n, t);
    }
  }
  
  // 4. generate fluorescence data
  
  for(int t = 0; t < T; t++){
    for(int n = 0; n < N; n++){
      double noise = R::rnorm(0, sn);
      Y(n, t) = b + c(n, t) + noise;
    }
  }
  
  return List::create(
    _["Y"] = Y,
    _["c"] = c,
    _["s"] = s
  );
}