#' (Internal) Lightweight Iterative Parameter Estimation
#'
#' Estimates parameters (g, lam) using a decimation strategy and RSS constraint.
#'
#' @param y Numeric vector. The raw fluorescence data.
#' @param decimate Integer. Downsampling factor.
#' @param max_iter Integer. Maximum iterations for lambda search.
#'
#' @return A list containing optimized \code{g} and \code{lam} for the original data scale.
#' 
#' @keywords internal
estimate_params_iterative <- function(y, decimate = 10, max_iter = 10) {
  
  T_full <- length(y)
  
  n_bins <- floor(T_full / decimate)
  y_small <- colMeans(matrix(y[1:(n_bins * decimate)], nrow = decimate))

  sn_small <- GetSnMAD(y_small)
  g_small <- GetDecay(y_small, sn = sn_small)
  target_rss <- (sn_small^2) * length(y_small)
  
  
  calc_rss_diff <- function(log_lam) {
    lam_curr <- exp(log_lam)
    fit <- oasis_ar1(y_small, g = g_small, lam = lam_curr)
    
    rss <- sum((y_small - fit$c)^2)
    return(rss - target_rss)
  }
  
  opt_res <- tryCatch({
    stats::uniroot(calc_rss_diff, interval = c(-5, 5), tol = 0.1)
  }, error = function(e) {
    # Fallback if root not found (e.g. noise constraint cannot be met)
    list(root = log(2 * sn_small)) 
  })
  
  lam_small_opt <- exp(opt_res$root)
  g_full <- g_small^(1/decimate)
  lam_full <- lam_small_opt * (1 - g_full) / (1 - g_small)
  
  return(list(g = g_full, lam = lam_full, sn = sn_small * sqrt(decimate)))
}