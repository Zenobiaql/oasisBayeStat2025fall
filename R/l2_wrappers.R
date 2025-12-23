# -------------------------------------------------------------------------
# L2 Component: Noise Estimation
# -------------------------------------------------------------------------

#' Estimate Noise Standard Deviation (PSD Method)
#'
#' Estimates the noise standard deviation from the high-frequency component
#' of the Power Spectral Density (PSD).
#'
#' @param y Numeric vector. The observed fluorescence trace.
#' @param range_ff Numeric vector of length 2. Normalized frequency range [0, 0.5].
#' @param method Character string. Method to average the spectrum: "mean" or "median".
#' 
#' @examples
#' 
#' # 1. Simulate data with known noise (sn = 0.5)
#' dat <- generate_data(N = 1, sn = 0.5, seed = 123)
#' y <- dat$Y
#' 
#' # 2. Estimate noise using the PSD method
#' sn_est <- GetSn(y)
#' 
#' # 3. Check accuracy
#' cat(sprintf("True sn: 0.5, Estimated: %.4f\n", sn_est))
#' 
#' # 4. Compare method variants
#' sn_mean <- GetSn(y, method = "mean")
#' sn_median <- GetSn(y, method = "median")
#' 
#' cat(sprintf("sn_mean: %.4f\n, sn_median: %.4f\n", sn_mean, sn_median))
#'
#' @return A numeric scalar (estimated sn).
#' @export
GetSn <- function(y, range_ff = c(0.25, 0.5), method = c("mean", "median")) {

  if (!is.numeric(y)) {
    stop("Input 'y' must be a numeric vector.")
  }
  
  if (length(range_ff) != 2 || any(range_ff < 0) || any(range_ff > 0.5)) {
    stop("'range_ff' must be a vector of length 2 within [0, 0.5].")
  }
  
  method <- match.arg(method)
  
  return(GetSn_cpp(y, range_ff, method))
}


# -------------------------------------------------------------------------
# L2 Component: Decay Factor Estimation
# -------------------------------------------------------------------------

#' Estimate AR(1) Decay Factor (Bias Corrected)
#'
#' Estimates the decay factor 'g' using autocorrelation, corrected for
#' noise variance to avoid underestimation.
#'
#' @param y Numeric vector. The observed fluorescence trace.
#' @param sn Numeric scalar. The noise standard deviation.
#'   If \code{NULL} (default), 
#'   it will be estimated automatically using \code{GetSn(y)}.
#'   
#' @examples
#' # 1. Simulate data with known decay (g = 0.95)
#' dat <- generate_data(N = 1, g = 0.95, seed = 123)
#' y <- dat$Y
#' 
#' # 2. Estimate decay (automatically estimates sn internally)
#' g_auto <- GetDecay(y)
#' 
#' # 3. Estimate decay with explicit noise (recommended if sn is known)
#' #    This avoids re-calculating sn and allows for manual tuning
#' sn_known <- 0.3
#' g_explicit <- GetDecay(y, sn = sn_known)
#' 
#' cat(sprintf("True g: 0.95, Est (Auto): %.4f, Est (Explicit): %.4f\n", 
#'             g_auto, g_explicit))
#'
#' @return A numeric scalar (estimated g).
#' @export
GetDecay <- function(y, sn = NULL) {
  if (!is.numeric(y)) {
    stop("Input 'y' must be a numeric vector.")
  }
  
  if (is.null(sn)) {
    sn <- GetSn(y)
  }
  
  return(estimate_g_corrected_cpp(y, sn))
}


# -------------------------------------------------------------------------
# L2 Component: Core Deconvolution Wrapper
# -------------------------------------------------------------------------

#' OASIS AR(1) Deconvolution (Core Algorithm)
#'
#' A low-level wrapper for the OASIS deconvolution algorithm.
#'
#' @param y Numeric vector. Fluorescence trace.
#' @param g Numeric scalar. Decay factor. Must be in [0, 1).
#' @param lam Numeric scalar. Sparsity penalty. 
#' @param s_min Numeric scalar. Minimum spike size.
#' 
#' @examples
#' # 1. Simulate data
#' dat <- generate_data(N = 1, g = 0.95, sn = 0.3, seed = 42)
#' y <- dat$Y
#' 
#' # 2. Step-by-step Parameter Estimation (Manual Pipeline)
#' #    A. Estimate noise
#' sn_est <- GetSn(y)
#' #    B. Estimate decay (using estimated noise)
#' g_est <- GetDecay(y, sn = sn_est)
#' #    C. Set sparsity penalty (heuristic)
#' lam_est <- 2 * sn_est
#' 
#' # 3. Run Core Deconvolution
#' #    Note: 'g' is mandatory here!
#' fit <- oasis_ar1(y, g = g_est, lam = lam_est)
#' 
#' # 4. Inspect results (fit is a raw list here, not an S3 object)
#' 
#' cat("Number of spikes detected:", sum(fit$s > 0), "\n")
#' cat("Max calcium value:", max(fit$c), "\n")
#' 
#' # 5. Compare with Ground Truth
#' cor_c <- cor(fit$c, dat$c)
#' cat(sprintf("Correlation with truth: %.4f\n", cor_c))
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{c}: Denoised calcium trace.
#'   \item \code{s}: Deconvolved spike train.
#' }
#' @export
oasis_ar1 <- function(y, g, lam = 0, s_min = 0) {
  
  if (missing(g) || is.null(g)) {
    stop("Parameter 'g' is required for the core algorithm.")
  }
  
  
  if (g < 0 || g >= 1) {
    stop("Parameter 'g' must be >= 0 and < 1.")
  }
  
  oasisAR1_cpp(y, g, lam, s_min)
}