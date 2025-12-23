#' Generate Synthetic Calcium Imaging Data
#'
#' Simulates spike trains, underlying calcium concentrations, and noisy fluorescence
#' traces based on a homogeneous Poisson process and an Auto-Regressive (AR) model.
#'
#' @param nf Integer. The total duration of the simulation (number of frames). 
#' @param N Integer. The number of neurons to simulate.
#' @param g Numeric vector. The parameters of the AR process. 
#'   If length is 1, an AR(1) process is used.
#'   If length is 2, an AR(2) process is used.
#' @param sn Double. The standard deviation of the Gaussian observation noise.
#' @param b Double. The baseline fluorescence value. 
#' @param firerate Double. The target firing rate (Hz). 
#' @param framerate Integer. The recording frame rate (Hz). 
#' @param seed Integer. Optional. Random seed for reproducibility. 
#'
#' @return A list containing three matrices (each with dimension N x T):
#' \describe{
#'   \item{Y}{The simulated noisy fluorescence data.}
#'   \item{c}{The ground truth calcium concentration traces.}
#'   \item{s}{The binary ground truth spike trains.}
#' }
#' 
#' @examples
#' # 1. Simulate a single neuron (returns vectors)
#' # -------------------------------------------
#' sim_one <- generate_data(N = 1, nf = 100, seed = 123)
#' 
#' # Check the summary
#' print(sim_one)
#' 
#' # Inspect the first 5 values of the noisy fluorescence
#' head(sim_one$Y, 5)
#' 
#' 
#' # 2. Simulate multiple neurons (returns matrices)
#' # ---------------------------------------------
#' sim_multi <- generate_data(N = 5, nf = 100, firerate = 1.0, seed = 42)
#' 
#' # Check dimensions
#' print(sim_multi)
#' 
#' # Inspect the raw data structure (5 neurons x 100 frames)
#' dim(sim_multi$Y)
#'
#' @export
#' 
generate_data <- function(nf = 3000, 
                          N = 20, 
                          g = 0.95, 
                          sn = 0.3, 
                          b = 0, 
                          firerate = 0.5, 
                          framerate = 30,
                          seed = NULL) {
  
  # --- 1. Set seeds and input validation ---
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  N <- as.integer(N)
  nf <- as.integer(nf)
  
  if (length(g) > 2 || length(g) < 1) {
    stop("Parameter 'g' must be a vector of length 1 (AR1) or 2 (AR2).")
  }
  
  if (firerate <= 0) {
    warning("Firing rate should be positive.")
  }
  
  # --- 2. Invoke cpp backend ---
  
  raw_res <- generate_data_impl(
    T = as.integer(nf),
    N = as.integer(N),
    g = as.numeric(g),
    sn = as.numeric(sn),
    b = as.numeric(b),
    firerate = as.numeric(firerate),
    framerate = as.integer(framerate)
  )
  
  # --- 3. Post-process output ---
  
  if (N == 1) {
    raw_res$Y <- as.vector(raw_res$Y)
    raw_res$c <- as.vector(raw_res$c)
    raw_res$s <- as.vector(raw_res$s)
  }
  
  attr(raw_res, "params") <- list(
    nf = nf,
    N = N,
    g = g,
    sn = sn,
    b = b,
    firerate = firerate,
    framerate = framerate
  )
  
  class(raw_res) <- c("calcium_sim", "list")
  return(raw_res)
}

# Print simulated calcium traces
#' 
#' @export
#' 
print.calcium_sim <- function(x, ...) {
  params <- attr(x, "params")
  
  cat("=== Calcium Imaging Simulation Data ===\n")
  cat(sprintf("Dimensions   : %d Neurons x %d Frames\n", params$N, params$nf))
  cat(sprintf("Length       : %d \n", params$nf))
  cat(sprintf("Neuron Number: %d \n", params$N))
  cat(sprintf("Model        : AR(%d) process\n", length(params$g)))
  cat(sprintf("Components   : $Y (Noisy Fluorescence), $c (Calcium Trace), $s (Spike Train)\n"))
  
  invisible(x)
}