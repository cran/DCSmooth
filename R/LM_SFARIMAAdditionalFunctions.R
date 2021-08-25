################################################################################
#                                                                              #
#                DCSmooth Package: additional Functions for SFARIMA              #
#                                                                              #
################################################################################

#----------------------------Simulation Function-------------------------------#

#' Simulation of a \eqn{SFARIMA(p, q, d)}-process
#' 
#' @description \code{sfarima.sim} simulates a specified SFARIMA-model
#'  on a lattice with normally distributed innovations.
#' 
#' @section Details:
#' Simulation of a separable spatial fractionally ARIMA process (SFARIMA). This 
#' function returns an object of class \code{"sfarima"}. The simulated
#' innovations are created from a normal distribution with specified variance
#' \eqn{\sigma^2}{sigma^2}.
#' 
#' @param n_x Number of simulated observation rows.
#' @param n_t Number of simulated observation columns.
#' @param model A list containing the coefficient matrices \code{ar} and 
#'  \code{ma} of the QARMA model, the long memory parameter vector \code{d} as
#'  well as the standard deviation of innovations \code{sigma}.
#' 
#' @return The function returns an object of class \code{"sfarima"}, consisting
#' of
#' 
#'  \tabular{ll}{
#'   \code{Y} \tab A \eqn{n_x \times n_t}{n_x x n_t}-matrix of simulated values
#'   of the specified SFARIMA process.\cr
#'   \code{innov} \tab The innovations used for simulation, iid. drawn from a 
#'    normal distribution with zero mean and variance 
#'    \eqn{\sigma^2}{(sigma)^2}.\cr
#'   \code{model} \tab The model used for simulation, inherited from input.\cr
#'   \code{stnry} \tab An logical variable indicating whether the simulated 
#'   model is stationary.\cr
#' }
#' 
#' @section Details: see the vignette for further details.
#' 
#' @seealso \code{\link{qarma.est}}
#' 
#' @examples
#' # See vignette("DCSmooth") for examples and explanation
#'  
#' ma = matrix(c(1, 0.2, 0.4, 0.1), nrow = 2, ncol = 2)
#' ar = matrix(c(1, 0.5, -0.1, 0.1), nrow = 2, ncol = 2)
#' d = c(0.1, 0.1)
#' sigma = 0.5
#' sfarima_model = list(ar = ar, ma = ma, d = d, sigma = sigma)
#' 
#' sfarima_sim = sfarima.sim(100, 100, model = sfarima_model)
#' surface.dcs(sfarima_sim$Y)
#' 
#' @export

sfarima.sim <- function(n_x, n_t, model)
{
  ar_mat = as.matrix(model$ar); ma_mat = as.matrix(model$ma)
  ar_x = -ar_mat[-1, 1]; ar_t = -ar_mat[1, -1]
  ma_x = ma_mat[-1, 1]; ma_t = ma_mat[1, -1]
  
  # check if provided model is correctly specified
  if (isFALSE(all.equal(as.matrix(ar_mat[-1, -1]), ar_x %*% ar_t))  ||
      isFALSE(all.equal(as.matrix(ma_mat[-1, -1]), ma_x %*% ma_t)))
  {
    warning("Provided coefficient matrices do not specify a separable process.")
  }
  if (!any(is.numeric(model$d)) || any(abs(model$d) > 0.5))
  {
    stop("Long memory parameter \"d\" incorrectly specified.")
  }
  

  # meta options
  nstart = max(floor(1.5 * c(n_x, n_t)), 150)
  k_x = min(50, n_x); k_t = min(50, n_t)
  
  n_x = n_x + nstart
  n_t = n_t + nstart
  eps_mat = matrix(stats::rnorm(n_x * n_t), n_x, n_t) * model$sigma
  
  ma_inf_x = c(1, stats::ARMAtoMA(ar = ar_x, ma = ma_x, lag.max = k_x))
  d_x = choose(-model$d[1], 0:k_x) * ((-1)^(0:k_x))
  coef_x = cumsum_part_reverse(d_x, ma_inf_x)
  
  ma_inf_t = t(c(1, stats::ARMAtoMA(ar = ar_t, ma = ma_t, lag.max = k_t)))
  d_t = choose(-model$d[2], 0:k_t) * ((-1)^(0:k_t))
  coef_t = cumsum_part_reverse(d_t, ma_inf_t)
  
  X1.sim = X2.sim = matrix(0, n_x, n_t)

  for(j in 1:n_t) {
    if (j <= k_t) {
      X2.sim[, j] = eps_mat[, j:1, drop = FALSE] %*% coef_t[1:j, drop = FALSE] 
    }
    else {
      X2.sim[, j] = eps_mat[, j:(j - k_t)] %*% coef_t
    }
  }
  for(i in 1:n_x) {
    if (i <= k_x) {
      X1.sim[i, ] = coef_x[1:i, drop = FALSE] %*% X2.sim[i:1, , drop = FALSE]
    }
    else {
      X1.sim[i, ] = t(coef_x) %*% X2.sim[i:(i - k_x), ]
    }
  }
    
  sfarima_out = X1.sim[(nstart + 1):n_x, (nstart + 1):n_t]
  error_out = eps_mat[(nstart + 1):n_x, (nstart + 1):n_t]
  coef_out = list(Y = sfarima_out, innov = error_out, model = model,
                  stnry = TRUE)
  class(coef_out) = "sfarima"
  attr(coef_out, "subclass") = "sim"
    
  return(coef_out)
}

#----------------------------------------------------------------#

macoef <- function(ar = 0, ma = 0, d = 0, k = 50) {
  p = length(ar[ar != 0])
  q = length(ma[ma != 0])
  if (p == 0) {
    ar = 0
  }
  if (q == 0) {
    ma = 0
  }
  ma.coef = c(1, ma, rep(0, k - q))
  arma.coef = (1:(k + 1)) * 0
  arma.coef[1] = 1
  
  if (p > 0 | q > 0) {
    for (i in 2:(k + 1)) {
      if ((i - p) < 1) {
        arma.coef[i] = sum(ar[1:(p - abs(i - p) - 1)] * arma.coef[(i - 1):1]) - ma.coef[i]
      } else {
        arma.coef[i] = sum(ar[1:p] * arma.coef[(i - 1):(i - p)]) - ma.coef[i]
      }
    }
  } else {
    arma.coef + 1
  }
  
  d.coef = choose(-d, 0:k) * ((-1)^(0:k))
  
  coef.all = (1:(k + 1)) * 0
  for (j in 1:(k + 1)) {
    coef.all[j] = sum(d.coef[1:j] * arma.coef[j:1])
  }
  return(coef.all)
}



