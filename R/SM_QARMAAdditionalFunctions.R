################################################################################
#                                                                              #
#                DCSmooth Package: additional Functions for QARMA              #
#                                                                              #
################################################################################

# Simulation of QARMA-processes, estimation of the spectral density


#----------------------------Simulation Function-------------------------------#

#' Simulation of a \eqn{QARMA(p, q)}-process
#' 
#' @description \code{qarma.sim} is simulates a specified QARMA-model
#'  on a lattice, with normally distributed innovations.
#' 
#' @section Details:
#' Simulation of a top-left dependent spatial ARMA process (QARMA). This 
#' function returns an object of class \code{"qarma"}. The simulated innovations
#' are created from a normal distribution with specified variance
#' \eqn{\sigma^2}{sigma^2}.
#' 
#' @param n_x Number of simulated observation rows.
#' @param n_t Number of simulated observation columns.
#' @param model A list containing the coefficient matrices \code{ar} and 
#'  \code{ma} of the QARMA model as well as the standard deviation of 
#'  innovations \code{sigma}.
#' 
#' @return The function returns an object of class \code{"qarma"}, consisting of
#' 
#'  \tabular{ll}{
#'   \code{Y} \tab A \eqn{n_x \times n_t}{n_x x n_t}-matrix of simulated values
#'   of the specified QARMA process.\cr
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
#' sigma = 0.5
#' q_model = list(ar = ar, ma = ma, sigma = sigma)
#' 
#' q_sim = qarma.sim(100, 100, model = q_model)
#' surface.dcs(q_sim$Y)
#' 
#' @export

qarma.sim = function(n_x, n_t, model)
{
  ar_mat = as.matrix(model$ar); ma_mat = as.matrix(model$ma)
  ar_x = dim(ar_mat)[1] - 1; ar_t = dim(ar_mat)[2] - 1
  ma_x = dim(ma_mat)[1] - 1; ma_t = dim(ma_mat)[2] - 1
  x_init = max(ar_x, ma_x) + 1
  t_init = max(ar_t, ma_t) + 1
  
  # set coefficients for zero-lags
  ma_mat[1, 1] = 1 # MA-coefficients
  ar_mat[1, 1] = 0 # AR-coefficients
  
  n_mat = floor(1.25 * c(n_x, n_t))
  error_mat = matrix(stats::rnorm(prod(n_mat)), nrow = n_mat[1], 
                     ncol = n_mat[2]) * model$sigma
  arma_mat = error_mat
  
  for (i in x_init:n_mat[1])
  {
    for (j in t_init:n_mat[2])
    {
      arma_mat[i, j] = sum(-ar_mat * arma_mat[i:(i - ar_x), j:(j - ar_t)]) +
        sum(ma_mat * error_mat[i:(i - ma_x), j:(j - ma_t)])
    }
  }
  
  arma_out = arma_mat[(n_mat[1] - n_x + 1):
                        n_mat[1], (n_mat[2] - n_t + 1):n_mat[2]]
  error_out = error_mat[(n_mat[1] - n_x + 1):
                          n_mat[1], (n_mat[2] - n_t + 1):n_mat[2]]
  
  coef_out = list(Y = arma_out, innov = error_out, model = model,
                  stnry = TRUE)
  class(coef_out) = "qarma"
  attr(coef_out, "subclass") = "sim"
  
  return(coef_out)
}

#----------------------Calculation of spectral density-------------------------#

qarma.ssde = function(Y, ar, ma, st_dev)
{
  X = T = seq(from = -pi, to = pi, length.out = 100)
  y_spectrum = matrix(NA, nrow = 100, ncol = 100)
  
  for (i in 1:100)
  {
    for (j in 1:100)
    {
      omega = c(X[i], T[j])
      y_spectrum[i, j] = qarma.spectral.density(Y, ar, ma, st_dev, omega)
    }
  }
  return(y_spectrum)
}

qarma.spectral.density = function(Y, ar, ma, st_dev, omega)
{
  lag_ar = dim(ar) - 1; lag_ma = dim(ma) - 1
  z1 = complex(argument = omega[1]); z2 = complex(argument = omega[2])
  
  g00 = stats::sd(Y)^2
  
  ar_sum = Re(sum(ar * (z1^(lag_ar[1]:0)) %*% t(z2^(lag_ar[2]:0))) * 
               sum(ar * ((1/z1)^(lag_ar[1]:0)) %*% t((1/z2)^(lag_ar[2]:0))))
  ma_sum = Re(sum(ma * (z1^(lag_ma[1]:0)) %*% t(z2^(lag_ma[2]:0))) * 
               sum(ma * ((1/z1)^(lag_ma[1]:0)) %*% t((1/z2)^(lag_ma[2]:0))))
  
  spectral_dens_out = st_dev^2/g00 * ma_sum/ar_sum
  
  return(spectral_dens_out)
}
