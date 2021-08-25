################################################################################
#                                                                              #
#                     DCSmooth Package: SFARIMA Estimation                     #
#                                                                              #
################################################################################

#' Estimation of a SFARIMA-process
#'
#' @description Parametric Estimation of a \eqn{SFARIMA(p, q, d)}-process on a 
#'  lattice.
#' 
#' @section Details:
#' The MA- and AR-parameters as well as the long-memory parameters \deqn{d} of a
#' SFARIMA process are estimated by minimization of the residual sum of squares
#' RSS. Lag-orders of \eqn{SFARIMA(p, q, d)} are given by \eqn{p = (p_1, p_2), 
#' q = (q_1, q_2)}{p = (p1, p2), q = (q1, q2)}, where \eqn{p_1, q_1}{p1, q1} are
#' the lags over the rows and \eqn{p_2, q_2}{p2, q2} are the lags over the 
#' columns. The estimated process is based on the (separable) model 
#' \deqn{\varepsilon_{ij} = \Psi_1(B) \Psi_2(B) \eta_{ij}}, where \deqn{\Psi_i =
#' (1 - B_i)^{-d_i}\phi^{-1}_i(B_i)\psi_i(B_i), i = 1,2}.

#' 
#' @param Y A numeric matrix that contains the demeaned observations of the
#'   random field or functional time-series.
#' @param model_order A list containing the orders of the SFARIMA model in the
#'   form \code{model_order = list(ar = c(p1, p2), ma = c(q1, q2))}. Default
#'   value is a \eqn{SFARIMA((1, 1), (1, 1), d)} model.
#' 
#' @return The function returns an object of class \code{"sfarima"} including
#'  \tabular{ll}{
#'   \code{Y} \tab The matrix of observations, inherited from input.\cr
#'   \code{innov} The estimated innovations.\cr
#'   \code{model} \tab The estimated model consisting of the coefficient 
#'   matrices \code{ar} and \code{ma}, the estimated long memory parameters
#'   \code{d} and standard deviation of innovations \code{sigma}.\cr
#'   \code{stnry} \tab An logical variable indicating whether the estimated
#'   model is stationary.\cr
#' }
#' 
#' @seealso \code{\link{qarma.est}, \link{sfarima.sim}}
#' 
#' @examples
#' # See vignette("DCSmooth") for examples and explanation
#' 
#' ## simulation of SFARIMA process
#' ma = matrix(c(1, 0.2, 0.4, 0.1), nrow = 2, ncol = 2)
#' ar = matrix(c(1, 0.5, -0.1, 0.1), nrow = 2, ncol = 2)
#' d = c(0.1, 0.1)
#' sigma = 0.5
#' sfarima_model = list(ar = ar, ma = ma, d = d, sigma = sigma)
#' sfarima_sim = sfarima.sim(50, 50, model = sfarima_model)
#' 
#' ## estimation of SFARIMA process
#' sfarima.est(sfarima_sim$Y)$model
#' sfarima.est(sfarima_sim$Y, 
#'            model_order = list(ar = c(1, 1), ma = c(0, 0)))$model
#' 
#' @export

sfarima.est = function(Y, model_order = list(ar = c(1, 1), ma = c(1, 1)))
{
  n_x = dim(Y)[1]; n_t = dim(Y)[2]
  
  theta_init = rep(0, times = sum(unlist(model_order)) + 2)
  theta_opt  = stats::optim(theta_init, sfarima.rss, R_mat = Y, 
                            model_order = model_order,
                            method = "Nelder-Mead")
  
  # put coefficients into matrices
  d_vec = theta_opt$par[1:2]
  ar_x = c(1, -theta_opt$par[2 + seq_len(model_order$ar[1])])
  ar_t = c(1, -theta_opt$par[2 + model_order$ar[1] + 
                             seq_len(model_order$ar[2])])
  ma_x = c(1, theta_opt$par[2 + sum(model_order$ar) +
                            seq_len(model_order$ma[1])])
  ma_t = c(1, theta_opt$par[2 + sum(model_order$ar) + model_order$ma[1] + 
                            seq_len(model_order$ma[2])])
  
  # prepare results for output
  ar_mat = ar_x %*% t(ar_t)
  ma_mat = ma_x %*% t(ma_t)
  stdev = sqrt(theta_opt$value/(n_x * n_t))
  model = list(ar = ar_mat, ma = ma_mat, d = d_vec, sigma = stdev)
  innov = sfarima.residuals(R_mat = Y, model = model)
  
  # check stationarity
  statTest = qarma.statTest(ar_mat)
  if (!statTest)
  {
    warning("SFARIMA model not stationary, try another order for the AR-parts.")
  }
  
  coef_out = list(Y = Y, innov = innov, model = list(ar = ar_mat, ma = ma_mat,
                                  d = d_vec, sigma = stdev), stnry = statTest)
  class(coef_out) = "qarma"
  attr(coef_out, "subclass") = "est"
  return(coef_out)
}

sfarima.rss = function(theta, R_mat,
                       model_order = list(ar = c(1, 1), ma = c(1, 1)))
{
  n_x = dim(R_mat)[1]; n_t = dim(R_mat)[2]
  k_x = min(50, n_x); k_t = min(100, n_t)
  
  # get coefficients from theta
  # theta = (d1, d2, ar_x, ar_t, ma_x, ma_t)
  d_vec = theta[1:2]
  ar_x = theta[2 + seq_len(model_order$ar[1])]
  ar_t = theta[2 + model_order$ar[1] + seq_len(model_order$ar[2])]
  ma_x = theta[2 + sum(model_order$ar) + seq_len(model_order$ma[1])]
  ma_t = theta[2 + sum(model_order$ar) + model_order$ma[1] + 
                 seq_len(model_order$ma[1])]
  
  # result matrices
  E_itm = R_mat * 0   # intermediate results
  E_fnl = R_mat * 0   # final results
  
  # calculate AR(inf) coefficients with long-memory
  ar_inf_x = c(1, astsa::ARMAtoAR(ar = ar_x, ma = ma_x, lag.max = k_x))
  d_x = choose(d_vec[1], 0:k_x) * ((-1)^(0:k_x))
  coef_x = cumsum_part_reverse(d_x, ar_inf_x)
  
  ar_inf_t = t(c(1, astsa::ARMAtoAR(ar = ar_t, ma = ma_t, lag.max = k_t)))
  d_t = choose(d_vec[2], 0:k_t) * ((-1)^(0:k_t))
  coef_t = cumsum_part_reverse(d_t, ar_inf_t)
  
  for (j in 1:n_t)
  {
    E_itm[, j] = R_mat[, j:max(1, j - k_t + 1), drop = FALSE] %*%
      coef_t[1:min(j, k_t), drop = FALSE]
  }
  
  for (i in 1:n_x)
  {
    E_fnl[i, ] = coef_x[1:min(i, k_x), drop = FALSE] %*%
      E_itm[i:max(1, i - k_x + 1), , drop = FALSE]
  }
  
  RSS = sum(E_fnl^2)
  
  return(RSS)
}

sfarima.residuals = function(R_mat, model)
{
  n_x = dim(R_mat)[1]; n_t = dim(R_mat)[2]
  k_x = min(50, n_x); k_t = min(100, n_t)
  
  # result matrices
  E_itm = R_mat * 0   # intermediate results
  E_fnl = R_mat * 0   # final results
  
  # calculate AR(inf) coefficients with long-memory
  ar_inf_x = c(1, astsa::ARMAtoAR(ar = model$ar[-1, 1],
                                  ma = model$ma[-1, 1], lag.max = k_x))
  d_x = choose(-model$d[1], 0:k_x) * ((-1)^(0:k_x))
  coef_x = cumsum_part_reverse(d_x, ar_inf_x)
  
  ar_inf_t = t(c(1, astsa::ARMAtoAR(ar = model$ar[1, -1],
                                    ma = model$ma[1, -1], lag.max = k_t)))
  d_t = choose(-model$d[2], 0:k_t) * ((-1)^(0:k_t))
  coef_t = cumsum_part_reverse(d_t, ar_inf_t)
  
  for (j in 1:n_t)
  {
    E_itm[, j] = R_mat[, j:max(1, j - k_t + 1), drop = FALSE] %*%
      ar_inf_t[1:min(j, k_t), drop = FALSE]
  }
  
  for (i in 1:n_x)
  {
    E_fnl[i, ] = ar_inf_x[1:min(i, k_x), drop = FALSE] %*%
      E_itm[i:max(1, i - k_x + 1), , drop = FALSE]
  }
  
  return(E_fnl)
}