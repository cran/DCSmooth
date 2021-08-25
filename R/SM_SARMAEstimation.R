################################################################################
#                                                                              #
#                DCSmooth Package: estimation of SARMA for cf                  #
#                                                                              #
################################################################################

# This file includes all functions related to the estimation of sSARMA-processes
# for the cf coefficient of the bandwidth selection procedure

#-----------------------Calculation of cf coefficient--------------------------#

sarma.cf = function(Y, model_order = list(ar = c(1, 1), ma = c(1, 1)))
{
  # estimation of qarma model
  sarma_est = sarma.est(Y, model_order = model_order)
  
  # get fractions for spectral density
  cf = sum(sarma_est$model$ma)^2/sum(sarma_est$model$ar)^2 *
    sarma_est$model$sigma^2
  sarma_model = sarma_est$model
  return_list = list(cf = cf, model = sarma_model)
  
  return(return_list)
}

#------------------------RSS Estimation of SARMA-------------------------------#


sarma.est = function(Y, model_order = list(ar = c(1, 1), ma = c(1, 1)))
{
  n_x = dim(Y)[1]; n_t = dim(Y)[2]
  theta_init = rep(0, times = sum(unlist(model_order)))
  theta_opt  = stats::optim(theta_init, sarma.rss, R_mat = Y, 
                              method = "Nelder-Mead")
  
  # put coefficients into matrices
  ar_x = c(1, -theta_opt$par[seq_len(model_order$ar[1])])
  ar_t = c(1, -theta_opt$par[model_order$ar[1] + seq_len(model_order$ar[2])])
  ma_x = c(1, theta_opt$par[sum(model_order$ar) + seq_len(model_order$ma[1])])
  ma_t = c(1, theta_opt$par[sum(model_order$ar) + model_order$ma[1] + 
                              seq_len(model_order$ma[2])])
  
  # prepare results for output
  ar_mat = ar_x %*% t(ar_t)
  ma_mat = ma_x %*% t(ma_t)
  stdev = sqrt(theta_opt$value/(n_x * n_t))
  model = list(ar = ar_mat, ma = ma_mat, sigma = stdev)
  innov = sarma.residuals(R_mat = Y, model = model)
  
  # check stationarity
  statTest = qarma.statTest(ar_mat)
  if (!statTest)
  {
    warning("QARMA model not stationary, try another order for the AR-parts.")
  }
  
  coef_out = list(Y = Y, innov = innov, model = list(ar = ar_mat, ma = ma_mat,
                  sigma = stdev), stnry = statTest)
  class(coef_out) = "qarma"
  attr(coef_out, "subclass") = "est"
  return(coef_out
         )
}

sarma.rss = function(theta, R_mat,
                     model_order = list(ar = c(1, 1), ma = c(1, 1)))
{
  n_x = dim(R_mat)[1]; n_t = dim(R_mat)[2]
  k_x = min(50, n_x); k_t = min(50, n_t)
  
  # get coefficients from theta
  # theta = (ar_x, ar_t, ma_x, ma_t)
  ar_x = theta[seq_len(model_order$ar[1])]
  ar_t = theta[model_order$ar[1] + seq_len(model_order$ar[2])]
  ma_x = theta[sum(model_order$ar) + seq_len(model_order$ma[1])]
  ma_t = theta[sum(model_order$ar) + model_order$ma[1] + 
                 seq_len(model_order$ma[2])]
  
  # result matrices
  E_itm = R_mat * 0   # intermediate results
  E_fnl = R_mat * 0   # final results
  
  # Two-step estimation of e_ij
  ar_inf_x = c(1, astsa::ARMAtoAR(ar = ar_x, ma = ma_x, lag.max = k_x))
  ar_inf_t = t(c(1, astsa::ARMAtoAR(ar = ar_t, ma = ma_t, lag.max = k_t)))
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
  
  RSS = sum(E_fnl^2)
  
  return(RSS)
}

sarma.residuals = function(R_mat, model)
{
  n_x = dim(R_mat)[1]; n_t = dim(R_mat)[2]
  k_x = min(50, n_x); k_t = min(50, n_t)
  
  # result matrices
  E_itm = R_mat * 0   # intermediate results
  E_fnl = R_mat * 0   # final results
  
  ar_inf_x = c(1, astsa::ARMAtoAR(ar = -model$ar[-1, 1], ma = model$ma[-1, 1],
                                  lag.max = k_x))
  ar_inf_t = t(c(1, astsa::ARMAtoAR(ar = -model$ar[1, -1], ma = model$ma[1, -1],
                                    lag.max = k_t)))
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