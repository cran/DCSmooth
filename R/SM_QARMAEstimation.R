################################################################################
#                                                                              #
#                DCSmooth Package: estimation of QARMA for cf                  #
#                                                                              #
################################################################################

# This file includes all functions related to the estimation of QARMA-processes
# for the cf coefficient of the bandwidth selection procedure


### Part A: Estimation Functions

#-----------------------Calculation of cf coefficient--------------------------#

qarma.cf = function(Y, model_order = list(ar = c(1, 1), ma = c(1, 1)))
{
  # estimation of qarma model
  qarma_model = qarma.est(Y, model_order = model_order)

  # get fractions for spectral density
  cf = sum(qarma_model$model$ma)^2/sum(qarma_model$model$ar)^2 *
           qarma_model$model$sigma^2
  return_list = list(cf = cf, qarma_model = qarma_model)
  
  return(return_list)
}

#-------------Hannan-Rissanen Estimation Function for QARMA--------------------#

#' Estimation of a QARMA-process
#'
#' @description Parametric Estimation of a \eqn{QARMA(p, q)}-process on a 
#'  lattice using the Hannan-Rissanen algorithm.
#' 
#' @section Details:
#' The MA- and AR-parameters of a top-left quadrant ARMA process are estimated
#' using the Hannan-Rissanen Algorithm. The lag-orders of
#' the \eqn{QARMA(p, q)} are given by \eqn{p = (p_1, p_2), q = (q_1, q_2)}{p =
#' (p1, p2), q = (q1, q2)}, where \eqn{p_1, q_1}{p1, q1} are the lags over the
#' rows and \eqn{p_2, q_2}{p2, q2} are the lags over the columns. The estimation
#' process is based on the model
#' \deqn{\phi(B_{1}B_{2})X_{i,j} = \theta(B_{1}B_{2})u_{i,j}}{\phi(B1 B2)
#' X[i,j] = \theta(B1 B2)u[i,j]}
#' 
#' @param Y A numeric matrix that contains the demeaned observations of the
#'   random field or functional time-series.
#' @param model_order A list containing the orders of the QARMA model in the
#'   form \code{model_order = list(ar = c(p1, p2), ma = c(q1, q2))}. Default
#'   value is a \eqn{QARMA((1, 1), (1, 1))} model.
#' 
#' @return The function returns an object of class \code{"qarma"} including
#'  \tabular{ll}{
#'   \code{Y} \tab The matrix of observations, inherited from input.\cr
#'   \code{innov} The estimated innovations.\cr
#'   \code{model} \tab The estimated model consisting of the coefficient 
#'   matrices \code{ar} and \code{ma} and standard deviation of innovations
#'   \code{sigma}.\cr
#'   \code{stnry} \tab An logical variable indicating whether the estimated
#'   model is stationary.\cr
#' }
#' 
#' @seealso \code{\link{qarma.sim}}
#' 
#' @examples
#' # See vignette("DCSmooth") for examples and explanation
#' 
#' ## simulation of QARMA process
#' ma = matrix(c(1, 0.2, 0.4, 0.1), nrow = 2, ncol = 2)
#' ar = matrix(c(1, 0.5, -0.1, 0.1), nrow = 2, ncol = 2)
#' sigma = 0.5
#' q_model = list(ar = ar, ma = ma, sigma = sigma)
#' qarma_simulated = qarma.sim(100, 100, model = q_model)
#' qarma_simulated$model
#' 
#' ## estimation of QARMA process
#' qarma.est(qarma_simulated$Y)$model
#' qarma.est(qarma_simulated$Y, 
#'            model_order = list(ar = c(1, 1), ma = c(0, 0)))$model
#' 
#' @export

qarma.est = function(Y, model_order = list(ar = c(1, 1), ma = c(1, 1)))
{
  nX = dim(Y)[1]; nT = dim(Y)[2]
  
  # calculate ARMA-orders for different purposes
  max_lag_x = max(model_order$ar[1], model_order$ma[1])
  max_lag_t = max(model_order$ar[2], model_order$ma[2])
  total_lag_ar = (model_order$ar[1] + 1) * (model_order$ar[2] + 1) - 1
  total_lag_ma = (model_order$ma[1] + 1) * (model_order$ma[2] + 1) - 1
  max_lag_ar = max(model_order$ar)
  max_lag_ma = max(model_order$ma)
  
  ### AR-only auxiliary estimation model by YW-estimation ###
  m1_ar = max(max_lag_x, 1) + ifelse(total_lag_ma > 0, 2, 0) 
                                          # increase order of aux AR model here
  m2_ar = max(max_lag_t, 1) + ifelse(total_lag_ma > 0, 2, 0)
  ar_aux = qarma.yw_matrix(Y, ar_order = c(m1_ar, m2_ar))
  
  # calculate residuals from auxiliary AR model
  residuals_mat = matrix(0, nrow = nX - m1_ar, ncol = nT - m2_ar)
  for (i in 1:(nX - m1_ar))
  {
    for (j in 1:(nT - m2_ar))
    {
      residuals_mat[i, j] = sum(ar_aux * Y[(i + m1_ar):i, (j + m2_ar):j])
    }
  }
  
  if (total_lag_ma > 0)
  {
    # set up submatrices of Y, res for use in fill_matrix procedur
    Y_mat = Y[(m1_ar + model_order$ma[1] - model_order$ar[1] + 1):nX,
              (m2_ar + model_order$ma[2] - model_order$ar[2] + 1):nT]
    
    # set up matrix for complete arma (dimension (nX - ar_x)*(nT - ar_t))
    matrix_arma = cbind(.qarma.fill_lag_matrix(Y_mat, lag_x = model_order$ar[1],
                            lag_t = model_order$ar[2], include_00 = TRUE,
                            name_prefix = "ar"),
                .qarma.fill_lag_matrix(residuals_mat, lag_x = model_order$ma[1],
                            lag_t = model_order$ma[2], include_00 = FALSE,
                            name_prefix = "ma"))
    
    # regression for QARMA model
    arma_reg = stats::lm(ar_00 ~ . + 0, data = matrix_arma)
    
    # fill ar and ma matrices
    # updated to lag-order (phi_00/psi_00 in upper left corner)
    # ar_mat is lhs of QARMA-equation (phi_00 = 1)
    ar_mat = matrix(c(1, -arma_reg$coef[1:total_lag_ar]), byrow = TRUE,
                    nrow = (model_order$ar[1] + 1), 
                    ncol = (model_order$ar[2] + 1))
    ma_mat = matrix(c(1, arma_reg$coef[(total_lag_ar + 1):
                    (total_lag_ar + total_lag_ma)]), byrow = TRUE,
                    nrow = (model_order$ma[1] + 1),
                    ncol = (model_order$ma[2] + 1))
    if (total_lag_ar == 0)
    {
      ar_mat[1, 1] = 1
    }
  } else {
    arma_reg = list(ar_aux = as.vector(ar_aux)[-1])
    
    # fill ar matrix (and ma matrix)
    ar_mat = matrix(c(1, arma_reg$ar_aux[1:total_lag_ar]),
                nrow = (model_order$ar[1] + 1), ncol = (model_order$ar[2] + 1))
    ma_mat = matrix(1)
    
    arma_reg$residuals = residuals_mat
  }
  
  innov = arma_reg$residuals #residuals_mat
  
  improve = FALSE
  if (improve == TRUE)
  {
    stdev = sqrt(sum(innov^2)/((nX - max_lag_x - model_order$ar[1]) * 
                                 (nT - max_lag_t - model_order$ar[2])))
    
    arma_reg$coef = qarma.est_3rdstep(Y, ar_mat, ma_mat, model_order) +
                  arma_reg$coef
    
    ar_mat = matrix(c(1, -arma_reg$coef[1:total_lag_ar]), byrow = TRUE,
                    nrow = (model_order$ar[1] + 1), 
                    ncol = (model_order$ar[2] + 1))
    ma_mat = matrix(c(1, arma_reg$coef[(total_lag_ar + 1):
                    (total_lag_ar + total_lag_ma)]), byrow = TRUE,
                    nrow = (model_order$ma[1] + 1),
                    ncol = (model_order$ma[2] + 1))
  }
  
  # check stationarity
  statTest = qarma.statTest(ar_mat)
  if (!statTest)
  {
    warning("QARMA model not stationary, try another order for the AR-parts.")
  }
  
  # preparation of output
  rownames(ar_mat) = paste0("lag ", 0:model_order$ar[1])
  colnames(ar_mat) = paste0("lag ", 0:model_order$ar[2])
  rownames(ma_mat) = paste0("lag ", 0:model_order$ma[1])
  colnames(ma_mat) = paste0("lag ", 0:model_order$ma[2])
  
  stdev = sqrt(sum(innov^2)/((nX - max_lag_x - model_order$ar[1]) * 
                               (nT - max_lag_t - model_order$ar[2])))
  coef_out = list(Y = Y, innov = innov, model = list(ar = ar_mat, ma = ma_mat,
                  sigma = stdev), stnry = statTest)
  class(coef_out) = "qarma"
  attr(coef_out, "subclass") = "est"
  return(coef_out)
}

# Function for refinement of QARMA estimation
# (does not seem to provide a significant improvement)


# TODO tidy up function, include .qarma.fill_lag matrix instead of loops
qarma.est_3rdstep = function(Y, ar_mat, ma_mat, model_order)
{
  nX = dim(Y)[1]; nT = dim(Y)[2]
  
  # value (x,t) of ma_mat has to be zero
  ma_mat_0 = ma_mat; ma_mat_0[1, 1] = 0
  ar_mat_0 = ar_mat; ar_mat_0[1, 1] = 0
  
  # calculate ARMA-orders for different purposes
  max_lag_x = max(model_order$ar[1], model_order$ma[1])
  max_lag_t = max(model_order$ar[2], model_order$ma[2])
  total_lag_ar = (model_order$ar[1] + 1) * (model_order$ar[2] + 1) - 1
  total_lag_ma = (model_order$ma[1] + 1) * (model_order$ma[2] + 1) - 1
  
  z_matrix = matrix(0, nrow = nX, ncol = nT)
  v_matrix = w_matrix = z_matrix
  
  for (i in (max_lag_x + 1):nX)
  {
    for (j in (max_lag_t + 1):nT)
    {
      ind_ar_x = i:(i - model_order$ar[1])
      ind_ar_t = j:(j - model_order$ar[2])
      ind_ma_x = i:(i - model_order$ma[1])
      ind_ma_t = j:(j - model_order$ma[2])
      
      # Strange Part here
      z_matrix[i, j] = sum(ar_mat * Y[ind_ar_x, ind_ar_t]) -
        sum(ma_mat_0 * z_matrix[ind_ma_x, ind_ma_t])
      v_matrix[i, j] = sum(-ar_mat_0 * v_matrix[ind_ar_x, ind_ar_t]) +
        z_matrix[i, j]
      w_matrix[i, j] = sum(-ma_mat_0 * w_matrix[ind_ma_x, ind_ma_t]) +
        z_matrix[i, j]
    }
  }
  
  zvw_lm_matrix = data.frame(matrix(NA,
                                    nrow = (nX - max_lag_x) * (nT - max_lag_t),
                                    ncol = total_lag_ar + total_lag_ma + 1))
  zvw_lm_matrix[, 1] =
    as.vector(z_matrix[(max_lag_x + 1):nX, (max_lag_t + 1):nT])
  names(zvw_lm_matrix)[1] = "Z"
  
  for (i in 0:model_order$ar[1])
  {
    for (j in 0:model_order$ar[2])
    {
      if (!(i == 0 && j == 0))
      {
        zvw_lm_matrix[, i*(model_order$ar[2] + 1) + j + 1] =
          as.vector(v_matrix[(max_lag_x + 1):nX - i, (max_lag_t + 1):nT - j])
        names(zvw_lm_matrix)[i*(model_order$ar[2] + 1) + j + 1] =
          paste0("V_", i, j)
      }
    }
  }
  for (i in 0:model_order$ma[1])
  {
    for (j in 0:model_order$ma[2])
    {
      if (!(i == 0 && j == 0))
      {
        zvw_lm_matrix[, total_lag_ar + i*(model_order$ma[2] + 1) + j + 1] =
          as.vector(w_matrix[(max_lag_x + 1):nX - i, (max_lag_t + 1):nT - j])
        names(zvw_lm_matrix)[total_lag_ar + i*(model_order$ma[2] + 1) + j + 1] =
          paste0("W_", i, j)
      }
    }
  }
  
  # linear regression
  zvw_lm = stats::lm(Z ~ . + 0, data = zvw_lm_matrix)
  
  return(zvw_lm$coef)
}

#-----------------------------Auxiliary Functions------------------------------#

.qarma.fill_lag_matrix = function(Y, lag_x, lag_t, include_00 = TRUE,
                                  name_prefix)
{
  nX = dim(Y)[1]; nT = dim(Y)[2]
  lag_obs = (nX - lag_x) * (nT - lag_t)
  lag_n = (lag_x + 1) * (lag_t + 1)
  
  matrix_out = data.frame(matrix(NA, nrow = lag_obs, ncol = lag_n))
  
  for (i in 0:lag_x) # use lag orders as loop indices
  {
    for (j in 0:lag_t)
    {
      index_x = (lag_x - i + 1):(nX - i)
      index_t = (lag_t - j + 1):(nT - j)
      matrix_out[, (j + 1) + i * (lag_t + 1)] = as.vector(Y[index_x, index_t])
      colnames(matrix_out)[(j + 1) + i * (lag_t + 1)] = 
        paste0(name_prefix, "_", i, j)
    }
  }
  
  if (include_00 == FALSE && lag_n > 1)
  {
    matrix_out = matrix_out[2:lag_n]
  } else if (include_00 == FALSE && lag_n == 1) {
    matrix_out = matrix_out * 0
  }
  
  return(matrix_out)
}

#-------------------Yule-Walker-Estimation Function for AR---------------------#

qarma.yw_matrix = function(Y, ar_order = c(1, 1), ma_order = c(0, 0))
{
  # set up parameters, result matrices and vectors etc.
  nX = dim(Y)[1]; nT = dim(Y)[2]; n = nX * nT
  p1_ar = ar_order[1]
  p2_ar = ar_order[2]
  d_ar = (p1_ar + 1) * (p2_ar + 1) - 1

  acf_matrix_ar = matrix(NA, nrow = d_ar, ncol = d_ar)
  acf_vector_ar = vector(mode = "numeric", length = d_ar)
  # st_mat is over rows first (i, j) -> (i - 1, j), which is normal R order
  st_mat = expand.grid(p1 = c(p1_ar:0), p2 = c(p2_ar:0))
  st_mat = st_mat[order(st_mat[, 1], st_mat[, 2]), ]

  # calculate acf_matrix and acf_vector
  for (k in 1:d_ar)
  {
    for (l in 1:d_ar)
    {
      st_vec = st_mat[k + 1, ] - st_mat[l + 1, ] + ma_order
      acf_matrix_ar[k, l] = .matrix_acf(Y, st_vec)
    }
    acf_vector_ar[k] = .matrix_acf(Y, st_mat[k + 1, ] + ma_order)
  }
  
  # solving Yule-Walker equations
  yw.estimators = solve(acf_matrix_ar) %*% acf_vector_ar
  phi_ar_matrix = matrix(c(1, -yw.estimators), nrow = (p1_ar + 1),
                         ncol = (p2_ar + 1), byrow = TRUE)

  return(phi_ar_matrix)
}

# auxiliary functions for Yule-Walker estimation
.matrix_acf = function(Y, st_vec)
{
  s = as.numeric(st_vec[1])
  t = as.numeric(st_vec[2])
  n = dim(Y)
  
  if ((s >= 0 && t >= 0) || (s < 0 && t < 0))
  {
    acf_out = .matrix_acf_positive(Y, abs(s), abs(t), n)
  } else if (s < 0 && t >= 0) {
    acf_out = .matrix_acf_negative(Y, s = -s, t, n)
  } else {
    acf_out = .matrix_acf_negative(Y, s, t = -t, n)
  }
  
  # unbiased estimator from Ha/Newton(1993)
  unbiased_factor = (n[1] - s)*(n[2] - t)/prod(n)
  return(acf_out/unbiased_factor)
}

.matrix_acf_positive = function(Y, s, t, n)
{
  Y0 = Y[1:(n[1] - s), 1:(n[2] - t)]
  Y1 = Y[(s + 1):n[1], (t + 1):n[2]]
  acf_out = sum(Y0 * Y1)
  return(acf_out)
}

.matrix_acf_negative = function(Y, s, t, n)
{
  Y0 = Y[1:(n[1] - s), (t + 1):n[2]]
  Y1 = Y[(s + 1):n[1], 1:(n[2] - t)]
  acf_out = sum(Y0 * Y1)
  return(acf_out)
}

#-------------------------Test for QARMA stationarity--------------------------#

# Test Function
qarma.statTest = function(ar)
{
  # make set.seed locally
  old_state = get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit(assign(".Random.seed", old_state, envir = .GlobalEnv, 
                 inherits = FALSE))
  
  ar = as.matrix(ar)
  ar[1, 1] = 1
  outValue = TRUE
  
  # compute reference sign at center (0, 0)
  signRef = sign(qarma.auxF(0, 0, ar))
  
  # check for random numbers inside unit circle
  
  set.seed(314)
  if (outValue == TRUE)
  {
    for (k in 1:1000)
    {
      # draw random point inside unit circle
      phi = stats::runif(2) * 2*pi
      rad = sqrt(stats::runif(2)) * 1.01
      x = rad * cos(phi)
      y = rad * sin(phi)
      z1c = complex(real = x[1], imaginary = y[1])
      z2c = complex(real = x[2], imaginary = y[2])
      
      signTest = sign(Re(qarma.auxF(z1c, z2c, ar)))
      if (signTest != signRef)
      {
        outValue = FALSE
        break()
      }
    }
  }
  return(outValue)
}

# Compute characteristic function of AR part
qarma.auxF = function(z1c, z2c, ar)
{
  ar_x = dim(ar)[1] - 1; ar_t = dim(ar)[2] - 1
  
  # set up vectors for z^i
  zx_vec = z1c^(0:ar_x)
  zt_vec = z2c^(0:ar_t)
  
  out = t(zx_vec) %*% ar %*% zt_vec
  
  return(out)
}


### Part B: Order Selection

#-------------------------------Model Selection--------------------------------#

# Order selection by GPAC (see Illig/Truong-Van (2006))
qarma.order.gpac = function(Y, order_max = list(ar = c(1, 1), ma = c(1, 1)))
{
  # set up vectors for ar and ma
  ar_vec_test = expand.grid(0:order_max$ar[1], 0:order_max$ar[2])
  ar_vec_test = ar_vec_test[2:dim(ar_vec_test)[1], ]
  names(ar_vec_test) = NULL
  ar_vec = expand.grid(0:(order_max$ar[1] + 1), 0:(order_max$ar[2] + 1))
  ar_vec = ar_vec[2:dim(ar_vec)[1], ]
  ar_vec = ar_vec[order(ar_vec[, 1], ar_vec[, 2]), ]
  names(ar_vec) = NULL
  #ar_vec = rbind(ar_vec, ar_vec[dim(ar_vec)[1], ] + 1)
  ma_vec = expand.grid(0:(order_max$ma[1]), 0:(order_max$ma[2]))
  ma_vec = ma_vec[order(ma_vec[, 1], ma_vec[, 2]), ]
  names(ma_vec) = NULL
  #ma_vec = rbind(ma_vec, ma_vec[dim(ma_vec)[1], ] + 1)
  
  # set up matrix for results of yw estimation
  gpac_mat = matrix(NA, nrow = dim(ma_vec)[1], ncol = dim(ar_vec)[1])
  colnames(gpac_mat) = paste0("(", ar_vec[, 1], ", ", ar_vec[, 2], ")")
  rownames(gpac_mat) = paste0("(", ma_vec[, 1], ", ", ma_vec[, 2], ")")
  ar_test_names = paste0("(", ar_vec_test[, 1], ", ", ar_vec_test[, 2], ")")
  
  # calculate GPAC
  for (i in 1:dim(ar_vec)[1])
  {
    for (j in 1:dim(ma_vec)[1])
    {
      ar = ar_vec[i, ]
      ma = ma_vec[j, ]
      gpac_mat[j, i] = qarma.yw_matrix(Y, ar_order = as.numeric(ar), 
                               ma_order = as.numeric(ma))[unlist(ar)[1] + 1,
                                                         unlist(ar)[2] + 1]
    }
  }
  
  # find zeros of GPAC
  nRow = dim(gpac_mat)[1]; nCol = dim(gpac_mat)[2]
  gpac_count = gpac_mat*0 + 1
  
  # values chosen by hand, might be improved  
  gpac_count[which(abs(gpac_mat) < max(stats::quantile(abs(gpac_mat), 0.35), 0.1)
                   , arr.ind = TRUE)] = 0
  
  count_zeros = 1:nRow*0
  
  ar_out = vector(mode = "numeric")
  ma_out = vector(mode = "numeric")
  # score = vector(mode = "numeric")
  
  for (i in 1:nRow)
  {
    for (j in 1:(prod(order_max$ar + 1) - 1))
    {
      ar = as.vector(ar_vec_test[j, ], mode = "numeric")
      jCol = which(ar_test_names[j] == colnames(gpac_count))
      colIndex = which(as.logical(vapply(ar[1], '<=', 
                                  logical(length(ar_vec[, 1])), ar_vec[, 1]) *
                                  vapply(ar[2], '<=', 
                                  logical(length(ar_vec[, 2])), ar_vec[, 2])))
      gpac_count_sub = gpac_count[i, colIndex]
      if ((gpac_count_sub[1] == 1) & (sum(gpac_count_sub) == 1))
      {
        ar_out = rbind(ar_out, ar)
        ma_out = rbind(ma_out, ma = as.vector(ma_vec[i, ], mode = "numeric"))
        # score = rbind(score, sum(gpac_mat[i, colIndex]^2))
      }
      
    }
  }
  
  if (length(ar_out) != 0)
  {
    # ar_out = ar_out[which.min(score), ]
    # ma_out = ma_out[which.min(score), ]
    ar_out = ar_out[dim(ar_out)[1], ]
    ma_out = ma_out[dim(ma_out)[1], ]
  } else if (length(ar_out) == 0) {
    warning("No order selection by GPAC possible, iid. case is assumed")
    ar_out = c(0, 0)
    ma_out = c(0, 0)
  }
  
  if (length(ar_out) == 0) {
    warning("No order selection by GPAC possible, iid. case is assumed")
    ar_out = c(0, 0)
    ma_out = c(0, 0)
  }
  
  return(list(ar = ar_out, ma = ma_out))
}

# Order selection by BIC
qarma.order.bic = function(Y, order_max = list(ar = c(1, 1), ma = c(1, 1)))
{
  n = prod(dim(Y))
  ar_matrix = expand.grid(0:order_max$ar[1], 0:order_max$ar[2])
  names(ar_matrix) = NULL
  ma_matrix = expand.grid(0:order_max$ma[1], 0:order_max$ma[2])
  names(ma_matrix) = NULL
  bic_matrix = matrix(NA, nrow = dim(ar_matrix)[1], ncol = dim(ma_matrix)[1])
  
  for (i in 1:dim(ar_matrix)[1])
  {
    for (j in 1:dim(ma_matrix)[1])
    {
      model_order = list(ar = as.numeric(ar_matrix[i, ]),
                         ma = as.numeric(ma_matrix[j, ]))
      qarma_model = qarma.est(Y, model_order = model_order)
      log_L = -n/2 * log(4*pi^2*qarma_model$model$sigma^2) -
        sum(qarma_model$innov^2)/(2*qarma_model$model$sigma^2)
      bic_matrix[i, j] = -2*log_L + sum(unlist(model_order)) * log(n)
    }
  }
  
  opt_index = which(bic_matrix == min(bic_matrix, na.rm = TRUE), arr.ind = TRUE)
  model_order_opt = list(ar = as.numeric(ar_matrix[opt_index[, 1], ]),
                         ma = as.numeric(ma_matrix[opt_index[, 2], ]))
  
  return(model_order_opt)
}
