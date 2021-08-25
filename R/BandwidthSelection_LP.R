################################################################################
#                                                                              #
#          DCSmooth Package: R-Functions for LP Bandwidth Selection            #
#                                                                              #
################################################################################

LP.bndw = function(Y, dcs_options, add_options)
{
  n_x = dim(Y)[1]; n_t = dim(Y)[2]
  n = n_x * n_t                   # total number of observations
  p_order = dcs_options$p_order
  drv_vec = dcs_options$drv
  
  # set variables for weight type
  kern_type_vec = sub("^([[:alpha:]]*).*", "\\1", dcs_options$kern)
                                  # extract weighting type
  mu_vec = as.numeric(substr(dcs_options$kern, nchar(dcs_options$kern) - 1,
                  nchar(dcs_options$kern) - 1))
                                  # extract kernel parameter mu
  weight_x = weight_fcn_assign(kern_type_vec[1])
  weight_t = weight_fcn_assign(kern_type_vec[2])
                                  # weight functions for LP-regression
  kernel_x = kernel_fcn_assign(dcs_options$kerns[1])
  kernel_t = kernel_fcn_assign(dcs_options$kerns[2])
                                  # (equivalent) kernel function for calculation
                                  # of kernel parameters
  
  h_opt = pmax(c(0.1, 0.1), (dcs_options$p_order + 2)/(c(n_x, n_t) - 1))
                                  # initial (arbitrary) values for h_0

  iterate = TRUE                  # iteration indicator
  iteration_count   = 0
  while(iterate)                  # loop for IPI
  {
    iteration_count = iteration_count + 1
    h_opt_temp = h_opt[1:2]       # store old bandwidths for breaking condition
    h_infl = inflation.LP(h_opt_temp, dcs_options, n_x, n_t)  
                                  # inflation of bandwidths for drv estimation
    h_defl = deflation.LP(h_opt_temp, dcs_options, n_x, n_t)
                                  # deflation if derivatives are estimated from
                                  # function
   
    # constant window width only reasonable for estimation of derivatives
    if (dcs_options$IPI_options$const_window == TRUE)
    {
      # smoothing of Y for variance factor estimation
      Y_smth = LP_dcs_const0_BMod(yMat = Y, hVec = h_defl, polyOrderVec
                                  = c(1, 1), drvVec = c(0, 0), muVec = mu_vec,
                                  weightFcnPtr_x = weight_x,
                                  weightFcnPtr_t = weight_t)
      # smoothing of derivatives m(2,0) and m(0,2)
      mxx = LP_dcs_const1_BMod(yMat = Y, hVec = h_infl$h_xx, polyOrderVec =
                          c(2*p_order[1] - drv_vec[1] + 1, p_order[2]), drvVec =
                          c(p_order[1] + 1, drv_vec[2]), muVec = mu_vec,
                          weightFcnPtr_x = weight_x, weightFcnPtr_t = weight_t)
      mtt = LP_dcs_const1_BMod(yMat = Y, hVec = h_infl$h_tt, polyOrderVec =
                          c(p_order[1], 2*p_order[2] - drv_vec[2] + 1), drvVec =
                          c(drv_vec[1], p_order[2] + 1), muVec = mu_vec,
                          weightFcnPtr_x = weight_x, weightFcnPtr_t = weight_t)
    } else {
      # smoothing of Y for variance factor estimation
      Y_smth = LP_dcs_const0_BMod(yMat = Y, hVec = h_defl, polyOrderVec
                                  = c(1, 1), drvVec = c(0, 0), muVec = mu_vec,
                                  weightFcnPtr_x = weight_x,
                                  weightFcnPtr_t = weight_t)

      # smoothing of derivatives m(2,0) and m(0,2)
      mxx = LP_dcs_const0_BMod(yMat = Y, hVec = h_infl$h_xx, polyOrderVec =
                          c(2*p_order[1] - drv_vec[1] + 1, p_order[2]), drvVec =
                          c(p_order[1] + 1, drv_vec[2]),  muVec = mu_vec,
                          weightFcnPtr_x = weight_x, weightFcnPtr_t = weight_t)
      mtt = LP_dcs_const0_BMod(yMat = Y, hVec = h_infl$h_tt, polyOrderVec =
                         c(p_order[1], 2*p_order[2] - drv_vec[2] + 1), drvVec =
                         c(drv_vec[1], p_order[2] + 1), muVec = mu_vec,
                         weightFcnPtr_x = weight_x, weightFcnPtr_t = weight_t)
    }
      
    # shrink mxx, mtt from boundaries if delta > 0
    if (dcs_options$IPI_options$delta[1] != 0 ||
        dcs_options$IPI_options$delta[2] != 0)
    {
      shrink_x = ceiling(dcs_options$IPI_options$delta[1] * n_x):
                         (n_x - floor(dcs_options$IPI_options$delta[1] * n_x))
      shrink_t = ceiling(dcs_options$IPI_options$delta[2] * n_t):
                         (n_t - floor(dcs_options$IPI_options$delta[2] * n_t))
      
      mxx = mxx[shrink_x, shrink_t]
      mtt = mtt[shrink_x, shrink_t]
      n_sub = dim(mxx)[1] * dim(mxx)[2]   # number of used observations
    } else {
      n_sub = n                           # all observations are used
    }
    
    if (dcs_options$var_est == "lm") ### Long-memory estimation
    {
      # calculate variance factor
      var_est = suppressWarnings(cf.estimation.LM(Y - Y_smth,
                                                  add_options$model_order))
      var_coef = var_est$cf_est
      var_model = var_est$var_model

      # calculate optimal bandwidths for next step
      h_opt = h.opt.LM(mxx, mtt, var_coef, var_model, n_sub, dcs_options, 
                       n_x, n_t)
    } else {                         ### Short-memory or iid. estimation
      # calculate variance factor
      var_est   = suppressWarnings(cf.estimation(Y - Y_smth, dcs_options,
                                                 add_options))
      var_coef  = var_est$cf_est
      var_model = var_est$var_model
      
      # calculate optimal bandwidths for next step
      h_opt = h.opt.LP(mxx, mtt, var_coef, n_sub, p_order, drv_vec, kernel_x,
                       kernel_t)
    }

    # break condition
    if( all(!is.nan(h_opt)) &&
        ((h_opt[1]/h_opt_temp[1] - 1 < 0.001) &&
         (h_opt[2]/h_opt_temp[2] - 1 < 0.001) &&
         (iteration_count > 2)) || iteration_count > 15)
    {
      iterate = FALSE
    }
  }
  return(list(h_opt = h_opt, iterations = iteration_count, var_coef = var_coef,
              var_model = var_model))
}