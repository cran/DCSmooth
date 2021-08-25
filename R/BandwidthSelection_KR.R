###############################################################################
#                                                                             #
#         DCSmooth Package: R-Functions for KR Bandwidth Selection            #
#                                                                             #
###############################################################################

#------------------Function for the optimal bandwidth via IPI-----------------#

KR.bndw = function(Y, dcs_options, add_options)
{
  n_x = dim(Y)[1]; n_t = dim(Y)[2]
  n  = n_x * n_t                            # total number of observations
  
  # set variables for weight type
  kern_type_vec = sub("^([[:alpha:]]*).*", "\\1", dcs_options$kern)
                                                # extract weighting type
  mu_vec = as.numeric(substr(dcs_options$kern, nchar(dcs_options$kern) - 1,
                             nchar(dcs_options$kern) - 1))
                                                # extract kernel parameter mu
  
  # set kernel Function to use in optimization
  kernel_x = kernel_fcn_assign(dcs_options$kerns[1])
  kernel_t = kernel_fcn_assign(dcs_options$kerns[2])

  kernel_prop_x = kernel.prop.KR(kernel_x)  # kernel properties R and mu_2
  kernel_prop_t = kernel.prop.KR(kernel_t)
  
  # TODO: more flexibility here
  kern_fcn_0 = kernel_fcn_assign("MW_220")  # kernel for regression surface
  kern_fcn_2 = kernel_fcn_assign("MW_422")  # kernel for 2nd derivative
  
  h_opt = c(0.1, 0.1)                       # initial (arbitrary) values for h_0
  
  iterate = TRUE                            # iteration indicator
  iteration_count = 0
  while(iterate)                            # loop for IPI
  {
    iteration_count = iteration_count + 1
    h_opt_temp   = pmin(h_opt[1:2], c(0.45, 0.45)) 
                          # KR can't handle too large bandwidths
    h_infl  = inflation.KR(h_opt_temp, c(n_x, n_t), dcs_options)
                          # inflation of bndws for estimation of derivatives
    
    # constant bandwidth only reasonable for estimation of derivatives
    if (dcs_options$IPI_options$const_window == TRUE)
    {
      # pre-smoothing of the surface function m(0,0) for estimation of variance
      Y_smth = KR_dcs_const0(yMat = Y, hVec = h_opt_temp, 
                             drvVec = c(0, 0), kernFcnPtrX = kern_fcn_0,
                             kernFcnPtrT = kern_fcn_0)
      # smoothing of derivatives m(2,0) and m(0,2)
      mxx = KR_dcs_const1(yMat = Y, hVec = h_infl$h_xx, drvVec = c(2, 0),
                        kernFcnPtrX = kern_fcn_2, kernFcnPtrT = kern_fcn_0)
      mtt = KR_dcs_const1(yMat = Y, hVec = h_infl$h_tt, drvVec = c(0, 2),
                    kernFcnPtrX = kern_fcn_0, kernFcnPtrT = kern_fcn_2)
    } else if (dcs_options$IPI_options$const_window == FALSE) {
      # pre-smoothing of the surface function m(0,0) for estimation of variance
      Y_smth = KR_dcs_const0(yMat = Y, hVec = h_opt_temp, drvVec = c(0, 0),
                             kernFcnPtrX = kern_fcn_0, kernFcnPtrT = kern_fcn_0)
      # smoothing of derivatives m(2,0) and m(0,2)
      mxx = KR_dcs_const0(yMat = Y, hVec = h_infl$h_xx, drvVec = c(2, 0),
                          kernFcnPtrX = kern_fcn_2, kernFcnPtrT = kern_fcn_0)
      
      mtt = KR_dcs_const0(yMat = Y, hVec = h_infl$h_tt, drvVec = c(0, 2),
                          kernFcnPtrX = kern_fcn_0, kernFcnPtrT = kern_fcn_2)
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
      n_sub = dim(mxx)[1]*dim(mxx)[2]   # number of used observations
    } else {
      n_sub = n                         # all observations are used
    }
    
    ### Estimation of Variance Factor and Model ###
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
      h_opt = h.opt.KR(mxx, mtt, var_coef, n, n_sub, kernel_prop_x,
                       kernel_prop_t)
    }
    
    # break condition
    if( ((h_opt[1]/h_opt_temp[1] - 1 < 0.001) && (h_opt[2]/h_opt_temp[2] - 1 
        < 0.001) && (iteration_count > 3)) || (iteration_count > 15) )
    {
      iterate = FALSE
    }
  }
  return(list(h_opt = h_opt, iterations = iteration_count, var_coef = var_coef,
              var_model = var_model))
}