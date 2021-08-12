################################################################################
#                                                                              #
#                DCSmooth Package: estimation of cf coefficients               #
#                                                                              #
################################################################################

# Different methods for estimation of the cf coefficient for bandwidth
# selection are joined in this function

# Y should be the residuals, e.g. Y - YSmth in the estimation codes
cf.estimation = function(Y, dcs_options, add_options)
{
  if (dcs_options$var_est == "iid")
  {
    cf_est = stats::sd(Y)^2
    var_model = list(sigma = stats::sd(Y), stnry = TRUE)
  } else if (dcs_options$var_est == "qarma") {
    qarma = qarma.cf(Y, model_order = add_options$model_order)
    cf_est = qarma$cf
    var_model = qarma$qarma_model$model
    var_model$stnry = qarma$qarma_model$stnry
  } else if (dcs_options$var_est == "qarma_gpac") {
    model_order = qarma.order.gpac(Y, order_max = add_options$order_max)
    qarma = qarma.cf(Y, model_order = model_order)
    cf_est = qarma$cf
    var_model = qarma$qarma_model$model
    var_model$stnry = qarma$qarma_model$stnry
  } else if (dcs_options$var_est == "qarma_bic") {
    model_order = qarma.order.bic(Y, order_max = add_options$order_max)
    qarma = qarma.cf(Y, model_order = model_order)
    cf_est = qarma$cf
    var_model = qarma$qarma_model$model
    var_model$stnry = qarma$qarma_model$stnry
  } else if (dcs_options$var_est == "np") {
    cf_est = specDens(Y, omega = c(0, 0))$cf
    var_model = NA
    var_model$stnry = TRUE
  } else if (dcs_options$var_est == "sarma") {
    sarma = sarma.cf(Y, model_order = add_options$model_order)
    cf_est = sarma$cf
    var_model = sarma$model
    var_model$stnry = TRUE
  }
  
  if (is.na(cf_est))
  {
    stop("Non-finite variance estimated. Check input matrix Y.")
  }
  
  return(list(cf_est = cf_est, var_model = var_model))
}