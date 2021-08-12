################################################################################
#                                                                              #
#       DCSmooth Package: Auxiliary Funcions for Bandwidth Selection           #
#                                                                              #
################################################################################

#----------------------Formula for optimal bandwidths--------------------------#

# Local Polynomial Regression
h.opt.LP = function(mxx, mtt, var_coef, n_sub, p_order, drv_vec, kern_fcn_x,
                    kern_fcn_t)
{
  # calculation of integrals
  i11 = sum(mxx^2)/n_sub; i22 = sum(mtt^2)/n_sub; i12 = sum(mxx * mtt)/n_sub
  
  # kernel constants (kernel Functions may also depend on p, drv)
  kernel_prop_1 = kernel.prop.LP(kern_fcn_x, p_order[1], drv_vec[1])
  kernel_prop_2 = kernel.prop.LP(kern_fcn_t, p_order[2], drv_vec[2])
  
  # relation factor gamma_21 (h_1 = gamma_21 * h_2)
  delta = (p_order - drv_vec)[1] # should be the same for both entries
  gamma_21 = (kernel_prop_1$mu/kernel_prop_2$mu)^(1/(delta + 1)) *
          ( i12/i11 * (drv_vec[1] - drv_vec[2])/(2*drv_vec[2] + 1) +
          sqrt(i12^2/i11^2 * (drv_vec[1] - drv_vec[2])^2/(2*drv_vec[2] + 1)^2 +
          i22/i11 * (2*drv_vec[1] + 1)/(2*drv_vec[2] + 1)) )^(1/(delta + 1))
  gamma_12 = 1/gamma_21
  
  # optimal bandwidths
  I1 = kernel_prop_1$mu^2 * i11 + kernel_prop_1$mu * kernel_prop_2$mu * 
            i12 * gamma_12^(delta + 1)
  hx_opt = (2*drv_vec[1] + 1)/(2*(delta + 1)) * (kernel_prop_1$R *
           kernel_prop_2$R * var_coef)/
           (n_sub * gamma_12^(2*drv_vec[1] + 1) * I1) 
  I2 = kernel_prop_2$mu^2 * i22 + kernel_prop_1$mu * kernel_prop_2$mu * 
           i12 * gamma_21^(delta + 1)
  ht_opt = (2*drv_vec[2] + 1)/(2*(delta + 1)) * (kernel_prop_1$R *
           kernel_prop_2$R * var_coef) /
           (n_sub * gamma_21^(2*drv_vec[2] + 1) * I2) 
  
  if(hx_opt < 0) { hx_opt = -hx_opt } # ensure positive bandwidth
  if(ht_opt < 0) { ht_opt = -ht_opt }
  
  hx_opt = hx_opt^(1/(2*(delta + sum(drv_vec) + 2)))
  ht_opt = ht_opt^(1/(2*(delta + sum(drv_vec) + 2)))
  
  return(c(hx_opt, ht_opt))
}

# Kernel Regression
h.opt.KR = function(mxx, mtt, var_coef, n, n_sub, kernel_prop_x, kernel_prop_t)
{
  i0x = integral.calc.KR(mxx, mtt, n_sub)[1]
  i0t = integral.calc.KR(mtt, mxx, n_sub)[1]
  
  R_tot = kernel_prop_x$R * kernel_prop_t$R
  mu_tot = kernel_prop_x$mu * kernel_prop_t$mu
  
  hx_opt = (R_tot * var_coef)/(n * mu_tot * i0x)
  ht_opt = (R_tot * var_coef)/(n * mu_tot * i0t)
  
  hx_opt = hx_opt^(1/6)
  ht_opt = ht_opt^(1/6)

  return(c(hx_opt, ht_opt))
}

#----------------------Integrals over mxx^2, mtt^2-----------------------------#

# calculation for LP done in h.opt.LP function

integral.calc.KR = function(m11, m22, n_sub)
{
  i11 = sum(m11 * m11)/n_sub
  i22 = sum(m22 * m22)/n_sub
  i12 = sum(m11 * m22)/n_sub

  i_out = (i11/i22)^0.75 * (sqrt(i11 * i22) + i12)
  
  return(c(i_out, i11/i22))   # second element of return actually not needed
}

#-------------------------Kernel property calculation--------------------------#

kernel.prop.LP = function(kernel_fcn, p, drv, n_int = 5000)
{
  u_seq  = seq(from = -1, to = 1, length.out = (2 * n_int + 1))
  np_matrix = np_matrix(kernel_fcn, p, n_int)
  add_weights = m_weights(np_matrix, u_seq, drv)
  
  val_R  = sum((add_weights^2 * 
                  kernel_fcn_use(u_seq, q = 1, kernel_fcn))^2)/n_int
  val_mu = sum((add_weights * kernel_fcn_use(u_seq, q = 1, kernel_fcn)) *
                u_seq^(p + 1)) / (n_int * factorial(p + 1))
  
  return(list(R = val_R, mu = val_mu))
}

kernel.prop.KR = function(kernel_fcn, n_int = 5000)
{
  u_seq  = seq(from = -1, to = 1, length.out = (2 * n_int + 1))
  val_R  = sum((kernel_fcn_use(u_seq, q = 1, kernel_fcn))^2)/n_int
  val_mu = sum((kernel_fcn_use(u_seq, q = 1, kernel_fcn)) * u_seq^2)/n_int
  
  return(list(R = val_R, mu = val_mu))
}

#-----------------bandwidth inflation for derivative estimations---------------#

inflation.LP = function(h, dcs_options, n_x, n_t)
{
  if (dcs_options$IPI_options$infl_exp[1] == "auto")
  {
    infl_exp = c(0, 0)
    infl_exp[1] = (dcs_options$p_order[1] + dcs_options$drv[2] + 2)/
                  (sum(dcs_options$p_order) + 3)
    infl_exp[2] = (dcs_options$p_order[2] + dcs_options$drv[1] + 2)/
                  (sum(dcs_options$p_order) + 3)
  } else {
    infl_exp = dcs_options$IPI_options$infl_exp
  }
  
  infl_par = dcs_options$IPI_options$infl_par
  h_infl_xx = c(infl_par[1] * h[1]^infl_exp[1], infl_par[2] * h[2]^infl_exp[1])
  h_infl_tt = c(infl_par[2] * h[1]^infl_exp[2], infl_par[1] * h[2]^infl_exp[2])
  
  # minimum values for LPR
  h_xx = pmax(h_infl_xx, (dcs_options$p_order + 2)/(c(n_x, n_t) - 1))
  h_tt = pmax(h_infl_tt, (dcs_options$p_order + 2)/(c(n_x, n_t) - 1))
  
  return(list(h_xx = h_xx, h_tt = h_tt))
}

deflation.LP = function(h, dcs_options, n_x, n_t)
{
  defl_exp = (dcs_options$p_order[1] + dcs_options$drv[2] + 2)/
              (dcs_options$p_order[1] - dcs_options$drv[1] + 2)
  
  h_defl = c(h[1]^defl_exp, h[2]^defl_exp)
  
  # minimum values for LPR
  h_defl = pmax(h_defl, c((dcs_options$p_order[1] + 1)/(n_x - 1), 
                      (dcs_options$p_order[2] + 1)/(n_t - 1)))
  
  return(h_defl)
}

inflation.KR = function(h, n, IPI_options)
{
  h_infl_xx = c(IPI_options$infl_par[1] * h[1]^IPI_options$infl_exp[1],
              IPI_options$infl_par[2] * h[2]^IPI_options$infl_exp[2])
  h_infl_tt = c(IPI_options$infl_par[2] * h[1]^IPI_options$infl_exp[2],
              IPI_options$infl_par[1] * h[2]^IPI_options$infl_exp[1])
  
  # ensure no bandwidth > 0.5 (or 0.45) as KR cannot handle estimation windows
  # > 1
  h_infl_xx = pmin(h_infl_xx, c(0.45, 0.45))
  h_infl_tt = pmin(h_infl_tt, c(0.45, 0.45))  
  
  return(list(h_xx = h_infl_xx, h_tt = h_infl_tt))
}