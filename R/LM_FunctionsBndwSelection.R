################################################################################
#                                                                              #
#       DCSmooth Package: Auxiliary Functions for Long Memory Bandwidths       #
#                                                                              #
################################################################################

#----------------------Formula for optimal bandwidths--------------------------#

h.opt.LM = function(mxx, mtt, var_coef, var_model, n_sub, p_order, drv_vec,
                    n_x, n_t, kernel_x, kernel_t)
{
  d_vec = var_model$d_vec
  # calculation of integrals
  i11 = sum(mxx^2)/n_sub; i22 = sum(mtt^2)/n_sub; i12 = sum(mxx * mtt)/n_sub

  # kernel constants (kernel Functions may also depend on p, drv)
  kernel_prop_x = kernel.prop.LM(kernel_x, p_order[1], drv_vec[1], d_vec[1])
  kernel_prop_t = kernel.prop.LM(kernel_t, p_order[2], drv_vec[2], d_vec[2])

  # compute additional values
  cb = h.coef.cb(i11, i22, i12, kernel_prop_x, kernel_prop_t, drv_vec, p_order,
                 d_vec)
  cn = n_t/n_x
  C1 = 4 * var_coef * var_model$sigma^2 *
       kernel_prop_x$V * kernel_prop_t$V
  delta = p_order[1] + 1 - drv_vec[1]

  C1A = (1 - 2*d_vec[1] + 2*drv_vec[1]) * C1 * cn^(diff(d_vec)) /
        (2 * delta * kernel_prop_x$mu * cb^(1 - 2*d_vec[2] + 2*drv_vec[2]) *
        (kernel_prop_x$mu * i11 + kernel_prop_t$mu * i12 * cb^delta))
  C2A = (1 - 2*d_vec[2] + 2*drv_vec[2]) * C1 * cn^(diff(d_vec)) /
        (2 * delta * kernel_prop_t$mu * cb^(1 - 2*d_vec[2] + 2*drv_vec[2]) *
        (kernel_prop_x$mu * i12 * cb^(-delta) + kernel_prop_t$mu * i22))

  b1A = (C1A / n_sub^(1 - sum(d_vec)))^(1/
                        (2*(1 - sum(d_vec) + delta + sum(drv_vec))))
  b2A = (C2A / n_sub^(1 - sum(d_vec)))^(1/
                        (2*(1 - sum(d_vec) + delta + sum(drv_vec))))

  return(c(b1A, b2A))
}

cf.estimation.LM = function(R_mat, model_order =
                              list(ar = c(0, 0), ma = c(0, 0)))
{
  sfarima = sfarima.est(R_mat, model_order = model_order)

  cf_est = sum(sfarima$ma)^2/sum(sfarima$ar)^2 * sfarima$sigma^2
  var_model = list(d_vec = sfarima$d_vec, ar = sfarima$ar, ma = sfarima$ma,
                   sigma = sfarima$sigma, stnry = TRUE)

  return(list(cf_est = cf_est, var_model = var_model))
}

#------------------------Formula for coefficient cb----------------------------#

h.coef.cb = function(i11, i22, i12, kernel_prop_1, kernel_prop_2, drv_vec,
                     p_order, d_vec)
{
  denom_value = (1 - 2*d_vec[1] + 2*drv_vec[1]) * kernel_prop_2$mu^2 * i22
  sec_term = (diff(d_vec) - diff(drv_vec)) * kernel_prop_1$mu *
              kernel_prop_2$mu * i12
  sqrt_value = sec_term^2 + denom_value *
                (1 - 2*d_vec[2] + 2*drv_vec[2]) * kernel_prop_1$mu^2 * i11

  delta = p_order[1] + 1 - drv_vec[1]
  cb_return = 1/denom_value * (sqrt(sqrt_value) - sec_term)

  return(cb_return^(1/delta))
}

#------------------------------Long-Memory KDF---------------------------------#

kernel.prop.LM = function(kernel_fcn, p, drv, d, n_int = 5000)
{
  mu = 2 # change later
  u_seq  = seq(from = -1, to = 1, length.out = n_int)

  val_V  = gamma(1 - 2*d) * sin(pi * d) * lookup$p1p3_lookup[mu + 1, 2][[1]](d)

  val_mu = sum((kernel_fcn_use(u_seq, q = 1, kernel_fcn)) *
                 u_seq^(p + 1)) / (n_int * factorial(p + 1))

  return(list(V = val_V, mu = val_mu))
}

# kernel.tlm = function(l, m, delta)
# {
#   i = 0:m
#   j = matrix(0, nrow = m + 1, ncol = l + m + 1)
#   for (k in 1:(m + 1))
#   {
#     j[k, 1:(l + m - i[k] + 1)] = 0:(l + m - i[k])
#   }
#
#   sum_j = apply(j, 1, kernel.tlm.sum2, l = l, m = m, delta = delta)
#   sum_i = 2*sum(choose(m, i) * 1/(2*delta + i) * sum_j)
# }
#
# kernel.tlm.sum2 = function(j, l, m, delta)
# {
#   i = length(j) - max(j) - 1
#   sum_terms = (-1)^j * choose(l + m - i, j) * 2^(2*delta + j + i + 1) /
#                   (2*delta + j + i + 1)
#   sum_return = sum(sum_terms[1:(length(j) - i)])
#
#   return(sum_return)
# }


#----------------------------SFARIMA estimation--------------------------------#

# sfarima.est = function(R_mat, model_order = list(ar = c(1, 1), ma = c(1, 1)))
# {
#   # coefficient estimation
#   n_x = dim(R_mat)[1]; n_t = dim(R_mat)[2]
#   R_t = as.vector(t(R_mat))
#   sfarima_t = fracdiff::fracdiff(R_t, nar = model_order$ar[2],
#                        nma = model_order$ma[2])
#   R_x = as.vector(matrix(sfarima_t$residuals, nrow = n_x, ncol = n_t,
#                          byrow = TRUE))
#   sfarima_x = fracdiff::fracdiff(R_x, nar = model_order$ar[1],
#                        nma = model_order$ma[1])
#
#   # computing matrices
#   ar_mat = c(1, sfarima_x$ar) %*% t(c(1, sfarima_t$ar))
#   ma_mat = c(1, sfarima_x$ma) %*% t(c(1, sfarima_t$ma))
#
#   # prepare output
#   sfarima_return = list(d_vec = c(sfarima_x$d, sfarima_t$d),
#                         ar = ar_mat, ma = ma_mat,
#                         sigma = sfarima_x$sigma)
#
#   return(sfarima_return)
# }


#----------------------------SFARIMA estimation alternativ (SL)--------------------------------#

sfarima.est = function(R_mat, model_order = list(ar = c(1, 1), ma = c(1, 1)),
                       K1 = 200, K2 = 200)
{
  #browser()
  n_x = dim(R_mat)[1]; n_t = dim(R_mat)[2]
  K1 = min(n_x, K1) ; K2 = min(n_t, K2)
  nar1 = model_order$ar[1]
  nar2 = model_order$ar[2]
  nma1 = model_order$ma[1]
  nma2 = model_order$ma[2]
  narma1 = model_order$ar[1] + model_order$ma[1]
  narma2 = model_order$ar[2] + model_order$ma[2]

  theta0 = c(0, rep(0, narma1), 0, rep(0, narma2))
  k.max = 3
  for(k in 1:k.max){
    result = stats::optim(theta0, sfarima_rss, R_mat = R_mat, model_order = model_order,
                   K1 = K1, K2 = K2)
    theta = result$par
    rss = result$value
    if (max(abs(theta - theta0)) < 5 / 1000) {
      k = k
      break
    }
    else {
      theta0 = theta
    }
  }

  sig_eta = sqrt(rss / (n_x * n_t))
  #sig_eta = sqrt(rss / (K1 * K2))

  ar_mat = c(1, theta[-1][0:nar1]) %*% t(c(1, theta[-(1:(2 + narma1))][0:nar2]))

  ma_mat = c(1, theta[-(1 :(1 + nar1))][0:nma1]) %*%
    t(c(1, theta[-(1:(2 + narma1 + nar2))][0:nma2]))

  sfarima_return = list(d_vec = c(theta[1], theta[narma1 + 2]),
                        ar = ar_mat, ma = ma_mat,
                        sigma = sig_eta)
  return(sfarima_return)
}

# #####################################################################################
sfarima_rss = function(theta,
                       R_mat,
                       model_order = list(ar = c(1, 1), ma = c(1, 1)),
                       K1 = 200, K2 = 200) {
  nar1 = model_order$ar[1]
  nar2 = model_order$ar[2]
  nma1 = model_order$ma[1]
  nma2 = model_order$ma[2]
  narma1 = model_order$ar[1] + model_order$ma[1]
  narma2 = model_order$ar[2] + model_order$ma[2]

  d1 = theta[1]
  phi1 = theta[-1][0:nar1]
  psi1 = theta[-(1 :(1 + nar1))][0:nma1]
  d2 = theta[narma1 + 2]
  phi2 = theta[-(1:(2 + narma1))][0:nar2]
  psi2 = theta[-(1:(2 + narma1 + nar2))][0:nma2]

  n1 = nrow(R_mat)
  n2 = ncol(R_mat)
  K1 = min(K1, n1)
  K2 = min(K2, n2)

  bk1 = arcoef(ar = phi1, ma = psi1, d = d1, k = K1)

  bk2 = arcoef(ar = phi2, ma = psi2, d = d2, k = K2)


  ######### calculating the RSS

  RSS = 0
  #e.est = xi.est = matrix(0, K1, K2)
  e.est = xi.est = matrix(0, n1, n2)

  ### 1st stage

  bk2 = t(bk2)
  R_mat = t(R_mat)

  # for(j in 1:K2) {
  #   xi.est[, j] = bk2[1:j] %*% R_mat[j:1, 1:K1]
  # }

  for(j in 1:n2) {
    if (j <= K2) {
      xi.est[, j] = bk2[1:j] %*% R_mat[j:1, ]
    }
    else {
      xi.est[, j] = bk2 %*% R_mat[j:(j - K2), ]
    }
  }

  ### 2nd stage

  bk1 = as.matrix(t(bk1))

  # for(i in 1:K1) {
  #   e.est[i, ] = bk1[1:i] %*% xi.est[i:1, ]
  # }

  for(i in 1:n1) {
    if (i <= K1) {
      e.est[i, ] = bk1[1:i] %*% xi.est[i:1, ]
    }
    else {
      e.est[i, ] = bk1 %*% xi.est[i:(i - K1), ]
    }

  }

  rss = sum(e.est^2)
  return(rss)

}

# ############################################################################

arcoef <- function(ar = 0, ma = 0, d = 0, k = 50) {
  p = length(ar[ar != 0])
  q = length(ma[ma != 0])
  if (p == 0) {
    ar = 0
  }
  if (q == 0) {
    ma = 0
  }

  ar.coef = c(1, -ar, rep(0, k - p))
  arma.coef = (1:(k + 1)) * 0
  arma.coef[1] = 1

  if (p | q > 0) {
    for (i in 2:(k + 1)) {
      if ((i - q) < 1) {
        arma.coef[i] = -sum(ma[1:(q - abs(i - q) - 1)] * arma.coef[(i - 1):1]) - (-ar.coef[i])
      } else {
        arma.coef[i] = -sum(ma[1:q] * arma.coef[(i - 1):(i - q)]) - (-ar.coef[i])
      }
    }
  } else {
    arma.coef**0
  }

  d.coef = choose(d, 0:k) * ((-1)^(0:k))

  coef.all = (1:(k + 1)) * 0
  for (j in 1:(k + 1)) {
    coef.all[j] = sum(d.coef[1:j] * arma.coef[j:1])
  }
  drop(coef.all)
}


