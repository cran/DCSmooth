sarma.cf = function(Y, model_order = list(ar = c(1, 1), ma = c(1, 1)))
{
  n_x = dim(Y)[1]; n_t = dim(Y)[2]
  Y_t = as.vector(t(Y))
  arma_t = stats::arima(Y_t, order = c(model_order$ar[2], 0, model_order$ma[2]),
                 include.mean = FALSE)
  Y_x = as.vector(matrix(arma_t$residuals, nrow = n_x, ncol = n_t,
                         byrow = TRUE))
  arma_x = stats::arima(Y_x, order = c(model_order$ar[1], 0, model_order$ma[1]),
                 include.mean = FALSE)
  
  ar_coefs_x = arma_x$coef[1:model_order$ar[1]]
  ar_coefs_t = arma_t$coef[1:model_order$ar[2]]
  ma_coefs_x = arma_x$coef[(model_order$ar[1] + 1):length(arma_x$coef)]
  ma_coefs_t = arma_t$coef[(model_order$ar[2] + 1):length(arma_t$coef)]
  
  ar_mat = c(1, -ar_coefs_x) %*% t(c(1, -ar_coefs_t))
  ar_mat = ar_mat[1:(model_order$ar[1] + 1), 1:(model_order$ar[2] + 1)]
  ma_mat = c(1, ma_coefs_x) %*% t(c(1, ma_coefs_t))
  ma_mat = ma_mat[1:(model_order$ma[1] + 1), 1:(model_order$ma[2] + 1)]
  sigma = sqrt(arma_x$sigma2)
  
  cf = sum(ma_mat)^2/sum(ar_mat)^2 * sigma^2
  model = list(ar = ar_mat, ma = ma_mat, sigma = sigma)
  
  return(list(cf = cf, model = model))
}