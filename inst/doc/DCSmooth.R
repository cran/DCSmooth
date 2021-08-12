## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- echo = FALSE, warning = FALSE-------------------------------------------
library(DCSmooth)
set.seed(123)
summary(set.options(type = "KR"))

## ---- echo = FALSE, warning = FALSE-------------------------------------------
summary(set.options(type = "KR"))

## -----------------------------------------------------------------------------
# surface.dcs(y.norm1)
# surface.dcs(y.norm2)
# surface.dcs(y.norm3)

## -----------------------------------------------------------------------------
# surface.dcs(y.norm1)

y_iid = y.norm1 + rnorm(101^2)
# surface.dcs(y_iid)

## -----------------------------------------------------------------------------
opt_iid_KR = set.options(type = "KR")
dcs_iid_KR = dcs(y_iid, opt_iid_KR)

# print results
dcs_iid_KR
# print options used for DCS procedure
dcs_iid_KR$dcs_options

# plot regression surface
# surface.dcs(dcs_iid_KR, plot_choice = 2)

## -----------------------------------------------------------------------------
summary(dcs_iid_KR)

## -----------------------------------------------------------------------------
dcs_LP_iid = dcs(y_iid)
dcs_LP_iid
summary(dcs_LP_iid)

# plot regression surface
# surface.dcs(dcs_LP_iid, plot_choice = 2)

## -----------------------------------------------------------------------------
ar_mat = matrix(c(1, 0.4, -0.3, 0.2), nrow = 2, ncol = 2)
ma_mat = matrix(c(1, 0.2, 0.2, -0.5), nrow = 2, ncol = 2)
sigma  = 0.5
model_list = list(ar = ar_mat, ma = ma_mat, sigma = sigma)
sim_qarma = qarma.sim(n_x = 101, n_t = 101, model = model_list)

# QARMA observations
y_qarma = y.norm1 + sim_qarma$Y
# surface.dcs(y_qarma)

## -----------------------------------------------------------------------------
estim_qarma = qarma.est(sim_qarma$Y, 
                        model_order = list(ar = c(1, 1), ma = c(1, 1)))

estim_qarma$model

## -----------------------------------------------------------------------------
# QARMA((1, 1), (1, 1))
opt_qarma_1 = set.options(var_est = "qarma")
dcs_qarma_1 = dcs(y_qarma, opt_qarma_1)
dcs_qarma_1
dcs_qarma_1$var_model

# QARMA((1, 1), (0, 0))
order_list = list(ar = c(1, 1), ma = c(0, 0))
dcs_qarma_2 = dcs(y_qarma, opt_qarma_1, model_order = order_list)
dcs_qarma_2
dcs_qarma_2$var_model

## -----------------------------------------------------------------------------
opt_lm = set.options(var_est = "lm")
dcs_lm = dcs(y_iid, opt_lm)

dcs_lm
dcs_lm$var_model

## -----------------------------------------------------------------------------
opt_drv_1 = set.options(drv = c(1, 0))
dcs_drv_1 = dcs(y_iid, opt_drv_1)

dcs_drv_1
# surface.dcs(dcs_drv_1, plot_choice = 2)

## -----------------------------------------------------------------------------
opt_drv_2 = set.options(drv = c(1, 2))
dcs_drv_2 = dcs(y_iid, opt_drv_2)

dcs_drv_2
# surface.dcs(dcs_drv_2, plot_choice = 2)

