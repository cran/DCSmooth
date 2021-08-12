################################################################################
#                                                                              #
#                           Test for set.options()                             #
#                                                                              #
################################################################################

### Test exception handling (positives)
context("set.options() exception handling (positives)")

test_that("set.options() gets \"type\" correctly.", {
  expect_equal(set.options()$type, "LP")
  expect_equal(set.options(type = "LP")$type, "LP")
  expect_equal(set.options(type = "KR")$type, "KR")
  expect_error(set.options(type = "test"), "Unsupported values in argument")
})

test_that("set.options() gets \"kerns\" correctly.",{
  for (test in seq_along(length(dcs_list_kernels)))
  {
    expect_equal(set.options(kerns = c(dcs_list_kernels[test], 
                                       dcs_list_kernels[test]))$kerns,
                 c(dcs_list_kernels[test], dcs_list_kernels[test]))
    expect_equal(set.options(kerns = dcs_list_kernels[test])$kerns,
                 c(dcs_list_kernels[test], dcs_list_kernels[test]))
  }
})

test_that("set.options() gets \"drv\" correctly.",{
  expect_equal(set.options(drv = c(0, 0))$drv, c(0, 0))
  expect_equal(set.options(drv = c(2, 1))$drv, c(2, 1))
})

test_that("set.options() gets \"var_est\" correctly.",{
  expect_equal(set.options()$var_est, "iid")
  expect_equal(set.options(var_est = "qarma")$var_est, "qarma")
  expect_equal(set.options(var_est = "qarma_gpac")$var_est, "qarma_gpac")
  expect_equal(set.options(var_est = "qarma_bic")$var_est, "qarma_bic")
  expect_equal(set.options(var_est = "sarma")$var_est, "sarma")
  expect_equal(set.options(var_est = "lm")$var_est, "lm")
})

test_that("set.options() gets \"IPI_options\" correctly.",{
  IPI_LP = list(infl_exp = c("auto", " "), infl_par = c(1, 1), 
                delta = c(0.05, 0.05), const_window = FALSE)
  IPI_KR = list(infl_exp = c(0.5, 0.5), infl_par = c(2, 1), 
                delta = c(0.05, 0.05), const_window = FALSE)
  expect_equal(set.options(type = "LP")$IPI_options, IPI_LP)
  expect_equal(set.options(type = "KR")$IPI_options, IPI_KR)
})

### Test exception handling (negatives)
context("set.options() exception handling (positives)")

test_that("set.options() gives correct error messages for problems", {
  expect_error(set.options(type = "test"), "Unsupported values in argument")
  expect_error(set.options(type = c("KR", "LP")),
               "Unsupported values in argument")
})
