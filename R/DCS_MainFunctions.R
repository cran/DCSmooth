################################################################################
#                                                                              #
#                      DCSmooth Package: User Functions                        #
#                                                                              #
################################################################################

#------------------------Set Options via Function------------------------------#

#' Set Options for the DCS procedure
#' 
#' @param type either local polynomial regression (\code{"LP"}, the default) or 
#'  kernel regression (\code{"KR"}).
#' @param kerns a character vector of length 2 containing the identifier for the
#'  kernels to be used in kernel regression. Weighting functions in local
#'  polynomial regression are computed according to the identifier. Default value
#'  is \code{MW_220}, the Mueller-Wang kernel of order \eqn{(2, 2, 0)}.
#' @param drv A non-negative vector of length 2, containing the derivative
#'  orders to be estimated from the given data. The default is \code{c(0, 0)}. For
#'  LP-regression, polynomial order is selected as \eqn{(\nu_1 + 1, \nu_2 + 1)}.
#' @param var_est the method of estimating the variance coefficient \eqn{c_f}. 
#'  Currently available are \code{var_est = "iid"} (the default), where an iid. 
#'  error structure is modeled and \code{var_est = "qarma"}, where a spatial
#'  ARMA model is applied. Use \code{"qarma_gpac"} for automated order 
#'  selection.
#' @param ... Additional arguments passed to \code{set.options()}. This includes
#'  \code{IPI_options}, a list containing further options used by the iterative
#'  plug-in algorithm.
#' 
#' @section Details:
#' This function is used to set the options for bandwidth selection in the 
#' \code{dcs} function. Detailed information can be found in the vignette.
#' 
#' @return  An object of class \code{"dcs_options"}.
#' 
#' @seealso \code{\link{dcs}}
#' 
#' @export
#' 
#' @examples
#' # See vignette("DCSmooth") for examples and explanation
#' 
#' set.options()
#' 
#' myOpt <- set.options(type = "KR", var_est = "qarma")
#' y <- y.norm1 + matrix(rnorm(101^2), nrow = 101, ncol = 101)
#' dcs(y, dcs_options = myOpt)

set.options <- function(    # inside function with default values in arguments
  # standard options
  type      = "LP",        # either "LP" for local polynomial regression or
  # "KR" for kernel regression
  kerns     = c("MW_220", "MW_220"), # choose a kernel function
  drv       = c(0, 0),
  var_est   = "iid",
  
  # advanced options in ellipsis
  ...
)
{
  # get ellipsis
  IPI_options <- list(...)
  
  # check if inputs are vectors
  if (length(kerns) == 1) { kerns <- c(kerns, kerns) }
  if (length(drv) == 1) { drv <- c(drv, drv) }
  
  exception.check.options.input(type, kerns, drv, var_est, IPI_options)
  
  if (var_est == "lm" && type == "KR")
  {
    warning("Long-Memory bandwidth selection only supported for local ", 
            "polynomial regression. Type set to \"LP\".")
    type <- "LP"
  }
  if (var_est == "lm")
  {
    message("Estimation under long-memory errors (SFARIMA) is currently in ", 
            "experimental state.")
  }
  
  # Select options according to type ("LP", "KR")
  if (type == "LP")
  {
    p_order <- drv + 1
    if (!exists("infl_exp", IPI_options))
    {
      IPI_options$infl_exp <- c("auto", " ")
    }
    if (!exists("infl_par", IPI_options)) { IPI_options$infl_par <- c(1, 1) }
    if (!exists("delta", IPI_options)) { IPI_options$delta <- c(0.05, 0.05) }
    if (!exists("const_window", IPI_options))
    {
      IPI_options$const_window <- FALSE
    }
    
    options_list <- list(type = type, kerns = kerns, p_order = p_order,
                        drv = drv, var_est = var_est, IPI_options = IPI_options)
  } else if (type == "KR") {
    p_order <- NA
    if (!exists("infl_exp", IPI_options))
    {
      IPI_options$infl_exp <- c(0.5, 0.5)
    }
    if (!exists("infl_par", IPI_options)) { IPI_options$infl_par <- c(2, 1) }
    if (!exists("delta", IPI_options)) { IPI_options$delta <- c(0.05, 0.05) }
    if (!exists("const_window", IPI_options))
    {
      IPI_options$const_window <- FALSE
    }
    
    options_list <- list(type = type, kerns = kerns, p_order = p_order,
                        drv = drv, var_est = var_est, IPI_options = IPI_options)
  } else {
    stop("Unknown type \"", type, "\"")
  }
  
  # apply class to output object
  class(options_list) <- "dcs_options"
  
  exception.check.options(options_list)
  
  return(options_list)
}

#--------------Function for smoothing and bandwidth estimation----------------#

#' Nonparametric Double Conditional Smoothing for 2D Surfaces
#' 
#' \code{DCSmooth} provides a double conditional nonparametric smoothing of the
#' expectation surface of a functional time series or a random field on a
#' lattice. Bandwidth selection is done via an iterative plug-in method.
#' 
#' @param Y A numeric matrix that contains the observations of the random field
#'   or functional time-series.
#' @param dcs_options An object of class \code{"dcs_options"}, specifying the
#'  parameters for the smoothing and bandwidth selection procedure.
#' @param h Bandwidth for smoothing the observations in \code{Y}. Can be a
#'  two-valued numerical vector with bandwidths in row- and column-direction.
#'  If the value is \code{"auto"} (the default), bandwidth selection will be 
#'  carried out by the iterative plug-in algorithm.
#' @param ... Additional arguments passed to \code{dcs}. These might include
#'  numerical vectors \code{X} and/or \code{T} containing the exogenous
#'  covariates with respect to the rows and columns. If \code{var_est} is any of
#'  \code{"qarma"}, \code{"sarma"} or \code{"lm"}, \code{model_order} is an 
#'  optional list containing the two-dimensional model orders in the form 
#'  \code{list(ar = c(1, 1), ma = c(1, 1)}. If \code{var_est} is any of
#'  \code{"qarma_gpac"} or \code{"qarma_bic"}, \code{order_max} an optional list
#'  containing the two dimensional maximum orders for the order selection of the
#'  QARMA estimation in the form \code{list(ar = c(1, 1), ma = c(1, 1)}.
#' 
#' @return \code{DCSmooth} returns an object of class "dcs", including
#'  \tabular{ll}{
#'  \code{Y} \tab matrix of original observations. \cr
#'  \code{X, T} \tab vectors of covariates over rows (\code{X}) and columns 
#'   (\code{T}). \cr
#'  \code{M} \tab resulting matrix of smoothed values. \cr
#'  \code{R} \tab matrix of residuals of estimation, \eqn{Y - M}. \cr
#'  \code{h} \tab optimized or given bandwidths. \cr
#'  \code{c_f} \tab estimated variance coefficient. \cr
#'  \code{dcs_options} \tab an object of class \code{cds_options} containing the
#'   initial options of the dcs procedure. \cr
#'  \code{iterations} \tab number of iterations of the IPI-procedure. \cr
#'  \code{time_used} \tab time spend searching for optimal bandwidths (not
#'   overall runtime of the function). \cr
#'  \code{qarma} \tab optional return, if method \code{"qarma"} is chosen for 
#'   estimation of the variance factor. Omitted, if \code{"iid"} is used. \cr
#'  \code{dcs_options} \tab an object of class \code{cds_options} containing the
#'   initial options of the dcs procedure. \cr
#' }
#' 
#' @section Details:
#' See the vignette for a more detailed description of the function.
#' 
#' @seealso \code{\link{set.options}}
#' 
#' @examples 
#' # See vignette("DCSmooth") for examples and explanation
#' 
#' y <- y.norm1 + matrix(rnorm(101^2), nrow = 101, ncol = 101)
#' dcs(y)
#' 
#' @export
#' 

dcs <- function(Y, dcs_options = set.options(), h = "auto", ...)
{
  #-------Front Matter-------#
  
  # check for correct inputs of data and options
  exception.check.Y(Y)
  exception.check.bndw(h, dcs_options)
  exception.check.options(dcs_options)
  n_x = dim(Y)[1]; n_t = dim(Y)[2]
  
  # process ellipsis arguments from (...)
  {
    # set up vectors for X and T if neccessary
    args_list <- list(...)
    add_options <- list() # additional options depending on type and var_est
    
    if (!exists("X", where = args_list))
    {
      X <- 0:(n_x - 1)/(n_x - 1)
    } else {
      X <- args_list$X
    }
    if (!exists("T", where = args_list))
    {
      T <- 0:(n_t - 1)/(n_t - 1)
    } else {
      T <- args_list$T
    }
    exception.check.XT(Y, X, T)
    
    # qarma order if var_est == "qarma"
    if (exists("model_order", where = args_list) &&
        dcs_options$var_est %in% c("qarma", "lm", "sarma"))
    {
      exception.check.model_order(args_list$model_order, dcs_options)
      add_options$model_order <- args_list$model_order
    } else if (!exists("model_order", where = args_list) &&
               dcs_options$var_est %in% c("qarma", "lm", "sarma")) {
      add_options$model_order <- list(ar = c(1, 1), ma = c(1, 1))
    }
    # maximum order of order selection is applied
    if (exists("order_max", where = args_list) &&
        ((dcs_options$var_est == "qarma_gpac") ||
        (dcs_options$var_est == "qarma_bic")))
    {
      exception.check.model_order(args_list$order_max, dcs_options)
      add_options$order_max <- args_list$order_max
    } else if (!exists("order_max", where = args_list) &&
               ((dcs_options$var_est == "qarma_gpac" ||
               dcs_options$var_est == "qarma_bic"))) {
      add_options$order_max <- list(ar = c(1, 1), ma = c(1, 1))
    }
  }

  #-------Bandwidth Selection-------#
  
  # check if given bandwidths exist
  if (h[1] == "auto")
  {
    h_select_auto <- TRUE
    
    time_start <- Sys.time() # get starting time for performance measuring
  
    # bandwidth selection process
    if (dcs_options$type == "KR")
    {
      ### kernel regression
      # TODO: allow for derivatives in kernel regression
      h_select_obj <- KR.bndw(Y, dcs_options, add_options)
      h_opt <- pmin(h_select_obj$h_opt, c(0.45, 0.45)) 
                                          # KR cannot handle larger bandwidths
    } else if (dcs_options$type == "LP") {
      ### local polynomial regression
      h_select_obj <- LP.bndw(Y, dcs_options, add_options)
      h_opt <- h_select_obj$h_opt
    }
    
    if (h_select_obj$var_model$stnry == FALSE)
    {
      warning("QARMA model not stationary, a change of the option ",
              "\"qarma_order\" is advised.")
    }
    
    time_end <- Sys.time()
    
  } else {
    h_opt <- h
    h_select_auto <- FALSE
    time_end <- time_start <- Sys.time()
  }
  
  #-------Estimation of Result-------#
  
  if (dcs_options$type == "KR") # kernel regression
  {
    ### prepare kernel functions
      kernel_x <- kernel_fcn_assign(dcs_options$kerns[1])
      kernel_t <- kernel_fcn_assign(dcs_options$kerns[2])
      
    ### check bandwidth
      if (any(h_opt > c(0.45, 0.45)))
      {
        h_opt = pmin(h_opt, c(0.45, 0.45))
        warning("Bandwidth h too large for \"KR\", changed to largest ",
                "working value.")
      }
    
    dcs_out <- KR_dcs_const0(yMat = Y, hVec = h_opt, drvVec = c(0, 0), 
                            kernFcnPtrX = kernel_x,
                            kernFcnPtrT = kernel_t)
  } else if (dcs_options$type == "LP") {     # local polynomial regression
    ### prepare weight functions
      # set variables for weight type
      kern_type_vec = sub("^([[:alpha:]]*).*", "\\1", dcs_options$kern)
      # extract weighting type
      mu_vec = as.numeric(substr(dcs_options$kern, 
                                 nchar(dcs_options$kern) - 1,
                                 nchar(dcs_options$kern) - 1))
      # extract kernel parameter mu
      weight_x = weight_fcn_assign(kern_type_vec[1])
      weight_t = weight_fcn_assign(kern_type_vec[2])
    
    ### check bandwidth
      if (any(h_opt < (dcs_options$p_order + 2)/(c(n_x, n_t) - 1)))
      {
        h_opt = pmax(h_opt, (dcs_options$p_order + 2)/(c(n_x, n_t) - 1))
        warning("Bandwidth h too small for \"LP\", changed to smallest ",
                "working value.")
      }
      
    dcs_out <- LP_dcs_const0_BMod(yMat = Y, hVec = h_opt, polyOrderVec = 
                            dcs_options$p_order, drvVec = dcs_options$drv,
                            muVec = mu_vec, weightFcnPtr_x = weight_x,
                            weightFcnPtr_t = weight_t)
  }
  
  # calculate residuals
  R <- Y - dcs_out
  
  if (h_select_auto == TRUE)
  {
    dcs_out <- list(Y = Y, X = X, T = T, M = dcs_out, R = R, h = h_opt,
                   c_f = h_select_obj$var_coef, 
                   var_model = h_select_obj$var_model,
                   dcs_options = dcs_options,
                   iterations = h_select_obj$iterations,
                   time_used = difftime(time_end, time_start, units = "secs"))
    attr(dcs_out, "h_select_auto") <- h_select_auto
  } else if (h_select_auto == FALSE) {
    # probably unadvised, that the dcs object differs according to the type
    # of h selection
    dcs_out <- list(X = X, T = T, Y = Y, M = dcs_out, R = R, h = h_opt,
                   c_f = NA, var_model = NA, dcs_options = dcs_options,
                   iterations = NA, time_used = NA)
    attr(dcs_out, "h_select_auto") <- h_select_auto
  }
  
  # apply class to output object
  class(dcs_out) <- "dcs"

  return(dcs_out)
}

#-------------------Function for plotting smoothed surface--------------------#

#' 3D Surface Plot of "dcs"-object or numeric matrix
#' 
#' @section Details:
#' \code{surface.dcs} uses the plotly device to plot the 3D surface of the given
#' \code{"dcs"}-object or matrix. If a "dcs"-object is passed to the function, 
#'  it can be chosen between plots of the original data (1), smoothed surface 
#'  (2) and residuals (3).
#' @param Y an object of class \code{"dcs"} or a numeric matrix that contains the
#'   values to be plotted.
#' @param plot_choice override the prompt to specify a plot, can be 
#'  \code{c(1, 2, 3)}.
#' @param ... optional arguments passed to the plot function.
#' 
#' @return \code{dcs.3d} returns an object of class "plotly" and "htmlwidget".
#' 
#' @seealso \code{\link{plot.dcs}}
#' 
#' @examples
#' # See vignette("DCSmooth") for examples and explanation
#' 
#' smth =  dcs(y.norm1 + rnorm(101^2))
#' surface.dcs(smth, plot_choice = 2)
#' 
#' @export
#' 

surface.dcs <- function(Y, plot_choice = "choice", ...)
{
  if (class(Y)[1] == "dcs")
  {
    fcn_arg <- list(...)
    
    if (plot_choice == "choice")
    {
      cat("Plot choices for dcs object:", fill = TRUE)
      choices <- c(1, 2, 3)
      choice_names <- c("original observations", "smoothed surface",
                        "residuals")
      choices_df <- data.frame(choices)
      colnames(choices_df) <- ""
      rownames(choices_df) <- choice_names
      print.data.frame(choices_df)
      plot_choice <- readline("Please enter the corresponding number: ")
      plot_choice <- as.numeric(plot_choice)
    } else if (!(plot_choice %in% 1:3)) {
      stop("Invalid value in argument \"plot_choice\". Use c(1, 2, 3).")
    }
    
    if (plot_choice == 1) {
      Y_mat <- Y$Y
    } else if (plot_choice == 2) {
      Y_mat <- Y$M
    } else if (plot_choice == 3) {
      Y_mat <- Y$R
    } else {
      stop(plot_choice, " is not a valid plot-type.")
    }
    
    .plotly.3d(Y = Y_mat, X = Y$X, T = Y$T, 
              color = c("#444C5C", "#78A5A3", "#E1B16A", "#CE5A57"), ...)
  } else {
    .plotly.3d(Y = Y, ...)
  }
}