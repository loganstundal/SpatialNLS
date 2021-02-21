#' Estimate spatial lag models using non-linear least squares
#'
#'

#-----------------------------------------------------------------------------#
#
# Author:        Logan Stundal
# Date:          February 02, 2021
# Purpose:       Function to estimate a spatial lag model via NLS
#
#
# Copyright (c): Logan Stundal, 2021
# Email:         stund005@umn.edu
#
#-----------------------------------------------------------------------------#
#
# "To do list:
#   1. sp_sim - add panel data functionality. Include setting temporal lag
#               and temporal length as well as cross-section size options."
#-----------------------------------------------------------------------------#


sp_lag <- function(lags,W, y, X){
  l<-lapply(1:lags, function(x){
    as.numeric(W^x %*% y)
  })
  l <- do.call(cbind, lapply(l, as.data.frame))
  colnames(l) <- paste("SpLag",1:lags, sep ="")

  return(cbind(X,l))
}

# w_subset = function(W, X){
#   # to do
# }

lagsarnls <- function(formula    = NULL,
                      W          = NULL,
                      data       = NULL,
                      subset_W   = FALSE,
                      lags       = 6,
                      intercept  = TRUE,
                      start_vals = NULL,
                      rho_start  = 0.2,
                      rho_bounds = NULL){

  # ----------------------------------- #
  # Initial setup based on code above
  mod_df <- model.frame(formula = formula, data = data)
  X <- model.matrix(formula, mod_df)
  y <- model.response(mod_df)
  # ----------------------------------- #


  # ----------------------------------- #
  # Check if W and model data incompatible
  if((nrow(W) != nrow(X)) & subset_W == FALSE){stop("Neighbors do not match data.")}

  # TO DO
  # if((nrow(W) != nrow(X)) & subset_W == TRUE){
  #   W = w_subset(W,X)
  # }
  # ----------------------------------- #


  # ----------------------------------- #
  # Set names
  y_name = names(mod_df)[1]
  colnames(X)[1] = "Intercept"
  colnames(X) = gsub(pattern = ":", replacement = "_", x = colnames(X))


  # Extract covariate names for nls formula construction
  if(intercept == FALSE){
    X = X[,-1]
    covariates = colnames(X)[1:ncol(X)]
  } else{
    covariates = colnames(X)[1:ncol(X)]
  }
  # ----------------------------------- #


  # ----------------------------------- #
  # Set starting values for nls optimizer
  if(!is.null(start_vals)){
    if(length(start_vals) != ncol(X)){
      stop("Starting values do not equal the number of covariates")
    } else{
      b_start = start_vals
    }
  } else{
    # Default - uses OLS estimates as starting values
    bstart = c(as.vector(solve(t(X) %*% X) %*% t(X) %*% y), rho_start)
  }
  names(bstart) = c(sprintf("b_%s",covariates), "rho")
  # ----------------------------------- #


  # ----------------------------------- #
  # Generate spatial lag approximations
  X <- sp_lag(lags = lags, W = W, y = y, X = X)
  # ----------------------------------- #


  # ----------------------------------- #
  # Build NLS model formula ---------------- This could be made a standalone function in the package for generic flexibility (i.e., non-spatial users)
  el_sp = paste(sprintf("rho^%s*SpLag%s", 1:lags, 1:lags), collapse = " + ")
  el_fn = paste(sprintf("b_%s * %s",covariates, covariates), collapse = " + ")
  nls_form = as.formula(sprintf("%s ~ %s + %s", y_name, el_fn, el_sp))
  # ----------------------------------- #


  # ----------------------------------- #
  # Bind DV to X and save as dataframe for nls data
  dat = cbind(X, mod_df[,1])
  names(dat)[ncol(dat)] = y_name
  # ----------------------------------- #


  # ----------------------------------- #
  # Estimate model
  if(is.null(rho_bounds)){
    m = nls(formula = nls_form,
            start   = bstart,
            data    = dat)
  } else{
    m = nls(formula   = nls_form,
            start     = bstart,
            algorithm = "port",
            lower     = c(rep(-Inf,length(covariates)), rho_bounds[1]),
            upper     = c(rep(Inf,length(covariates)),  rho_bounds[2]),
            data      = dat)
  }
  # ----------------------------------- #


  return(m)
}



