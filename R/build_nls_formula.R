#' Helper function to build nls spatial formulas
#'

# Note - in future versions I can explore alternative functional forms here
# for temporal, nonlinear, or error process models.

build_nls_formula <- function(lags, covariates, y_name,
                           model = NULL){

  # Can include a "temporal approximation" of the lag process here in the
  # future.


  # Buld spatial approximation
  f_sp = paste(sprintf("rho^%s*SpLag%s", 1:lags, 1:lags), collapse = " + ")

  # Build functional form
  f_ff = paste(sprintf("b_%s * %s",covariates, covariates), collapse = " + ")

  # Build full formula
  return(as.formula(sprintf("%s ~ %s + %s", y_name, f_ff, f_sp)))
}

