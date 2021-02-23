#' Estimate spatial and spatiotemporal effects
#'
#' @description \code{impacts} estimates dynamic marginal effects for spatial and
#'   spatiotemporal dynamic regression models estimated with \code{\link{lagsarnls}}.
#'   For spatial models (i.e., no temporal lag), the function returns the spatial
#'   impacts (direct, indirect, and total) for unit-changes in exogenous variables.
#'   For spatiotemporal models (i.e., those with temporal lags), the function returns
#'   the long-run steady state (LRSS) response to unit-changes in exogenous variables.
#'
#' @param model Model object estimated with \code{lagsarnls}
#' @param variable Vector of one or more characters matching variable names in model
#'   formula. If NULL (default), the function returns impact estimates for all
#'   variables.
#' @param W Spatial weights matrix.
#' @param lag_dv Character string naming lagged dependent variable in model.
#' @param groups Vector of IDs for cross-sectional units in data. NOTE: the
#'   number of unique groups must match the dimensions of W
#'   (i.e, nrow(W) == length(unique(groups))), or the function will return
#'   an error. In this case, subset your W matrix for ONE temporal
#'   cross-section.
#' @return Returnd a data frame with estimated impacts
#'
#' @importFrom stringr str_detect str_to_lower str_remove_all
#' @importFrom magrittr "%>%"
#' @importFrom dplyr select everything
#'
#' @export
#'
#' @examples
#'  # Simulate spatiotemporal data using defaults from \code{spacetime_sim2}
#'  x = spacetime_sim2()
#'
#'  # Estimate model
#'  m = lagsarnls(formula = y ~ ylag + x1 + x2,
#'                data    = x$data,
#'                W       = x$W_matrix)
#'
#'  # Subset w for ONE cross-section
#'  n  = length(unique(x$data$Group))
#'  W2 = x$W_matrix[1:n,1:n]
#'
#'  impacts(model  = m,
#'          W      = W2,
#'          lag_dv = "ylag",
#'          groups = x$data$Groups)
#'

impacts <- function(model,
                    variable = NULL,
                    W,
                    lag_dv = NULL,
                    groups = NULL){

  if(!is.null(lag_dv) & is.null(groups)){
    stop("Please supply unit IDs in 'groups' parameter to estimate LRSS.")
  }

  if(is.null(lag_dv)){
    # cat("\14");
    cat(sprintf("Spatial (Marginal) Effects\n"))
  } else{
    # cat("\14");
    cat(sprintf("LRSS Spatiotemproal Response\n"))
  }

  rho <- coef(model)["rho"]
  if(is.null(variable)){
    variable <- names(coef(model))
    variable <- variable[!str_detect(str_to_lower(variable),
                                     pattern = c("int|rho|ylag"))]
  } else{
    if(!str_detect(variable, "b_")){
      variable = paste0("b_", variable)
    }
  }
  n  <- nrow(W)
  Id <- diag(1, n, n)

  if(n > length(unique(groups))){
    stop(paste0("W > groups. LRSS relies on ONE cross-section for impact estimation.\n",
                "Please supply a W-matrix where nrow(W) == one cross-section in your data."))
  }

  if(!is.null(lag_dv)){
    phi <- coef(model)[paste0("b_",lag_dv)]
    M   <- solve(Id - rho * W - phi * Id)
    variable <- variable[!str_detect(str_to_lower(variable),
                                     pattern = lag_dv)]
  } else{
    M   <- solve(Id - rho * W)
  }

  im <- sapply(variable, function(x){
    b  <- coef(model)[x]
    bM <- b * M

    dir <- mean(diag(bM))
    ind <- (sum(bM) - sum(diag(bM))) / n
    tot <- sum(bM) / n

    return(cbind("Direct"   = dir,
                 "Indirect" = ind,
                 "Total"    = tot))
  }, simplify = FALSE)

  im <- do.call(rbind, im) %>%
    as.data.frame()
  im$Variable <- str_remove_all(variable, "b_")
  im <- im %>%
    select(Variable, everything())
  return(im)
}
