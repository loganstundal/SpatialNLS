#' Estimate spatial (marginal) effects
#'
#' @description impacts_e uses parametric simulation to estimate spatial effects
#'  for parameters in a model containing an endogenous spatial lag.
#'
#' @param model      an nls model object estimated with lagsarnls()
#' @param eigen_vals pre-computed eigen values from a spatial weights matrix
#' @param w          not used yet
#' @param sims       number of simulation iterations
#' @param cred_int   credibility interval to calculate
#' @param vcv        option to supply a variance-covariance matrix. Useful to
#'  estimate spatial effects with robust or clustered standard errors. Default
#'  uses model variance-covariance matrix extracted using generic vcv(model)
#' @param params     vector of parameters to estimate spatial effects
#' @param ignore_params vector of parameters to ignore. Note - if "params" option
#'  is specified, this option is ignored and only results for those variables
#'  in "params" are returned. If "params" is empty, the function will estimate
#'  spatial effects for all variables in a model except those listed in this
#'  option. If both params and ignore_params are null, the function estimates
#'  spatial effects for all variables in the model.
#'
#' @return Returns a tidied dataframe containing spatial effects for all or
#'  selected variable subsets from a spatial lag model.
#' @export
#'
#' @importFrom dplyr mutate relocate bind_rows select
#' @importFrom magrittr "%>%"
#' @importFrom stringr str_detect
#' @importFrom MASS mvrnorm
#' @importFrom coda as.mcmc HPDinterval
#' @importFrom tibble rownames_to_column
#'
#' @examples
#' x <- spacetime_sim2(phi_true = 0)
#' d <- x$data
#' w <- x$W_matrix
#' e <- eigen(w, only.values = TRUE)$values
#' m <- lagsarnls(y ~ x1 + x2, data = d, W = w)
#' impacts_e(model      = m,
#'           eigen_vals = e)
#'


impacts_e <- function(model,
                      eigen_vals,
                      w        = NULL,
                      sims     = 1000,
                      cred_int = 0.95,
                      vcv      = NULL,
                      params   = NULL,
                      ignore_params = NULL){

  if(!is.null(params)){
    params <- paste(params, collapse = "|")
    params <- coef(model)[str_detect(names(coef(model)), params)]
  } else{
    if(!is.null(ignore_params)){
      ignore_params <- paste("Int|rho",
                             paste(ignore_params, sep = "", collapse = "|"), sep = "|")
    } else{
      ignore_params <- "Int|rho"
    }
    params <- coef(model)[!str_detect(names(coef(model)), ignore_params)]
  }

  rho    <- coef(model)["rho"]

  nobs      <- length(eigen_vals)
  tmp_names <- apply(expand.grid(names(params),
                                 c("Direct","Indirect","Total")), 1, paste, collapse=".")

  param_names <- c(names(params),"rho")

  if(is.null(vcv)){
    vcv <- vcov(model)[param_names, param_names]
  } else if(dim(vcv)[1] != length(param_names)){
    vcv <- vcv[param_names, param_names]
  }

  # DF to store simulation results
  tmp       <- data.frame()

  # Simulation
  for(i in 1:sims){

    # Sample coefs
    param_sim <- mvrnorm(n         = 1,
                         mu        = c(params,rho),
                         Sigma     = vcv)
    sim_rho   <- param_sim["rho"]
    sim_betas <- param_sim[!str_detect(names(param_sim), "rho")]

    # Calculate effects
    direct   <- as.numeric((sim_betas * sum(1 / (1 - sim_rho * eigen_vals))) / nobs)
    total    <- sim_betas/(1-sim_rho)
    indirect <- total - direct

    # Gather results
    effs        <- c(direct, indirect, total)
    names(effs) <- tmp_names

    tmp <- bind_rows(tmp, as.data.frame(t(effs)))
  }

  # Note - use the credibility interval and median due to potential
  # for explosive draws on rho (rho>abs(1)). Median and quantile
  # mitigate impact of such draws on inference.

  tmp <- as.mcmc(tmp)
  res <- as.data.frame(HPDinterval(tmp, prob = cred_int))

  res <- res %>% mutate(median = apply(tmp, 2, median)) %>%
    relocate("median", .after = "lower")


  # Tidy results data frame
  res <- res %>%
    rownames_to_column(var = "var") %>%
    mutate(Effect = str_split(var, pattern = "[.]") %>% purrr::map(2) %>% unlist,
           Variable = str_split(var, pattern = "[.]") %>% purrr::map(1) %>% unlist) %>%
    mutate(Variable = str_remove(Variable, "b_")) %>%
    dplyr::select(Variable, Effect, lower, median, upper) %>%
    arrange(Variable, Effect)
  return(res)
}
# ----------------------------------- #












