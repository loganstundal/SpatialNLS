#' Simulate spatial or spatiotemporal data
#'
#' @description spacetime_sim - easy simulation of spatial or spatiotemporal data.
#'
#' @param sample_size   sample size
#' @param k             number of vars
#' @param b_true        true beta values
#' @param rho_true      true rho (spatial) correlation
#' @param phi_true      true phi (temporal) correlation (not working)
#' @param prob_connect  density of connections in W
#' @param e_sd          noise in the dgp
#' @param probit        boolean simulate binary outcome data
#' @param temporal      boolean simulate temporal data
#' @param groups        number of groups if temporal is true
#' @param time_periods  number of time periods if temporal is true
#'
#' @returns A list consisting of 4 elements: data as a tibble, simulated
#'   parameters (user set and observed values), w, w as a list object
#'
#' @importFrom dplyr as_tibble
#' @importFrom spdep mat2listw listw2mat
#'
#' @export


spacetime_sim <- function(sample_size   = 100,
                          k             = 2,
                          b_true        = c(1, 2.5, -1),
                          rho_true      = 0.6,
                          phi_true      = 0.3,
                          prob_connect  = 0.3,
                          e_sd          = 0.2,
                          probit        = FALSE,
                          temporal      = FALSE,
                          groups        = NULL,
                          time_periods  = NULL){

  # NOTE : SAR Reduced Form:
  # y = (I - rho*W)^(-1) * (XB + e)
  if(!is.null(groups) & !is.null(time_periods)){
    # Add 1 time period to user-defined in order to create a 1-unit temporal lag
    time_periods = time_periods + 1

    sample_size = groups * time_periods
    warning("Sample size provided as well as group and time periods. Sample size updated to reflect groups*time_periods.")
  }

  # Exogenous variables:
  dat = X = replicate(k, rnorm(sample_size, runif(1,-1,1), sd = runif(1,0.5,2)))
  xnames = colnames(dat) = colnames(X) = sprintf("x%s",1:k)

  # Spatial Weights (row-standardized): W
  Id = diag(1, sample_size, sample_size)
  W  = matrix(data = rbinom(n    = sample_size*sample_size,
                            size = 1,
                            prob = prob_connect),
              nrow = sample_size,
              ncol = sample_size)

  # Clean-up W and row-standardize
  diag(W) = 0
  W[lower.tri(W)] <- t(W)[lower.tri(W)]

  # Record actual simulated connections:
  prob_connect_actual = as.numeric(prop.table(table(W))[2])

  ## ??? I'm not sure why I have this part, but I think it works better than the row stand. for accuracy...?
  nbl = spdep::mat2listw(x = W, style = "W")
  W   = spdep::listw2mat(nbl)


  # Solve - Spatial Multiplier: M
  if(temporal){

    # Do all the time stuff
    L <- diag(sample_size - time_periods)
    L <- rbind(matrix(0,time_periods,
                      sample_size - time_periods),L)
    L <- cbind(L,matrix(0,sample_size,
                        time_periods))
    M <- solve(Id - rho_true*W - phi_true*L)
  } else{
    M = solve(Id - rho_true * W)
  }


  # Row-standardized (n by n) W
  # W = apply(W, 1, function(x){x/sum(x)})


  # Homoskedastic error term - FUTURE CAN BUILD OUT FOR HETEROSKEDASTIC ERRORS HERE.
  e = rnorm(n    = sample_size,
            mean = 0,
            sd   = e_sd)

  # Matrix of Exogenous regressors and constant:
  X = cbind(1, X)

  # OUTCOME VARIABLE:
  if(probit){
    z = (M %*% (X %*% b_true)) + e
    z = pnorm(z)

    y = rbinom(n = sample_size, size = 1, prob = z)
  } else{
    y = (M %*% (X %*% b_true)) + e
  }


  # Return relevant simulation data and parameters
  dat = as_tibble(cbind(dat, y))
  colnames(dat) = c(xnames, "y")

  if(temporal){

    "THIS BLOC OF CODE [MAYBE L BLOC ABOVE] PRESENTS THE CURRENT PROBLEM - how to create the temporal lag with 'L' while
     preserving the original group / time structure defined by the fn call.

     If I can not figure out how badly I messed this function up, rectify by resolving any
     if(temporal){} chunks."

    # Generate temporal lag on y and drop NAs
    dat$ylag = as.vector(L %*% y)
    bob <- dat[dat$ylag>0,]


    # Add temporal and group variables if temporal true - NOTE ORDER IS IMPORTANT HERE (TIMES FIRST!)
    time_periods = time_periods - 1
    sample_size  = time_periods * groups

    dat$times <- rep(1:(sample_size/groups), times = groups)
    dat$group <- rep(1:(sample_size/time_periods), each  = time_periods)
  }

  d = list("data" = dat,
           "sim_params" = list("betas"               = b_true,
                               "rho"                 = rho_true,
                               "nbr_connections_set" = prob_connect,
                               "nbr_connections_sim" = prob_connect_actual,
                               "sample_size"         = sample_size),
           "W_matrix" = W,
           "W_list"   = nbl)

  return(d)
}





#-----------------------------------------------------------------------------#
# FOR TESTING PURPOSES - DELETE LATER
        # rm(list=ls())
        # sample_size   = 100
        # k             = 2
        # b_true        = c(1, 2.5, -1)
        # rho_true      = 0.6
        # phi_true      = 0.3
        # prob_connect  = 0.3
        # e_sd          = 0.2
        # temporal      = T
        # groups        = 10
        # time_periods  = 12
        # probit        = F
#-----------------------------------------------------------------------------#

