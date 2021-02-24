#' Estimate credibility intervals for spatial and spatiotemporal impacts using
#'   parametric simulation.
#'
#' @description This function uses a parametric simulation to estimate
#'  credibility intervals for impacts from a spatial or spatiotemporal
#'  regression model estimate with \code{\link{lagsarnls}}.
#'
#' @param model Model object of class nls estimated with
#'  \code{\link{lagsarnls}}
#' @param W Spatial weights matrix SUBSET for one cross-section if
#'  spatiotemporal
#' @param beta Character string indicating variable to estimate impacts
#' @param phi Character string indicating name of temporally lagged dependent
#'  variable
#' @param rho Character string indicating name of spatial lag (set to "rho" by
#'  default)
#' @param unit_ids Vector of unit IDs from panel data
#' @param time_ids Vector of time IDs from panel data
#' @param parametric Not used. Relevant for delta method implementation in
#'  the future
#' @param truncated_normal Boolean. Use truncated normal for sampling? Default
#'  True. See details for more.
#' @param cred_interval Credibility interval to estimate as a decimal
#'  (0.95 ~ 95% by default)
#' @param sims Simulation iterations
#' @param seed Set seed (relevant for replication purposes)
#'
#' @details Note: Sampling from a truncated normal distribution
#'  (\code{truncated_normal = TRUE}, i.e., the default) is highly recommended.
#'  Sampling from a standard multivariate normal may result in sampled values of
#'  Rho or Phi which exceed the theoretical bounds of [-1,1] resulting in an
#'  explosive system and invalid quantites of interest. To avoid this possibility,
#'  the truncated normal sampling option limits draws for Rho and Phi between
#'  -1 and 1.
#'
#' @return Returns a data frame of Direct, Indirect, and Total Impact estimates
#'  along with simulated credibility intervals. "Estimate" reflects the median
#'  response from simulations while the Lower and Uppwer bounds are the quantile
#'  at probabilities specified by the user (95% default). Quantiles help to
#'  further mitigate the problems of explosive parameter draws.
#' @export
#'
#' @examples
#' m <- lagsarnls(formula = y ~ ylag + x1 + x2,
#'                data    = x$data,
#'                W       = x$W_matrix)
#' n  <- length(unique(x$data$Group))
#' W2 <- x$W_matrix[1:n,1:n]
#'
#' effs <- impacts_ci(model            = m,
#'                    W                = W2,
#'                    beta             = "x1",
#'                    phi              = "ylag",
#'                    rho              = "rho",
#'                    unit_ids         = x$data$Group,
#'                    time_ids         = x$data$Time,
#'                    parametric       = TRUE,
#'                    truncated_normal = TRUE,
#'                    cred_interval    = 0.95,
#'                    sims             = 100,
#'                    seed             = 123)
#'
#' print(effs)


impacts_ci <- function(model,
                       W,
                       beta,
                       phi        = NULL,
                       rho        = "rho",
                       unit_ids   = NULL,
                       time_ids   = NULL,
                       parametric = TRUE,
                       truncated_normal = TRUE,
                       cred_interval    = 0.95, # 95% credibility interval
                       sims       = 100,
                       seed       = 1){

  # IMPORTS LIST [WORKING] - ADD TO ROXYGEN2 HEADER AND DELETE LATER
  # MASS mvrnorm
  # tmvtnorm rtmvnorm

  # unit_ids <- head(dat$statenm, n)

  # ----------------------------------- #
  # Setup
  # ----------------------------------- #
  if(class(model) != "nls"){
    stop("Model object not of nls class.")
  }

  #####
  # will need error handling here for unbalanced panels


  # sub unit-ids for one time slice
  n <- unique(table(time_ids))
  if(length(n) > 1){
    stop("Panel unbalanced, I'm sure solving this will require staring at my screen while frowning.")
  }
  unit_ids <- unit_ids[1:n]


  # Create Identity matrix
  Id <- diag(1, n, n)

  # Verify W-dims == n
  if(any(dim(W) != c(n,n))){
    stop("W matrix dims do not match observations, n")
  }
  #####


  # Ensure parameters have relevant "b_" prevfix
  if(!str_starts(beta,"b_")){
    beta <- paste0("b_", beta)
  }
  if(!str_starts(phi,"b_")){
    phi <- paste0("b_", phi)
  }

  # Extract variance-covariance from model
  vcvm <- vcov(model)[c(rho, phi, beta),
                      c(rho, phi, beta)]

  if(!all(dim(vcvm) == c(3,3))){
    stop("Problems extracting variance-covariance matrix from model. God speed.")
  }
  if(any(is.na(vcvm))){
    stop("NAs in vcov... how did you mangage to get this far?")
  }


  # Extract relevant model parameters
  beta_param <- as.numeric(coef(model)[beta])
  phi_param  <- as.numeric(coef(model)[phi])
  rho_param  <- as.numeric(coef(model)[rho])

  # # Handle phi (lag dv - temporal)
  # if(is.null(phi)){
  #   stop("Dependent varible missing, please specify phi option to coefficient name")
  # } else{
  #   if(!str_starts(phi, "b_")){
  #     phi <- paste0("b_", phi)
  #     phi <- coef(model)[phi]
  #   }
  # }
  #
  # # Handle rho (lag dv - spatial)
  # if(is.null(rho)){
  #   rho <- coef(model)["rho"]
  #   if(is.na(rho)){
  #     stop("Unable to find spatial lag - rho - in model. Please set rho parameter to coefficient name")
  #   }
  # } else if(is.character(rho)){
  #   rho <- coef(model)[rho]
  #   if(is.na(rho)){
  #     stop("Unable to find spatial lag - rho - in model. Please set rho parameter to coefficient name")
  #   }
  # }
  # ----------------------------------- #

  # Compute user-desired significance bounds
  cred_interval <- 0.95
  a  <- (1-cred_interval) / 2
  lb <- a
  ub <- cred_interval + a
  cred_interval <- cred_interval + a

  # ----------------------------------- #
  # Parametric simulation
  # ----------------------------------- #
  # r_draws = c()
  # p_draws = c()
  # b_draws = c()

  # sims = 1e3
  pb = txtProgressBar(min = 0, max = sims, initial = 0)

  # Create array to store n*n matrices over n sims
  sim.effs <- array(data     = 0,
                    dim      = c(n,n,sims),
                    dimnames = list(unit_ids,
                                    unit_ids))

  set.seed(seed)

  cat("\14");for(i in 1:sims){
    setTxtProgressBar(pb,i)

    if(truncated_normal){
      draws <- rtmvnorm(n     = 1,
                        mean  = c(rho_param, phi_param, beta_param),
                        sigma = vcvm,
                        lower = c(-1, -1, -Inf),
                        upper = c(1, 1, Inf))
      draws <- as.vector(draws)
      names(draws) <- c(rho, phi, beta)
    } else{
      draws <- mvrnorm(n         = 1,
                       mu        = c(rho_param, phi_param, beta_param),
                       Sigma     = vcvm)
    }

    sim.rho  = draws[rho]
    sim.phi  = draws[phi]
    sim.beta = draws[beta]

    # r_draws = c(r_draws, sim.rho)
    # p_draws = c(p_draws, sim.phi)
    # b_draws = c(b_draws, sim.beta)

    sim.m  = solve(Id - sim.rho * W - sim.phi * Id)
    sim.bm = sim.beta * sim.m

    rownames(sim.bm) = colnames(sim.bm) = unit_ids

    sim.effs[,,i] = sim.bm
  }
  close(pb)
  # ----------------------------------- #


  # ----------------------------------- #
  # Calculate impact estimates and CIs
  # ----------------------------------- #
  # TO NOTE IN DOCUMENTATION - here going to use median absolute deviation (mad)
  # to construct confidence intervals on sampled values due to mad being more
  # robust to outliers than sd. Outliers are a problem, particularly when
  # sampling from a standard (i.e., non-truncated) mvn.

  # DIRECT EFFECTS
  all_directs      <- apply(sim.effs, 3, function(x){median(diag(x))})

  direct           <- median(all_directs)
  # all_directs_mad  <- mad(all_directs)
  # direct_lb        <- direct - abs(qt(cred_interval, Inf) * direct)
  # direct_ub        <- direct + abs(qt(cred_interval, Inf) * direct)
  direct_lb        <- quantile(all_directs, probs = lb)
  direct_ub        <- quantile(all_directs, probs = ub)


  # c(direct_lb, direct, direct_ub)

  # # AVERAGE DIRECT EFFECT (ADE) - across all simulations
  # ade.mean <- mean(apply(sim.effs, 3, function(x) median(diag(x))))
  # hist(apply(sim.effs, 3, function(x) median(diag(x))))
  #
  # # 95% credibility interval
  # ade.ci   <-   quantile(apply(sim.effs, 3, function(x) median(diag(x))),
  #                        probs = c(0.025, 0.975))

  # INDIRECT EFFECTS
  all_indirects    <- apply(sim.effs, 3, function(x){
    (sum(x) - sum(diag(x))) / n
  })
  indirect           <- median(all_indirects)
  # all_indirects_mad  <- mad(all_indirects)
  # indirect_lb        <- indirect - abs(qt(cred_interval, Inf) * indirect)
  # indirect_ub        <- indirect + abs(qt(cred_interval, Inf) * indirect)
  indirect_lb        <- quantile(all_indirects, probs = lb)
  indirect_ub        <- quantile(all_indirects, probs = ub)

  # c(indirect_lb, indirect, indirect_ub)


  # TOTAL EFFECTS
  all_totals <- apply(sim.effs, 3, function(x){
    sum(x) / n
  })
  total       <- median(all_totals)
  # totals_mad  <- mad(all_totals)
  # total_lb    <- total - abs(qt(cred_interval, Inf) * total)
  # total_ub    <- total + abs(qt(cred_interval, Inf) * total)
  total_lb        <- quantile(all_totals, probs = lb)
  total_ub        <- quantile(all_totals, probs = ub)

  c(total_lb, total, total_ub)


  # ----------------------------------- #
  # EFFECTS DATA FRAME
  # ----------------------------------- #

  # Collect Effects into data frame for return statement
  d = data.frame("Effect"     = c("Direct","Indirect","Total"),
                 "LowerBound" = c(direct_lb, indirect_lb, total_lb),
                 "Estimate"   = c(direct, indirect, total),
                 "UpperBound" = c(direct_ub, indirect_ub, total_ub))

  # ----------------------------------- #


  # ----------------------------------- #
  # Return statement
  # ----------------------------------- #
  return(d)
  if(!truncated_normal){
    cat(paste0("Sampling done using multivariate normal. Use cation interpreting results\n",
               "as there is a possibility of explosive parameter draws on phi or rho which\n",
               "may have been sampled outside of -1, 1 bounds."))
  }
  # ----------------------------------- #


  # ----------------------------------- #
  # Testing = delete later
  # ----------------------------------- #
  # rm(draws, sim.rho, sim.beta, sim.m, sim.bm)
  # ----------------------------------- #
}
