#' Simulate spatiotemporal data
#'
#' @description spacetime_sim2 simulates spatiotemporal data with the option of
#'   unit and time effects in addtion to space and time dynamic processes.
#'
#' @param rho_true Measure of spatial dependence (-1 to 1)
#' @param phi_true Measure of temporal dependence (-1 to 1)
#' @param b_true True paramaeter values for intercept (if included) and Xs (must
#'   equal length of k+1)
#' @param k Number of exogenous independent variables to include
#' @param e_sd Disturbance term standard deviation (assumed homoskedastic)
#' @param intercept Boolean to include intercept
#' @param panel_units Number of units in panel
#' @param panel_times Number of time periods in panel
#' @param forfeit_lag Boolean to forfeit temporal lag. FALSE by default. Mostly
#'   to experiment with consequences of dropping first lagged unit on model
#'   performance.
#' @param unit_effects Boolean to include unit-fixed-effects
#' @param time_effects Boolean to include time-fixed-effects
#' @param unit_betas,time_betas
#'   Vector for unit or time FE betas if included in DGP. Must equal length
#'   (unit_effects - 1) or (time_effects - 1) to avoid perfect collinearity.
#'   If NULL (default), random values are generated and returned.
#' @param prob_connect Percentage of connections in matrix W, 0.3 by default.
#'   Higher values can help to avoid "island" cases in W where a unit has no
#'   neighbors.
#'
#' @importFrom magic adiag
#' @importFrom spdep mat2listw listw2mat
#' @importFrom magrittr "%>%"
#' @importFrom dplyr mutate group_by ungroup lag bind_cols
#' @importFrom tidyr drop_na
#'
#' @return Returns a list containing a dataframe with simulated spatiotemporal
#'   data, input parameters, and spatial weights information in matrix and
#'   spdep list formats.
#' @export
#'
#' @examples
#' # Rapid simulation using defualt values
#' x = spacetime_sim2()
#' head(x$data)
#'


spacetime_sim2 <- function(rho_true = 0.3,
                           phi_true = 0.8,
                           b_true   = c(1, 2.5, -1),
                           k        = 2,
                           e_sd     = 1,
                           intercept= TRUE,
                           panel_units = 20,
                           panel_times = 10,
                           forfeit_lag = FALSE,
                           unit_effects = T,
                           time_effects = T,
                           unit_betas   = NULL,
                           time_betas   = NULL,
                           prob_connect = 0.3){

  #-----------------------------------------------------------------------------#
  # NOTES                                                                   ----
  #-----------------------------------------------------------------------------#
  # To dos:
  #
  # 1. add "newdata" parameter. Thoughts: can run spactime_sim2 once saving
  #    spatially depdent y and associated W.
  #    - can then use this spatially dependent y (AND W) in a second run of
  #      spacetime_sim2 using the same W and y_orig to generate a new new y
  #      which is both spatially dynamic and depends on an exogenous spatially
  #      dependent "x" (y_orig)
  # 2. add a set.seed option
  #-----------------------------------------------------------------------------#


  #-----------------------------------------------------------------------------#
  # SETUP                                                                   ----
  #-----------------------------------------------------------------------------#
  # ----------------------------------- #
  # USER INPUT MODIFICATION
  # ----------------------------------- #
  if(!forfeit_lag){
    panel_times = panel_times + 1
  }

  nobs <- panel_units * panel_times


  # ----------------------------------- #
  # EXOGENOUS VARIABLES
  # ----------------------------------- #
  X <- replicate(k, rnorm(nobs, runif(1,-1,1), sd = runif(1,0.5,2)))
  colnames(X) <- sprintf("x%s",1:k)
  #-----------------------------------------------------------------------------#



  #-----------------------------------------------------------------------------#
  # CREATE W                                                                ----
  #-----------------------------------------------------------------------------#
  W  = matrix(data = rbinom(n    = panel_units,
                            size = 1,
                            prob = prob_connect),
              nrow = panel_units,
              ncol = panel_units)

  # Clean-up W and row-standardize
  diag(W) = 0
  W[lower.tri(W)] <- t(W)[lower.tri(W)]

  # Record actual simulated connections:
  prob_connect_actual = as.numeric(prop.table(table(W))[2])

  # W <- adiag(W, time_periods)
  W <- do.call(adiag, replicate(panel_times, W, simplify = FALSE))

  # Create spdep weigh objects
  nbl = suppressWarnings({spdep::mat2listw(x = W, style = "W")})
  W   = spdep::listw2mat(nbl)
  #-----------------------------------------------------------------------------#



  #-----------------------------------------------------------------------------#
  # GENERATE TEMPORAL LAG MATRIX (L) AND SPATIO-TEMPORAL MULTIPLIER (M)     ----
  #-----------------------------------------------------------------------------#
  Id <- diag(1, nobs, nobs)

  # Build L, working off conversations with Rob here re: Jude's idea for L
  L <- diag(nobs - panel_units)
  L <- rbind(matrix(0, panel_units, nobs - panel_units), L)
  L <- cbind(L, matrix(0, nobs, panel_units))

  M <- solve(Id - rho_true*W - phi_true*L)
  #-----------------------------------------------------------------------------#



  #-----------------------------------------------------------------------------#
  # GENERATE ERROR TERM                                                     ----
  #-----------------------------------------------------------------------------#
  e = rnorm(n    = nobs,
            mean = 0,
            sd   = e_sd)
  # Homoskedastic error term - FUTURE CAN BUILD OUT FOR HETEROSKEDASTIC ERRORS HERE.
  #-----------------------------------------------------------------------------#



  #-----------------------------------------------------------------------------#
  # GENERATE UNIT AND TIME IDs                                              ----
  #-----------------------------------------------------------------------------#
  X <- X %>%
    as.data.frame() %>%
    bind_cols(.,data.frame("Group" = factor(rep(LETTERS[1:panel_units], panel_times),
                                            levels = LETTERS[1:panel_units]),
                           "Time"  = factor(rep(1:panel_times, each = panel_units),
                                            levels = 1:panel_times)))

  # Create a data frame here to merge generated "y" onto for return
  # Necessary to do before overwriting X as a "model.matrix" below.
  d <- X
  #-----------------------------------------------------------------------------#



  #-----------------------------------------------------------------------------#
  # CREATE MODEL FRAME (UNIT & TIME EFFECTS IF SELECTED)                    ----
  #-----------------------------------------------------------------------------#
  if(unit_effects & is.null(unit_betas)){
    unit_betas = runif(n = (panel_units-1), min = -1, max = 1)
  }
  if(time_effects & is.null(time_betas)){
    time_betas = runif(n = (panel_times-1), min = -1, max = 1)
  }


  if(unit_effects & !time_effects){
    print("ONLY UNITS")
    form   <- ~ x1 + x2 + Group
    mod_df <- model.frame(formula = form, data = X)
    X      <- model.matrix(object = form, data = mod_df)
    all_bs <- c(b_true, unit_betas)

  } else if(time_effects & !unit_effects){
    print("ONLY TIME")
    form   <- ~ x1 + x2 + Time
    mod_df <- model.frame(formula = form, data = X)
    X      <- model.matrix(object = form, data = mod_df)
    all_bs <- c(b_true, time_betas)

  } else if(unit_effects & time_effects){
    print("UNITS AND TIME")
    form   <- ~ x1 + x2 + Group + Time
    mod_df <- model.frame(formula = form, data = X)
    X      <- model.matrix(object = form, data = mod_df)
    all_bs <- c(b_true, unit_betas, time_betas)

  } else{
    print("NEITHER UNITS NOR TIME")
    form   <- ~ x1 + x2
    mod_df <- model.frame(formula = form, data = X)
    X      <- model.matrix(object = form, data = mod_df)
    all_bs <- c(b_true)
  }



  #-----------------------------------------------------------------------------#
  # DEPENDENT VARIABLE CREATION                                             ----
  #-----------------------------------------------------------------------------#
  y <- as.vector((M %*% (X %*% all_bs)) + e)
  #-----------------------------------------------------------------------------#



  #-----------------------------------------------------------------------------#
  # TIDY DATA FOR RETURN                                                    ----
  #-----------------------------------------------------------------------------#
  d <- d %>%
    mutate(y = y) %>%
    group_by(Group) %>%
    mutate(ylag = lag(y, 1)) %>%
    ungroup() %>%
    tidyr::drop_na(ylag) %>%
    mutate(Time = factor(rep(1:(panel_times-1), each = panel_units),
                         levels = 1:(panel_times-1)))

  # Rectify W to drop lag if necessary
  if(!forfeit_lag){
    W <- W[(nobs + 1 - (panel_units * (panel_times - 1))):nobs,
           (nobs + 1 - (panel_units * (panel_times - 1))):nobs]
  }

  # Recreate spdep weigh objects for revised W
  nbl = suppressWarnings({spdep::mat2listw(x = W, style = "W")})
  W   = spdep::listw2mat(nbl)

  #-----------------------------------------------------------------------------#



  #-----------------------------------------------------------------------------#
  # REUTRN STATEMENT
  #-----------------------------------------------------------------------------#
  # Return panel_times and nobs to original user-input values when lag not dropped
  if(!forfeit_lag){
    panel_times = panel_times - 1
    nobs        = panel_units * panel_times
  }

  inputs = list(
    # Parameters
    "rho_true" = rho_true,
    "phi_true" = phi_true,
    "b_true"   = b_true,
    "k"        = k,
    "e_sd"     = e_sd,
    # "intercept"= TRUE,

    # Panel characteristics
    "panel_units" = panel_units,
    "panel_times" = panel_times,
    "panel_nobs"  = nobs,
    "forfeit_lag" = forfeit_lag,

    "unit_effects" = unit_effects,
    "time_effects" = time_effects,
    "unit_betas"   = unit_betas,
    "time_betas"   = time_betas,

    # W Connection density
    "prob_connect"     = prob_connect,
    "prob_connect_sim" = prob_connect_actual)


  rt = list("data"     = d,
            "params"   = inputs,
            "W_matrix" = W,
            "W_list"   = nbl)

  return(rt)
  #-----------------------------------------------------------------------------#
}
