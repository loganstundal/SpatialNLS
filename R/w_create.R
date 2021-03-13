#' w create
#'
#' @description w_create automates spatial weights construction based on a
#'  supplied model formula. The user supplies a formula class object, their
#'  data, and a simple features object containing the geomety they with to
#'  use for identifying neighbors. The function will return a cross-section
#'  or, if a unique time_id is supplied, a panel spatial weights object.
#'
#'  NOTE - at present, the function will not work wtih unbalanced panels.
#'
#' @param formula A formula class object which will be used for subsetting data
#' @param data    The data frame containing rows which will be subset
#' @param sp_data Simple feature object containing geometry to identify neighbors
#' @param unit_id A character string identifying units in the data
#' @param time_id A character string identifying times in the data. If supplied,
#'   the function returns a panel spatial weights matrix. Note, for now only
#'   balanced panels are supported - i.e., units must be identical from year-
#'   to-year.
#' @param style   Definition of neighbors passed to poly2nb() from spedep. Default
#'   is "W" for row-standardized spatial weights.
#' @param nb_list Boolean True or False (default) indicating whether the function
#'   should return a spatial weights list object in addition to a spatial weights
#'   matrix
#' @param silent Boolean True or False (default) indicating whether the function
#'  should print progress messages to the console.
#'
#' @importFrom magrittr "%>%"
#' @importFrom sf st_drop_geometry as_Spatial
#' @importFrom tibble rowid_to_column
#' @importFrom dplyr select filter
#'
#' @return Returns a sparse spatial weights matrix and, if requested, a
#'  neighbor list which can be used with models in spatialreg
#' @export


w_create <- function(formula,                 # A formula class object to extact relevant variables from data
                     data,                    # Full data frame to subset in model
                     sp_data,                 # Single cross section of data to rapidly copy spatial weights
                     unit_id = NULL,          # String for unique unit-ids
                     time_id = NULL,          # String for unique time-ids
                     style   = "W",           # Spatial weights style (row-standardized by default)
                     nb_list = FALSE,
                     silent  = FALSE){
  # ----------------------------------- #
  # This function takes a formula and dataframe and returns a
  # spatial weights subset with no missing values ready for
  # modeling in lagsarnls()
  # ----------------------------------- #

  # ----------------------------------- #
  # Setup - Subset data based on formala
  # ----------------------------------- #
  if(!silent){cat("Extracting observations based on formula...\n")}

  if(any(class(data) == "sf")){
    data = as.data.frame(data %>% st_drop_geometry())
  }

  datf <- model.frame(formula = formula, data = data)

  # Variales to preserve: formula variables, unit and time IDS, and
  # a row "id" which will be created
  cs   <- c(colnames(datf), unit_id, time_id, "id")

  # Create row IDs from unique row-ids extracted from the
  # model frame - these are the observations from the
  # data which have all observations based on the formula.
  rs   <- as.numeric(rownames(datf))


  # subset data for spatial weights creation here retaining only
  # rows in the data corresponding to the supplied formual
  data <- data %>%
    rownames_to_column(var = "id") %>%
    select(cs)

  data       <- data[rs,]
  unique_obs <- unique(data[[unit_id]])
  # ----------------------------------- #


  # ----------------------------------- #
  # FOR FUTURE - unbalanced panels code possibility
  # ----------------------------------- #
  # units <- unique(groupN)
  # N <- length(units)
  # time <- unique(groupT)
  # T <- length(time)
  # brows <- c()
  # for (i in 1:T) {
  #   br <- which(groupT == time[i])
  #   check <- length(br) == N
  #   if (check) {
  #     brows <- c(brows, br)
  #   }
  # }

  # Source: Bailey and Katz in their PCSE package which requires
  # similar subsetting
  # https://cran.r-project.org/web/packages/pcse/vignettes/pcse.pdf
  # ----------------------------------- #


  # ----------------------------------- #
  # Subset supplied spatial data
  # ----------------------------------- #
  sp_data <- sp_data %>%
    filter(!!sym(unit_id) %in% unique_obs) %>%
    as_Spatial(.)
  # ----------------------------------- #


  # ----------------------------------- #
  # Create spatial neighbors and matrices
  # ----------------------------------- #
  if(!silent){cat("Creating single spatial weights cross-section...\n")}

  nbs <- spdep::poly2nb(pl        = sp_data,
                        row.names = sp_data[[unit_id]])

  nbw <- spdep::nb2mat(neighbours = nbs,
                       style      = style)
  colnames(nbw) <- rownames(nbw)

  # Create neighbor list object if requested
  if(nb_list & is.null(time_id)){
    if(!silent){cat("Working on panel list...\n")}
    nbl <- spdep::nb2listw(nbs, style = "W")
  }


  # ----------------------------------- #
  # Create panel spatial weights if time_id supplied
  # ----------------------------------- #
  if(!is.null(time_id)){
    if(!silent){cat("Creating panel spatial weights...\n")}
    nbw <- do.call(what = magic::adiag,
                   args = replicate(n        = length(unique(data[[time_id]])),
                                    expr     = nbw,
                                    simplify = FALSE))

    if(nb_list){
      if(!silent){cat("Working on panel list...\n")}
      nbl <- spdep::mat2listw(nbw, style = "W")
    }
  }
  # ----------------------------------- #


  # ----------------------------------- #
  # Tidy up
  # ----------------------------------- #
  # Convert matrices to sparse matrices to reduce file size and increase
  # computation speed during modeling
  nbw <- as(nbw, "sparseMatrix")

  # Assign column name ID's to spatial weights
  colnames(nbw) <- rownames(nbw)
  # ----------------------------------- #


  # ----------------------------------- #
  # Return statement
  # ----------------------------------- #
  if(nb_list){
    return(list("w_matrix" = nbw,
                "w_list"   = nbl))
  } else{
    return(nbw)
  }
  # ----------------------------------- #
}

