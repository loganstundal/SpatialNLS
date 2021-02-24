#' Helper function to generate W in spacetime_sim
#'
#' @description Not exported. Function borrowed (directly) from Davide Martinetti
#'   and Ghislain Geniaux's package "ProbitSpatial" which is no longer (2021-02-23)
#'   hosted on CRAN. All code credit belongs to them.
#'
#' @param n      number of observations in W
#' @param nneigh number of nearerst neighbors
#' @param seed   simulation seed for setting random draws
#'
#' @return Returns sparse weights matrix

simulate_W <- function (n, nneigh, seed = 123)
{
  set.seed(seed)
  coord <- cbind(runif(n, 0, 1), runif(n, 0, 1))
  k = nneigh + 1
  nb1 <- RANN::nn2(as.matrix(coord), k = k, treetype = c("bd"))
  W <- Matrix::sparseMatrix(i = rep(seq_along(rep(k, n)),
                                    rep(k, n)), j = t(nb1$nn.idx), x = 1/k)
  diag(W) <- 0
  if (class(W) == "matrix")
    W <- Matrix::Matrix(W)
  ret <- Matrix::summary(as(W, "dgCMatrix"))
  ti <- tapply(ret$x, ret$i, function(x) sum(x, na.rm = TRUE))
  ret$x <- as.numeric(ret$x/ti[match(ret$i, as.numeric(names(ti)))])
  W <- Matrix::sparseMatrix(i = ret$i, j = ret$j, x = ret$x,
                            dims = dim(W))
  Matrix::drop0(W)
  W
}
