% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/impacts_e.R
\name{impacts_e}
\alias{impacts_e}
\title{Estimate spatial (marginal) effects}
\usage{
impacts_e(
  model,
  eigen_vals,
  w = NULL,
  sims = 1000,
  cred_int = 0.95,
  vcv = NULL,
  params = NULL,
  ignore_params = NULL
)
}
\arguments{
\item{model}{an nls model object estimated with lagsarnls()}

\item{eigen_vals}{pre-computed eigen values from a spatial weights matrix}

\item{w}{not used yet}

\item{sims}{number of simulation iterations}

\item{cred_int}{credibility interval to calculate}

\item{vcv}{option to supply a variance-covariance matrix. Useful to
estimate spatial effects with robust or clustered standard errors. Default
uses model variance-covariance matrix extracted using generic vcv(model)}

\item{params}{vector of parameters to estimate spatial effects}

\item{ignore_params}{vector of parameters to ignore. Note - if "params" option
is specified, this option is ignored and only results for those variables
in "params" are returned. If "params" is empty, the function will estimate
spatial effects for all variables in a model except those listed in this
option. If both params and ignore_params are null, the function estimates
spatial effects for all variables in the model.}
}
\value{
Returns a tidied dataframe containing spatial effects for all or
selected variable subsets from a spatial lag model.
}
\description{
impacts_e uses parametric simulation to estimate spatial effects
for parameters in a model containing an endogenous spatial lag.
}
\examples{
x <- spacetime_sim2(phi_true = 0)
d <- x$data
w <- x$W_matrix
e <- eigen(w, only.values = TRUE)$values
m <- lagsarnls(y ~ x1 + x2, data = d, W = w)
impacts_e(model      = m,
          eigen_vals = e)

}
