\name{fit_cope}
\alias{fit_cope}
\title{Fit Coherent Marginal OS and PFS Curves}
\usage{
fit_cope(data_os, data_pfs, population, degree = 3, eta = 2,
         tail = c("free", "no_cure"),
         control = list(maxit = 2000, reltol = 1e-10))
}
\arguments{
\item{data_os, data_pfs}{Non-empty data frames with numeric \code{time} and binary
\code{status} columns (1 for an event, 0 for right censoring). PFS events include
progression or death. Times must be finite and non-negative. Rows are not linked
between endpoints. Fit each treatment group separately.}
\item{population}{A list with vectorised \code{survival(t)} and \code{hazard(t)}
functions representing the same population mortality process and the same
target cohort for both endpoints. Survival must start at one.}
\item{degree}{Positive integer Bernstein degree.}
\item{eta}{Positive time-transformation scale, in the same units as the data.
The time map is \eqn{t/(t+\eta)}.}
\item{tail}{\code{"free"} estimates the last relative-survival coefficients;
\code{"no_cure"} fixes both to zero and requires degree at least two.
These are relative-survival limits, not limits of all-cause survival.}
\item{control}{Control list passed to \code{\link[stats]{constrOptim}}.}
}
\value{An object of class \code{cope_fit}. Includes coefficients, degree, eta,
tail specification, composite log-likelihood, optimisation diagnostics,
endpoint sample and event counts, and the population functions for prediction.
The stored Hessian and covariance are working-likelihood diagnostics, not
validated sampling covariance estimates for linked endpoints.}
\description{Fits ordered Bernstein coefficients with a common population
survival component. Coefficient restrictions imply
\eqn{0 \le S_{\mathrm{PFS}}(t) \le S_{\mathrm{OS}}(t) \le 1}.}
\details{Independent right censoring is assumed for each endpoint. The objective
adds two marginal log-likelihoods without recovering patient-level dependence.
The implementation is unpenalised and uses linear constraints with a strictly
feasible start and one fallback attempt.
The caller must ensure that supplied population hazard and survival are
mathematically consistent; basic range and monotonicity checks do not prove
that relationship. Use \code{\link{population_mortality}} for a matched pair.
No automatic model selection or inferential intervals are provided.}
\examples{
d <- cope_example()
fit <- fit_cope(d$os, d$pfs, d$population)
summary(fit)
predict(fit, times = c(0, 1, 5, 10), type = "survival")
}
\seealso{\code{\link{cope_fit_methods}}, \code{\link{state_outcomes}}}
