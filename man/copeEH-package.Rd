\name{copeEH-package}
\alias{copeEH}
\alias{copeEH-package}
\docType{package}
\title{Coherent Population-Mortality-Adjusted Survival Extrapolation}
\description{Jointly constrained marginal overall and progression-free survival
extrapolation using Bernstein relative-survival polynomials and a common
population mortality model.}
\details{The method fits marginal curves using a working-independence composite
likelihood. It is not a patient-level joint event-history model and does not
identify transition hazards or cross-endpoint covariance from unlinked data.
The public interface provides point estimates, not confidence intervals.
Start with \code{\link{cope_example}} and \code{\link{fit_cope}}.}
\author{Hao Chen}
\seealso{\code{\link{population_mortality}}, \code{\link{state_outcomes}}}
