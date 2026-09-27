\name{state_outcomes}
\alias{state_outcomes}
\title{Integrate Health-State Times, QALYs and State Costs}
\usage{
state_outcomes(object, horizon, utility_pf = 0.8, utility_pd = 0.65,
               annual_cost_pf = 0, annual_cost_pd = 0,
               discount = 0.03, step = 0.02)
}
\arguments{
\item{object}{A fitted \code{cope_fit}.}
\item{horizon}{Positive horizon in years. Endpoint times and eta must therefore
also have been expressed in years.}
\item{utility_pf, utility_pd}{Non-negative finite constant state utility weights.
Defaults are illustrative, not clinical recommendations.}
\item{annual_cost_pf, annual_cost_pd}{Non-negative annual state costs in a
common currency.}
\item{discount}{Non-negative effective annual discount rate for QALYs and costs.}
\item{step}{Positive integration step in years; the exact horizon is always
included. Reduce it to assess numerical sensitivity.}
}
\value{A one-row data frame with horizon, discount rate, undiscounted PF and PD
years, undiscounted life-years, discounted QALYs, and discounted state costs.}
\description{Uses trapezoidal integration from zero to the supplied horizon.}
\details{Discounting is \eqn{(1+r)^{-t}}. These are restricted horizon outcomes,
not necessarily lifetime outcomes. Drug acquisition, administration,
adverse-event, subsequent-treatment and terminal costs are not automatically
included. Add independently calculated discounted incremental costs through
\code{\link{compare_outcomes}}. Death has zero utility and state cost.}
\examples{
d <- cope_example()
fit <- fit_cope(d$os, d$pfs, d$population)
state_outcomes(fit, horizon = 10, annual_cost_pf = 20000,
               annual_cost_pd = 100000)
}
