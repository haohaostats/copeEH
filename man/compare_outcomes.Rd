\name{compare_outcomes}
\alias{compare_outcomes}
\title{Compare Deterministic Outcomes Between Treatment Groups}
\usage{compare_outcomes(treatment, control, wtp, incremental_cost = 0)}
\arguments{
\item{treatment, control}{Single-row outputs from \code{\link{state_outcomes}},
using the same horizon, discount rate, currency and utility scale.}
\item{wtp}{Positive willingness-to-pay per QALY in the same currency as costs.}
\item{incremental_cost}{Additional discounted treatment-minus-control costs
not included in state costs. May be negative.}
}
\value{One row containing incremental life-years, QALYs, costs, ICER, INMB,
and a descriptive cost-effectiveness quadrant classification.}
\description{Computes treatment-minus-control differences and
\eqn{\mathrm{INMB}=\lambda\Delta Q-\Delta C}.}
\details{The ICER is returned only when incremental costs and QALYs have the
same non-zero sign. Otherwise it is \code{NA}, with the classification providing
the relevant interpretation. An ICER with both increments negative represents
savings per QALY lost. INMB is a deterministic value, not a probability.
Exactly zero effects or costs require interpretation alongside the increments.}
\examples{
d <- cope_example()
fit <- fit_cope(d$os, d$pfs, d$population)
x <- state_outcomes(fit, horizon = 10)
compare_outcomes(x, x, wtp = 100000, incremental_cost = 5000)
}
