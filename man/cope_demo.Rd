\name{cope_demo}
\alias{cope_demo}
\docType{data}
\title{Fixed Synthetic Two-Arm Teaching Dataset with Known Truth}
\usage{data(cope_demo)}
\format{A list with six components:
\describe{
\item{os}{4,000 rows with arm, time in years, and status (1 = death, 0 = censored).}
\item{pfs}{4,000 rows with arm, time and status (1 = progression or death).}
\item{truth}{OS, PFS and PF/PD/death probabilities on a 0--20 year grid.}
\item{population}{Interval breaks and hazards for population_mortality().}
\item{parameters}{True OS/PFS Bernstein coefficient vectors by arm.}
\item{settings}{Sample size, seed, cutoff, horizon, degree, eta, illustrative
utilities, costs, discount rate and willingness-to-pay.}
}}
\description{A fixed teaching dataset, generated within the COPE-EH model
family to illustrate fitting, prediction and comparison to a known truth.
This is not clinical trial data and not evidence of comparative superiority.}
\details{There are 2,000 patients per arm with administrative censoring at
10 years. Degree is three and eta is two years. True control coefficients are
OS (1, 0.60, 0.28, 0.06) and PFS (1, 0.32, 0.08, 0.01); treatment coefficients
are OS (1, 0.72, 0.42, 0.18) and PFS (1, 0.48, 0.18, 0.055).
Population hazards are 0.012, 0.02 and 0.035 per year on intervals starting at
0, 10 and 20. Independent latent death and progression times were sampled by
inverse survival. Death survival was population survival times relative OS;
progression survival was relative PFS divided by relative OS. PFS was the
minimum of progression and death, giving coherent marginal times.
The fixed seed was 2026; it was not selected by prediction error.
The dataset intentionally has substantial sample size and follow-up and a
correctly specified model. Its performance must not be generalised to small,
immature or misspecified clinical datasets. Generator code is not distributed
in the package.}
\source{Artificial teaching data constructed for this package.}
\examples{
data(cope_demo)
head(cope_demo$os)
cope_demo$settings
}
\keyword{datasets}
