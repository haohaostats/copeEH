\name{run_cope_demo}
\alias{run_cope_demo}
\alias{plot.cope_demo_result}
\title{Analyse the Fixed Teaching Dataset and Plot Truth Comparisons}
\usage{
run_cope_demo()
\method{plot}{cope_demo_result}(x, type = c("survival", "state"), ...)
}
\arguments{
\item{x}{The result of run_cope_demo().}
\item{type}{Plot marginal survival or PF/PD state probabilities.}
\item{...}{Reserved; additional plot arguments are rejected.}
}
\value{run_cope_demo() returns a cope_demo_result containing the fixed data,
fitted models, predictions, observed/extrapolated curve errors, state and
economic outcome errors, and true/estimated incremental economic outcomes.
Plot returns its input invisibly.}
\description{Fits each treatment arm without using the stored truth to select
or estimate parameters. Evaluates predictions against the truth only after
fitting. No data generation is performed.}
\details{The plotting method uses a four-panel layout with a shared bottom
legend. It resets graphical parameters on exit but is intended for a dedicated
graphics device rather than an existing layout. Shading indicates extrapolation
beyond ten years; dashed curves show the known truth. Observed KM curves are
included in survival panels. Economic settings are illustrative and fixed;
see \code{\link{cope_demo}} for the matched-model teaching design.}
\examples{
result <- run_cope_demo()
result$accuracy
result$outcomes
result$economic
plot(result)
}
