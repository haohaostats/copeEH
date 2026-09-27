\name{cope_fit_methods}
\alias{cope_fit_methods}
\alias{coef.cope_fit}
\alias{predict.cope_fit}
\alias{print.cope_fit}
\alias{summary.cope_fit}
\alias{print.summary.cope_fit}
\alias{plot.cope_fit}
\title{Summarise, Predict and Plot a COPE-EH Fit}
\usage{
\method{coef}{cope_fit}(object, free = FALSE, ...)
\method{predict}{cope_fit}(object, times,
  type = c("all", "survival", "state"), ...)
\method{print}{cope_fit}(x, ...)
\method{summary}{cope_fit}(object, ...)
\method{print}{summary.cope_fit}(x, ...)
\method{plot}{cope_fit}(x,
  times = seq(0, max(5, 2 * x$observed_end), length.out = 301),
  type = c("survival", "state"), xlab = "Time", ylab = "Probability",
  col = c("#0072B2", "#D55E00", "#666666"), lty = c(1, 2, 3),
  lwd = 2, legend = TRUE, ...)
}
\arguments{
\item{object, x}{A \code{cope_fit}; the summary print method accepts its summary.}
\item{free}{Return only estimated coefficients, excluding fixed endpoints.}
\item{times}{Finite non-negative times in the units used for fitting. Prediction
preserves their order. Plotting sorts and deduplicates them.}
\item{type}{Select survival, health-state probabilities, or all prediction
columns. The plotting method accepts survival or state.}
\item{xlab, ylab}{Axis labels.}
\item{col, lty, lwd}{Line colours, types and widths.}
\item{legend}{Whether to show the legend.}
\item{...}{Additional graphical arguments for plotting. Extra prediction
arguments are rejected; other methods reserve this argument for S3 compatibility.}
}
\value{Coefficients are a named vector. Summary returns endpoint counts and a
coefficient table with fit information. Print returns its input invisibly.
Prediction returns a data frame. \code{type = "all"} includes relative and
all-cause OS/PFS, PF/PD/death probabilities, excess hazards and total hazards.
Plot returns the plotted predictions invisibly.}
\description{Standard R methods for the public fitted-model class.}
\details{State probabilities are PFS, OS minus PFS, and one minus OS.
Only point predictions are exposed. No \code{vcov}, \code{confint} or likelihood
ratio inference is provided. A non-converged fit cannot be used for prediction.}
\examples{
d <- cope_example()
fit <- fit_cope(d$os, d$pfs, d$population)
coef(fit)
plot(fit, times = seq(0, 20, by = 0.1))
}
