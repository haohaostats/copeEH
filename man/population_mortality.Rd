\name{population_mortality}
\alias{population_mortality}
\title{Matched Piecewise-Constant Population Mortality Functions}
\usage{population_mortality(breaks = 0, hazard)}
\arguments{
\item{breaks}{Strictly increasing interval start times, starting at zero.}
\item{hazard}{One finite non-negative hazard rate per interval.}
}
\value{A list containing vectorised \code{survival}, \code{hazard} and
\code{cumhazard} functions, and the supplied breaks and rates.}
\description{Constructs a matched hazard, integrated hazard and survival pair.}
\details{Intervals are left closed and right open. The final hazard is held
constant indefinitely; this is an explicit extrapolation assumption.
Time units must agree with those of the endpoint data.
For annual death probabilities, transform them using \code{-log1p(-qx)}.
This helper is not a life-table downloader or automatic cohort standardiser.
For fractional baseline ages, the caller must supply the corresponding
attained-age interval boundaries rather than rounding age.}
\examples{
pop <- population_mortality(c(0, 1, 2), c(0.01, 0.02, 0.03))
pop$survival(c(0, 1, 2, 3))
pop$hazard(c(0, 1, 2, 3))
}
