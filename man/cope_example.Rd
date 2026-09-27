\name{cope_example}
\alias{cope_example}
\title{Load One Arm of the Fixed Synthetic Teaching Dataset}
\usage{cope_example(arm = c("control", "treatment"))}
\arguments{
\item{arm}{Treatment arm to load.}
}
\value{A list of OS and PFS endpoint data frames and matched population
mortality functions.}
\description{Loads fixed bundled synthetic teaching data. Does not simulate
patients or modify the random-number state. This is not clinical trial data.}
\details{See \code{\link{cope_demo}} for the design, provenance and known truth.}
\examples{
d <- cope_example()
head(d$os)
head(d$pfs)
}
