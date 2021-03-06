% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mr_ivw-methods.R
\name{r.weights}
\alias{r.weights}
\title{Calculates p-values for penalization of weights}
\usage{
r.weights(.Byse, .pen.weights)
}
\arguments{
\item{.Byse}{Standard error of genetic association with the outcome.}

\item{.pen.weights}{Factors for penalizing weights.}
}
\value{
Penalized weights.
}
\description{
Internal function for calculating penalized weights in conjunction with \code{penalised.weights}.

These weights are used in either the \code{mr_ivw} or \code{mr_egger} functions when \code{penalized = TRUE}.
}
\details{
None.
}
\examples{
r.weights(chdloddsse, penalised.weights(ldlc, ldlcse, chdlodds, chdloddsse))

}
\keyword{internal}
