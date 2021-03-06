% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R
\docType{class}
\name{Egger-class}
\alias{Egger-class}
\title{Egger Class}
\description{
An object containing the estimate produced using the MR-Egger method as well as various statistics.

The MR-Egger model uses a random-effects model; a fixed-effect model does not make sense as pleiotropy leads to heterogeneity between the causal estimates targeted by the genetic variants. The (multiplicative) random-effects model allows over-dispersion in the regression model. Under-dispersion is not permitted (in case of under-dispersion, the residual standard error is set to 1).
}
\section{Slots}{

\describe{
\item{\code{Model}}{Model always takes the value \code{random}, as only random-effects analyses are permitted.}

\item{\code{Exposure}}{The name of the exposure variable.}

\item{\code{Outcome}}{The name of the outcome variable.}

\item{\code{Correlation}}{The matrix of correlations between genetic variants.}

\item{\code{Robust}}{Whether robust regression was used in the regression model relating the genetic associations with the outcome and those with the exposure.}

\item{\code{Penalized}}{Whether weights in the regression model were penalized for variants with heterogeneous causal estimates.}

\item{\code{Estimate}}{The causal point estimate from the MR-Egger method.}

\item{\code{StdError.Est}}{The standard error associated with \code{Estimate}.}

\item{\code{Pvalue.Est}}{P-value associated with the causal estimate from the Wald method.}

\item{\code{CILower.Est}}{The lower bound of the confidence interval for \code{Estimate} based on \code{StdError.Est}.}

\item{\code{CIUpper.Est}}{The upper bound of the confidence interval for \code{Estimate} based on \code{StdError.Est}.}

\item{\code{Intercept}}{The intercept estimate from the MR-Egger method. Under the InSIDE assumption, the intercept represents the average pleiotropic effect (average direct effect on the outcome) of a genetic variant. If the intercept differs from zero, this is evidence that the genetic variants are not all valid instruments; specifically, there is directional pleiotropy.}

\item{\code{StdError.Int}}{The standard error associated with \code{Intercept}.}

\item{\code{Pvalue.Int}}{P-value associated with the intercept from the Wald method.}

\item{\code{CILower.Int}}{The lower bound of the confidence interval for \code{Intercept} based on \code{StdError.Int}.}

\item{\code{CIUpper.Int}}{The upper bound of the confidence interval for \code{Estimate} based on \code{StdError.Int}.}

\item{\code{Alpha}}{The significance level used in constructing the confidence interval (default is 0.05).}

\item{\code{SNPs}}{The number of SNPs that were used in the calculation.}

\item{\code{Causal.pval}}{P-value associated with the causal estimate.}

\item{\code{Pleio.pval}}{P-value associated with the intercept (p-value for the MR-Egger intercept test of directional pleiotropy).}

\item{\code{RSE}}{The estimated residual standard error from the regression model.}

\item{\code{Heter.Stat}}{Heterogeneity statistic (Cochran's Q statistic) and associated p-value: the null hypothesis is that the MR-Egger regression model describes the associations with the outcome with no excess heterogeneity.}

\item{\code{I.sq}}{A measure of heterogeneity between the genetic associations with the exposure (see Bowden IJE 2016: "Assessing the suitability of summary data for Mendelian randomization analyses using MR-Egger regression: The role of the I2 statistic."). Low values of \code{I.sq} relate both to large differences in precision between MR-Egger and IVW estimates, and to more weak instrument bias (in a two-sample setting, this is attenuation of MR-Egger estimate towards the null).}
}}

