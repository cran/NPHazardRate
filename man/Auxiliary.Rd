\name{RdistSwitch, PdfSwitch, CdfSwitch, HazardRate}
\title{User driven input for random number generation and pdf, survival and hazard rate function calculation}
\alias{RdistSwitch}
\alias{PdfSwitch}
\alias{CdfSwitch}
\alias{HazardRate}

\description{
Auxiliary functions that help automate the process of random number generation or pdf, survival function or hazard rate functions
}


\usage{
RdistSwitch(dist, SampleSize, par1, par2)
PdfSwitch(xout, dist, par1, par2)
CdfSwitch(xout, dist, par1, par2)
HazardRate(xout, dist, par1, par2)
}

\arguments{
\item{dist}{A string. Corresponds to one of weibull, lognorm, chisquare, exponential, binomial, geometric, poisson, negativebinomial, uniform}
\item{SampleSize}{The size of the random sample to be drawn}
\item{xout}{Grid points where the function will be evaluated}
\item{par1}{parameter 1 of the distirbution}
\item{par2}{parameter 2 of the distirbution}
}

\details{Implements random number generation and density, survival and hazard rate estimates for several distributions. These functions are mainly used when simulating the mean square error etc from known distributions.}

\value{
A vector with the values of the function at the designated points xout or the random numbers drawn.
}

%\references{
%Bagkavos, Dimitrios and Patil, Prakash (2009) Variable Bandwidths for Nonparametric Hazard Rate Estimation, Communications in Statistics - %Theory and Methods, 38, 1055-1078, http://dx.doi.org/10.1080/03610920802364088
%}

%\author{
%Dimitrios Bagkavos
%
%R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com> }

%\note{
%%  ~~further notes~~
%}

%\seealso{ \code{\link{HazardRateEst}, \link{TransHazRateEst} } }


