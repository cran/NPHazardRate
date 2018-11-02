\name{KMest}
\title{Kaplan-Meier Estimate}
\alias{KMest}
\description{
Custom implementation of the Kaplan Meier estimate. The major difference with existing implementations is that the user can specify exactly the grid points where the estimate is calculated. The implementation corresponds to \eqn{1-\hat H(x)} of \href{https://onlinelibrary.wiley.com/doi/full/10.1111/anzs.12224}{Hua, Patil and Bagkavos (2018)}, and is used mainly for estimation of the censoring distribution.
}

\usage{KMest(xin, cens, xout) }

\arguments{
 \item{xin}{ A vector of data points}
 \item{xout}{ The point at which the estimates should be calculated.}
\item{cens}{ Censoring Indicators.}
}

\details{
Calculates the well known Kaplan-Meier estimate
\deqn{1-\hat H(x) =  1,   0 \leq x \leq X_{(1)}} or
\deqn{
1-\hat H(x) = \prod_{i=1}^{k-1} \left ( \frac{n-i+1}{n-i+2} \right )^{1-\delta_{(i)}},   X_{(k-1)} <x \leq X_{(k)}, k=2,\dots,n } or
\deqn{
1-\hat H(x) = \prod_{i=1}^{n} \left ( \frac{n-i+1}{n-i+2} \right )^{1-\delta_{(i)}},  X_{(n)}<x.}
The implementation is mainly for estimating the censoring distribution of the available sample.
}

\value{
A vector with the Kaplan-Meier estimate at xout.
}

\references{
\href{https://www.jstor.org/stable/2281868}{Kaplan, E. L., and Paul Meier. Nonparametric Estimation from Incomplete Observations., J. of the American Statist. Association 53, (1958): 457-81.}
}

%\author{
%Dimitrios Bagkavos
%
%R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com>
%}



%\seealso{ \code{\link{HazardRateEst}, \link{LLHRPlugInBand} } }

\examples{
x<-seq(0, 5,length=100) #design points where the estimate will be calculated
SampleSize<-100 #amount of data to be generated
ti<- rweibull(SampleSize, .6, 1) # draw a random sample
ui<-rexp(SampleSize, .2)         # censoring sample
cat("\n AMOUNT OF CENSORING: ", length(which(ti>ui))/length(ti)*100, "\n")
x1<-pmin(ti,ui)                  # observed data
cen<-rep.int(1, SampleSize)      # initialize censoring indicators
cen[which(ti>ui)]<-0             # 0's correspond to censored indicators

arg1<- KMest(x1, cen, x)
plot(x, arg1, type="l")
}
