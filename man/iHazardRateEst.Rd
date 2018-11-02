\name{iHazardRateEst}
\title{Kernel Integrated Hazard Rate Estimation}
\alias{iHazardRateEst}
\description{
Implements the integrated kernel hazard rate estimator for right censored data, i.e. a kernel estimate of the cummulative hazard function.
}

\usage{
iHazardRateEst(xin, xout, ikfun, h, ci)
}

\arguments{
\item{xin}{A vector of data points.  Missing values not allowed.}
 \item{xout}{A vector of grid points at which the estimates will be calculated. }
  \item{ikfun}{Integrated kernel function to use}
 \item{h}{A scalar, the bandwidth to use in the estimate. }
 \item{ci}{A vector of censoring indicators: 1's indicate uncensored observations, 0's correspond to censored obs.}
}

\details{
The function \code{\link{iHazardRateEst}} implements the cummulative hazard rate estimator \eqn{\hat \Lambda(x; h_1)} given by
\deqn{\hat \Lambda(x; h_1) =   \sum_{i=1}^n \frac{k\left \{(x-X_{(i)})h_1^{-1}\right \}\delta_{(i)}}{n-i+1}}{\hat \Lambda(x; h_1) = \int_{-\infty}^x \lambda(t;h_1)dt}
 where
\deqn{k(x) = \int_{-\infty}^x K(y)\,dy}

Note that  \code{\link{iHazardRateEst}} is used in the implementation of the transformed hazard rate estimate  \code{\link{TransHazRateEst}}.
}

\value{
A vector with the cummulative hazard rate estimates at the designated points xout.
}

\references{
\href{https://projecteuclid.org/download/pdf_1/euclid.aos/1176346265}{Tanner and Wong (1983), The Estimation Of The Hazard Function From Randomly Censored Data By The Kernel Method, Annals of Statistics,  3, pp. 989-993.}
}

%\author{
%Dimitrios Bagkavos
%
%R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com>
%}

\seealso{ \code{\link{VarBandHazEst}, \link{TransHazRateEst}, \link{PlugInBand}} }

\examples{
x<-seq(0, 5,length=100) #design points where the estimate will be calculated

SampleSize <- 100
ti<- rweibull(SampleSize, .6, 1)  #draw a random sample from the actual distribution
ui<-rexp(SampleSize, .2)  #draw a random sample from the censoring distribution
cat("\n AMOUNT OF CENSORING: ", length(which(ti>ui))/length(ti)*100, "\n")
x1<-pmin(ti,ui)             #this is the observed sample
cen<-rep.int(1, SampleSize) #censoring indicators
cen[which(ti>ui)]<-0        #censored values correspond to zero
huse<-PlugInBand(x1, x,   cen, Biweight)
arg2<-iHazardRateEst(x1, x, IntEpanechnikov, huse, cen) #Calculate the estimate
}
