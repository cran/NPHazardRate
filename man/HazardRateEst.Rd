\name{HazardRateEst}
\title{Kernel Hazard Rate Estimation}
\alias{HazardRateEst}
\description{
Implements the (classical) kernel hazard rate estimator for right censored data defined in \href{https://projecteuclid.org/download/pdf_1/euclid.aos/1176346265}{Tanner and Wong (1983)}.
}

\usage{
HazardRateEst(xin, xout, kfun, h, ci)
}

\arguments{
\item{xin}{A vector of data points.  Missing values not allowed.}
 \item{xout}{A vector of grid points at which the estimate will be calculated. }
 \item{kfun}{Kernel function to use. Supported kernels: Epanechnikov, Biweight, Gaussian, Rectangular, Triangular, HigherOrder. }
 \item{h}{A scalar, the bandwidth to use in the estimate. }
 \item{ci}{A vector of censoring indicators: 1's indicate uncensored observations, 0's correspond to censored obs.}
}

\details{
The kernel hazard rate estimator of \href{https://projecteuclid.org/download/pdf_1/euclid.aos/1176346265}{Tanner and Wong (1983)} is given by
\deqn{\hat \lambda(x;h) = \sum_{i=1}^n \frac{K_h(x-X_{(i)})\delta_{(i)}}{n-i+1}}{\lambda(x;h) = sum(K_h(x-X_i)*\delta_i*(n-i+1)^{-1})}
\eqn{h} is determined by  a bandwidth rule such as \code{\link{PlugInBand}}. \code{\link{HazardRateEst}} is also used as a pilot estimate in the implementation of both the variable bandwidth estimate \code{\link{VarBandHazEst}} and the transformed hazard rate estimate  \code{\link{TransHazRateEst}}.
}

\value{
A vector with the hazard rate estimates at the designated points xout.
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
plot(x, HazardRate(x, "weibull", .6, 1),  type="l", xlab = "x",
                   ylab="Hazard rate") #plot true hazard rate function
SampleSize <- 100
ti<- rweibull(SampleSize, .6, 1)  #draw a random sample from the actual distribution
ui<-rexp(SampleSize, .2)  #draw a random sample from the censoring distribution
cat("\n AMOUNT OF CENSORING: ", length(which(ti>ui))/length(ti)*100, "\n")
x1<-pmin(ti,ui)             #this is the observed sample
cen<-rep.int(1, SampleSize) #censoring indicators
cen[which(ti>ui)]<-0        #censored values correspond to zero
huse<-PlugInBand(x1, x,   cen, Biweight)
arg2<-HazardRateEst(x1, x, Epanechnikov, huse, cen) #Calculate the estimate
lines(x, arg2, lty=2)       #draw the result on the graphics device.
}
