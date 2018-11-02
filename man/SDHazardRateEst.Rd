\name{SDHazardRateEst}
\title{Kernel Second Derivative Hazard Rate Estimation}
\alias{SDHazardRateEst}
\description{
Implements the kernel estimate of the second derivative of the hazard rate for right censored data defined - based on the estimate of \href{https://projecteuclid.org/download/pdf_1/euclid.aos/1176346265}{Tanner and Wong (1983)}. The implementation is based on the second derivative of the Biweight Kernel.
}

\usage{
SDHazardRateEst(xin, xout, h, ci)
}

\arguments{
\item{xin}{A vector of data points.  Missing values not allowed.}
 \item{xout}{A vector of grid points at which the estimates will be calculated. }
 \item{h}{A scalar, the bandwidth to use in the estimate. }
 \item{ci}{A vector of censoring indicators: 1's indicate uncensored observations, 0's correspond to censored obs.}
}

\details{
The function \code{\link{SDHazardRateEst}} implements the kernel estimate of the second derivative of the hazard rate estimator, given by
\deqn{\hat \lambda_2(x;h) = \sum_{i=1}^n \frac{K_h''(x-X_{(i)})\delta_{(i)}}{n-i+1}}
where \eqn{K} is taken to be the \code{\link{Biweight}} kernel. The function is used for estimation of the functional \eqn{R(\lambda'')} in \code{\link{PlugInBand}} so a default bandwidth rule is used for \eqn{h} provided in (16), \href{https://onlinelibrary.wiley.com/doi/full/10.1111/anzs.12224}{Hua, Patil and Bagkavos (2018)}.
}


\value{
A vector with the second derivative of the hazard rate at the designated points xout.
}

\references{
\enumerate{
\item \href{https://projecteuclid.org/download/pdf_1/euclid.aos/1176346265}{Tanner and Wong (1983), The Estimation Of The Hazard Function From Randomly Censored Data By The Kernel Method, Annals of Statistics,  3, pp. 989-993.}
\item \href{https://onlinelibrary.wiley.com/doi/full/10.1111/anzs.12224}{Hua, Patil and Bagkavos, An $L_1$ analysis of a kernel-based hazard rate estimator, Australian and New Zealand J. Statist., (60), 43-64, (2018).}
}
}

%\author{
%Dimitrios Bagkavos
%
%R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com>
%}

%\seealso{ \code{\link{VarBandHazEst}, \link{TransHazRateEst}, \link{PlugInBand}} }

%\examples{
%x<-seq(0, 5,length=100) #design points where the estimate will be calculated
%plot(x, HazardRate(x, "weibull", .6, 1),  type="l", xlab = "x",
%                   ylab="Hazard rate") #plot true hazard rate function
%SampleSize <- 100
%ti<- rweibull(SampleSize, .6, 1)  #draw a random sample from the actual distribution
%ui<-rexp(SampleSize, .2)  #draw a random sample from the censoring distribution
%cat("\n AMOUNT OF CENSORING: ", length(which(ti>ui))/length(ti)*100, "\n")
%x1<-pmin(ti,ui)             #this is the observed sample
%cen<-rep.int(1, SampleSize) #censoring indicators
%cen[which(ti>ui)]<-0        #censored values correspond to zero
%huse<-PlugInBand(x1, x,   cen, Biweight)
%arg2<-HazardRateEst(x1, x, Epanechnikov, huse, cen) #Calculate the estimate
%lines(x, arg2, lty=2)       #draw the result on the graphics device.
%}
