\name{LocLinEst}
\title{Local Linear Hazard Rate Estimator}
\alias{LocLinEst}
\description{Implements the local linear kernel hazard rate estimate of \href{https://ieeexplore.ieee.org/document/4385743}{Bagkavos and Patil (2008)} and \href{https://link.springer.com/article/10.1007/s10463-010-0277-6}{Bagkavos (2011)}. The estimate assumes binned data (fixed design), of the form \eqn{(x_i, y_i)} where \eqn{x_i} are the bin centers and \eqn{y_i} are empirircal hazard rate estimates at each \eqn{x_i}. These are calculated  via  the \code{\link{DiscretizeData}} function. The estimate then smooths the empircal hazard rate estimates and achieves  automatic  boundary adjustments through approrpiately defined kernel weights. The user is able to supply their own bandwidth values through the \eqn{h} argument.

Currently only the \code{\link{LLHRPlugInBand}} bandwidth selector is provided which itself it depends on the \code{\link{bw.nrd}} distribution function default bandwidth rule of \href{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9469.2005.00472.x}{Swanepoel and Van Graan (2005)} for the constant estimate.
\itemize{
\item TO DO: In future implementations the EBBS (empirical bias bandwidth) and AIC based bandwidth methods (see \href{https://link.springer.com/article/10.1007/s10463-010-0277-6}{Bagkavos (2011)}) will be added to the package}
}

\usage{
LocLinEst(BinCenters, xout, h, kfun, ci)
}

\arguments{
\item{BinCenters}{ A vector with the bin centers of the discretized data.}
 \item{xout}{ A vector of points at which the hazard rate function will be estimated. }
 \item{h}{ A scalar, the bandwidth to use in the estimate. }
 \item{kfun}{ Kernel function to use. Supported kernels: Epanechnikov, Biweight, Gaussian, Rectangular, Triangular}
 \item{ci}{ Empirical hazard rate estimates.}
}

\details{
The estimate in both \href{https://ieeexplore.ieee.org/document/4385743}{Bagkavos and Patil (2008)} and \href{https://link.springer.com/article/10.1007/s10463-010-0277-6}{Bagkavos (2011)} is given by
\deqn{
 \hat \lambda_L(x)= \frac{T_{n,1}(x) S_{n,1}(x) - T_{n,0}(x) S_{n,2}(x)}{S_{n,1}(x)S_{n,1}(x)-S_{n,0}(x)S_{n,2}(x)}.
}
The difference between the censored and the uncensored cased is only on the calculation of the empirical hazard rate estimates.
}

\value{
A vector with the values of the function at the designated points xout.
}

\references{
\enumerate{
\item \href{https://ieeexplore.ieee.org/document/4385743}{Bagkavos and Patil, Local Polynomial Fitting in Failure Rate Estimation, IEEE Transactions on Reliability, 57, (2008)},
\item \href{https://link.springer.com/article/10.1007/s10463-010-0277-6}{Bagkavos (2011),  Annals of the Institute of Statistical Mathematics, 63(5), 1019-1046},
}
}

%\author{
%Dimitrios Bagkavos
%
%R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com>
%}

%\note{
%%  ~~further notes~~
%}

\seealso{ \code{\link{HazardRateEst}, \link{LLHRPlugInBand} }}

\examples{
x<-seq(0.05, 5,length=80) #grid points to calculate the estimates
plot(x, HazardRate(x,"weibull", .6, 1),type="l", xlab = "x",ylab="Hazard rate")

SampleSize = 100                 #select sample size
ti<- rweibull(SampleSize, .6, 1) # draw a random sample
ui<-rexp(SampleSize, .2)         # censoring sample
cat("\n AMOUNT OF CENSORING: ", length(which(ti>ui))/length(ti)*100, "\n")
x1<-pmin(ti,ui)                  # observed data
cen<-rep.int(1, SampleSize)      # initialize censoring indicators
cen[which(ti>ui)]<-0             # 0's correspond to censored indicators

a.use<-DiscretizeData(ti, x)     # discretize the data
BinCenters<-a.use$BinCenters     # get the data centers
ci<-a.use$ci                     # get empircal hazard rate estimates
Delta=a.use$Delta                # Binning range
h2<-bw.nrd(ti)                   # Bandwidth to use in constant est. of the plug in rule
h.use<-h2                        # the first element is the band to use

# Calcaculate the plug-in bandwidth:
huse1<- LLHRPlugInBand(BinCenters,h.use,Epanechnikov,Delta,ti,x,IntEpanechnikov,ci, cen)
arg2<-HazardRateEst(x1,x,Epanechnikov, huse1, cen)      # Tanner-Wong Estimate
lines(x, arg2, lty=2) # draw the Tanner-Wong   estimate # Draw TW estimate
arg5<-HazardRateEst(x1,x,BoundaryBiweight,huse1,cen)    # Boundary adjusted TW est
lines(x, arg5, lty=2, col=4) # draw the variable bandwidth # Draw the estimate
arg6<-LocLinEst(BinCenters ,x, huse1, Epanechnikov, ci) # Local linear est.
lines(x, arg6, lty=5, col=5)                             # Draw the estimate
legend("topright", c("Tanner-Wong",  "TW - Boundary Corrected", "Local Linear"),
          lty=c(2,2, 5), col=c(1,4, 5)) # add legend
}
