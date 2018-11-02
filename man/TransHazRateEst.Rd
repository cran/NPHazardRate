\name{TransHazRateEst}
\alias{TransHazRateEst}
\title{Transformation Based Hazard Rate Estimator}

\description{
Implements the transformated kernel hazard rate estimator of \href{http://dx.doi.org/10.1080/10485250802440184}{Bagkavos (2008)}. The estimate is expected to have less bias compared to the ordinary kernel estimate \code{\link{HazardRateEst}}. The estimate results by first transforming the data to a sample from the exponential distribution through the integrated hazard rate function, estimated by \code{\link{iHazardRateEst}} and uses the result as input to the classical kernel hazard rate estimate  \code{\link{HazardRateEst}}. An inverse transform turn the estimate to a hazard rate estimate of the original sample. See section "Details" below.
}

\usage{
TransHazRateEst(xin, xout, kfun, ikfun, h1, h2, ci)
}

\arguments{
\item{xin}{A vector of data points.  Missing values not allowed.}
 \item{xout}{A vector of points at which the hazard rate function will be estimated. }
\item{kfun}{Kernel function to use. Supported kernels: Epanechnikov, Biweight, Gaussian, Rectangular, Triangular, HigherOrder.}
\item{ikfun}{An integrated kernel function to use. Supported kernels: Epanechnikov, Biweight, Gaussian, Rectangular, Triangular, HigherOrder.}
\item{h1}{A scalar, pilot bandwidth.}
\item{h2}{A scalar, transformed kernel bandwidth.}
 \item{ci}{A vector of censoring indicators: 1's indicate uncensored observations, 0's correspond to censored obs. }
}

\details{
The transformed kernel hazard rate estimate of \href{http://dx.doi.org/10.1080/10485250802440184}{Bagkavos (2008)} is given by
\deqn{\hat \lambda_t(x;h_1, h_2) = \sum_{i=1}^n  \frac{K_{h_2}\left \{ (\hat \Lambda(x;h_1 ) - \hat \Lambda(X_{(i)};h_1 ) ) \right \}\delta_{(i)}}{n-i+1}\hat \lambda(x;h_1 ).}{\lambda_t(x;h_1, h_2).}
The estimate uses the classical kernel hazard rate estimate \eqn{\lambda(x; h_1)} implemented in \code{\link{HazardRateEst}} and its integrated version
\deqn{\hat \Lambda(x; h_1) =  \sum_{i=1}^n \frac{k\left \{(x-X_{(i)})h_1^{-1}\right \}\delta_{(i)}}{n-i+1}}{\Lambda(x; h_1) = \int_{-\infty}^x \lambda(t;h_1)dt}
 where
\eqn{k(x) = \int_{-\infty}^x K(y)\,dy}
implemented in \code{\link{iHazardRateEst}}. The pilot bandwidth \eqn{h_1} is determined by an optimal bandwidth rule such as \code{\link{PlugInBand}}.
\itemize{
\item TO DO: Insert a  rule for the adaptive bandwidth  \eqn{h_2}. }
}

\value{
A vector with the values of the function at the designated points xout.
}

\references{
\href{http://dx.doi.org/10.1080/10485250802440184}{Bagkavos (2008), Transformations in hazard rate estimation, J. Nonparam. Statist., 20, 721-738}
}

%\author{
%Dimitrios Bagkavos\cr
%R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com>
%}

%\note{
%%  ~~further notes~~
%}

\seealso{ \code{\link{VarBandHazEst}, \link{HazardRateEst}, \link{PlugInBand} } }

\examples{
x<-seq(0, 5,length=100) #design points where the estimate will be calculated
plot(x, HazardRate(x, "weibull", .6, 1),  type="l",
                   xlab = "x", ylab="Hazard rate") #plot true hazard rate function
SampleSize <- 100
mat<-matrix(nrow=SampleSize, ncol=20)
for(i in 1:20)
{ #Calculate the average of 20 estimates and draw on the screen
 ti<- rweibull(SampleSize, .6, 1)  #draw a random sample from the actual distribution
 ui<-rexp(SampleSize, .05)          #draw a random sample from the censoring distribution
 cat("\n AMOUNT OF CENSORING: ", length(which(ti>ui))/length(ti)*100, "\n")
 x1<-pmin(ti,ui)                   #this is the observed sample
 cen<-rep.int(1, SampleSize)       #censoring indicators
 cen[which(ti>ui)]<-0              #censored values correspond to zero

 h2<-DefVarBandRule(ti, cen)     #Deafult Band. Rule - Weibull Reference
 huse1<- PlugInBand(x1, x,   cen, Biweight) #
 mat[,i]<-TransHazRateEst(x1,x,Epanechnikov,IntEpanechnikov,huse1,h2,cen)
}
lines(x, rowMeans(mat) , lty=2) #draw the average transformed estimate
}
