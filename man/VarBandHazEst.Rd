\name{VarBandHazEst}
\alias{VarBandHazEst}
\title{Variable Bandwidth Hazard Rate Estimator}

\description{
Implements the adaptive variable bandwidth hazard rate estimator of \href{http://dx.doi.org/10.1080/03610920802364088}{Bagkavos and Patil (2009)}. The estimate itself is an extension of the classical kernel hazard rate estimator of \href{https://projecteuclid.org/download/pdf_1/euclid.aos/1176346265}{Tanner and Wong (1983)} implemented in \code{\link{HazardRateEst}}. The difference is that instead of \eqn{h}, the variable bandwidth estimate uses bandwidth \eqn{h \lambda(X_i)^{-1/2}}. This particular choice cancels the second order term in the bias expansion of the hazard rate estimate and thus it is expected to result in a more precise estimation compared to \code{\link{HazardRateEst}}.
}

\usage{VarBandHazEst(xin, xout, kfun, h1, h2, ci)}

\arguments{
\item{xin}{A vector of data points.  Missing values not allowed.}
 \item{xout}{A vector of points at which the hazard rate function will be estimated. }
\item{kfun}{Kernel function to use. Supported kernels: Epanechnikov, Biweight, Gaussian, Rectangular, Triangular, HigherOrder}
 \item{h1}{A scalar, pilot bandwidth. }
 \item{h2}{A scalar, variable kernel (adaptive) bandwidth.}
 \item{ci}{A vector of censoring indicators: 1's  indicate uncensored observations, 0's correspond to censored obs.}
}

\details{
Implements the adaptive variable bandwidth hazard rate estimator of Bagkavos and Patil (2009), Comm. Statist. Theory and Methods.
\deqn{\hat \lambda_v(x;h_1, h_2) = \sum_{i=1}^n \hat \lambda^{-1/2}(x;h_1 ) \frac{K_{h_2}\left \{ (x-X_{(i)})\hat \lambda^{-1/2}(x;h_1 ) \right \}\delta_{(i)}}{n-i+1}}{\lambda_v(x;h_1, h_2).}
The pilot bandwidth \eqn{h_1} is determined by an optimal bandwidth rule such as \code{\link{PlugInBand}}. and used as input to the pilot kernel estimate, implemented by \code{\link{HazardRateEst}}.

\itemize{
\item TO DO: Insert a  rule for the adaptive bandwidth  \eqn{h_2}. }
}

\value{
A vector with the values of the function at the designated points xout.
}

\references{
\href{http://dx.doi.org/10.1080/03610920802364088}{Bagkavos and Patil (2009), Variable Bandwidths for Nonparametric Hazard Rate Estimation, Communications in Statistics - Theory and Methods, 38:7, 1055-1078}
}

%\author{
%Dimitrios Bagkavos
%
%R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com>
%}

%\note{
%%  ~~further notes~~
%}

\seealso{ \code{\link{HazardRateEst}, \link{TransHazRateEst}, \link{PlugInBand} } }

\examples{
x<-seq(0, 5,length=100) #design points where the estimate will be calculated
plot(x, HazardRate(x, "weibull", .6, 1),  type="l",
     xlab = "x", ylab="Hazard rate") #plot true hazard rate function
SampleSize <- 100
mat<-matrix(nrow=SampleSize, ncol=20)
for(i in 1:20)
{
  ti<- rweibull(SampleSize, .6, 1)#draw a random sample from the actual distribution
  ui<-rexp(SampleSize, .05)       #draw a random sample from the censoring distribution
  cat("\n AMOUNT OF CENSORING: ", length(which(ti>ui))/length(ti)*100, "\n")
  x1<-pmin(ti,ui)                 #this is the observed sample
  cen<-rep.int(1, SampleSize)     #censoring indicators
  cen[which(ti>ui)]<-0            #censored values correspond to zero

  h2<-DefVarBandRule(ti, cen)     #Deafult Band. Rule - Weibull Reference
  huse1<-  PlugInBand(x1, x,   cen, Biweight)
  mat[,i]<- VarBandHazEst(x1, x, Epanechnikov, huse1,h2, cen) #Var. bandwidth est.
}
lines(x, rowMeans(mat) , lty=2)   #draw the average  vb estimate
}
