\name{HRSurv}
\title{Estimate of the constant in the optimal AMISE expression}
\alias{HRSurv}
\description{
Calculation of the integrand of the contant term in the AMISE plugin bandwidth rule implemented in \code{\link{PlugInBand}}.
}

\usage{HRSurv(x, xin, cens, h, kfun) }

\arguments{
 \item{xin}{ A vector of data points}
 \item{x}{ The point at which the estimates should be calculated.}
\item{cens}{ Censoring Indicators.}
\item{h}{ bandwidth to use.}
\item{kfun}{ The kernel function to use.}
}

\details{
Calculates the term
\deqn{  \frac{\lambda_T(x)}{1-F(x)}\,dx }
which is passed then as argument to the function  \code{\link{NP.M.Estimate}} for numerical integtaion. Currrently the fraction is estimated by
\deqn{\frac{\hat \lambda(x;b)}{1-\hat F(x)}}
where \eqn{\hat \lambda(x;b)} is implemented by \code{\link{HazardRateEst}} using bandwidth \code{bw.nrd{xin}}. For \eqn{1-\hat F(x)} the Kaplan-Meier estimate \code{\link{KMest}} is used.
}

\value{
A vector with the value of the fraction.
}

\references{
\href{https://onlinelibrary.wiley.com/doi/full/10.1111/anzs.12224}{Hua, Patil and Bagkavos, An $L_1$ analysis of a kernel-based hazard rate estimator, Australian and New Zealand J. Statist., (60), 43-64, (2018).}
}

%\author{
%Dimitrios Bagkavos
%
%R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com>
%}



\seealso{ \code{\link{PlugInBand}, \link{NP.M.Estimate} } }

\examples{
x<-seq(0, 5,length=100) #design points where the estimate will be calculated
SampleSize<-100 #amount of data to be generated
ti<- rweibull(SampleSize, .6, 1) # draw a random sample
ui<-rexp(SampleSize, .2)         # censoring sample
cat("\n AMOUNT OF CENSORING: ", length(which(ti>ui))/length(ti)*100, "\n")
x1<-pmin(ti,ui)                  # observed data
cen<-rep.int(1, SampleSize)      # initialize censoring indicators
cen[which(ti>ui)]<-0             # 0's correspond to censored indicators

HRSurv(x, x1, cen, bw.nrd(x1), Biweight)
}
