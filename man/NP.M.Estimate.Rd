\name{NP.M.Estimate}
\title{Estimate of bandwidth constant}
\alias{NP.M.Estimate}
\description{
Calculation of the contant term in the AMISE plugin bandwidth rule \code{\link{PlugInBand}}.
}

\usage{NP.M.Estimate(xin, cens, xout) }

\arguments{
 \item{xin}{ A vector of data points}
 \item{xout}{ The point at which the estimates should be calculated.}
\item{cens}{ Censoring Indicators.}
}

\details{
Approximates the term
\deqn{
 M = \int_0^T  \frac{ \lambda_T(x)  }{1-F(x)} \,dx
}
 which is needed in the optimal AMISE bandwidth expression of \code{\link{PlugInBand}}. The integrand
\deqn{  \frac{\lambda_T(x)}{1-F(x)}\,dx }
is calculated by \code{\link{HRSurv}} and integration is performed via the extended Simpson's numerical integration rule (\code{\link{SimpsonInt}}).
}

\value{
A scalar with the value of the constant.
}

\references{
\href{https://onlinelibrary.wiley.com/doi/full/10.1111/anzs.12224}{Hua, Patil and Bagkavos, An $L_1$ analysis of a kernel-based hazard rate estimator, Australian and New Zealand J. Statist., (60), 43-64, (2018).}
}

%\author{
%Dimitrios Bagkavos
%
%R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com>
%}



%\seealso{ \code{\link{HazardRateEst}, \link{LLHRPlugInBand} } }

%\examples{
%x<-seq(0, 5,length=100) #design points where the estimate will be calculated
%SampleSize<-100 #amount of data to be generated
%ti<- rweibull(SampleSize, .6, 1) # draw a random sample
%ui<-rexp(SampleSize, .2)         # censoring sample
%cat("\n AMOUNT OF CENSORING: ", length(which(ti>ui))/length(ti)*100, "\n")
%x1<-pmin(ti,ui)                  # observed data
%cen<-rep.int(1, SampleSize)      # initialize censoring indicators
%cen[which(ti>ui)]<-0             # 0's correspond to censored indicators



%huse1<- PlugInBand(x1, x, cen, Biweight)
%huse1


%}
