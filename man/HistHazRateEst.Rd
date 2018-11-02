\name{HazardHistogram}
\title{Histogram Hazard Rate Estimator}
\alias{HazardHistogram}
\description{
Implements the histogram hazard rate estimator of \href{https://www.jstor.org/stable/42003749}{Patil and Bagkavos (2012)}
}

\usage{
HazardHistogram(xin, xout, cens, bin)
}

\arguments{
\item{xin}{ A vector of data points.  Missing values not allowed.}
\item{xout}{A vector of grid points at which the histogram will be calculated.}
 \item{cens}{ A vector of 1s and zeros, 1's indicate uncensored observations, 0's correspond to censored obs. }
\item{bin}{ Number of bins to use in construction of the histogram.}
}

\details{
The histogram hazard rate estimator is defined in (1), \href{https://www.jstor.org/stable/42003749}{Patil and Bagkavos (2012)} by
\deqn{
\hat \lambda (x) = h_n^{-1} C_{i_{(x)}} = h_n^{-1}f_{i_{(x)}}^0(\bar F_{i_{(x)}}+1)^{-1}.
}}

\value{
A vector with the values of the histogram estimate at each bin.
}

\references{
\href{https://www.jstor.org/stable/42003749}{Patil and Bagkavos (2012), Histogram for hazard rate estimation,  pp. 286-301, Sankhya, B.}

}

%\author{
%Dimitrios Bagkavos
%
%R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com>
%}

%\note{
%%  ~~further notes~~
%}

%\seealso{
%\code{\link{alfa.profile}, \link{alfa}, \link{alfainv}, \link{alpha.mle}
%}
%}

\examples{
SampleSize <-400
ti<-rweibull(SampleSize,0.5,0.8)
xout<-seq(0.02, 3.5, length=80)
true.hazard<-dweibull(xout,0.5, 0.8)/(1-pweibull(xout, 0.5, 0.8))
cen<-rep.int(1, SampleSize)
cen[sample(1:SampleSize, SampleSize/10)]<-0

band<-nlminb(start= 2, obj=cvfunction, control = list(iter.max = 100, x.tol = .001)
            ,xin=ti, xout= xout, cens = cen, 	lower=.01, upper=max(xout))
bin<- 3.49 * sd(ti)^2 * SampleSize^(-1/3) /50 #Scott 1979 Biometrika default rule
bin<-unlist(band[1])
histest<- HazardHistogram(ti,xout, cen, bin+0.013 )

plot(xout, true.hazard, type="l")
lines(histest[,1], histest[,2], col=2, type="s")
barplot( histest[,2], rep(bin, times=length(histest[,2])))
lines(xout, true.hazard, type="l", lwd=2, col=2)
}
