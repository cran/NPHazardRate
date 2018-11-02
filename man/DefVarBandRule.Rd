\name{DefVarBandRule}
\alias{DefVarBandRule}
\title{Default adaptive bandwidth rule}

\description{
Implements an adaptive variable bandwidth hazard rate rule for use with the \code{\link{VarBandHazEst}} based on the Weibull distribution, with parameters estimated by maximum likelihood
}

\usage{DefVarBandRule(xin, cens)}

\arguments{
\item{xin}{A vector of data points.  Missing values not allowed.}
 \item{cens}{A vector of censoring indicators: 1's  indicate uncensored observations, 0's correspond to censored obs.}
}

\details{
The adaptive AMISE optimal bandwidth for the variable bandwidth hazard rate estimator \code{\link{VarBandHazEst}} is given by
\deqn{ h_2 =  \left [ \frac{R(K) M_2}{8n\mu_4^2(K) R(g)} \right ]^{1/14}}
 where
 \deqn{ M_2 =   \int \frac{\lambda^{3/2}(x)}{1-F(x)} \,dx}
and
\deqn{ g(x)=\frac{1}{24\lambda(x)^5} \Bigl (24{\lambda'(x)}^4-36{\lambda'(x)}^2{\lambda''(x)}^2\lambda(x)+6{\lambda''(x)}^2\lambda^2(x)
+ 8\lambda'(x)\lambda'''(x)\lambda^2(x) -\lambda^{(4)}(x)\lambda^3(x)\Bigr ) }
}
\value{
the value of the adaptive bandwidth
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
library(survival)
x<-seq(0, 5,length=100) #design points where the estimate will be calculated

SampleSize <- 100

ti<- rweibull(SampleSize, .6, 1)#draw a random sample from the actual distribution
ui<-rexp(SampleSize, .05)       #draw a random sample from the censoring distribution
cat("\n AMOUNT OF CENSORING: ", length(which(ti>ui))/length(ti)*100, "\n")
x1<-pmin(ti,ui)                 #this is the observed sample
cen<-rep.int(1, SampleSize)     #censoring indicators
cen[which(ti>ui)]<-0            #censored values correspond to zero

h2<-DefVarBandRule(ti, cen)     #Deafult Band. Rule - Weibull Reference
}

