\name{TutzPritscher}
\title{Discrete non parametric kernel hazard rate estimator}
\alias{TutzPritscher}

\description{
Implementation of the kernel discrete hazard rate estimator of \href{https://doi.org/10.1007/BF00128979}{Tutz and Pritscher (1996)} based on the discrete \code{\link{Habbema}} kernel. The estimate is used for comparison with the semiparametric estimate deveoped in \href{https://doi.org/10.1007/BF00128979}{Tutz and Pritscher (1996)}.
}

\usage{
 TutzPritscher(xin, cens, xout)
}

\arguments{
\item{xin}{A vector of data points.  Missing values not allowed.}
\item{cens}{Censoring indicators as a vector of 1s and zeros, 1's indicate uncensored observations, 0's correspond to censored obs.}
\item{xout}{The grid points where the estimates will be calculated.}
}

\details{
The discrete kernel estimate of \href{https://doi.org/10.1007/BF00128979}{Tutz and Pritscher (1996)} is defined by
\deqn{\hat \lambda(t_m|v) = \sum_{s=1}^q \sum_{i=1}^{m_s} w_m \left ( (t,x), (s, x_{is}) \right )\tilde \lambda(s|x_{is})}{\lambda(t_m|v) =  }
where \eqn{w_m} is the discrete Habbema kernel.
}

\value{
Returns a vector with the values of the hazard rate estimates at \eqn{x=xout}.
}

\references{
\href{https://doi.org/10.1007/BF00128979}{Tutz, G. and Pritscher, L. Nonparametric Estimation of Discrete Hazard Functions, Lifetime Data Anal, 2, 291-308 (1996)}
}



%\author{
%Dimitrios Bagkavos
%
%R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com>
%}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{SemiparamEst}
}
}

\examples{
options(echo=FALSE)
xin<-c(7,34,42,63,64, 74, 83, 84, 91, 108, 112,129, 133,133,139,140,140,146,
      149,154,157,160,160,165,173,176,185, 218,225,241, 248,273,277,279,297,
      319,405,417,420,440, 523,523,583, 594, 1101, 1116, 1146, 1226, 1349,
      1412, 1417)
cens<-c(1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,
        0,1,0,1,1,1,1,1,0,1,1,1,0,1)
xin<-xin/30.438    #Adjust the data
storage.mode(xin)<-"integer"  # turn the data to integers
xout<-seq(1,47, by=1)  # define the grid points to evaluate the estimate
arg<-TutzPritscher(xin,cens,xout)   #Discrete kernel estimate
plot(xout, arg, type="l", ylim=c(0, .35), lty=2,  col=6) # plot the estimate
argSM<-lambdahat(xin, cens, xout)  #crude nonparametric estimate
lines(xout, argSM, lty=3, col=5) # plot the crude estimate
}
