\name{lambdahat}
\title{Discrete non parametric mle hazard rate estimator}
\alias{lambdahat}

\description{
Implementation of the purely nonparametric discrete hazard rate estimator \code{lambdahat} discussed among others in \href{https://onlinelibrary.wiley.com/doi/abs/10.1002/bimj.201100058}{Patil and Bagkavos (2012)}.  \code{lambdahat} is also used as the nonparametric component in the implementation of \code{\link{SemiparamEst}}.
}

\usage{
 lambdahat(xin, cens, xout)
}

\arguments{
\item{xin}{A vector of data points.  Missing values not allowed.}
\item{cens}{Censoring indicators as a vector of 1s and zeros, 1's indicate uncensored observations, 0's correspond to censored obs.}
\item{xout}{The grid points where the estimates will be calculated.}
}

\details{
The discrete - crude - hazard rate estimator (NPMLE) in \href{https://onlinelibrary.wiley.com/doi/abs/10.1002/bimj.201100058}{Patil and Bagkavos (2012)} is given by
\deqn{\hat \lambda(t_k) =   \frac{n^0_k}{m_k+1} }{\lambda(t_k) =n^0_k * (m_k+1)^{-1}}
}

\value{
Returns a vector with the values of the hazard rate estimates at \eqn{x=xout}.
}

\references{
 \href{https://onlinelibrary.wiley.com/doi/abs/10.1002/bimj.201100058}{Patil and Bagkavos (2012), Semiparametric smoothing of discrete failure time data,  Biometrical Journal, 54, (2012), 5--19.}
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
