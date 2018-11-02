\name{SemiparamEst}
\title{Discrete hazard rate estimator}
\alias{SemiparamEst}


\description{
Implements the semiparametric hazard rate estimator for discrete data developed in \href{https://onlinelibrary.wiley.com/doi/abs/10.1002/bimj.201100058}{Patil and Bagkavos (2012)}. The estimate is obtained by semiparametric smoothing of the (nonsmooth) nonparametric maximum likelihood estimator, which is achieved by repeated multiplication of a Markov chain transition-type matrix. This matrix is constructed with basis
a parametric discrete hazard rate model (vehicle model). }

\usage{
SemiparamEst(xin, cens, xout, Xdistr, Udistr, vehicledistr, Xpar1=1, Xpar2=0.5,
              Upar1=1, Upar2=0.5, vdparam1=1, vdparam2=0.5)
}

\arguments{
\item{xin}{A vector of data points.  Missing values not allowed.}
\item{cens}{Censoring indicators as a vector of 1s and zeros, 1's indicate uncensored observations, 0's correspond to censored obs. }
\item{xout}{Design points where the estimate will be calculated.}
\item{Xdistr}{The distribution where the data are coming from, currently ignored}
\item{Udistr}{Censoring distribution, currently ignored}
\item{vehicledistr}{String specifying the vehicle hazard rate (the assumed parametric model)}
\item{Xpar1}{Parameter 1 for the X distr, currently ignored}
\item{Xpar2}{Parameter 2 for the X distr, currently ignored}
\item{Upar1}{Parameter 1 for the Cens. distr., currently ignored}
\item{Upar2}{Parameter 2 for the Cens. distr., currently ignored}
\item{vdparam1}{Parameter 1 for the vehicle hazard rate. }
\item{vdparam2}{Parameter 2 for the vehicle hazard rate.}
}

\details{
The semiparmaetric estimator implemented is defined in (1) in \href{https://onlinelibrary.wiley.com/doi/abs/10.1002/bimj.201100058}{Patil and Bagkavos (2012)} by
\deqn{\tilde \lambda  = \hat \lambda \Gamma^S}{\tilde \lambda  = \hat \lambda \Gamma^S}
where \eqn{S} determines the number of repetions and hence the amount of smoothing applied to the estimate. For \eqn{S=0} the semiparametric estimate equals the nonparmaetric estimate \code{\link{lambdahat}}. On the other hand, if the true unknown underlying probability model is known (up to an unknown constant or constants) then, the greater the \eqn{S}, the closer the semiparmaetric estimate to the vehicle hazard rate model.

\itemize{
\item  TO DO: The extension to hazard rate estimation with covariates will be added in a future release.
\item TO DO: Also, the data driven estimation of the parameter \eqn{S} will be also added in a future release; this will inlcude the \eqn{SC} product and \eqn{C} and \eqn{\gamma} parameter calculations.
}
}


\value{
A vector with the values of the discrete hazard rate estimate, calculated at \eqn{x=xout}.
}

\references{
\href{https://onlinelibrary.wiley.com/doi/abs/10.1002/bimj.201100058}{Patil and Bagkavos (2012), Semiparametric smoothing of discrete failure time data,  Biometrical Journal, 54, (2012), 5-19}
}

%\author{
%Dimitrios Bagkavos
%
%R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com>
%}

%\note{
%%  ~~further notes~~
%}

\seealso{\code{\link{lambdahat}, \link{TutzPritscher}}}

\examples{
options(echo=FALSE)
xin<-c(7,34,42,63,64, 74, 83, 84, 91, 108, 112,129, 133,133,139,140,140,146,
      149,154,157,160,160,165,173,176,185, 218,225,241, 248,273,277,279,297,
      319,405,417,420,440, 523,523,583, 594, 1101, 1116, 1146, 1226, 1349,
      1412, 1417) #head and neck data set
cens<-c(1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
        1,0,1,1,1,1,1,1,0,1,0,1,1,1,1,1,0,1,1,1,0,1) #censoring indicators
xin<-xin/30.438 #mean adjust the data
storage.mode(xin)<-"integer"  # turn the data to integers
xout<-seq(1,47, by=1) #design points where to calculate the estimate
arg<-TutzPritscher(xin,cens,xout) #Kernel smooth estimate
plot(xout, arg, type="l", ylim=c(0, .35), lty=2,  col=6)
argSM<-SemiparamEst(xin, cens, xout, "geometric", "uniform",
                    "geometric", 0.2, .6, 0, 90, .25, .9) #semipar. est.
lines(xout, argSM[,2], lty=3, col=5) #add tilde lambda to the plot
}
