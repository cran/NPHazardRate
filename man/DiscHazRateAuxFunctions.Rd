\name{nsf, Tm, CparamCalculation, power.matrix, base, SmoothedEstimate}
\title{Auxiliary functions for discrete hazard rate estimators}
\alias{nsf}
\alias{Tm}
\alias{CparamCalculation}
\alias{power.matrix}
\alias{base}
\alias{SmoothedEstimate}

\description{
Auxiliary functions for discrete semiparametric and kernel smooth hazard rate estimation
}

\usage{
 nsf(xin, cens, xout)
 Tm(tk, xout, distribution, par1, par2)
 CparamCalculation(gamparam, VehHazard)
 power.matrix(M, n)
 base(m, b)
 SmoothedEstimate(NonParEst, VehHazard, gammapar, SCproduct, Cpar)
}

\arguments{
\item{xin}{A vector of data points.  Missing values not allowed.}
\item{cens}{A vector of 1s and zeros, 1's indicate uncensored observations, 0's correspond to censored obs. }
\item{xout}{The points where the estimate should be calculated.}
\item{tk}{desing points for the NPMLE estimate.}
\item{distribution}{which distribution to use?}
\item{par1}{distribution parameter 1}
\item{par2}{distribution parameter 2}
\item{gamparam}{gamma parameter}
\item{M}{a matrix to be raised to a power}
\item{n}{the power the matrix will be raised at}
\item{m}{express m as a power of b}
\item{b}{express m as a power of b}
\item{NonParEst}{The crude nonparametric hazard rate estimate.}
\item{VehHazard}{Vehicle hazard rate}
\item{gammapar}{gamma parameter }
\item{SCproduct}{SC product, the result of DetermineSCprod}
\item{Cpar}{C parameter, the result of CparamCalculation.}
}

\details{
Auxiliary functions for discrete hazard rate estimators. The function \code{nsf} is used for the kernel smooth estimate \code{\link{TutzPritscher}}.
\itemize{
\item{Tm}{ used to calculate \eqn{\max(t_k; 1-\sum_{l=0}^k \eta_l > \epsilon), \epsilon>0} in the implementation of the semiparametric estimate}
\item{CparamCalculation}{returns the C smoothing parameter calculated as \eqn{C= \gamma/\max_{k \geq 0} ( \lambda(t_{k-1}) + \lambda(t_k) + \lambda(t_{k+1}) )}}
\item{DetermineSCprod}{this finds \eqn{SC = \gamma((n+1) \hat B_1)^{-1} \hat V_1}  n = number of obs, gammapar = sum of vehicle haz at xout (computed elsewhere)}
}
}

\value{
A vector with the values of the hazard rate estimates.
}

\references{
\enumerate{
\item \href{https://onlinelibrary.wiley.com/doi/abs/10.1002/bimj.201100058}{Patil and Bagkavos (2012), Semiparametric smoothing of discrete failure time data,  Biometrical Journal, 54, (2012), 5--19.}
\item\href{https://doi.org/10.1007/BF00128979}{Tutz, G. and Pritscher, L. Nonparametric Estimation of Discrete Hazard Functions, Lifetime Data Anal, 2, 291-308 (1996)}
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

%\seealso{
%\code{\link{alfa.profile}, \link{alfa}, \link{alfainv}, \link{alpha.mle}
%}
%}

