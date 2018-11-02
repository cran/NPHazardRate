\name{sn.i, tn.i}
\title{Local kernel weights}
\alias{sn.0}
\alias{sn.1}
\alias{sn.2}
\alias{sn.3}
\alias{sn.4}
\alias{sn.5}
\alias{sn.6}
\alias{tn.0}
\alias{tn.1}
\alias{tn.2}
\alias{tn.3}
\description{
Implements the local kernel weights which are used in the implementation of \code{\link{LocLinEst}} and the second derivative estimate used in \code{\link{PlugInBand}}.
}

\usage{
sn.0(xin, xout, h, kfun)
sn.1(xin, xout, h, kfun)
sn.2(xin, xout, h, kfun)
sn.3(xin, xout, h, kfun)
sn.4(xin, xout, h, kfun)
sn.5(xin, xout, h, kfun)
sn.6(xin, xout, h, kfun)
tn.0(xin, xout, h, kfun, Y)
tn.1(xin, xout, h, kfun, Y)
tn.2(xin, xout, h, kfun, Y)
tn.3(xin, xout, h, kfun, Y)
}

\arguments{
\item{xin}{ A vector of data points, typicaly these are the bin centers.  Missing values not allowed.}
 \item{xout}{ A vector of data points where the estimate will be evaluated.}
\item{h}{A scalar. The bandwidth to use.}
\item{kfun}{The kernel function to use.}
\item{Y}{Empirical hazard rate estimates.}
}

\details{
The functions calculate the quantities
\deqn{S_{n,l}(x) = \sum_{i=1}^n K \left (\frac{x_i-x}{h}\right ) (x_i-x)^l, l=0,\dots,6 }
and
\deqn{T_{n,l}(x) = \sum_{i=1}^n K \left (\frac{x_i-x}{h}\right ) (x_i-x)^l Y_i, l=0,\dots,3}
These qunatities are used to adjust the hazard rate estimate and its second derivative in the boundary.
}

\value{
The weight of the functional at \eqn{x}
}

\references{
\enumerate{
\item \href{https://ieeexplore.ieee.org/document/4385743}{Bagkavos and Patil, Local Polynomial Fitting in Failure Rate Estimation, IEEE Transactions on Reliability, 57, (2008)},
\item \href{https://link.springer.com/article/10.1007/s10463-010-0277-6}{Bagkavos (2011),  Annals of the Institute of Statistical Mathematics, 63(5), 1019-1046}.
}
}

%\author{
%Dimitrios Bagkavos
%
%R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com>
%}


