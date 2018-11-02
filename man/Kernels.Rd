\name{Kernels}
\alias{Biweight}
\alias{Epanechnikov}
\alias{Triangular}
\alias{Gaussian}
\alias{HigherOrder}
\alias{Rectangular}
\alias{IntBiweight}
\alias{IntEpanechnikov}
\alias{IntRectangular}
\alias{IntTriangular}
\alias{IntGaussian}
\alias{SDBiweight}
\alias{a0}
\alias{a1}
\alias{a2}
\alias{BoundaryBiweight}
\alias{b0}
\alias{b1}
\alias{b2}
\alias{BoundaryEpanechnikov}
\alias{Habbema}

\title{Kernel functions}

\description{
Implements various kernel functions, including boundary, integrated and discrete kernels for use in the definition of the nonparametric estimates
}

\usage{
Biweight(x, ...)
Epanechnikov(x, ...)
Triangular(x, ...)
Gaussian(x, ...)
HigherOrder(x, ...)
Rectangular(x, ...)
IntBiweight(x)
IntEpanechnikov(x)
IntRectangular(x)
IntTriangular(x)
IntGaussian(x)
SDBiweight(x)
a0(x,h)
a1(x,h)
a2(x,h)
BoundaryBiweight(x, h)
b0(x,h)
b1(x,h)
b2(x,h)
BoundaryEpanechnikov(x, h)
Habbema(xin, x)
}

\arguments{
\item{x}{ A vector of data points where the kernel will be evaluated.}
\item{h}{A scalar.}
\item{xin}{Discrete data inputs especially for the Habbema discrete kernel.}
\item{...}{Further arguments.}
}

\details{
Implements the Biweight, Second Derivative Biweight,  Epanechnikov, Triangular, Guassian, Rectangular, the Boundary adjusted Biweight and Epanechnikov kernels. It also provides the kernel distribution functions for the Biweight, Epanechnikov, Rectangular, Triangular and Guassian kernels. Additionally it implements the discrete kernel Habbema.
}

\value{
The value of the kernel at \eqn{x}
}

\references{
\enumerate{
\item \href{https://ieeexplore.ieee.org/document/4385743}{Bagkavos and Patil, Local Polynomial Fitting in Failure Rate Estimation, IEEE Transactions on Reliability, 57, (2008)},
\item \href{https://link.springer.com/article/10.1007/s10463-010-0277-6}{Bagkavos (2011),  Annals of the Institute of Statistical Mathematics, 63(5), 1019-1046},
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

%\examples{
%
%}
