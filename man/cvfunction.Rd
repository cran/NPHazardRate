\name{cvfunction}
\title{Cross Validation for Histogram Hazard Rate Estimator}
\alias{cvfunction}
\description{
Implements the cross validation function for determining the optimal number of bins for the histogram hazard rate estimator of \href{https://www.jstor.org/stable/42003749}{ Patil and Bagkavos (2012)}. It is used as input in   \code{\link{HazardHistogram}}.
}

\usage{
cvfunction(h, xin, xout, cens)
}

\arguments{
\item{h}{Target number of bins.}
\item{xin}{ A vector of data points.  Missing values not allowed.}
\item{xout}{ A vector of grid points at which the histogram will be calculated.}
 \item{cens}{ A vector of 1s and zeros, 1's indicate uncensored observations, 0's correspond to censored obs. }
}

\details{
The least square cross validation criterion, defined in (12), \href{https://www.jstor.org/stable/42003749}{Patil and Bagkavos (2012)} is
\deqn{
CV (h) = \frac{1}{h}\sum_k \big\{(2 f^0_k - f^{0^2}_k)[\bar F_k (\bar F_k +1)]^{-1} - f^{0^2}_k[\bar F_k (\bar F_k +1)^2]^{-1}\big\}.
}
Optimization of the criterion is done through a nonlinear optimization function such as \code{\link{nlminb}} as illustrated also in the example of \code{\link{HazardHistogram}}.
}

\value{
Returns the optimal number of bins.
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

\seealso{
\code{\link{HazardHistogram}
}
}


