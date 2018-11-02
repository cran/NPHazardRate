\name{PlugInBand}
\title{Simple Plug in badnwidth selector}
\alias{PlugInBand}
\description{
Provides  the asymptotic MISE optimal plug-in bandwidth for the hazard rate estimator \code{\link{HazardRateEst}}, see \href{https://onlinelibrary.wiley.com/doi/full/10.1111/anzs.12224}{Hua, Patil and Bagkavos (2018).}  The bandwidth is also suitable for use   as a pilot bandwidth in  \code{\link{TransHazRateEst} } and \code{\link{VarBandHazEst}}.}


\usage{PlugInBand(xin, xout,   cens, kfun ) }

\arguments{
 \item{xin}{ A vector of data points}
 \item{xout}{ The point at which the estimates should be calculated.}
\item{cens}{ Censoring Indicators.}
\item{kfun}{ A kernel function. }
}

\details{
The asymptotic MISE optimal plug-in bandwidth selector for \code{\link{HazardRateEst}} is defined by
\deqn{h_{ opt} = \left[\frac{R(K)}{nR(\lambda_T'')\mu_{2,K}^2}\int \frac{\lambda_T(x)}{1-F(x)}\,dx \right]^{1/5} }
see (9) in \href{https://onlinelibrary.wiley.com/doi/full/10.1111/anzs.12224}{Hua, Patil and Bagkavos (2018)}. The   estimate of \eqn{R(\lambda_T'')} to be used in \eqn{h_{opt}} is
\deqn{
R(\hat \lambda_T'') = \int_0^\xi \left (\hat \lambda_T''(x|\hat b_n^\ast) \right )^2\,dx.
 }
Also,
\deqn{ \int_0^T \frac{\lambda_T(x)}{1-F(x)}\,dx }
is estimated by applying the extended Simpson's numerical integration rule, \code{\link{SimpsonInt}}, on
\deqn{
 \frac{\hat \lambda_T(x|\hat b_n^\ast)  }{1-F(x)}
 }
 where \eqn{1-F(x)} is estimated by \code{\link{KMest}}. The estimation is implemented in the \code{\link{NP.M.Estimate}} function.

Currently \eqn{b_n^\ast} is estimated by \code{\link{bw.nrd}}. However according to (11) in \href{https://onlinelibrary.wiley.com/doi/full/10.1111/anzs.12224}{Hua, Patil and Bagkavos (2018).}, in future versions this package will support
\deqn{
b_n^\ast = \left \{ \frac{5R(K'')}{n \mu_{2,K}^2 R(\lambda_T^{(4)})} \int \frac{\lambda_T(x)}{1-F(x)}\,dx  \right \}^{1/9}.}
where
\deqn{
R(\hat  \lambda_T^{(4)}) =  \frac{(\hat a(\hat a-1)(\hat a-2)(\hat a-3)(\hat a-4))^2}{(2\hat a-9){\hat{b}}^{2\hat a} } (\xi^{2\hat a-9} - {p_\alpha}^{2\hat a-9}), \hat a\neq 9/2
}
and \eqn{\hat M} is already estimated by \code{\link{NP.M.Estimate}} as expalined above (it will be much more stable than using a Weibull reference model).
}

\value{
A scalar with the value of the suggested bandwidth.
}

\references{
\href{https://onlinelibrary.wiley.com/doi/full/10.1111/anzs.12224}{Hua, Patil and Bagkavos, An $L_1$ analysis of a kernel-based hazard rate estimator, Australian and New Zealand J. Statist., (60), 43-64, (2018).}
}

%\author{
%Dimitrios Bagkavos
%
%R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com>
%}



\seealso{ \code{\link{HazardRateEst}, \link{LLHRPlugInBand} } }

\examples{
x<-seq(0, 5,length=100) #design points where the estimate will be calculated
SampleSize<-100 #amount of data to be generated
ti<- rweibull(SampleSize, .6, 1) # draw a random sample
ui<-rexp(SampleSize, .2)         # censoring sample
cat("\n AMOUNT OF CENSORING: ", length(which(ti>ui))/length(ti)*100, "\n")
x1<-pmin(ti,ui)                  # observed data
cen<-rep.int(1, SampleSize)      # initialize censoring indicators
cen[which(ti>ui)]<-0             # 0's correspond to censored indicators

huse1<- PlugInBand(x1, x, cen, Biweight)
huse1


}
