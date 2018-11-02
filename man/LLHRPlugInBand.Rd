\name{LLHRPlugInBand}
\title{Simple Plug in badnwidth selector}
\alias{LLHRPlugInBand}
\description{
Provides  the asymptotic MISE optimal plug-in bandwidth for the local linear hazard rate estimator \code{\link{LocLinEst}}, defined in (4),  \href{https://link.springer.com/article/10.1007/s10463-010-0277-6}{Bagkavos (2011)}. This is the binned data version of the \code{\link{PlugInBand}} AMISE optimal bandwidth rule.
}

\usage{
LLHRPlugInBand(BinCenters, h, kfun, Delta, xin, xout, IntKfun, ci, cens)
}

\arguments{
\item{BinCenters}{ A vector of data points, the centers of the bins resulting from the discretization of the data.}
 \item{h}{ Bandwidth for the estimate of the distribution function. }
 \item{kfun}{ A kernel function. }
 \item{Delta}{A scalar. The length of the bins. }
 \item{xin}{ A vector of data points}
 \item{xout}{ The point at which the estimates should be calculated.}
\item{IntKfun}{ The integrated kernel function.}
\item{ci}{ Crude hazard rate estimates.}
\item{cens}{ Censoring Indicators.}
}

\details{
The bandwidth selector requires binned data, i.e. data in the form \eqn{(x_i, y_i)} where \eqn{x_i} are the bin centers and \eqn{y_i} are empirircal hazard rate estimates at each \eqn{x_i}. This is achieved  via  the \code{\link{DiscretizeData}} function. As it can be seen from (4) in \href{https://link.springer.com/article/10.1007/s10463-010-0277-6}{Bagkavos (2011)}, the bandwidth selector also requires an estimate of the functional
\deqn{
\int \left \{ \lambda^{(2)}(x) \right \}^2 \,dx
}
which is readily implemented in \code{\link{PlugInBand}}. It also requires an estimate of the constant
\deqn{
\int \frac{\lambda(x)}{1-F(x)} \,dx
}
For this reason additionally the  plug in bandwidth rule is also used, as it is implemented in the  \code{\link{bw.nrd}} distribution function default bandwidth rule of \href{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9469.2005.00472.x}{Swanepoel and Van Graan (2005)}. The constants \eqn{R(K)} and \eqn{\mu_2^2(K)} are deterministic and specific to the kernel used in the implementation hence can be claculated precisely.
}

\value{
A scalar with the value of the suggested bandwidth.
}

\references{
\href{https://link.springer.com/article/10.1007/s10463-010-0277-6}{Bagkavos (2011), Annals of the Institute of Statistical Mathematics, 63(5), 1019-1046.}
}

%\author{
%Dimitrios Bagkavos
%
%R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com>
%}



\seealso{ \code{\link{PlugInBand} } }

\examples{
x<-seq(0, 5,length=100) #design points where the estimate will be calculated
SampleSize<-100 #amount of data to be generated
ti<- rweibull(SampleSize, .6, 1) # draw a random sample
ui<-rexp(SampleSize, .2)         # censoring sample
cat("\n AMOUNT OF CENSORING: ", length(which(ti>ui))/length(ti)*100, "\n")
x1<-pmin(ti,ui)                  # observed data
cen<-rep.int(1, SampleSize)      # initialize censoring indicators
cen[which(ti>ui)]<-0             # 0's correspond to censored indicators

a.use<-DiscretizeData(ti, x)     # discretize the data
BinCenters<-a.use$BinCenters     # get the data centers
ci<-a.use$ci                     # get empircal hazard rate estimates
Delta=a.use$Delta                # Binning range
h2<-bw.nrd(ti)                   # Bandwidth to use in constant est. of the plug in rule
h.use<-h2                        # the first element is the band to use

huse1<- LLHRPlugInBand(BinCenters,h.use,Epanechnikov,Delta,ti,x,IntEpanechnikov,ci,cen)
huse1


}
