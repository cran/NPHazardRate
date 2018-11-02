\name{DiscretizeData}
\title{Discretize the available data set}
\alias{DiscretizeData}

\description{
Defines equispaced disjoint intervals based on the range of the sample and calculates empirical hazard rate estimates at each interval center
}


\usage{
DiscretizeData(xin, xout)
}

\arguments{
\item{xin}{A vector of input values}
\item{xout}{Grid points where the function will be evaluated}
}

\details{The function defines the subinterval length \eqn{\Delta = (0.8\max(X_i) - \min(X_i))/N} where \eqn{N} is the sample size. Then at each bin (subinterval) center, the empirical hazard rate estimate is calculated by
\deqn{
c_i = \frac{f_i}{\Delta(N-F_i +1) }
}
where \eqn{f_i} is the frequency of observations in the ith bin and \eqn{F_i = \sum_{j\leq i} f_j} is the empirical cummulative distribution estimate.
}

\value{
A vector with the values of the function at the designated points xout or the random numbers drawn.
}

%\references{
%Bagkavos, Dimitrios and Patil, Prakash (2009) Variable Bandwidths for Nonparametric Hazard Rate Estimation, Communications in Statistics - %Theory and Methods, 38, 1055-1078, http://dx.doi.org/10.1080/03610920802364088
%}

%\author{
%Dimitrios Bagkavos
%
%R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com> }

%\note{
%%  ~~further notes~~
%}

%\seealso{ \code{\link{HazardRateEst}, \link{TransHazRateEst} } }

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


}
