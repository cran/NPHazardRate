\name{l1-l4, lw, lwF, gx}
\alias{l1}
\alias{l2}
\alias{l3}
\alias{l4}
\alias{lw}
\alias{lwF}
\alias{gx}
\title{Weibull hazard rate functionals}

\description{
Privides the various hazard rate function derivatives and related functionals with reference to the Weibull function
}

\usage{
l1(x,p,l)
l2(x,p,l)
l3(x,p,l)
l4(x,p,l)
lw(x,p,l)
lwF(x,p,l)
gx(x,p,l)
}

\arguments{
 \item{x}{A vector of points at which the hazard rate function will be estimated. }
\item{p}{MLE estimate of the shape parameter}
 \item{l}{MLE estimate of the scale parameter}
}

\details{
Implements the necessary functions for calculating the squared bias term of the variable bandwidth estimate.


}

\value{
A vector with the values of the function at the designated points x.
}

\references{
\href{http://dx.doi.org/10.1080/03610920802364088}{Bagkavos and Patil (2009), Variable Bandwidths for Nonparametric Hazard Rate Estimation, Communications in Statistics - Theory and Methods, 38:7, 1055-1078}
}

%\author{
%Dimitrios Bagkavos
%
%R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com>
%}

%\note{
%%  ~~further notes~~
%}

%\seealso{ \code{\link{HazardRateEst}, \link{TransHazRateEst}, \link{PlugInBand} } }

%\examples{
%x<-seq(0, 5,length=100) #design points where the estimate will be calculated
%plot(x, HazardRate(x, "weibull", .6, 1),  type="l",
%     xlab = "x", ylab="Hazard rate") #plot true hazard rate function
%SampleSize <- 100
%mat<-matrix(nrow=SampleSize, ncol=20)
%for(i in 1:20)
%{
%  ti<- rweibull(SampleSize, .6, 1)#draw a random sample from the actual distribution
%  ui<-rexp(SampleSize, .05)       #draw a random sample from the censoring distribution
%  cat("\n AMOUNT OF CENSORING: ", length(which(ti>ui))/length(ti)*100, "\n")
%  x1<-pmin(ti,ui)                 #this is the observed sample
%  cen<-rep.int(1, SampleSize)     #censoring indicators
%  cen[which(ti>ui)]<-0            #censored values correspond to zero
%
%  h2<-BandsIkde(ti)               #temporary default rule for var. bandwidth
%  h2<-h2[2]                       #extract the temp bandwidth
%  huse1<-  PlugInBand(x1, x,   cen, Biweight)
%  mat[,i]<- VarBandHazEst(x1, x, Epanechnikov, huse1,h2, cen) #Var. bandwidth est.
%}
%lines(x, rowMeans(mat) , lty=2)   #draw the average  vb estimate
%}
