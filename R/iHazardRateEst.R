"iHazardRateEst"<- function(xin, xout, ikfun, h, ci)
{
  n<- length(xin)
  n1<-1:n
  xin.use<-sort(xin)
  arg1<-(sapply(xout, "-", xin.use))/h
  arg2<-ikfun(arg1)* ci/(n-n1+1)
  arg3<-colSums(arg2)
  arg3
}
