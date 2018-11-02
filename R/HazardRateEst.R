

"HazardRateEst"<- function(xin, xout, kfun, h, ci)
  {
    n<- length(xin)
    n1<-1:n
    xin.use<-sort(xin)
    arg1<-(sapply(xout, "-", xin.use))/h
    arg2<-kfun(arg1, h)* ci/(n-n1+1)
    arg3<-colSums(arg2)/h
    arg3
  }


