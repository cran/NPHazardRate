"SDHazardRateEst"<- function(xin, xout, h, ci)
  {
    n<- length(xin)
    n1<-1:n
    xin.use<-sort(xin)
    arg1<-(sapply(xout, "-", xin.use))/h
    arg2<-SDBiweight(arg1)* ci/(n-n1+1)
    arg3<-colSums(arg2)/h^3
    arg3
  }


