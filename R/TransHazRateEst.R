"TransHazRateEst"<-
  #Transformed hazard rate estimator (2nd stage) using nonparametric transformation -log(1-\hat F(x))
  function(xin, xout, kfun, ikfun, h1, h2, ci)
  {
    n<-length(xin)
    n1<-1:n
    nn<-length(xout)
    xin.use<-sort(xin)
    GofX<- iHazardRateEst(xin.use, xout, kfun, h1, ci)
    Tsample<- iHazardRateEst(xin.use, xin.use, ikfun, h1, ci)

    GofXdashed<- HazardRateEst(xin.use, xout,  kfun, h1, ci)
    arg1<-(sapply(GofX, "-", Tsample))/h2
    arg2<-kfun(arg1)* ci/(n-n1+1)
    arg3<-colSums(arg2)/h2
    TransEstimate<-  arg3*GofXdashed

    TransEstimate
  }
