VarBandHazEst<-
  #Adaptive variable bandwidth hazard rate estimator
  function(xin, xout,  kfun, h1, h2, ci)
  {
    n<- length(xin)
    nn<-length(xout)
    n1<-1:n
    xin.use<-sort(xin)
    pilot.est<-HazardRateEst(xin.use, xin.use, kfun, h1, ci) #Pilot estimate (W-L)
    sqrt.pilot<-sqrt(pilot.est)
    vhre<-sapply(1:nn, function(i, xin.use, xout, h2, kfun, n, n1,ci, sqrt.pilot)
    {
      sum( (sqrt.pilot *  kfun( ((xout[i]-xin.use)/h2) * sqrt.pilot)	)*ci/(n-n1+1))
    }, xin.use, xout, h2, kfun, n, n1,ci, sqrt.pilot)
    vhre/h2
  }
