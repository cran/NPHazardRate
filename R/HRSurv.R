HRSurv<-function(x, xin, cens, h, kfun)
{
  HazardRateEst(xin, x, kfun, h, cens)/KMest(xin, cens, x) #denom
}

