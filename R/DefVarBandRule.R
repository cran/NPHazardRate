"DefVarBandRule"<-function(xin, cens)
{
  m4K<-0.0857 # for Epanechnikov kernel
  RofK<- 3/5
  n<-length(xin)
  par<- survreg(Surv(xin, cens) ~1, dist="weibull")
  pl<- 1/par$scale
  pp<- exp(par$coefficients)
  uplim<- max(rweibull(1000, shape=pp, scale=pl))
  Rg<- integrate(gx, lower=0.04, upper = uplim, p=pp, l=pl)$value
  xuse<-seq(0.04, uplim, length=100)
  hatmint<-  lwF(xuse, pp, pl)
  #cat(hatmint, "\n")
  hatM<-  integrate(lwF, lower=0.04, upper = uplim, p=pp, l=pl)$value
  ((RofK * hatM)/(8*n* m4K^2 * Rg))^(1/14)
}
