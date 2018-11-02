"PlugInBand"<-function( xin, xout,   cens, kfun )
{
  n<-length(xin)
  RofK <- 5/7 # for biweight - for epanech use 3/5
  MuofK <- 1/7 # for biweight - for epanech use 1/5

  hsd<- bw.nrd(xin) # RefBand.b.Sec.Der(xin, cens, xout) #
  theta<-SDHazardRateEst(xin, xout, hsd, cens)
  theta22<-SimpsonInt(theta^2, xout[2]-xout[1])

  M.funct<-HazardRateEst(xin, xout, kfun, bw.nrd(xin), cens)/KMest(xin, cens, xout)

  Muse<-  SimpsonInt(M.funct, xout[2]-xout[1])
  #cat("Muse PluginBand= ", Muse, " Sec. der =", theta22, "\n")
  huse1<-( RofK*Muse /(n * MuofK^2 * theta22))^0.2
  huse1
}
