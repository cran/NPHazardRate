"LLHRPlugInBand"<-function(BinCenters, h, kfun, Delta, xin, xout, IntKfun,   ci, cens)
{
  sn0<-sn.0(BinCenters, BinCenters, h, kfun)
  sn1<-sn.1(BinCenters, BinCenters, h, kfun)
  sn2<-sn.2(BinCenters, BinCenters, h, kfun)
  sn3<-sn.3(BinCenters, BinCenters, h, kfun)
  sn4<-sn.4(BinCenters, BinCenters, h, kfun)
  sn5<-sn.5(BinCenters, BinCenters, h, kfun)
  sn6<-sn.6(BinCenters, BinCenters, h, kfun)
  tn0<-tn.0(BinCenters, BinCenters, h, kfun, ci)
  tn1<-tn.1(BinCenters, BinCenters, h, kfun, ci)
  tn2<-tn.2(BinCenters, BinCenters, h, kfun, ci)
  tn3<-tn.3(BinCenters, BinCenters, h, kfun, ci)
  #COMPUTE THE LOCAL LINEAR DERIVATIVE ESTIMATE:

  den<-  sn0*sn2*sn4*sn6 - sn0*sn2*(sn5^2) - sn0*(sn3^2) *sn6 + 2*sn0*sn3*sn4*sn5-sn0*(sn4^3)-sn1^2*sn4*sn6 +( sn1^2)*(sn5^2 )+ 2*sn1*sn2*sn3*sn6
  - 2*sn1*(sn3^2) *sn5 - 2*sn1*sn4*sn2*sn5 + 2*sn1*sn3*(sn4^2) - (sn2^3) *sn6 + 2*(sn2^2)*sn3*sn5 + (sn2^2)*(sn4^2)-3*sn2*sn4*sn3^2 + sn3^4
  A1<-   sn1*sn3*sn6 -sn1*sn4*sn5 - (sn2^2) *sn6+sn2*sn3*sn5+sn2*sn4^2 -sn4*( sn3^2)

  A2<- sn0*sn3*sn6-sn0*sn4*sn5-sn2*sn1*sn6+sn3*sn2*sn4+sn1*sn3*sn5-(sn3^3)
  A3<- sn0*sn2*sn6-sn0*(sn4^2) -(sn1^2) *sn6+2*sn3*sn1*sn4-sn2*sn3^2
  A4<-  sn0*sn2*sn5 - sn0*sn3*sn4 - (sn1^2) *sn5+sn1*sn2*sn4+sn1*(sn3^2)-sn3*(sn2^2)
  est2<- 2* (A1*tn0 - A2*tn1 + A3*tn2 - A4*tn3)/den

  theta22<- Delta * sum(est2^2)
  RofK <- 3* sqrt(5)/25
  MuofK <- 3* sqrt(5)	/35
  Muse<-  NP.M.Estimate(xin, cens, xout) #1/(1-IntKde(xin, max(xout), h, IntKfun)) -1
  N<-length(xin)
  huse1<-( RofK*Muse /(N * MuofK^2 * theta22))^0.2
  huse1
}


