"sn.0"<-function(xin, xout, h, kfun) { sapply(1:length(xout), function(i, xin, xout, h, kfun) sum(kfun((xin-xout[i])/h)), xin, xout, h, kfun) }

"sn.1"<-function(xin, xout, h, kfun) { sapply(1:length(xout), function(i, xin, xout, h, kfun)  sum(kfun((xin-xout[i])/h)*(xin-xout[i])), xin, xout, h, kfun) }

"sn.2"<-function(xin, xout, h, kfun) { sapply(1:length(xout), function(i, xin, xout, h, kfun) sum(kfun((xin-xout[i])/h)*((xin-xout[i])^2)), xin, xout, h, kfun) }

"sn.3"<-function(xin, xout, h, kfun){sapply(1:length(xout), function(i, xin, xout, h, kfun)  sum(kfun((xin-xout[i])/h)*((xin-xout[i])^3)), xin, xout, h, kfun) }

"sn.4"<-function(xin, xout, h, kfun) {sapply(1:length(xout), function(i, xin, xout, h, kfun)  sum(kfun((xin-xout[i])/h)*((xin-xout[i])^4)), xin, xout, h, kfun) }

"sn.5"<-function(xin, xout, h, kfun){sapply(1:length(xout), function(i, xin, xout, h, kfun)  sum(kfun((xin-xout[i])/h)*((xin-xout[i])^5)), xin, xout, h, kfun) }

"sn.6"<-function(xin, xout, h, kfun) {sapply(1:length(xout), function(i, xin, xout, h, kfun)  sum(kfun((xin-xout[i])/h)*((xin-xout[i])^6)), xin, xout, h, kfun) }


"tn.0"<-function(xin, xout, h, kfun, Y) { sapply(1:length(xout), function(i, xin, xout, h, kfun,Y)  sum(kfun((xin-xout[i])/h)*Y), xin, xout, h, kfun, Y) }

"tn.1"<-function(xin, xout, h, kfun, Y){sapply(1:length(xout), function(i, xin, xout, h, kfun,Y)  sum(kfun((xin-xout[i])/h)*((xin-xout[i])*Y)),xin, xout, h, kfun, Y) }

"tn.2"<-function(xin, xout, h, kfun, Y){sapply(1:length(xout), function(i, xin, xout, h, kfun,Y) sum(kfun((xin-xout[i])/h)*((xin - xout[i])^2 )*Y),xin, xout, h, kfun, Y) }

"tn.3"<-function(xin, xout, h, kfun, Y){sapply(1:length(xout), function(i, xin, xout, h, kfun,Y) sum(kfun((xin-xout[i])/h)*((xin - xout[i])^3 )*Y),xin, xout, h, kfun, Y) }





"LocLinEst"<-function(BinCenters, xout, h, kfun, ci)
{
	B1<- sn.2(BinCenters, xout, h, kfun)*tn.0(BinCenters, xout, h, kfun,ci) - sn.1(BinCenters, xout, h, kfun)*tn.1(BinCenters, xout, h, kfun,ci)
	B2<- sn.2(BinCenters, xout, h, kfun)*sn.0(BinCenters, xout, h, kfun) - sn.1(BinCenters, xout, h, kfun)*sn.1(BinCenters, xout, h, kfun)
	llest1<- B1/B2
	llest1
}



