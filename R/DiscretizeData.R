"DiscretizeData"<-function(xin, xout)
{
  MinX<- min(xin)
  MaxX<- max(xin)-0.2*max(xin)
  N<-length(xin)
  n<-length(xout)
  Delta<-(MaxX-MinX)/n
  PartInt<-c(MinX, MinX+(1:n)*Delta) #Partion the interval (x[1], x[n]) into subintervals
  FirstBinCenter<-(PartInt[1]+PartInt[2])/2
  BinCenters<- c(FirstBinCenter, FirstBinCenter+(1:(n-1))*Delta)
  fi<-sapply(1:n, function(i, xin, PartInt) length(which(xin > PartInt[i] & xin <= PartInt[i+1])), xin, PartInt) #crete the bin counts
  citmp<-fi/(N-cumsum(fi)+1)
  ci<-citmp/Delta
  list(BinCenters=BinCenters, ci=ci, Delta=Delta)
}
