"HazardHistogram"<- function(xin, xout, cens, bin)
{
  MinX<-0.1#min(xin)
  MaxX<- 3#max(xin)#max(xin)
  N<-length(xin)
  nn<-length(xout)
  #n<-length(xout)
  #Delta<-(MaxX-MinX)/n
  #cat(Delta)
  n<-(MaxX - MinX)/bin
  PartInt<-c(MinX, MinX+(1:n)*bin) #Partion the interval (x[1], x[n]) into subintervals
  iRows = length(PartInt)
  MatOut<-matrix(ncol=2, nrow=iRows )
  MatOut[,1]<-PartInt
  ordxin<-sort(xin)
  ordcens<-cens[order(xin)]
  indexequals1<-which(ordcens == 1)
  xinnocens<- xin[indexequals1] # X_i for which \delta_i=1
  fizero<-sapply(1:iRows , function(i, xinnocens, PartInt) length(which(xinnocens >= PartInt[i] & xinnocens < PartInt[i+1])),xinnocens, PartInt) #crete the bin counts
  fi<-sapply(1:iRows , function(i, xin, PartInt) length(which(xin >= PartInt[i] & xin < PartInt[i+1])),xin, PartInt)
  MatOut[,2]<-fizero/(bin*cumsum(fi)+1)
  MatOut
}







