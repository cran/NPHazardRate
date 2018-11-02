"cvfunction"<- function(h, xin, xout, cens)
{
  MinX<- 0.01#min(xin)
  MaxX<- 3#max(xin)
  N<-length(xin)
  n<-(MaxX - MinX)/h
  PartInt<-c(MinX, MinX+(1:n)*h) #Partion the interval (x[1], x[n]) into subintervals
  ordxin<-sort(xin)
  ordcens<-cens[order(xin)]
  indexequals1<-which(ordcens == 1)
  xinnocens<- xin[indexequals1] # X_i for which \delta_i=1
  fizero<-sapply(1:n, function(i, xinnocens, PartInt) length(which(xinnocens >= PartInt[i] & xinnocens < PartInt[i+1])),xinnocens, PartInt) #crete the bin counts
  fi<-sapply(1:n, function(i, xin, PartInt) length(which(xin >= PartInt[i] & xin < PartInt[i+1])),xin, PartInt)
  Fi<-cumsum(fi)
  sum( (2*fizero - fizero^2)/(Fi*(Fi+1)) - fizero^2 /(Fi*((Fi+1)^2))) /h
}