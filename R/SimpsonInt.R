"SimpsonInt"<- function(xin, h)
{   #xin = data, h=distance of the points at which function is evaluated
  n<-length(xin)
  nn<-1:n
  int0<-xin[1]
  int2n<-xin[n]
  inteven<-xin[seq(3,n-1,by=2)]
  intodd<-xin[seq(2,n-2,by=2)]
  sumodd<-sum(intodd)
  sumeven<-sum(inteven)
  res<-(h/3)*(int0+4*sumodd+2*sumeven+int2n)
  res
}