"KMest"<-function(xin, cens, xout)
  #this is the actual KM estimate not the version that estimates the censoring distribution
{
  ordxin <- sort(xin) #order input data=Xi
  ordcens <- cens[order(xin)] #order censoring =Di
  n <- length(xin)
  nn <- length(xout)
  vKuse <- sapply(1:nn, function (i, ordxin, xout)
  {
    k.use<-c(max(which(ordxin <= xout[i])), min(which(xout[i] <= ordxin)))
    k.use.1<-ifelse(is.infinite(k.use[2]), c(k.use[1],k.use[1]) , k.use) #cat("k.use = ", k.use, "\n")
    l.use<- max(k.use.1)
    l.use
  }, ordxin, xout)
  vK <- vKuse -1

  Huse<-sapply(vK, function (i,n,ordcens) {
    i<-ifelse(i==-Inf, 0,i)
    index<-1:i
    expon.use<-ordcens[index]
    ((n-index+1)/(n-index+2))^expon.use
  }
  , n, ordcens)
  result<-sapply(Huse, prod)
  result
}
