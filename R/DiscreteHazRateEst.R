lambdahat<-function(xin, cens, xout)
  #non smooth estimate of the discrete hazard function
{
  n<-length(xout)
  mat<-sapply(1:n, function(i, xout, xin, cens)
  {
    xoutuse<-xout[i]
    whichtk<-which(xoutuse == xin) ## index of those T_j = t_k
    xinuse<-xin[whichtk]
    deltause<-cens[whichtk]
    num<- length(which(deltause !=	0))	#equivalent to \sum_{i=1}^n I(T_i = t_k)\Delta_i
    den<-length(xinuse)
    c(num, den)
  }, xout, xin, cens)
  mat<-matrix(mat, ncol=2, byrow=TRUE)

  numerator<-mat[,1]
  denuse<-mat[,2]
  m <- length(denuse)
  denominator<-sapply(1:m,  function(j, denuse, m) sum(denuse[j:m]), denuse=denuse, m=m) #denominator: \sum_{j \geq k} \sum_{k=1}^n I(T_j=t_k)
  numerator/(denominator+1)
}


nsf<-function(xin, cens, xout)
{
  n<-length(xout)
  mat<-sapply(1:n, function(i, xout, xin, cens)
  {
    xoutuse<-xout[i]
    whichtk<-which(xoutuse == xin) ## index of those T_j = t_k
    length(xin[whichtk]	)
  }, xout, xin, cens)
  #mat<-matrix(mat, ncol=2, byrow=T)
  m <- length(mat)
  denominator<-sapply(1:n,  function(j, mat, n) sum(mat[j:n]), mat=mat, n=n)
  denominator
}



"Habbema"<-function(xin, x)
{
  n<-length(x)
  nsuse<-nsf(xin, xin, x)
  #cat(nsuse)
  huse<- .8 * 0.08 * (nsuse/length(xin))^1.5
  #cat(huse)
  mat<-sapply(1:n, function(i, x, xin, huse)
  {
    xoutuse<-x[i]
    hh<-huse[i]
    (1-hh)^((xoutuse-xin)^2)
  }, x, xin, huse)
  mat
}

"TutzPritscher"<-function(xin, cens, xout)
{
  mle<-lambdahat(xin, cens, xout)
  nsuse<-nsf(xin, cens, xin)
  #nsuse<-matrix(nsuse, nrow=length(xin), ncol=length(xout))
  huse<-bw.nrd(xin)
  kernel<-Habbema(xin, xout)
  weightsuse<- nsuse * kernel
  weights<-weightsuse/colSums(weightsuse)
  colSums(mle * weights, na.rm=TRUE)
}

Tm<-function(tk, xout, distribution, par1, par2)
  # This is used to calculate max(t_k; 1-\sum_{l=0}^k \eta_l > \epsilon), \epseilon>0
{
  SumDistrTk<-PdfSwitch(xout, distribution, par1, par2)
  n<-length(SumDistrTk)
  SumTK<-sapply(1:n,  function(i, SumDistrTk, n) sum(SumDistrTk[i:n]), SumDistrTk=SumDistrTk, n=n) #\sum_{l=0}^n \eta_l
  Tmh<- 1-SumTK  #1-\sum_{l=0}^n \eta_l

  maxTmh<-max(which(Tmh == max(Tmh)) )#Tmh[Tmh>= 0 & Tmh<= min(Tmh[Tmh>=0])] # max{1-\sum_{l=0}^n \eta_l >0}
  tk[maxTmh] # find the gridpoint tk which corresponds to Tmh, return the result.
}



CparamCalculation<-function(gamparam, VehHazard)
  # returns the C smoothing parameter calculated as
  # C= gamma/max_{k \geq 0} ( \lambda(t_{k-1}) + \lambda(t_k) + \lambda(t_{k+1}) )
{
  VehHazard<-c(0, VehHazard) # \lambda(t_{-1}) is set to 0 here
  n<-length(VehHazard)-2
  Lsums<-sapply(1:n, function(i, VehHazard) sum(VehHazard[seq(i,i+2, by=1)]), VehHazard) #calculate \lambda(t_{k-1}) + \lambda(t_k) + \lambda(t_{k+1})

  gamparam/max(Lsums) #return C
}

##################################################################################################################
##############alternative DetermineSCprod<-function(NonParEst, VehHazard, Hvector, n, gammapar) function #########

# DetermineSCprod<-function(NonParEst, VehHazard, Hvector, n, gammapar)
#   #this finds SC = \gamma((n+1) \hat B_1)^{-1} \hat V_1
#   # n = number of obs, gammapar = sum of vehicle haz at xout (computed elsewhere)
# {
#   VehHazard<-c(0, VehHazard)
#   NonParEst<-c(0, NonParEst)
#   n<-length(VehHazard)-2
#   B1hatElements<-sapply(1:n, function(i, VehHazard, NonParEst)
#     ((NonParEst[i] + NonParEst[i+2])* VehHazard[i+1] - (VehHazard[i] + VehHazard[i+2])* NonParEst[i+1] )^2
#     , VehHazard, NonParEst) #calculate (\hat \lambda(t_{k-1}) + \hat \lambda(t_{k+1})) \lambda(t_k) - ...)^2
#
#   B1hat<-sum(B1hatElements)
#
#   V1hatElements<-sapply(1:n, function(i, VehHazard, NonParEst)
#     NonParEst[i+1]*(1-NonParEst[i+1]) *(VehHazard[i] + VehHazard[i+2])/Hvector[i]
#     , VehHazard, NonParEst)
#
#   V1hat<- sum(V1hatElements	)
#   gammapar *  V1hat / ((n+1)*B1hat)	#estimated optimal product for S and C
#
# }

#######################################################################################################################

# DetermineSCprod<-function(NonParEst, VehHazard, Hvector, nn, gammapar)
#   #this finds SC = \gamma((n+1) \hat B_1)^{-1} \hat V_1
#   # n = number of obs, gammapar = sum of vehicle haz at xout (computed elsewhere)
# {
#   VehHazard<-c(0, VehHazard)
#   NonParEst<-c(0, NonParEst)
#   n<-length(VehHazard)-2
#   B1hatElements<-sapply(1:n, function(i, VehHazard, NonParEst)
#     ( sum(NonParEst[seq(i,i+2, by=2)])* VehHazard[(i+i+2)/2] - sum(VehHazard[seq(i,i+2, by=2)])* NonParEst[(i+i+2)/2])^2
#     , VehHazard, NonParEst) #calculate (\hat \lambda(t_{k-1}) + \hat \lambda(t_{k+1})) \lambda(t_k) - ...)^2
#
#   B1hat<-sum(B1hatElements)
#   V1hatElements<-sapply(1:n, function(i, VehHazard, NonParEst, Hvector){
#     NonParEst[(i+i+2)/2]*(1-NonParEst[(i+i+2)/2]) *sum(VehHazard[seq(i,i+2, by=2)])/Hvector[i]}
#     , VehHazard, NonParEst, Hvector)
#   V1hat<- sum(V1hatElements)
#   (gammapar *  V1hat) / ((nn+1)*B1hat)	#estimated optimal product for S and C
# }

power.matrix <- function(M, n){
  #Raise matrix M to the n'th power
  powers <- base(n, 2)
  result <- diag( nrow(M) )
  for(i in 1:length(powers)) {
    if(powers[i])
      result <- result %*% M
    M <- M %*% M
  }
  result
}

base <- function(m, b){
  #Express m as powers of b; m may be vector
  maxi <- floor(log(max(m))/log(b))
  result <- matrix( 0, length(m), maxi+1 )
  for(i in 1:(maxi+1)) {
    result[,i] <- m %% b
    m <- m %/% b
  }
  result
}


SmoothedEstimate<-function(NonParEst, VehHazard, gammapar, SCproduct, Cpar)
{
  n<-length(VehHazard)-2 # maybe length(VehHazard)-2
  Sparam<- floor(SCproduct/Cpar) # S = S*C/C - this may need to change
  b0<- gammapar - Cpar*VehHazard[2] # b0 = \gamma - C \lambda_1 (counting starts at zero in the manuscript)
  bkplus1<- gammapar - Cpar * head(tail(VehHazard,2),1) # thisis the last value of the Bk vector (not specified in manuscript)
  ak <- Cpar * VehHazard
  bk <- sapply(1:n, function(i, VehHazard, Cpar, gammapar)
  { gammapar - Cpar * (VehHazard[i] + VehHazard[i+2])}, VehHazard, Cpar, gammapar)
  bk<- c(b0,bk, bkplus1)
  Wmatrix<-diag(bk)

  Wmatrix[col(Wmatrix) - row(Wmatrix)==1]<-head(ak, length(ak)-1)
  Wmatrix[ row(Wmatrix) - col(Wmatrix) ==1]<- tail(ak, length(ak)-1)
  LambdaMat<- (1/gammapar) * ( diag(1, nrow = length(VehHazard), ncol = length(VehHazard)) %*% Wmatrix)

  #################### test #################################
  #cat("\n Vehicle Hazard = ", VehHazard, "\n LAMBDA MAT = ", NonParEst %*% power.matrix(LambdaMat,200), "\ndiag= ", diag(power.matrix(LambdaMat,200)),  "\n")
  #aaaa<- NonParEst %*% power.matrix(LambdaMat,5)
  #plot(1:length(aaaa), VehHazard, type="l")
  #lines(1:length(aaaa), aaaa, lty=4)
  ###########################################################

  if(Sparam == 0) LambdaMatrixSpower <- diag(1, nrow=length(VehHazard), ncol= length(VehHazard))
  else LambdaMatrixSpower <- power.matrix(LambdaMat, Sparam)

  SmoothEst<- matrix(NonParEst, nrow=1, ncol=length(VehHazard)) %*% LambdaMatrixSpower
  SmoothEst

}

# VehicleHazardRate<-function(xout, vehicledistr, param1, param2)
# {
#   switch(vehicledistr, binomial=dbinom(xout, param1, param2)/(1-pbinom(xout-1, param1, param2))
#          , geometric = dgeom(xout, param1)/(1-pgeom(xout-1, param1))
#          , poisson = dpois(xout, param1)/(1-ppois(xout-1, param1))
#          , negativebinomial = dnbinom(xout, param1, param2)/(1-pnbinom(xout-1, param1, param2))
#
#   )
# }






SemiparamEst<-function(xin, cens, xout, Xdistr, Udistr, vehicledistr, Xpar1=1, Xpar2=0.5, Upar1=1, Upar2=0.5, vdparam1=1, vdparam2=0.5)
{


  TkCutOff<-min(Tm(xout, xout, Xdistr, Xpar1, Xpar2) , Tm(xout, xout, Udistr, Upar1, Upar2), Tm(xout, xout, vehicledistr, vdparam1, vdparam2))
  # cat(xout, " ", TkCutOff)
  xout<-head(xout, TkCutOff)

  result<-matrix(nrow=length(xout), ncol=2)

  Qvector <- 1-cumsum(PdfSwitch(xout-1, Xdistr, Xpar1, Xpar2))
  Pvector <- 1-cumsum(PdfSwitch(xout-1, Udistr, Upar1, Upar2))

  Hvector <-  Pvector * Qvector
  NonparamEst<-lambdahat(xin, cens, xout)

  VehicleHazard<-HazardRate(xout, vehicledistr, vdparam1, vdparam2)
  GammaParam <- 2 #sum(VehicleHazard)
  Cpar <-  0.5# CparamCalculation(GammaParam , VehicleHazard) #upper bound for C 0.5#
  SCproduct<- 32* Cpar # DetermineSCprod(NonparamEst, VehicleHazard, Hvector , length(xin), GammaParam)

  SmoothEst<- SmoothedEstimate(NonparamEst, VehicleHazard, GammaParam , SCproduct, Cpar )
  result[,1]<-xout
  result[,2]<-SmoothEst
  result
}


