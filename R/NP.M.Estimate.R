NP.M.Estimate<-function(xin, cens, xout)
  #purely non parametric M estimate
{
  Muse.np<-SimpsonInt(HRSurv(xout, xin, cens, bw.nrd(xin), Epanechnikov), xout[2]-xout[1])
  Muse.np
}
