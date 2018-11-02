RdistSwitch<-function(dist, SampleSize, par1, par2=.5)
{
  switch(dist,
         weibull = rweibull(SampleSize, par1, par2),
         lognorm = rlnorm(SampleSize),
         chisquare = rchisq(SampleSize, par1),
         exponential = rexp(SampleSize, par1),
         binomial = rbinom(SampleSize, par1, par2),
         geometric = rgeom(SampleSize, par1),
         poisson = rpois(SampleSize, par1),
         negativebinomial = rnbinom(SampleSize, par1, par2),
         uniform = runif(SampleSize, par1, par2)
         )
}

PdfSwitch<-function(xout, dist, par1, par2=.5)
{
  switch(dist,
         weibull = dweibull(xout, par1, par2),
         lognorm = dlnorm(xout, par1, par2),
         chisquare = dchisq(xout, par1),
         exponential = dexp(xout, par1),
         binomial = dbinom(xout, par1, par2),
         geometric = dgeom(xout, par1),
         poisson = dpois(xout, par1),
         negativebinomial = dnbinom(xout, par1, par2),
         uniform = dunif(xout, par1, par2)
         	)
}

CdfSwitch<-function(xout, dist, par1, par2=.5)
{
  switch(dist,
         weibull = 1 - pweibull(xout, par1, par2),
         lognorm = 1 - plnorm(xout, par1, par2),
         chisquare = 1 - pchisq(xout, par1),
         exponential = 1 - pexp(xout, par1),
         binomial = 1 -pbinom(xout, par1, par2),
         geometric = 1 - pgeom(xout, par1),
         poisson = 1 - ppois(xout, par1),
         negativebinomial = 1 - pnbinom(xout, par1, par2),
         uniform = 1 - punif(xout, par1, par2)
      )
}

HazardRate<-function(xout, dist, par1, par2)
{
  switch(dist,
         weibull=dweibull(xout, par1, par2)/(1-pweibull(xout, par1, par2)),
         lognorm = dlnorm(xout)/(1-plnorm(xout)),
         chisquare = dchisq(xout, par1)/(1-pchisq(xout, par1)),
         exponential = dexp(xout, par1)/(1-pexp(xout, par1)),
         normixt = ((2/3) * dnorm(xout,3,1)+(1/3) *dnorm(xout,2,.2))/ (1-((2/3) * pnorm(xout,3,1)+(1/3)*pnorm(xout,2,.2))),
         binomial=dbinom(xout, par1, par2)/(1-pbinom(xout-1, par1, par2)),
         geometric = dgeom(xout, par1)/(1-pgeom(xout-1, par1)),
         poisson = dpois(xout, par1)/(1-ppois(xout-1, par1)),
         negativebinomial = dnbinom(xout, par1, par2)/(1-pnbinom(xout-1, par1, par2))
  )
}

