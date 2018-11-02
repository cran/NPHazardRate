"Biweight"<-function(x, ...) { ifelse(abs(x) < 1, (15/16) * (  1 - x^2)^2, ifelse(abs(x) > 1, 0, 0)) }
"Epanechnikov"<-    function(x, ... ) { ifelse(abs(x)<1,  0.75 * (1-x^2) , ifelse(abs(x)>1, 0,0))  }
"HigherOrder"<-function(x, ...){ifelse(x < -4, 0, ifelse(x> 4, 0, dnorm(x, 0, 1) * (3-x^2)/2))}#, ifelse(abs(x) > 4, 0, 0)) }
"Rectangular"<- function(x, ...){ ifelse(abs(x)<1, .5, ifelse(abs(x)>1,0,0))}
"Triangular"<-function(x, ...){ ifelse(abs(x)<1, 1-abs(x), ifelse(abs(x)>1,0,0)) }
"Gaussian"<- function(x, ...){ dnorm(x) }

"SDBiweight"<-function(x){ifelse(abs(x) < 1, -15/4 + (45*x^2)/4, ifelse(abs(x) > 1, 0, 0))}

"IntBiweight"<- function(x){ ifelse(x< -1, 0,ifelse(x> 1, 1, .5 + ((15*x)+ 3*(x^5))/16 - (5*x^3)/8)) }
"IntEpanechnikov"<-function(x) { ifelse(x< -2.236068, 0, ifelse(x> 2.236068, 1, .5- (2.236068*x*(x^2-15))/100)) }

"IntRectangular"<-  function(x){ ifelse(x< -1, 0, ifelse(x>1, 1, (x+1)/2))}
"IntTriangular"<- function(x){  ifelse(x< -1, 0, ifelse(x>1, 1, x+ ((x^2)*sign(x) - 1)/2)) }
"IntGaussian"<-  function(x){ pnorm(x) }



"a0"<-function(x,h) { 3/(4*sqrt(5))  * (x/h - sqrt(5) - (1/15) * ((x/h)^3 - (-sqrt(5))^3)) }
"a1"<-function(x,h) { 3/(4*sqrt(5))  * (   ( (x/h)^2 - 5)/2 - ( (x/h)^4 - (sqrt(5))^4)/20 ) }
"a2"<-function(x,h) {  -3/500*sqrt(5)* (x/h)^5+1/2+(1/20)*sqrt(5)*(x/h)^3 }
"BoundaryBiweight"<-function(x, h){ifelse(abs(x)<2.236068, (a2(x,h)- a1(x,h)*x)*Biweight(x)/(a1(x,h)^2 - a0(x,h)*a2(x,h)), ifelse(abs(x)>2.236068, 0,0))}

"b2"<-function(x,h) { 0.75* (((x/h)^3)/3 - ((x/h)^5)/5) + 1/10 }
"b1"<-function(x,h) { 0.75* (((x/h)^2)/2 - ((x/h)^4)/4 - 1/4) }
"b0"<-function(x,h) { 0.75* ((x/h) - ((x/h)^3)/3 + 2/3) }
"BoundaryEpanechnikov"<- function(x, h){ ifelse(abs(x)<1, (a2(x,h) - a1(x,h)*x) * Epanechnikov(x)/(a1(x,h)^2 - a0(x,h)*a2(x,h)), ifelse(abs(x)>1, 0,0))}
