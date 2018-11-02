l1<-function(x,p,l){p*(p-1)* l^p * x^(p-2)}
l2<-function(x,p,l){p*(p-1)*(p-2)* l^p * x^(p-3)}
l3<-function(x,p,l){p*(p-1)*(p-2)*(p-3)* l^p * x^(p-4)}
l4<-function(x,p,l){p*(p-1)*(p-2)*(p-3)*(p-4) *l^p * x^(p-5)}
lw<-function(x,p,l){p* l^p * x^(p-1)}
gx<-function(x,p,l){ ((24* (l1(x,p,l))^4 -36*(l1(x,p,l))^2 * (l2(x,p,l))^2 * lw(x,p,l) + 6* (l2(x,p,l))^2 * (lw(x,p,l))^2
                       + 8* l1(x,p,l) *l3(x,p,l) * (lw(x,p,l))^2 -  l4(x,p,l) * (lw(x,p,l))^3)/(24*(lw(x,p,l))^5 ))^2}


lwF<-function(x,p,l){((p* l^p * x^(p-1))^{3/2})/(1-pweibull(x, scale=l, shape=p) )}
