library(survival)
library(ggplot2)

#Data simulation

N = 10000;
r = 2;
alpha = 10;
p = 0.35;
#censoring point > 1 (uniform for all subjects)
C=10;

#Inverse-CDF method to simulate from Lomax distribution with uniform censoring
rlomax_cen <- function(N, r, alpha, C) {
  
  f = runif(N);
  
  X = -alpha*( -1 + (1-f)**(1/r))*(1-f)**(-1/r) ; 
  
  X[which(X>=C)]=C;
  
  return ( X );
  
}

D = rbinom(N, 1, p);
t = (D==1)+(D==0)*rlomax_cen(N,r,alpha,C);

#Make a data frame for Kaplan-Meier
for_km <- data.frame( time = t,
                      status = (t<C) ) ;

#from Kaplan-Meier dataframe - can be modified for user's dataframe
t_uncen<-for_km$t[which(for_km$t<C)];
t_cen<-for_km$t[which(for_km$t>=C)];

#Create a logsumexp function to avoid computation error
logsumexp <- function (a,b) {
  
  c = pmax(a,b);
  return (c + log(exp(a-c)+exp(b-c))) ;
  
}

###Lomax + 1-inflated likelihood function
lomax_1f <- function(par,uncen,cen) {
  
  r=par[1];
  alpha=par[2]
  p=par[3];
  
  I0=which(uncen<=1);
  I1=which(uncen>1);
  
  #if uncen > 1  
  LL1=log(1-p)+log(r)-log(alpha)+(r+1)*log(alpha/(alpha+uncen[I1]));
  
  #if uncen <= 1  
  LL2=logsumexp(rep(log(p),length(I0)),
                log(1-p)+log(r)-log(alpha)+(r+1)*log(alpha/(alpha+uncen[I0])));
  
  #if cen  
  LL_cen = r*log(alpha/(alpha+cen)) + log(1-p);
  
  return (-sum(LL1)-sum(LL2)-sum(LL_cen));
}


solution=optim(par=c(1,1,.5),fn=lomax_1f,
               uncen=t_uncen, cen=t_cen,
               method="L-BFGS-B",
               lower=c(1e-5,1e-5,1e-5), upper=c(Inf,Inf,.999))


###Lomax + 1-inflated survival function
lomax1f_surv <- function(par,t) {
  
  r=par[1];
  alpha=par[2];
  p=par[3];
  
  cdf=t;
  surv=t;
  
  for (i in (1:length(t)) ) {
    
    if (t[i]>=1) {
      cdf[i] = 1-(1-p)*(alpha/(alpha+t[i]))^r;
      surv[i] = 1 - cdf[i];
    }
    else {
      cdf[i] = 1-(alpha/(alpha+t[i]))^r;
      surv[i] = 1-cdf[i];
    }
  }
  
  return ( data.frame(time=t,
                      surv=surv) );
  
}

lomax1f_df<-lomax1f_surv(solution$par,c(seq(0, 20, by = .1)))

#Kaplan-Meier curve
fit.km = survfit( Surv(for_km$time, for_km$status) ~ 1, conf.int=F)
km_df <- data.frame(fit.km$time,fit.km$surv)

#Plot them side by side
ggplot() + 
  geom_line(data = lomax1f_df, aes( x=time, y=surv, color = "Model")) +
  geom_line(data = km_df, aes( x=fit.km.time, y=fit.km.surv, color = "K-M")) + 
  scale_color_manual(name = "", 
                     values = c("Model"="blue", "K-M" = "red"))
#see the png file in repository

####Final Remark
#Simulation also reveals parameter recovery is less precise than one hopes. But this does not seem to affect the model fitting to the data.

