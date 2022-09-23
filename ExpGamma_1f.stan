// Exponential-Gamma with 1-month inflation


functions{

real gamma_Exp_ldensity(real t, real alpha, real r) {

  return log(r)-log(alpha) - (r+1)*log1p( t / alpha );

}

real gamma_Exp_lsurv(real t , real alpha, real r) {

  return   -r * log1p( t / (alpha) ) ;
    
  }

}

data{
  
  int Nuc;
  real yuc[Nuc];
  int Nc;
  real yc[Nc];
  int T2;

}

parameters{
  
  real<lower=0> inv_alpha;
  real<lower=0> r_over_a;
  real<lower=0,upper=1> p;

}
  
transformed parameters{
  real alpha = inv(inv_alpha);
  real r = r_over_a * alpha;
}
  
  
model{
  
  p ~ normal(.5,.25);
  
  for (i in 1:Nuc) {
    
    if (yuc[i]<=1) 
      { target += log_sum_exp(log(p), log(1-p)+gamma_Exp_ldensity(yuc[i],alpha,r)); } 
    else 
      { target += log(1-p)+gamma_Exp_ldensity(yuc[i],alpha,r); }
  }
  
  for (i in 1:Nc){
    target += log(1-p)+gamma_Exp_lsurv(yc[i],alpha,r);
  }
  
}

 
generated quantities{
   
   vector[T2*4+1] S = rep_vector(0,T2*4+1);

    {
      vector[Nuc + Nc] time;
  
      for (i in 1:(Nuc + Nc)) {
          real D = bernoulli_rng(p);
          if (D==1) {
            time[i] = 1; 
            }
          else {
            real f = uniform_rng(0,1);
            time[i] = -alpha*( -1 + pow( 1-f,inv(r) ) )* pow( 1-f, -inv(r) ) ;
            }; 
       }
  
      for (t in 0:T2*4){
        for (i in 1 :(Nuc+Nc) ) {
          S[t+1] += ( time[i]> t/4.0 ) ;
        }
      }
   
    }
    
    S=S/(Nuc+Nc);
}

