#Load libraries.  I use 'cmdstanr' for Stan. One may use 'rstan' instead.
library(cmdstanr)
library(posterior)
library(bayesplot)
library(tidyverse)
library(survival)

N = 10000;
r = .3;
alpha = 5;
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

#create data frame for Kaplan-Meier
for_km <- data.frame( time = t,
                      status = (t<C) ) ;

#compile Stan model.  The 'stanc_options =' can be omitted.
file <- file.path("~/R/ExpGamma_1f.stan")
mod <- cmdstan_model(file, stanc_options = list("O1"), quiet=TRUE)

#if your "time" or "t" is recorded in days, you should normalize it to month by dividing by 28
#if your "time" or "t" is recorded in weeks, you should normalize it to month by dividing by 4
#if your "time" or "t" is recorded in years, you should normalize it to month by multiplying it by 12
data_list <- list(
  Nuc = length(which(for_km$status==1)), 
  yuc = unlist(for_km$t[which(for_km$status==1)]),
  Nc = length(which(for_km$status==0)),
  yc = unlist(for_km$t[which(for_km$status==0)]),
  T2=20
)



fit_eg1f <- mod$sample(
  data = data_list,
  seed = 123,
  init = 1,
  chains = 4,
  parallel_chains = 4,
  refresh = 500,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = .8
)

fit_eg1f$summary()

posterior1 <- fit_eg1f$draws(format="df")

#check parameters convergence and diagnostics
mcmc_trace(posterior1, pars = c("r","alpha","p"))

mcmc_intervals(posterior1, pars =c("r","alpha","p"))

mcmc_areas(posterior1, prob = 0.90, prob_outer = 1, pars =c("r") )
mcmc_areas(posterior1, prob = 0.90, prob_outer = 1, pars =c("alpha") )
mcmc_areas(posterior1, prob = 0.90, prob_outer = 1, pars =c("p") )

mcmc_dens_overlay(posterior1, pars=c("r","alpha","p")) 

mcmc_pairs(posterior1,pars=c("r","alpha") )

mcmc_hex(posterior1, pars=c("r","alpha"))

mcmc_scatter(posterior1,pars=c("r","alpha") )

####################################################################

#summarize the posterior survival curves
model_draws_Surv <-
  posterior1 %>%
  select(.chain, .iteration, .draw, starts_with("S[")) %>%
  gather(key=key, value=Surv,starts_with("S[")) %>%
  mutate(timing = as.numeric(str_sub(key,3,-2))) %>%
  group_by(timing) %>%
  summarize(surv_mean = mean(Surv),
            surv_median = median(Surv),
            surv_p5 = quantile(Surv, probs = .05),
            surv_p95 = quantile(Surv, probs = .95)) %>%
  ungroup() %>%
  mutate(timing = (timing-1)/4)


fit.km<-survfit(Surv(time,  status) ~ 1, data = for_km, conf.int=F)  

km_df <- data.frame(fit.km$time,fit.km$surv)

names(km_df)

km_df <- km_df %>%
  rename(time = fit.km.time, surv = fit.km.surv) 


#Plot them side by side.  Output graph is omitted.
ggplot() + 
  geom_line(data = model_draws_Surv, aes( x=timing, y=surv_median, color = "posterior prediction")) +
  geom_ribbon(data = model_draws_Surv, aes(x=timing, ymin=surv_p5, ymax=surv_p95), alpha=.2) +
  geom_line(data = km_df, aes( x=time, y=surv, color = "K-M")) + 
  scale_color_manual(name = "", 
                     values = c("posterior prediction"="blue", "K-M" = "red"))

