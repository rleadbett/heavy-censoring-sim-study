// Independant-informative-prior 
//
// This is a Bayesian model for censored lifetime data 
// which followes a two-paramater Weibull distribution. 
// Informative indapendent marginal priors placed on 
// both the shape and scale paramaters.
//
// The informative priors are constructed from the domain 
// knowledge that:
// 1) The manufacturer claims that the mean lifetime is 5 
//    years with sd of 2 years.
// 2) The beta value of a steel rolling element bearing
//    ranges between 1.1 and 1.5.
//

data {
	int N_obs;
	int N_Rcens;
	int N_Icens;
	real lifetime_obs[N_obs];
	real lifetime_Rcens[N_Rcens];
	real lifetime_Icens_Upper[N_Icens];
	real lifetime_Icens_Lower[N_Icens];
}

parameters {
	real<lower = 0> beta;
	real<lower = 0> lambda;
}

transformed parameters {
  real log_lambda;
  real<lower = 0> eta;
  
  log_lambda = log(lambda);
  eta = exp((-1 / beta) * log_lambda);
}

model{

// Likelihood
// non-censored portion
for(i in 1:N_obs){
	target += weibull_lpdf(lifetime_obs[i]|beta, eta);
}
// censored portion
for(j in 1:N_Rcens){
	target += weibull_lccdf(lifetime_Rcens[j]|beta, eta);
}
// interval portion
for(k in 1:N_Icens){
	target += log_diff_exp(weibull_lcdf(lifetime_Icens_Upper[k]|beta, eta), weibull_lcdf(lifetime_Icens_Lower[k]|beta, eta));
}
  
// Prior models
// These two priors place 95% or the probability
lambda ~ gamma(6.25, 42.59);
beta ~ gamma(156.25, 120.19);

}