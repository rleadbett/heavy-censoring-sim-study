// Independant-uninformative-prior 
//
// This is a Bayesian model for censored lifetime data 
// which followes a two-paramater Weibull distribution. 
// Uninformative indapendent marginal priors are placed on 
// both the shape and scale paramaters.
//
// A gamma(0.001, 0.001) prior is used for both of the 
// paramaters. This prior is dominated by the likelihood as
// long as n >= 1.
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
// right censored portion
for(j in 1:N_Rcens){
	target += weibull_lccdf(lifetime_Rcens[j]|beta, eta);
}
// interval censored portion
for(k in 1:N_Icens){
	target += log_diff_exp(weibull_lcdf(lifetime_Icens_Upper[k]|beta, eta), weibull_lcdf(lifetime_Icens_Lower[k]|beta, eta));
}
  
// Prior models
lambda ~ gamma(0.0001, 0.0001);
beta ~ gamma(0.0001, 0.0001);
}