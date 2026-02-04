data {
  int<lower=0> N;
  vector[N] S;
  vector[N] d;
  vector[N] n;
  vector[N] th;
  array[N] int idx;  
  int n_idx;
  real sigma_prior;
  real alpha_prior;
  real beta_prior;
  real gamma_prior;
}
// 
parameters {
	// vector<lower=0>[n_idx] alpha; //OK
	// vector<lower=0>[N] alpha; //OK
  // vector[N] beta; //OK
  // real<lower=0> sigma;    //OK
  // 
  real kappa_0;
  real kappa_1;
  real kappa_2;
  real<lower=0> eta_0;
  real eta_1;
  real eta_2;
}
// 
transformed parameters{
	vector<lower=0>[N] alpha = exp(kappa_0) + exp(kappa_1) * log(n)+ exp(kappa_2) * th;
	vector[N] beta = eta_0 + eta_1 * log(n) - exp(eta_2) * th^2; //OK
  vector[N] mu = log(alpha) + beta .* log(d);
}
// 
model {
  // Priors
  kappa_0 ~ normal(-6,0.5);
  kappa_1 ~ normal(-8,0.5);
  kappa_2 ~ normal(-4,0.1);
  eta_0 ~ normal(2,0.5); 
  eta_1 ~ normal(-1,0.5);
  eta_2 ~ normal(-4,1);
  // 
  // Likelihood
  S ~ normal(exp(mu), 0.1);
}

// generated quantities {
//   array[N] real S_pred;         ////OK
// 
//   for (i in 1:N) {
//     S_pred[i] = normal_rng(alpha[idx[i]] * pow(d[i], beta[idx[i]]), sigma); //OK
//   }
// }
