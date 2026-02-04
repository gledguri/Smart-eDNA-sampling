data {
  int<lower=0> N;
  vector[N] S;
  vector[N] d;
  array[N] int idx;  
  int n_idx;
  real sigma_prior;
  real alpha_prior;
  real beta_prior;
  real intercept_prior;
  real gamma_prior;
}
// 
parameters {
	vector<upper=0>[n_idx] alpha; //OK
	vector[n_idx] intercept; //OK
  vector[n_idx] beta; //OK
  // vector[n_idx] gamma; //OK2
  real<lower=0> sigma;    //OK
}
// 
transformed parameters{
  vector[N] mu;
  // mu = alpha[idx] .* pow(d, beta[idx]); //OK
  // mu = intercept[idx]+((beta[idx].*d)/(alpha[idx]-d));
  mu = intercept[idx] + (beta[idx] .* d) ./ (alpha[idx] - d);
  // mu = alpha[idx] .* pow(d, beta[idx]) + exp(gamma[idx]) .* d; //OK2
}
// 
model {
  // Priors
  alpha ~ normal(0, alpha_prior); //OK
  beta ~ normal(0, beta_prior); //OK
  intercept ~ normal(0, intercept_prior); //OK
  // gamma ~ normal(-4, gamma_prior);
  sigma ~ normal(0,sigma_prior);
  // sigma ~ exponential(sigma_prior); //OK
  // 
  // Likelihood
  S ~ normal(mu, sigma); //OK
}

generated quantities {
  array[N] real S_pred;         ////OK

  for (i in 1:N) {
    S_pred[i] = normal_rng(alpha[idx[i]] * pow(d[i], beta[idx[i]]), sigma); //OK
  }
}
