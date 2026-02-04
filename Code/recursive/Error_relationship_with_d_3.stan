// Multiple Linear Regression Model for epsilon (dev) prediction
data {
  int<lower=0> N;           // number of observations
  vector[N] dev;            // dependent variable (epsilon) - raw scale
  vector[N] th;             // independent variable T
  vector[N] d;              // independent variable D
  vector[N] N_ini;          // independent variable N
}

parameters {
  real alpha;               // intercept
  real beta_th;             // coefficient for th (T)
  real beta_d;              // coefficient for d (D)
  real beta_N_ini;          // coefficient for N_ini (N)
  real<lower=0> sigma;      // error standard deviation
}

model {
  // Priors
  alpha ~ normal(0, 10);
  beta_th ~ normal(0, 10);
  beta_d ~ normal(0, 10);
  beta_N_ini ~ normal(0, 10);
  sigma ~ exponential(0.1);
  
  // Likelihood - modeling log(dev) as response
  log(dev) ~ normal(alpha + beta_th * th + beta_d * d + beta_N_ini * N_ini, sigma);
}

generated quantities {
  vector[N] y_pred_log;     // predicted values on log scale
  
  for (i in 1:N) {
    y_pred_log[i] = alpha + beta_th * th[i] + beta_d * d[i] + beta_N_ini * N_ini[i];
  }
}