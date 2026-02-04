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
  real beta_th;             // coefficient for log(th) (T)
  real beta_d;              // coefficient for d (D)
  real beta_N_ini;          // coefficient for N_ini (N)
  real beta_th_d;           // interaction: log(th) * d
  real beta_th_N;           // interaction: log(th) * N_ini
  real beta_d_N;            // interaction: d * N_ini
  real beta_th_d_N;         // three-way interaction: log(th) * d * N_ini
  real<lower=0> sigma;      // error standard deviation
}

model {
  // Priors
  alpha ~ normal(0, 10);
  beta_th ~ normal(0, 10);
  beta_d ~ normal(0, 10);
  beta_N_ini ~ normal(0, 10);
  beta_th_d ~ normal(0, 5);      // more conservative prior for interactions
  beta_th_N ~ normal(0, 5);
  beta_d_N ~ normal(0, 5);
  beta_th_d_N ~ normal(0, 2);    // even more conservative for 3-way interaction
  sigma ~ exponential(1);
  
  // Likelihood - modeling log(dev) with interactions
  log(dev) ~ normal(alpha + 
                    beta_th * th + 
                    beta_d * log(d) + 
                    beta_N_ini * log(sqrt(N_ini)) +
                    beta_th_d * th .* log(d) +
                    beta_th_N * th .* log(sqrt(N_ini)) +
                    beta_d_N * log(d) .* log(sqrt(N_ini)) +
                    beta_th_d_N * th .* log(d) .* log(sqrt(N_ini)), 
                    sigma);
///////////////////////////////////////////////////
//////////// Working very good FORMULA 1 //////////
///////////////////////////////////////////////////
 // log(dev) ~ normal(alpha + 
 //                    beta_th * th + 
 //                    beta_d * log(d) + 
 //                    beta_N_ini * log(N_ini) +
 //                    beta_th_d * th .* log(d) +
 //                    beta_th_N * th .* log(N_ini) +
 //                    beta_d_N * log(d) .* log(N_ini) +
 //                    beta_th_d_N * th .* log(d) .* log(N_ini), 
 //                    sigma);
// }
///////////////////////////////////////////////////
//////////// Working very good FORMULA 2 //////////
///////////////////////////////////////////////////
  // log(dev) ~ normal(alpha + 
  //                   beta_th * th + 
  //                   beta_d * log(d) + 
  //                   beta_N_ini * log(N_ini) +
  //                   beta_th_d * log(th) .* log(d) +
  //                   beta_th_N * log(th) .* log(N_ini) +
  //                   beta_d_N * log(d) .* log(N_ini) +
  //                   beta_th_d_N * log(th) .* log(d) .* log(N_ini), 
  //                   sigma);
}

generated quantities {
  vector[N] y_pred_log;     // predicted values on log scale
  vector[N] y_pred;         // predicted values on original scale
  
  for (i in 1:N) {
    y_pred_log[i] = alpha + 
                    beta_th * th[i] + 
                    beta_d * log(d)[i] + 
                    beta_N_ini * log(sqrt(N_ini))[i] +
                    beta_th_d * th[i] * log(d)[i] +
                    beta_th_N * th[i] * log(sqrt(N_ini))[i] +
                    beta_d_N * log(d)[i] * log(sqrt(N_ini))[i] +
                    beta_th_d_N * th[i] * log(d)[i] * log(sqrt(N_ini))[i];
    y_pred[i] = exp(y_pred_log[i]);  // back-transform to original scale
  }
}
