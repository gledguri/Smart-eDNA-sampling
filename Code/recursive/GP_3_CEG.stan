// Stan model using functions for cleaner GP implementation
functions {
  void gp_marginal_lp(array[] vector X, vector y, real alpha, real rho, 
                      real sigma, real mu_val) {
    int N = size(X);
    matrix[N, N] K = gp_exp_quad_cov(X, alpha, exp(rho));
    
    // Add noise to diagonal
    for (n in 1:N) {
      K[n, n] = K[n, n] + square(sigma);
    }
    
    matrix[N, N] L_K = cholesky_decompose(K);
    target += multi_normal_cholesky_lpdf(y | rep_vector(mu_val, N), L_K);
  }
  
  vector gp_predict_depth(int N,
                          int N_pred,
                          array[] vector X,
                          vector y,
                          array[] vector X_pred,
                          real alpha,
                          real rho,
                          real sigma,
                          real mu) {
    matrix[N, N] K = gp_exp_quad_cov(X, alpha, exp(rho));
    matrix[N, N_pred] k_pred = gp_exp_quad_cov(X, X_pred, alpha, exp(rho));
    matrix[N, N] L_K;
    vector[N] K_div_y;
    vector[N] y_centered = y - mu;
    
    for (n in 1:N)
      K[n, n] = K[n, n] + square(sigma);
    L_K = cholesky_decompose(K);
    K_div_y = mdivide_left_tri_low(L_K, y_centered);
    K_div_y = mdivide_right_tri_low(K_div_y', L_K)';
    return mu + (k_pred' * K_div_y);
  }
}
// 
data {
  int N_depths;
  int N_total;
  array[N_total] vector[2] X;
  vector[N_total] y;
  array[N_depths] int N_by_depth;
  array[N_depths] int start_idx;
  // Priors
  vector[2] alpha_prior;
  vector[2] sigma_prior;
  vector[2] mu_prior;
  vector[2] rho_prior;
  vector[2] rho_sd_prior;
  // Prediction data
  int N_pred;
  int N_depths_pred;
  array[N_pred] vector[2] X_pred;
  array[N_depths_pred] int N_by_depth_pred;
  array[N_depths_pred] int start_idx_pred;
}

parameters {
  vector<lower=0>[N_depths] alpha;
  real rho;
  real<lower=0> rho_sd;
  real<lower=0> sigma;
  vector[N_depths] mu;
}

model {
  // Priors
  alpha ~ normal(alpha_prior[1], alpha_prior[2]);
  // rho ~ normal(rho_prior[1], rho_prior[2]);
  rho ~ normal(rho_prior[1], rho_sd);
  rho_sd ~ normal(rho_sd_prior[1],rho_sd_prior[2]);
  sigma ~ normal(sigma_prior[1], sigma_prior[2]);
  mu ~ normal(mu_prior[1], mu_prior[2]);
  
  // Process each depth
  for (d in 1:N_depths) {
    int N_d = N_by_depth[d];
    int start = start_idx[d];
    int end = start + N_d - 1;
    
    array[N_d] vector[2] X_d = X[start:end];
    vector[N_d] y_d = y[start:end];
    
    gp_marginal_lp(X_d, y_d, alpha[d], rho, sigma, mu[d]);
  }
}

generated quantities {
  // Predictions for specified depth using your GP function
  vector[N_pred] y_pred;

  // Predictions
  for (d in 1:N_depths_pred) {
    // Changed: Use training data indices for this depth
    int N_d = N_by_depth[d];
    int start = start_idx[d];
    int end = start + N_d - 1;
    array[N_d] vector[2] X_d = X[start:end];
    vector[N_d] y_d = y[start:end];

    // Changed: Use prediction data indices for this depth
    int N_d_p = N_by_depth_pred[d];
    int start_p = start_idx_pred[d];
    int end_p = start_p + N_d_p - 1;
    array[N_d_p] vector[2] X_d_pred = X_pred[start_p:end_p];

    // Changed: Use correct dimensions
    y_pred[start_p:end_p] = gp_predict_depth(N_d, N_d_p, X_d, y_d, X_d_pred, alpha[d], rho, sigma, mu[d]);
  }
}
