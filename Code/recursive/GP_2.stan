// Stan model for 2D Gaussian Process regression of species abundance data
data {
  int N_0;                      // number of observations
  int N_50;                      // number of observations
  array[N_0] vector[2] X_0;                // observation locations (x,y coordinates)
  array[N_50] vector[2] X_50;                // observation locations (x,y coordinates)
  vector[N_0] y_0;                         // observed species abundances
  vector[N_50] y_50;                         // observed species abundances
  int N_pred;                 // number of prediction locations
  array[N_pred] vector[2] X_pred;
  // 
  // Priors
  vector[2] alpha_prior;
  vector[2] rho_prior;
  vector[2] sigma_prior;
  vector[2] mu_prior;
}
// 
parameters {
  vector<lower=0>[2] alpha;       // GP variance parameter
  real<lower=0> rho;              // GP length scale parameter
  real<lower=0> sigma;            // noise scale
  vector[2] mu;                   // mean parameter
}
// 
model {
  // Define the mean vector and covariance matrix
  matrix[N_0, N_0] K_0;
  matrix[N_50, N_50] K_50;
  matrix[N_0, N_0] L_K_0;
  matrix[N_50, N_50] L_K_50;
  
  // Priors for GP hyperparameters
  alpha ~ normal(alpha_prior[1], alpha_prior[2]);
  rho ~ normal(rho_prior[1], rho_prior[2]);
  sigma ~ normal(sigma_prior[1], sigma_prior[2]);
  mu ~ normal(mu_prior[1], mu_prior[2]);
  
    K_0 = gp_exp_quad_cov(X_0, alpha[1], exp(rho));
    K_50 = gp_exp_quad_cov(X_50, alpha[2], exp(rho));
    
    // Add noise to diagonal
    for (n in 1:N_0){
    	K_0[n, n] = K_0[n, n] + square(sigma);
    }
    for (n in 1:N_50){
    	K_50[n, n] = K_50[n, n] + square(sigma);
    }

  // Use Cholesky decomposition for numerical stability
  L_K_0 = cholesky_decompose(K_0);
  L_K_50 = cholesky_decompose(K_50);
  
  // GP likelihood
  y_0 ~ multi_normal_cholesky(rep_vector(mu[1], N_0), L_K_0);
  y_50 ~ multi_normal_cholesky(rep_vector(mu[2], N_50), L_K_50);
}

generated quantities {
  // Predict species abundance at new locations
  vector[N_pred] f_pred_0;
  vector[N_pred] f_pred_50;
  {
    matrix[N_0, N_0] K_0 = gp_exp_quad_cov(X_0, alpha[1], exp(rho));
    matrix[N_50, N_50] K_50 = gp_exp_quad_cov(X_50, alpha[2], exp(rho));
    matrix[N_0, N_pred] k_pred_0 = gp_exp_quad_cov(X_0, X_pred, alpha[1], exp(rho));
    matrix[N_50, N_pred] k_pred_50 = gp_exp_quad_cov(X_50, X_pred, alpha[2], exp(rho));
    matrix[N_0, N_0] L_K_0;
    matrix[N_50, N_50] L_K_50;
    vector[N_0] K_div_y_0;
    vector[N_50] K_div_y_50;

    // Add observation noise
    for (n in 1:N_0)
      K_0[n, n] = K_0[n, n] + square(sigma);
    for (n in 1:N_50)
      K_50[n, n] = K_50[n, n] + square(sigma);

    L_K_0 = cholesky_decompose(K_0);
    L_K_50 = cholesky_decompose(K_50);
    K_div_y_0 = mdivide_left_tri_low(L_K_0, y_0 - mu[1]);
    K_div_y_50 = mdivide_left_tri_low(L_K_50, y_50 - mu[2]);
    K_div_y_0 = mdivide_right_tri_low(K_div_y_0', L_K_0)';
    K_div_y_50 = mdivide_right_tri_low(K_div_y_50', L_K_50)';

    f_pred_0 = mu[1] + (k_pred_0' * K_div_y_0);
    f_pred_50 = mu[2] + (k_pred_50' * K_div_y_50);
  }
}
