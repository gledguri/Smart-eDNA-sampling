// Stan model for 2D Gaussian Process regression of species abundance data
data {
  int<lower=1> N;                      // number of observations
  int<lower=1> N_pred;                 // number of prediction locations
  array[N] vector[2] X;                // observation locations (x,y coordinates)
  vector[N] y;                         // observed species abundances
  array[N_pred] vector[2] X_pred;      // prediction locations
  
  // Priors
  vector[2] alpha_prior;
  vector[2] rho_prior;
  vector[2] sigma_prior;
  vector[2] mu_prior;
}

parameters {
  real<lower=0> alpha;            // GP variance parameter
  real<lower=0> rho;              // GP length scale parameter
  real<lower=0> sigma;            // noise scale
  real mu;                        // mean parameter
}

model {
  // Define the mean vector and covariance matrix
  vector[N] f;
  matrix[N, N] K;
  matrix[N, N] L_K;
  
  // Priors for GP hyperparameters
  alpha ~ normal(alpha_prior[1], alpha_prior[2]);
  rho ~ normal(rho_prior[1], rho_prior[2]);
  sigma ~ normal(sigma_prior[1], sigma_prior[2]);
  mu ~ normal(mu_prior[1], mu_prior[2]);
  
  // Construct covariance matrix
  K = gp_exp_quad_cov(X, alpha, exp(rho));
  
  // Add noise component
  for (n in 1:N)
    K[n, n] = K[n, n] + square(sigma);
  
  // Use Cholesky decomposition for numerical stability
  L_K = cholesky_decompose(K);
  
  // GP likelihood
  y ~ multi_normal_cholesky(rep_vector(mu, N), L_K);
}

generated quantities {
  // Predict species abundance at new locations
  vector[N_pred] f_pred;
  {
    matrix[N, N] K = gp_exp_quad_cov(X, alpha, exp(rho));
    matrix[N, N_pred] k_pred = gp_exp_quad_cov(X, X_pred, alpha, exp(rho));
    matrix[N, N] L_K;
    vector[N] K_div_y;

    // Add observation noise
    for (n in 1:N)
      K[n, n] = K[n, n] + square(sigma);

    L_K = cholesky_decompose(K);
    K_div_y = mdivide_left_tri_low(L_K, y - mu);
    K_div_y = mdivide_right_tri_low(K_div_y', L_K)';

    f_pred = mu + (k_pred' * K_div_y);
  }
}
