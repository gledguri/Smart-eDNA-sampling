// Stan model using functions for cleaner GP implementation
functions {
	void gp_marginal_lp(array[] vector X, vector y, real alpha, real rho, 
											real sigma, real mu_val) {
		int N = size(X);
		matrix[N, N] K = gp_exp_quad_cov(X, alpha, exp(rho));
		
		// Add noise to diagonal covariate matrix
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

data {
  int N_depths; // Number of depth levels (0,50,150,300, and 500m)
  int N_total; // Total number of observations
  array[N_total] vector[2] X; // 2D coordinates (lon, lat)
  vector[N_total] y; // eDNA observations
  // Depth indexing
  array[N_depths] int N_by_depth; // Number of observations per each depth
  array[N_depths] int start_idx; // Starting index for each depth
  // Priors
  vector[2] alpha_prior; // mean and sd (declared in R)
  vector[2] rho_prior; // mean and sd (declared in R)
  vector[2] sigma_prior; // mean and sd (declared in R)
  vector[2] mu_prior; // mean and sd (declared in R)
  // Prediction data
  int N_pred; // Total number of prediction points
  int N_depths_pred; // Number of depth levels for predictions (same as N_depths)
  array[N_pred] vector[2] X_pred; // 2D coordinates for predictions
  array[N_depths_pred] int N_by_depth_pred; // Number of points of predictions per each depth
  array[N_depths_pred] int start_idx_pred; // Starting index for each depth in predictions
}

parameters {
  real<lower=0> alpha;
  real<lower=-5> rho;
  real<lower=0> sigma;
  vector<lower=-5>[N_depths] mu;
}

model {
  // Priors
  alpha ~ normal(alpha_prior[1], alpha_prior[2]);
  rho ~ normal(rho_prior[1], rho_prior[2]);
  sigma ~ normal(sigma_prior[1], sigma_prior[2]);
  mu ~ normal(mu_prior[1], mu_prior[2]);
  
  // GP for each depth
  for (d in 1:N_depths) {
    int N_d = N_by_depth[d]; // Number of samples at depth d
    int start = start_idx[d]; // Starting index for depth d
    int end = start + N_d - 1; // Ending index for depth d
    
    array[N_d] vector[2] X_d = X[start:end]; // Get the coordinates of points for depth d
    vector[N_d] y_d = y[start:end]; // Get the 'observed' eDNA values at depth for X_d coordaintes
    
    gp_marginal_lp(X_d, y_d, alpha, rho, sigma, mu[d]); // Call the GP function 
  }
}

generated quantities {
  // Declare the vector where the predictions will be stored
  vector[N_pred] y_pred;

  // Predictions loop per depth
  for (d in 1:N_depths_pred) {
    // Get the index for the 'observed' data by depth
    int N_d = N_by_depth[d]; // Declare the number of samples at depth d
    int start = start_idx[d]; // Starting index for depth d
    int end = start + N_d - 1; // Ending index for depth d
    array[N_d] vector[2] X_d = X[start:end]; // Get the coordinates of points for depth d
    vector[N_d] y_d = y[start:end]; // Get the 'observed' eDNA values at depth for X_d coordaintes

    // Get the index for the predicted data by depth
    int N_d_p = N_by_depth_pred[d]; // Declare the number of samples at depth d (for predictions)
    int start_p = start_idx_pred[d]; // Starting index for depth d (for predictions)
    int end_p = start_p + N_d_p - 1; // Ending index for depth d (for predictions)
    array[N_d_p] vector[2] X_d_pred = X_pred[start_p:end_p]; // Get the coordinates of points for depth d (for predictions)

    // All predictions are made for specified depth using the GP function
    y_pred[start_p:end_p] = gp_predict_depth(N_d, N_d_p, X_d, y_d, X_d_pred, alpha, rho, sigma, mu[d]);
  }
}
