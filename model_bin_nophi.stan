// spatio-temporal binomial logit regression with explanatory variables
data {
  int<lower = 0> T; // number of time points
  int<lower = 0> N; // number of regions
  int<lower = 0> n_obs; // number of observations
  
  array[N * T] int<lower = 0> response; // whatever the response data is
  
  int<lower = 0> N_edges; // neighborhood
  array[N_edges] int<lower = 1, upper = N> node1;
  array[N_edges] int<lower = 1, upper = N> node2;
  
  array[N * T] int<lower = 0> N_all; // total counts

  array[n_obs] int<lower=0> ind; // indices for non-missing observations
  array[n_obs] int<lower=0> region_ind; // region indices for non-missing observations
  
  int<lower = 0> n_covariates_t; // number of temporal covariates
  matrix[N * T, n_covariates_t] covariates_t; // covariates

  array[N * T] int<lower=0> time; // number of current time point
}
parameters {
  // yearly levels
  vector[T] c;
  real<lower = 0> sigma_c;
  
  // regression coefficients, helper variables
  vector[n_covariates_t] b_1; // intercept for the first state
  vector[n_covariates_t] b_t; // intercept for the last state
  array[n_covariates_t] simplex[T - 1] b; // for the intercepts between
}
transformed parameters{
  // actual regression coefficients
  matrix[T, n_covariates_t] beta;
  for (i in 1:n_covariates_t) {
    beta[, i] = b_1[i] + (b_t[i] - b_1[i]) * cumulative_sum(append_row(0, b[i, ]));
  }
}
model {
  // yearly levels as a random walk
  c[1] ~ std_normal();
  c[2:T] ~ normal(c[1:(T - 1)], sigma_c);
  sigma_c ~ std_normal();

  // monotonic regression coefficients
  b_1 ~ std_normal();
  for (i in 1:n_covariates_t) {
    b_t[i] ~ std_normal();
  }

  // assign binomial distribution with logit link
  target += binomial_logit_lpmf(response[ind] | N_all[ind],
                              c[time[ind]] + rows_dot_product(covariates_t, beta[time, ])[ind]);
}
generated quantities {
  // get the log likelihood values
  vector[n_obs] log_lik;
  
  for(i in 1:n_obs) {
    log_lik[i] = binomial_logit_lpmf(response[ind[i]] | N_all[ind[i]],
    c[time[ind[i]]] + rows_dot_product(covariates_t, beta[time, ])[ind[i]]);
  }
}
