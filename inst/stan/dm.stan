functions {
  real dm_explicit_lpmf(int[] y, vector alpha) {
    real sum_alpha = sum(alpha);
    return lgamma(sum_alpha) - lgamma(sum(y) + sum_alpha)
           + lgamma(sum(y)+1) - sum(lgamma(to_vector(y)+1))
           + sum(lgamma(to_vector(y) + alpha)) - sum(lgamma(alpha));
  }
}

data {
  int K; // categories
  int N; // number of samples
  int y [N, K]; // counts matrix
}

parameters {
  real <lower=0> kappa;
  vector [K] alpha;
  vector<lower=0> [N] beta_sigma;
  vector [K] beta_z [N];
}

transformed parameters {
  vector [K] beta [N];
  for(i in 1:N) {
    beta[i] = 0 + beta_sigma[i]*beta_z[i];
  }
}

model {
  target += exponential_lpdf(2+kappa|0.01);
  target += normal_lpdf(alpha|0,3);
  for(i in 1:N) {
    target += cauchy_lpdf(beta_sigma[i]|0, 1);
    target += normal_lpdf(beta_z[i]|0,1);
    target += dm_explicit_lpmf(y[i]|kappa*softmax(alpha + beta[i]));
  }
}

generated quantities {
  int y_hat [N, K];
  simplex [K] p [N];
  for(i in 1:N) {
    p[i] = dirichlet_rng(kappa*softmax(alpha + beta[i]));
    y_hat[i] = multinomial_rng(p[i], sum(y[i]));
  }
}
