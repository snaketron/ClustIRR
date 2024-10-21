functions {
  real dm_lpmf(int[] y, vector alpha) {
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
  int x [2]; // not used
}

parameters {
  real <lower=0> kappa;
  vector [K] alpha;
  vector [K] beta [N];
}

model {
  target += exponential_lpdf(2+kappa|0.01);
  target += normal_lpdf(alpha|0,3);
  for(i in 1:N) {
    target += normal_lpdf(beta[i]|0, 1);
    target += dm_lpmf(y[i]|kappa*softmax(alpha + beta[i]));
  }
}

generated quantities {
  int y_hat [N, K];
  real log_lik [N];
  simplex [K] p [N];
  vector [K] delta [N*(N-1)/2];
  int k;
  
  k = 1;
  for(i in 1:N) {
    p[i] = dirichlet_rng(kappa*softmax(alpha + beta[i]));
    y_hat[i] = multinomial_rng(p[i], sum(y[i]));
    log_lik[i] = dm_lpmf(y[i]|kappa*softmax(alpha + beta[i]));
    if(i != N) {
      for(j in (i+1):N) {
        delta[k] = beta[i]-beta[j];
        k = k + 1;
      }
    }
  }
}
