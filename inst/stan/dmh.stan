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
    int <lower=1> G [N]; // groups
    int y [N, K]; // counts matrix
    int<lower=0, upper=1> compute_delta; // should we compute delta
}


transformed data {
    int N_delta;
    int K_delta;
    K_delta = 0;
    N_delta = 0;
    if(compute_delta==1) {
        K_delta = K;
        N_delta = max(G)*(max(G)-1)/2;
    }
}


parameters {
    vector <lower=0> [max(G)] kappa;
    vector [K] alpha;
    vector [K] beta_mu [max(G)];
    vector <lower=0> [max(G)] beta_sigma;
    vector [K] beta_z [N];
}


transformed parameters {
    vector [K] beta [N];
    for(i in 1:N) {
        beta[i] = beta_mu[G[i]] + beta_sigma[G[i]] * beta_z[i];
    }
}

model {
    target += cauchy_lpdf(beta_sigma|0, 1);
    target += normal_lpdf(alpha|0,3);
    for(j in 1:max(G)) {
        target += exponential_lpdf(2+kappa[j]|0.01);
        target += normal_lpdf(beta_mu[j]|0, 1);
    }
    for(i in 1:N) {
        target += std_normal_lpdf(beta_z[i]);
        target += dm_lpmf(y[i]|kappa[G[i]]*softmax(alpha + beta[i]));
    }
}

generated quantities {
    int y_hat [N, K];
    real log_lik [N];
    simplex [K] p;
    vector [K_delta] delta [N_delta];
    vector [K_delta] epsilon [N_delta];
    int k = 1;
    
    for(i in 1:N) {
        p = dirichlet_rng(kappa[G[i]]*softmax(alpha + beta[i]));
        y_hat[i] = multinomial_rng(p, sum(y[i]));
        log_lik[i] = dm_lpmf(y[i]|kappa[G[i]]*softmax(alpha + beta[i]));
    }
    
    if(compute_delta==1) {
        for(i in 1:max(G)) {
            for(j in (i+1):max(G)) {
                delta[k] = beta_mu[i] - beta_mu[j];
                epsilon[k] = softmax(alpha+beta_mu[i]) - softmax(alpha+beta_mu[j]);
                k += 1;
            }
        }
    }
}
