functions {
    real dm_lpmf(array[] int y, vector alpha) {
        real sum_alpha = sum(alpha);
        return lgamma(sum_alpha) - lgamma(sum(y) + sum_alpha)
        + lgamma(sum(y)+1) - sum(lgamma(to_vector(y)+1))
        + sum(lgamma(to_vector(y) + alpha)) - sum(lgamma(alpha));
    }
}


data {
    int K; // categories
    int N; // number of samples
    array[N, K] int y; // counts matrix
    int<lower=0, upper=1> compute_delta; // should we compute delta
}


transformed data {
    int N_delta;
    int K_delta;
    K_delta = 0;
    N_delta = 0;
    if(compute_delta==1) {
        K_delta = K;
        N_delta = N*(N-1)/2;
    }
}


parameters {
    real <lower=0> kappa;
    vector [K] alpha;
    vector [K] beta [N];
}

model {
    target += exponential_lpdf(kappa|0.01);
    target += normal_lpdf(alpha|0,1);
    
    for(i in 1:N) {
        target += normal_lpdf(beta[i]|0,1);
        target += dm_lpmf(y[i]|kappa*softmax(alpha + beta[i]));
    }
}

generated quantities {
    array[N, K] int y_hat;
    array[N] real log_lik;
    simplex[K] p;
    array[N_delta] vector[K_delta] epsilon;
    int k = 1;
    
    for(i in 1:N) {
        p = dirichlet_rng(kappa*softmax(alpha + beta[i]));
        y_hat[i] = multinomial_rng(p, sum(y[i]));
        log_lik[i] = dm_lpmf(y[i]|kappa*softmax(alpha + beta[i]));
    }
    
    if(compute_delta==1) {
        for(i in 1:(N-1)) {
            for(j in (i+1):N) {
                epsilon[k] = softmax(alpha+beta[i]) - softmax(alpha+beta[j]);
                k += 1;
            }
        }
    }
}
