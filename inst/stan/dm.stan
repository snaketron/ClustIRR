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
    int<lower=0, upper=1> compute_delta; // should we compute delta
}


transformed data {
    int N_beta;
    int N_delta;
    int K_delta;
    K_delta = 0;
    N_delta = 0;
    if(compute_delta==1) {
        K_delta = K;
        N_delta = N*(N-1)/2;
    }
    if(N==2) {
        N_beta = 1;
    } else {
        N_beta = N;
    }
}


parameters {
    real <lower=0> kappa;
    vector [K] alpha;
    vector [K] beta [N_beta];
}

model {
    target += exponential_lpdf(2+kappa|0.01);
    target += normal_lpdf(alpha|0,3);
    
    if(N==2) {
        target += normal_lpdf(beta[1]|0, 1);
        target += dm_lpmf(y[1]|kappa*softmax(alpha + beta[1]*1));
        target += dm_lpmf(y[2]|kappa*softmax(alpha + beta[1]*-1));
    } 
    else {
        for(i in 1:N) {
            target += normal_lpdf(beta[i]|0, 1);
            target += dm_lpmf(y[i]|kappa*softmax(alpha + beta[i]));
        }
    }
}

generated quantities {
    int y_hat [N, K];
    real log_lik [N];
    simplex [K] p [N];
    vector [K_delta] delta [N_delta];
    vector [K_delta] epsilon [N_delta];
    int k;
    
    if(N==2) {
        p[1] = dirichlet_rng(kappa*softmax(alpha + beta[1]*1));
        p[2] = dirichlet_rng(kappa*softmax(alpha + beta[1]*-1));
        delta[1] = beta[1];
        epsilon[1] = softmax(alpha+beta[1])-softmax(alpha-beta[1]);
        log_lik[1] = dm_lpmf(y[1]|kappa*softmax(alpha + beta[1]));
        log_lik[2] = dm_lpmf(y[2]|kappa*softmax(alpha - beta[1]));
        y_hat[1] = multinomial_rng(p[1], sum(y[1]));
        y_hat[2] = multinomial_rng(p[2], sum(y[2]));
    } 
    else {
        for(i in 1:N) {
            p[i] = dirichlet_rng(kappa*softmax(alpha + beta[i]));
            y_hat[i] = multinomial_rng(p[i], sum(y[i]));
            log_lik[i] = dm_lpmf(y[i]|kappa*softmax(alpha + beta[i]));
        }
        
        k = 1;
        for(i in 1:N) {
            if(compute_delta==1) {
                if(i != N) {
                    for(j in (i+1):N) {
                        delta[k] = beta[i]-beta[j];
                        epsilon[k] = softmax(alpha+beta[i])-softmax(alpha+beta[j]);
                        k = k + 1;
                    }
                }
            }
        }
    }
}
