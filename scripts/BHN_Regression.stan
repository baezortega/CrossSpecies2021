/*
* Bayesian hierarchical normal regression model (with offset)
* Adrian Baez-Ortega, 2020
* Adapted from
* https://datascienceplus.com/bayesian-regression-with-stan-beyond-normality
* https://stats.stackexchange.com/questions/11182/when-to-use-an-offset-in-a-poisson-regression
*/

data {
    int N;           // number of observations
    int G;           // number of groups (observation categories)
    int K;           // number of parameters (columns in the model matrix)
    int group[N];    // sample group indices
    real y[N];       // response values
    matrix[N, K] X;  // model matrix
    vector[N] t;     // offset (untransformed)
}

parameters {
    vector[K] beta;         // global regression coefficients
    matrix[G, K] beta_g;    // group-specific regression coefficients
    real<lower=0> sigma_g;  // std deviation of the beta_g around beta
    real<lower=0> sigma;    // std deviation of response
}

transformed parameters {
    vector<lower=0>[N] mu;  // linear predictor with offset
    for (n in 1:N) {
        mu[n] = t[n] * (X[n] * beta_g[group[n]]');
    }
}

model {
    // Priors for model parameters
    // Best to use Stan's default improper uniform priors
    // as noninformative priors for beta, sigma, sigma_g
    for (g in 1:G) {
        beta_g[g] ~ normal(beta, sigma_g);
    }
    
    // Normal likelihood for response
    y ~ normal(mu, sigma);
}
