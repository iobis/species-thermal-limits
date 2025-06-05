// basic occupancy model with thermal tolerance function
// flexible maximum occupancy probability for each species
data{
    int N;
    int N_spp;
    array[N] int sid; // species id
    array[N] real sst; // temps at surveys
    array[N] int y; // 0/1 obs
}
parameters{
    vector<lower=0>[N_spp] tmu;  // temp mean for each species
    vector<lower=0>[N_spp] tsd;  // temp sd for each species
    vector<lower=0,upper=1>[N_spp] tomax; // max occupancy prob
    real<lower=0,upper=1> p;    // detection probability

    // Hyperparameters for partial pooling
    real mu_tmu;
    real<lower=0> sigma_tmu;

    real mu_tsd;
    real<lower=0> sigma_tsd;
}
model{
    // Hyperpriors
    mu_tmu ~ normal(20, 5);
    sigma_tmu ~ normal(0, 2.5);

    mu_tsd ~ normal(5, 2);
    sigma_tsd ~ normal(0, 2.5);

    // Species-level priors with partial pooling
    tmu ~ normal(mu_tmu, sigma_tmu);
    tsd ~ normal(mu_tsd, sigma_tsd);

    p ~ beta(2,2);
    //tmu ~ normal(20,3);
    //tsd ~ normal(5,1);
    tomax ~ beta(5,1);

    // loop over surveys
    for ( i in 1:N ) {
        // prob of occupancy given suitability
        //real q = exp( -0.5*(( sst[i] - tmu[sid[i]] ) / tsd[sid[i]])^2 );
        int s = sid[i];
        real log_q = -0.5 * square((sst[i] - tmu[s]) / tsd[s]) + log(tomax[s]);
        //real q = 0.5;
        // check survey
        if ( y[i]==1 )
            // prob occupied and detected
            target += log(p) + log_q;
        if ( y[i]==0 ) {
            // absence, so use mixture
            // Prob(0|q,p) = (1-q) + q*(1-p)
            target += log_sum_exp(
                log1m_exp(log_q) ,
                log1m(p) + log_q
            );
        }
    }
}
