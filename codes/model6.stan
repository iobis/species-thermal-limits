// basic occupancy model with thermal tolerance function
// flexible maximum occupancy probability for each species
// add partial pooling across species
data{
    int N;
    int N_spp;
    array[N] int sid; // species id
    array[N] real sst; // temps at surveys
    array[N] int y; // 0/1 obs
}
parameters{
    vector[N_spp] log_tmu;  // temp mean for each species
    vector[N_spp] log_tsd;  // temp sd for each species
    vector<lower=0,upper=1>[N_spp] tomax; // max occupancy prob
    real<lower=0,upper=1> p;    // detection probability
    // pooling parameters - starting simple with no covariance matrix
    vector[2] spp_means; // mean values for each species parameter
    vector<lower=0>[2] spp_sds; // standard deviations for each species parameter
}
model{
    p ~ beta(2,2);
    tomax ~ beta(5,1);
    //tmu ~ normal(20,3);
    //tsd ~ normal(5,1);
    spp_means[1] ~ lognormal(20,3);
    spp_means[2] ~ lognormal(5,1);
    spp_sds ~ exponential(1);
    log_tmu ~ normal( spp_means[1] , spp_sds[1] );
    log_tsd ~ normal( spp_means[2] , spp_sds[2] );

    // loop over surveys
    for ( i in 1:N ) {
        // prob of occupancy given suitability
        //real q = exp( -0.5*(( sst[i] - tmu[sid[i]] ) / tsd[sid[i]])^2 );
        real log_q =  -0.5*(( sst[i] - exp(log_tmu[sid[i]]) ) / exp(log_tsd[sid[i]]))^2 + log(tomax[sid[i]]) ;
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
