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
}
model{
    p ~ beta(2,2);
    tmu ~ normal(20,3);
    tsd ~ normal(5,1);
    tomax ~ beta(5,1);

    // loop over surveys
    for ( i in 1:N ) {
        // prob of occupancy given suitability
        //real q = exp( -0.5*(( sst[i] - tmu[sid[i]] ) / tsd[sid[i]])^2 );
        real log_q =  -0.5*(( sst[i] - tmu[sid[i]] ) / tsd[sid[i]])^2 + log(tomax[sid[i]]) ;
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
                log(1-p) + log_q
            );
        }
    }
}
