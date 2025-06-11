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
    //vector[N_spp] log_tmu;  // temp mean for each species
    //vector[N_spp] log_tsd;  // temp sd for each species
    //vector<lower=0,upper=1>[N_spp] tomax; // max occupancy prob
    matrix[3,N_spp] z;
    real<lower=0,upper=1> p;    // detection probability
    // pooling parameters - starting simple with no covariance matrix
    vector<lower=0>[3] sigma;
    cholesky_factor_corr[3] L_Rho;
    vector[3] spp_means; // mean values for each species parameter
}
transformed parameters {
   matrix[N_spp,3] log_spp_pars;
   vector[N_spp] log_tmu;  // temp mean for each species
   vector[N_spp] log_tsd;  // temp sd for each species
   vector<lower=0,upper=1>[N_spp] tomax; // max occupancy prob
   log_spp_pars = (diag_pre_multiply(sigma, L_Rho) * z)';
   // convert to ordinary parameter vectors for convenience
   // implies these priors:
   //log_tmu ~ normal( spp_means[1] , sigma[1] );
   //log_tsd ~ normal( spp_means[2] , sigma[2] );
   log_tmu = spp_means[1] + log_spp_pars[,1];
   log_tsd = spp_means[2] + log_spp_pars[,2];
   tomax = inv_logit( spp_means[3] + log_spp_pars[,3] );
}
model{
    p ~ beta(2,2);
    //tomax ~ beta(5,1);
    //tmu ~ normal(20,3);
    //tsd ~ normal(5,1);
    spp_means[1] ~ normal(20,3);
    spp_means[2] ~ normal(5,1);
    spp_means[3] ~ normal(0,0.5); // will be converted to prob [0,1] scale
    sigma ~ exponential(1);
    L_Rho ~ lkj_corr_cholesky(4);
    to_vector(z) ~ std_normal();

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
