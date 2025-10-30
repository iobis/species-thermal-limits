functions {
  matrix cov_GPL1(matrix x, real sq_alpha, real sq_rho, real delta) {
        int N = dims(x)[1];
        matrix[N, N] K;
        for (i in 1:(N-1)) {
          K[i, i] = sq_alpha + delta;
          for (j in (i + 1):N) {
            K[i, j] = sq_alpha * exp(-sq_rho * x[i,j] );
            K[j, i] = K[i, j];
          }
        }
        K[N, N] = sq_alpha + delta;
        return K;
    }
}
data {
  int N;
  int<lower=1> D; // number of pars per species
  int N_spp;
  array[N] int sid; // species id
  array[N] real sst; // temps at surveys
  array[N] int y; // 0/1 obs
  matrix[N_spp,N_spp] M; // distance matrix
}
transformed data {
  real delta = 1e-9;
}
parameters {
  real<lower=0> rho; // GP max cov
  real<lower=0> sigma; // GP scale
  vector<lower=0>[D] tau; // scale of each effect
  vector[D] alpha; // mean of each effect on latent scale
  cholesky_factor_corr[D] L_Omega; // correlation structure of effects
  matrix[N_spp, D] eta; // unscaled effects (z scores)
  real<lower=0,upper=1> p;    // detection probability
}
transformed parameters {
  matrix[N_spp, D] f; // scaled effects
  {
    //matrix[N, N] K = cov_exp_quad(x, 1.0, rho);
    matrix[N,N] K = cov_GPL1(M,sigma,rho,delta);
    matrix[N,N] L_K = cholesky_decompose(K);
    f = L_K * eta * diag_pre_multiply(tau, L_Omega)'; // offsets
    for (i in 1:D) f[,i] = f[,i] + alpha[i]; // alpha = means for each param
  }
}
model {
  
  // priors
  p ~ beta(2,2);
  rho ~ inv_gamma(5, 5);

  tau ~ std_normal();
  alpha[1] ~ normal(3,1.5);
  alpha[2] ~ normal(1.5,0.5);
  alpha[3] ~ normal(0,0.5); // will be converted to prob [0,1] scale

  sigma ~ std_normal();
  L_Omega ~ lkj_corr_cholesky(4);
  to_vector(eta) ~ std_normal();

  // loop over surveys
    for ( i in 1:N ) {
        // prob of occupancy given suitability
        //real q = exp( -0.5*(( sst[i] - tmu[sid[i]] ) / tsd[sid[i]])^2 );
        real log_q =  -0.5*(( sst[i] - exp(f[sid[i],1]) ) / exp(f[sid[i],2]))^2 + log_inv_logit(f[sid[i],3]) ;
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
generated quantities {
  matrix[D, D] Omega;
  Omega = L_Omega * L_Omega';
}
