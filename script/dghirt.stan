// Stan model syntax combining the Modified DG-IRT for response accuracy and
// DG-LNRT for response time

  
data{
    int <lower=1> J;                       // number of examinees          
    int <lower=1> I;                       // number of items
    int <lower=1> n_obs;                   // number of observations (J x I - missing responses)
    array[n_obs] int<lower=1> p_loc;       // person indicator   
    array[n_obs] int<lower=1> i_loc;       // item indicator
    array[n_obs] real RT;                  // vector of log of responses
    array[n_obs] int<lower=0,upper=1> Y;   // vector of item responses
}

parameters {
  
// Parameters common for both response time and response accuracy component  
  
  vector<lower=0,upper=1>[I] pC; // vector of length I for the probability of item compromise status
  
  vector<lower=0,upper=1>[J] pH; // vector of length J for the probability of examinee item peknowledge 

// Parameters for the item parameters component  

  real mu_beta;                     // mean for time intensity parameters
  real mu_alpha;                    // mean for log of time discrimination parameters
  real mu_b;                        // mean for item difficulty parameters
  vector<lower=0>[3] sigma_I;       // vector of standard deviations for item parameters
  cholesky_factor_corr[3] Lomega_I; // Cholesky factors of 3x3 correlation matrix for item parameters
  
// Parameters for the person parameters component  
  
  real mu_thetat;                   // mean for theta_t
                                    // we will fix mu_taut to 0 in the the transformed parameters block
  vector<lower=0>[2] sigma_P;       // vector of standard deviations for latent trait parameters
  cholesky_factor_corr[2] Lomega_P; // Cholesky factors of 2x2 correlation matrix for person parameters


  real<lower=0> mu_delta_tau;         // mean for change from taut to tau_c
  real<lower=0> mu_delta_theta;       // mean for change from thetat to thetac
  vector<lower=0>[2] sigma_D;        // vector of standard deviations for change in person parameters
  cholesky_factor_corr[2] Lomega_D;  // Cholesky factors of 2x2 correlation matrix for person parameters
 
  array[I] vector[3] item;                 // I x 3 item parameter matrix
  array[J] vector[2] person;               // J x 2 person parameter matrixx
  array[J] vector<lower=0>[2] delta;       // J x 2 delta (change in person parameters) matrixx
}


transformed parameters{
  vector[2] mu_P = [0,mu_thetat]';                 // vector for mean of person parameters
                                                   // mu_taut is fixed to 0
                                                   
  vector[3] mu_I = [mu_alpha,mu_beta,mu_b]';       // vector for mean of item parameters
  
  vector[2] mu_D = [mu_delta_tau,mu_delta_theta]'; // vector for mean of change parameters

  vector[I] alpha;                        // a dedicated vector for time-discrimination
  vector[I] beta;                         // a dedicated vector for time-intensity
  vector[I] b_star;                       // a dedicated vector for item difficulty
  
  // below, we extract the item parameters from the item parameter matrix
  // and store them in their own specific vector
  // This is not necessary, but just for convenience, 
  // so the model block below is less cluttered and more readable
  
  for(i in 1:I){
   alpha[i] = 1/exp(item[i][1]);   # variance of the lognormal distribution
   beta[i]  = item[i][2];
   b_star[i]= item[i][3];
  }
  
  vector[I] b = b_star;
   b[1] = -sum(b_star[2:I]);    # we fix the sum of b parameters to zero by 
                                # realigning the item difficulty for the first item
                                # as the sum of difficulty values for the remaining items
  
  vector[J] taut;               # true latent speed parameter
  vector[J] tauc;               # cheating latent speed parameter
  vector[J] thetat;             # true latent trait parameter
  vector[J] thetac;             # cheating latent trait parameter
  
  // below, we extract the person parameters from the person parameter matrix
  // also calculate the theta_c and tau_c
  // and store them all in their own specific vector
  // This is not necessary, but just for convenience, 
  // so the model block below is less cluttered and more readable
  
  for(j in 1:J){
    taut[j] = person[j][1];                
    tauc[j] = person[j][1]+delta[j][1];
    thetat[j] = person[j][2];
    thetac[j] = person[j][2]+delta[j][2];
  }
}

model{

  # Prior distributions

  sigma_P   ~ exponential(1);
  mu_thetat ~ normal(0,1);
  Lomega_P  ~ lkj_corr_cholesky(1);
    
  sigma_I   ~ exponential(1);
  mu_beta   ~ normal(4,1);
  mu_alpha  ~ lognormal(0,0.5);
  mu_b      ~ normal(0,1);
  Lomega_I  ~ lkj_corr_cholesky(1);
   
  sigma_D        ~ exponential(1);
  mu_delta_tau   ~ normal(0,1);
  mu_delta_theta ~ normal(0,1);
  Lomega_D       ~ lkj_corr_cholesky(1);
 
  person     ~  multi_normal_cholesky(mu_P,diag_pre_multiply(sigma_P, Lomega_P));
  item       ~  multi_normal_cholesky(mu_I,diag_pre_multiply(sigma_I, Lomega_I));
  delta      ~  multi_normal_cholesky(mu_D,diag_pre_multiply(sigma_D, Lomega_D));

  pC ~ beta(1,1);
  pH ~ beta(1,1);
  
// Joint log_density of response time and response accuracy  

  vector[I] log1m_pC = log1m(pC) ;
  vector[I] log_pC   = log(pC) ;
  vector[J] log1m_pH = log1m(pH) ;
  vector[J] log_pH   = log(pH) ;
  
  for (i in 1:n_obs) {
    real p_taut = beta[i_loc[i]] - taut[p_loc[i]];   
    real p_tauc = beta[i_loc[i]] - tauc[p_loc[i]];  
    real norm_p_t    = normal_lpdf(RT[i] | p_taut, alpha[i_loc[i]]);  
    real norm_p_c    = normal_lpdf(RT[i] | p_tauc, alpha[i_loc[i]]);  
      
    real p_thetat = thetat[p_loc[i]] - b[i_loc[i]];  
    real p_thetac = thetac[p_loc[i]] - b[i_loc[i]];
    real bern_p_t      = bernoulli_logit_lpmf(Y[i] | p_thetat); 
    real bern_p_c      = bernoulli_logit_lpmf(Y[i] | p_thetac); 
    
    vector[4] log_prob;
    log_prob[1] = log1m_pC[i_loc[i]] + log1m_pH[p_loc[i]] + norm_p_t + bern_p_t;
    log_prob[2] = log1m_pC[i_loc[i]] + log_pH[p_loc[i]]   + norm_p_t + bern_p_t;
    log_prob[3] = log_pC[i_loc[i]]   + log1m_pH[p_loc[i]] + norm_p_t + bern_p_t;
    log_prob[4] = log_pC[i_loc[i]]   + log_pH[p_loc[i]]   + norm_p_c + bern_p_c;
    target += log_sum_exp(log_prob);
  
  }
}


