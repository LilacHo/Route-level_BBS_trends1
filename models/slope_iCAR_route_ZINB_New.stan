// iCAR route-level trend model — zero-inflated negative binomial (ZINB)
// Adds a single, GLOBAL zero-inflation probability theta on top of
// models/slope_iCAR_route_NB_New.stan: with probability theta an
// observation is a structural (non-sampling) zero, and with probability
// (1-theta) it is drawn from the same NB2(mu = exp(E), phi) process as the
// plain NB model. This is the iCAR analogue of
// models/slope_habitat_route_ZINB_New.stan (no habitat covariates here).
//
// Uses the same sum_to_zero_vector / target += ICAR reparameterization as
// slope_iCAR_route_NB_New.stan (see that file for notes).
//
// Paired CV model: models/slope_iCAR_route_ZINB_cv.stan
// Used by: 3_compare_NB_ZINB.R (NB vs ZINB model comparison)

data {
  int<lower=1> nroutes;
  int<lower=1> ncounts;
  int<lower=1> nyears;
  int<lower=0> nobservers;
  
  real<lower=0> sd_alpha_prior;
 
  array [ncounts] int<lower=0> count;              // count observations
  array [ncounts] int<lower=1> year; // year index
  array [ncounts] int<lower=1> route; // route index
  array [ncounts] int<lower=0> firstyr; // first year index
  array [ncounts] int<lower=1> observer;              // observer indicators

  int<lower=1> fixedyear; // centering value for years
 
 // spatial neighbourhood information
  int<lower=1> N_edges;
  array [N_edges] int<lower=1, upper=nroutes> node1;  // node1[i] adjacent to node2[i]
  array [N_edges] int<lower=1, upper=nroutes> node2;  // and node1[i] < node2[i]



}


parameters {

  
  sum_to_zero_vector[nroutes] beta_raw_space; // spatial slope term with improved sum to zero constraint
  real BETA; 

  sum_to_zero_vector[nroutes] alpha_raw; // spatial abundance term with improved sum to zero constraint
  real ALPHA; 

  real eta; //first-year intercept
  
  sum_to_zero_vector[nobservers] obs_raw; //observer effects improved sum to zero constraint

  real<lower=0> sdnoise;    // inverse of sd of over-dispersion
  real<lower=0> sdobs;    // sd of observer effects
  real<lower=0> sdbeta_space;    // sd of slopes 
  real<lower=0> sdalpha;    // sd of intercepts

  real<lower=0,upper=1> theta; // global zero-inflation (structural-zero) probability

}


model {

  vector[ncounts] E;           // log_scale additive likelihood
  vector[nroutes] beta_space;
  vector[nroutes] beta;
  vector[nroutes] alpha;
  vector[nobservers] obs;
  real phi;

// covariate effect on intercepts and slopes
   beta_space = (sdbeta_space*beta_raw_space);
   
   beta = beta_space + BETA;
   alpha = (sdalpha*alpha_raw) + ALPHA;
   obs = sdobs*obs_raw;

  for(i in 1:ncounts){
    E[i] =  beta[route[i]] * (year[i]-fixedyear) + alpha[route[i]] + obs[observer[i]] + eta*firstyr[i];
  }
  

  sdnoise ~ student_t(3,0,1); //prior on scale of inverse squared dispersion of NBinomial distribution phi = 1/sqrt(sdnoise)

  phi = 1/sqrt(sdnoise); //as recommended to avoid prior that places most prior mass at very high overdispersion by https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations

  sdobs ~ normal(0,0.3); //prior on sd of observer effects
 
  obs_raw ~ std_normal();//observer effects

  theta ~ beta(1,1); // flat prior on global zero-inflation probability

  // zero-inflated NB2 count likelihood (mixture: cannot vectorize because
  // the y==0 and y>0 cases use different mixture components)
  for(i in 1:ncounts){
    if(count[i] == 0){
      target += log_sum_exp(
        bernoulli_lpmf(1 | theta),
        bernoulli_lpmf(0 | theta) + neg_binomial_2_log_lpmf(0 | E[i], phi)
      );
    } else {
      target += bernoulli_lpmf(0 | theta) + neg_binomial_2_log_lpmf(count[i] | E[i], phi);
    }
  }
 
  BETA ~ normal(0,0.1);// prior on fixed effect mean slope
  ALPHA ~ std_normal();// prior on fixed effect mean intercept
  eta ~ std_normal();// prior on first-year observer effect
 
  
  //spatial iCAR intercepts and slopes by strata
 sdalpha ~ normal(0,sd_alpha_prior); //prior on sd of intercept variation

  sdbeta_space ~ gamma(3,30);//zero-avoiding prior on sd of slope spatial variation w mean = 0.1 and 99% < 0.3

   target += -0.5 * dot_self(beta_raw_space[node1] - beta_raw_space[node2]); // ICAR prior
   target += -0.5 * dot_self(alpha_raw[node1] - alpha_raw[node2]); // ICAR prior


}

 generated quantities {

  vector[nroutes] beta_space;
  vector[nroutes] beta;
  vector[nroutes] alpha;

// intercepts and slopes
   beta_space = (sdbeta_space*beta_raw_space);
   
   beta = beta_space + BETA;
    alpha = (sdalpha*alpha_raw) + ALPHA;
  

 }
