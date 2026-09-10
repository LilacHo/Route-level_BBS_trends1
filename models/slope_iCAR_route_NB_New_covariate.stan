// iCAR route-level slope model + a single land-cover covariate
// Single-covariate variant of slope_iCAR_route_NB_New_2covariates.stan, for
// when only one route-year covariate (proportion, range [0,1]) is available
// or wanted in the log-linear predictor:
//
//   log(lambda[r,t]) = alpha[r] + beta[r]*t + gamma1*Covariate[r,t]
//
// All other structure (iCAR spatial priors on alpha/beta, observer effects,
// first-year effect, NB2 overdispersion) is unchanged from
// slope_iCAR_route_NB_New.stan. "covariate" is generic on purpose (grassland,
// developed, or any other [0,1] route-year proportion can be passed in as
// this single column).


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

  // Route-year covariate, proportion of route in the given cover type
  // (range [0,1]), aligned 1:1 with count/year/route/etc. above.
  array [ncounts] real<lower=0, upper=1> covariate; // e.g. proportion grassland or developed cover

  int<lower=1> fixedyear; // centering value for years

 // spatial neighbourhood information
  int<lower=1> N_edges;
  array [N_edges] int<lower=1, upper=nroutes> node1;  // node1[i] adjacent to node2[i]
  array [N_edges] int<lower=1, upper=nroutes> node2;  // and node1[i] < node2[i]



}


parameters {


  sum_to_zero_vector[nroutes] beta_raw_space; // spatial slope term with improved sum to zero constraint
  //vector[nroutes] beta_raw_rand;
  real BETA;

  sum_to_zero_vector[nroutes] alpha_raw; // spatial abundance term with improved sum to zero constraint
  real ALPHA;

  real eta; //first-year intercept

  sum_to_zero_vector[nobservers] obs_raw; //observer effects improved sum to zero constraint

  real gamma1; // effect of the covariate proportion on log(lambda)

  real<lower=0> sdnoise;    // inverse of sd of over-dispersion
 //real<lower=1> nu;  //optional heavy-tail df for t-distribution
  real<lower=0> sdobs;    // sd of observer effects
  real<lower=0> sdbeta_space;    // sd of slopes
// real<lower=0> sdbeta_rand;    // sd of slopes
  real<lower=0> sdalpha;    // sd of intercepts


}


model {


  vector[ncounts] E;           // log_scale additive likelihood
   //vector[nroutes] beta_rand;
  vector[nroutes] beta_space;
  vector[nroutes] beta;
  vector[nroutes] alpha;
  vector[nobservers] obs;
  real phi;

// covariate effect on intercepts and slopes
   beta_space = (sdbeta_space*beta_raw_space);
   //beta_rand = (sdbeta_rand*beta_raw_rand);

   beta = beta_space + BETA;
   alpha = (sdalpha*alpha_raw) + ALPHA;
 //  noise = sdnoise*noise_raw;
   obs = sdobs*obs_raw;

  for(i in 1:ncounts){
    E[i] =  beta[route[i]] * (year[i]-fixedyear) + alpha[route[i]] + obs[observer[i]] + eta*firstyr[i]
            + gamma1*covariate[i];
  }


  // beta_raw_rand ~ normal(0,1);//random slope effects
  // sum(beta_raw_rand) ~ normal(0,0.001*nroutes);


  sdnoise ~ student_t(3,0,1); //prior on scale of inverse squared dispersion of NBinomial distribution phi = 1/sqrt(sdnoise)

  phi = 1/sqrt(sdnoise); //as recommended to avoid prior that places most prior mass at very high overdispersion by https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations

  sdobs ~ normal(0,0.3); //prior on sd of observer effects

  obs_raw ~ std_normal();//observer effects
  sum(obs_raw) ~ normal(0,0.001*nobservers);

  count ~ neg_binomial_2_log(E,phi); //vectorized count likelihood with log-transformation

  BETA ~ normal(0,0.1);// prior on fixed effect mean slope
  ALPHA ~ std_normal();// prior on fixed effect mean intercept
  eta ~ std_normal();// prior on first-year observer effect

  gamma1 ~ normal(0,1);// prior on the covariate effect


  //spatial iCAR intercepts and slopes by strata
 // sdalpha ~ gamma(2,2); // alternate zero-avoiding prior on sd of intercept variation with similar limits to normal
 sdalpha ~ normal(0,sd_alpha_prior); //prior on sd of intercept variation

  sdbeta_space ~ gamma(3,30);//zero-avoiding prior on sd of slope spatial variation w mean = 0.1 and 99% < 0.3
  //sdbeta_space ~ normal(0,0.1);//alternative prior

   target += -0.5 * dot_self(beta_raw_space[node1] - beta_raw_space[node2]); // ICAR prior
   target += -0.5 * dot_self(alpha_raw[node1] - alpha_raw[node2]); // ICAR prior

  // beta_raw_space ~ icar_normal(nroutes, node1, node2);
  // alpha_raw ~ icar_normal(nroutes, node1, node2);


}

 generated quantities {

   //vector[nroutes] beta_rand;
  vector[nroutes] beta_space;
  vector[nroutes] beta;
  vector[nroutes] alpha;
  real gamma1_out;

// intercepts and slopes
   beta_space = (sdbeta_space*beta_raw_space);
   //beta_rand = (sdbeta_rand*beta_raw_rand);

   beta = beta_space + BETA;
    alpha = (sdalpha*alpha_raw) + ALPHA;

// covariate effect, passed through so it's grouped with the other
// derived/reported quantities in the output summary
   gamma1_out = gamma1;

 }
