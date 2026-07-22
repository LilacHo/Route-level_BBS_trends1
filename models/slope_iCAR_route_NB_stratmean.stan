// iCAR route-level trend model — STRATUM-MEAN variant
//
// Identical to slope_iCAR_route_NB_New.stan (route-level slope model with
// spatial iCAR prior on alpha/beta, no random year effects) EXCEPT: the single
// global-mean scalars ALPHA/BETA are replaced by stratum-indexed vectors
// ALPHA_strat[nstrata]/BETA_strat[nstrata], one value per BCR/stratum, each
// route mapped to its stratum via route_strat[nroutes].
//
// The stratum means are themselves partially pooled toward a grand mean
// (MU_ALPHA/MU_BETA) with an estimated among-stratum sd (sd_ALPHA_strat/
// sd_BETA_strat). This nests the original global-mean model as the special
// case sd_ALPHA_strat, sd_BETA_strat -> 0 (all strata collapse to the grand
// mean), so the posteriors for sd_ALPHA_strat/sd_BETA_strat are themselves a
// direct test of whether a stratum mean is distinguishable from one global
// mean: if they concentrate near zero, using per-stratum means makes little
// difference; if they are clearly greater than zero, strata have materially
// different average intercepts/slopes.
//
// route-level alpha_r/beta_r still get their own iCAR spatial deviation on
// top of their stratum mean, exactly as before.

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

  // stratum (BCR) membership, one entry per route
  int<lower=1> nstrata;
  array [nroutes] int<lower=1, upper=nstrata> route_strat;

}


parameters {


  sum_to_zero_vector[nroutes] beta_raw_space; // spatial slope term with improved sum to zero constraint
  //vector[nroutes] beta_raw_rand;

  sum_to_zero_vector[nroutes] alpha_raw; // spatial abundance term with improved sum to zero constraint

  // stratum-level means (replace the single global ALPHA/BETA scalars)
  vector[nstrata] BETA_strat;
  vector[nstrata] ALPHA_strat;
  real MU_BETA;    // grand mean slope across strata
  real MU_ALPHA;   // grand mean intercept across strata
  real<lower=0> sd_BETA_strat;   // among-stratum sd of mean slope
  real<lower=0> sd_ALPHA_strat;  // among-stratum sd of mean intercept

  real eta; //first-year intercept

  sum_to_zero_vector[nobservers] obs_raw; //observer effects improved sum to zero constraint

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

   beta = beta_space + BETA_strat[route_strat];
   alpha = (sdalpha*alpha_raw) + ALPHA_strat[route_strat];
 //  noise = sdnoise*noise_raw;
   obs = sdobs*obs_raw;

  for(i in 1:ncounts){
    E[i] =  beta[route[i]] * (year[i]-fixedyear) + alpha[route[i]] + obs[observer[i]] + eta*firstyr[i];
  }


  // beta_raw_rand ~ normal(0,1);//random slope effects
  // sum(beta_raw_rand) ~ normal(0,0.001*nroutes);


  sdnoise ~ student_t(3,0,1); //prior on scale of inverse squared dispersion of NBinomial distribution phi = 1/sqrt(sdnoise)

  phi = 1/sqrt(sdnoise); //as recommended to avoid prior that places most prior mass at very high overdispersion by https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations

  sdobs ~ normal(0,0.3); //prior on sd of observer effects

  obs_raw ~ std_normal();//observer effects
  sum(obs_raw) ~ normal(0,0.001*nobservers);

  count ~ neg_binomial_2_log(E,phi); //vectorized count likelihood with log-transformation

  // stratum means, partially pooled toward a grand mean (nests the
  // global-mean model as sd_BETA_strat, sd_ALPHA_strat -> 0)
  MU_BETA ~ normal(0,0.1);// prior on grand mean slope (same scale as original BETA prior)
  MU_ALPHA ~ std_normal();// prior on grand mean intercept (same scale as original ALPHA prior)
  sd_BETA_strat ~ normal(0,0.1);// weakly-informative prior on among-stratum sd of slope means
  sd_ALPHA_strat ~ normal(0,1);// weakly-informative prior on among-stratum sd of intercept means
  BETA_strat ~ normal(MU_BETA, sd_BETA_strat);
  ALPHA_strat ~ normal(MU_ALPHA, sd_ALPHA_strat);

  eta ~ std_normal();// prior on first-year observer effect


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

// intercepts and slopes
   beta_space = (sdbeta_space*beta_raw_space);
   //beta_rand = (sdbeta_rand*beta_raw_rand);

   beta = beta_space + BETA_strat[route_strat];
   alpha = (sdalpha*alpha_raw) + ALPHA_strat[route_strat];


 }
