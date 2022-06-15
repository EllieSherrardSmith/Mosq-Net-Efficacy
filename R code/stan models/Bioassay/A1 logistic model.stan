functions{
  real logistic(real x, real a, real c){
    return(1 / (1 + exp(-(x - c) * a)));
  }
}

data{
  int S;
  int N_b[S];
  int N_h[S];
  int X_b[S];
  int X_h[S];
  
  // random effects
  int nsite;
  int site[S];
  
  int S_test;
  vector[S_test] theta_b_test;
  
  // hut
  //int is_ifa[S];
}

parameters{
  real<lower=0, upper=1> theta_b[S];
  real theta_h_logit[S];
  real<lower=0> a;
  real c;
  real<lower=0> sigma;
  vector[nsite] r_b;
  vector[nsite] r_h;
  real<lower=0> sigma_b;
  real<lower=0> sigma_h;
}

transformed parameters{
 real theta_h_bar_logit[S];
 for(i in 1:S)
  theta_h_bar_logit[i] = logit(logistic(theta_b[i], a, c));
}

model{
  // likelihood for bioassay
  X_b ~ binomial_logit(N_b, to_vector(logit(theta_b)) + r_b[site]);
  
  // likelihood for huts
  X_h ~ binomial_logit(N_h, to_vector(theta_h_logit) + r_h[site]);
  
  // priors for shape of logistic
  a ~ cauchy(0, 1);
  c ~ normal(0, 2);
  
  // prior for hut survival
  sigma ~ normal(0, 2);
  for(i in 1:S)
    theta_h_logit[i] ~ normal(theta_h_bar_logit[i], sigma);
    
  // priors for random effects
  r_b ~ normal(0, sigma_b);
  r_h ~ normal(0, sigma_h);
  sigma_b ~ normal(0, 1);
  sigma_h ~ normal(0, 1);
}

generated quantities{
  vector[S_test] theta_h_test;
  vector[S] logLikelihood;
  for(i in 1:S_test)
    theta_h_test[i] = logistic(theta_b_test[i], a, c);
  for(i in 1:S)
    logLikelihood[i] = binomial_logit_lpmf(X_b[i]|N_b[i], logit(theta_b[i]) + r_b[site[i]]) + binomial_logit_lpmf(X_h[i]|N_h[i], theta_h_logit[i] + r_h[site[i]]);
    //logLikelihood[i] = binomial_logit_lpmf(X_b[i]|N_b[i], logit(theta_b[i])) + binomial_logit_lpmf(X_h[i]|N_h[i], theta_h_logit[i]);
}
