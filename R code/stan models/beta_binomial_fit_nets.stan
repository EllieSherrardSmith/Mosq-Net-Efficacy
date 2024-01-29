// bernoulli_logistic transformed data function
data {
  
  int<lower=1> N;                  // rows of data
  
  int<lower=0> n_t[N];             // Total number of mosquitoes entering net huts
  int<lower=0> d_t[N];             // Number mosquites dead net hut

  vector<lower=0>[N] x;       // predictor (mortality in LLIN huts)
  
  int<lower = 1> R;                      //number of resistance points for plotting
  vector<lower = 0>[R] tau;              //time series for plotting
  int<lower = 1> N_M;                    //simulated individuals
  
}

parameters {
  //Consider death. This is the proportion of mosquitoes dying (d_t) in treated huts (n_t)
  real alpha1;
  real alpha2;
  
  real a1;
  real a2;
  
  //  real<lower=0,upper=10> sigma;
  real<lower = 0> k;                     //overdispersion parameter
}

transformed parameters {
  //mean mortality
  vector[N] m_t = (1.0 ./ (1.0 + exp(-(a1 + a2 * x))));
}

model {
  real sp[N];
  alpha1 ~ normal(0,100);
  alpha2 ~ normal(0,100);

  a1 ~ normal(0, 100);  // tried lognormal priors (1,1) and (0,1)
  a2 ~ normal(0, 100);

  //  study_a ~ normal(0,sigma);
  k ~ pareto(1.0, 1.5);
  
  for (n in 1:N) {
    sp[n] = alpha1  + alpha2 * x[n];
  }
  
  d_t ~ binomial_logit(n_t, sp);
  
  //binomial likelihood
  d_t[] ~ binomial(n_t, m_t);
  //beta-binomial likelihood
  d_t[] ~ beta_binomial(n_t, k * m_t, k * (1.0 - m_t));

}

generated quantities{
  real sp_ppc[100];

    for(xx in 1:100){
      sp_ppc[xx] = binomial_rng(100, inv_logit(alpha1 + alpha2 * xx)) / 100.0;
  }
  
  //predictions over fine-scaled time series
  vector[R] p = (1.0 ./ (1.0 + exp(-(a1 + a2 * tau)))); 
  int<lower = 0> M[R]; //number of individuals at each fine-scaled time point
  M = rep_array(N_M, R);
  int<lower = 0> rho[R]; //resistance counts for new fine-scaled time series
  rho = beta_binomial_rng(M, k * p, k * (1.0 - p));
}
