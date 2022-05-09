################################
##
## Data: from Mosha et al. 2022 The Lancet
##       Appendix Table S1


# Net use through the trial to estimate adherence to using the nets

################################
##
## Net retention
##
##################################

time_obs = c(3/12,12/12,18/12,24/12)
pyr_net_use = c(0.768, 0.616, 0.522, 0.495)
pyrrole_net_use = c(0.684,0.653,0.515,0.464)
pyrpbo_net_use = c(0.737,0.59,0.407,0.296)


y_standard_net_usage = log(pyr_net_use)
y_pbo_net_usage = log(pyrpbo_net_use)
y_pp_net_usage = log(pyrrole_net_use)

time_m = seq(0,3,0.01)
# 

stanmodelcode <- "
data {
int<lower=0> N;
int<lower=0> N2; //the size of the new_X matrix
vector[N] y;
vector[N] x;
vector[N2] New_x;

}
parameters {
real beta0;
real beta1;
real sigma;
}
transformed parameters {
vector[N] m;
m = beta0 + beta1 * x;
} 
model {
// priors
beta0 ~ cauchy(0, 10); 
beta1 ~ cauchy(0, 10); 

// likelihood
y ~ normal(m, sigma);   
}
generated quantities {
vector[N2] y_pred;
y_pred = beta0 + beta1 * New_x; //the y values predicted by the model
}
"
stanDso <- stan_model(model_code = stanmodelcode) 

dat_standard <- list(N = length(time_obs), 
                     N2 = length(seq(0,3,0.01)),
                     y = y_standard_net_usage, 
                     x = time_obs,
                     New_x = seq(0,3,0.01)); 
fit <- sampling(stanDso, data = dat_standard, iter = 10000, warmup=2000) 
fit

#plotting the posterior distribution for the parameters
post_beta<-As.mcmc.list(fit,pars="beta0")
plot(post_beta)

## gradient is fit to the data for alpha
## standard_net_usage ~ exp(-alpha*time_obs)
b0 <- rstan::extract(fit, 'beta0')
b0<- unlist(b0, use.names=FALSE)
b1 <- rstan::extract(fit, 'beta1')
b1<- unlist(b1, use.names=FALSE)

y_predicted_stn_exp = exp(mean(b0)) * exp(mean(b1)*time_m)
y_predicted_stn_exp_min = exp(quantile(b0,0.45)) * exp(quantile(b1,0.25)*time_m)
y_predicted_stn_exp_max = exp(quantile(b0,0.55)) * exp(quantile(b1,0.75)*time_m)

par(mfrow=c(1,1))
plot(y_predicted_stn_exp ~ time_m,type="l",
     ylab = "Proportion of community using nets (%)",
     xlab = "Time in months",
     ylim=c(0,1),xlim=c(0,3),xaxt="n")
axis(1, at=c(0,1,2,3),labels = c("0","12","24","36"))

polygon(c(time_m,rev(time_m)),
        c(y_predicted_stn_exp_min,rev(y_predicted_stn_exp_max)),border=NA,col = adegenet::transp("grey",0.6))
points(pyr_net_use~time_obs,pch=19,col="grey")



##########################
##
## pyrethroid-PBO nets
dat_pbo <- list(N = length(time_obs), 
                N2 = length(seq(0,3,0.01)),
                y = y_pbo_net_usage, 
                x = time_obs,
                New_x = seq(0,3,0.01)); 
fit2 <- sampling(stanDso, data = dat_pbo, iter = 10000, warmup=2000) 
fit2
#
# #plotting the posterior distribution for the parameters
# post_beta<-As.mcmc.list(fit,pars="beta0")
# plot(post_beta)

## gradient is fit to the data for alpha
## standard_net_usage ~ exp(-alpha*time_obs)
b0_pbo <- rstan::extract(fit2, 'beta0')
b0_pbo<- unlist(b0_pbo, use.names=FALSE)
b1_pbo <- rstan::extract(fit2, 'beta1')
b1_pbo<- unlist(b1_pbo, use.names=FALSE)

y_predicted_pbo_exp = exp(mean(b0_pbo)) * exp(mean(b1_pbo)*time_m)
y_predicted_pbo_exp_min = exp(quantile(b0_pbo,0.45)) * exp(quantile(b1_pbo,0.25)*time_m)
y_predicted_pbo_exp_max = exp(quantile(b0_pbo,0.55)) * exp(quantile(b1_pbo,0.75)*time_m)

polygon(c(time_m,rev(time_m)),
        c(y_predicted_pbo_exp_min,rev(y_predicted_pbo_exp_max)),border=NA,
        col = adegenet::transp("darkorange",0.6))
lines(y_predicted_pbo_exp ~ time_m,col="darkorange",lwd=2)
points(pyrpbo_net_use~time_obs,pch=19,col="darkorange")


##########################
##
## pyrethroid-pyrrole nets
dat_pp <- list(N = length(time_obs), 
               N2 = length(seq(0,3,0.01)),
               y = y_pp_net_usage, 
               x = time_obs,
               New_x = seq(0,3,0.01)); 
fit3 <- sampling(stanDso, data = dat_pp, iter = 10000, warmup=2000) 
# fit3
#
# #plotting the posterior distribution for the parameters
# post_beta<-As.mcmc.list(fit,pars="beta0")
# plot(post_beta)

## gradient is fit to the data for alpha
## standard_net_usage ~ exp(-alpha*time_obs)
b0_pp <- rstan::extract(fit3, 'beta0')
b0_pp<- unlist(b0_pp, use.names=FALSE)
b1_pp <- rstan::extract(fit3, 'beta1')
b1_pp<- unlist(b1_pp, use.names=FALSE)

y_predicted_pp_exp = exp(mean(b0_pp)) * exp(mean(b1_pp)*time_m)
y_predicted_pp_exp_min = exp(quantile(b0_pp,0.45)) * exp(quantile(b1_pp,0.25)*time_m)
y_predicted_pp_exp_max = exp(quantile(b0_pp,0.55)) * exp(quantile(b1_pp,0.75)*time_m)

polygon(c(time_m,rev(time_m)),
        c(y_predicted_pp_exp_min,rev(y_predicted_pp_exp_max)),border=NA,
        col = adegenet::transp("darkgreen",0.6))
points(pyrrole_net_use~time_obs,pch=19,col="darkgreen")
lines(y_predicted_pp_exp ~ time_m,col="darkgreen",lwd=2)

legend("topright",title="Net type",legend=c("Pyrethoid-only nets",
                                           "Pyrethroid-PBO nets",
                                           "Pyrethroid-pyrrole nets"),
       col=c("grey","orange","darkgreen"),pch=19,bty="n",lty=1)
######################
##
## Translate these to half-life in years
##
######################

parms_usage = data.frame(itn_leave_dur_pyr = sample(b1[b1 >= quantile(b1,0.25) & b1 <= quantile(b1,0.75)],1000,replace=FALSE),
                         itn_leave_dur_pbo = sample(b1_pbo[b1_pbo >= quantile(b1_pbo,0.25) & b1_pbo <= quantile(b1_pbo,0.75)],1000,replace=FALSE),
                         itn_leave_dur_pp = sample(b1_pp[b1_pp >= quantile(b1_pp,0.25) & b1_pp <= quantile(b1_pp,0.75)],1000,replace=FALSE)) ##

parms_usage$itn_leave_dur_pyr = -1/parms_usage$itn_leave_dur_pyr
parms_usage$itn_leave_dur_pbo = -1/parms_usage$itn_leave_dur_pbo
parms_usage$itn_leave_dur_pp = -1/parms_usage$itn_leave_dur_pp

quantile(parms_usage$itn_leave_dur_pyr,c(0.5,0.025,0.975))
quantile(parms_usage$itn_leave_dur_pbo,c(0.5,0.025,0.975))
quantile(parms_usage$itn_leave_dur_pp,c(0.5,0.025,0.975))



##############################
##
## Working out top up to stay at 80% if new nets each 6 months
##
##############################
## pyrethroid nets
(0.88+0.854+0.786+0.826)/4
1 - y_predicted_stn_exp[50]/0.8365
y_predicted_stn_exp_top_up1 = c(rep(0,50),y_predicted_stn_exp)
# lines(y_predicted_stn_exp_top_up1[1:301]~time_m)
1 - y_predicted_stn_exp_top_up1[100]/0.8365
# i.e. 0.119

## pyr-PBO nets
(0.858+0.868+0.744+0.765)/4
1 - y_predicted_pbo_exp[50]/0.80875
y_predicted_pbo_exp_top_up1 = c(rep(0,50),y_predicted_pbo_exp)
# lines(y_predicted_pbo_exp_top_up1[1:301]~time_m)
1 - y_predicted_pbo_exp_top_up1[100]/0.80875
# i.e. 0.124

## pyr-pyrrole nets
(0.844+0.849+0.748+0.77)/4
1 - y_predicted_pp_exp[50]/0.80275
y_predicted_pp_exp_top_up1 = c(rep(0,50),y_predicted_pp_exp)
# lines(y_predicted_pp_exp_top_up1[1:301]~time_m)
1 - y_predicted_pp_exp_top_up1[100]/0.80275
# i.e. 0.142

#############################
##
## How to stay at 56% prior to this
##
#############################
y_predicted_stn_exp[1] - 0.56 
historic_cover_estimates = y_predicted_stn_exp - 0.2379052
lines(historic_cover_estimates ~ time_m)
historic_cover_estimates2 = c(rep(0,50),historic_cover_estimates)
lines(historic_cover_estimates2[1:301]~time_m)
1 - historic_cover_estimates2[100]/0.56


###############################
##
## Mosquito resistance
##
###############################

fun = c(0.943,0.928,0.954)
gam = 1-fun

bio_mort_fun = c(0.325,0.422,0.325)
bio_mort_gam = c(0.333,0.460,0.333)

weighted_means = round(fun*bio_mort_fun + gam*bio_mort_gam,2)
resistance = 1 - weighted_means

############################################
##
## Simulating the Misungwe trial 
##
############################################


# Set up
library(malariasimulation)
library(malariaEquilibrium)
library(reshape2)
library(ggplot2)
library(cali)

set.seed(24345)
## Read in seasonality data
## ASA, Imo
## MORO, Kwara
## Using NNP pilots in Burkina as example sites to contrast net impacts.


dat_res_pyr = read.csv("C:/Users/esherrar/Documents/Rprojects/Mosq-Net-Efficacy/parameters/pyrethroid_only_nets.csv",header=TRUE) 
dat_res_pbo = read.csv("C:/Users/esherrar/Documents/Rprojects/Mosq-Net-Efficacy/parameters/pyrethroid_pbo_nets.csv",header=TRUE) 
dat_res_pp = read.csv("C:/Users/esherrar/Documents/Rprojects/Mosq-Net-Efficacy/parameters/pyrethroid_pyrrole_nets.csv",header=TRUE) 

dat_res_pyr[which(dat_res_pyr$resistance == data1$RESISTANCE[row_drawn]),c(3,2,4)]
dat_res_pbo[which(dat_res_pbo$resistance == data1$RESISTANCE[row_drawn]),c(3,2,4)]
dat_res_pp[which(dat_res_pp$resistance == data1$RESISTANCE[row_drawn]),c(3,2,4)]

dat_res_pyr[which(dat_res_pyr$resistance == data1$RESISTANCE[row_drawn]),c(6,5,7)]
dat_res_pbo[which(dat_res_pbo$resistance == data1$RESISTANCE[row_drawn]),c(6,5,7)]
dat_res_pp[which(dat_res_pp$resistance == data1$RESISTANCE[row_drawn]),c(6,5,7)]

dat_res_pyr[which(dat_res_pyr$resistance == data1$RESISTANCE[row_drawn]),c(9,8,10)]
dat_res_pbo[which(dat_res_pbo$resistance == data1$RESISTANCE[row_drawn]),c(9,8,10)]
dat_res_pp[which(dat_res_pp$resistance == data1$RESISTANCE[row_drawn]),c(9,8,10)]


sites = read.csv("Figures/data RCT Misungwi/site_file.csv",header=TRUE)
data1 = read.csv("Figures/data RCT Misungwi/parameter_estimates_Mosha2022.csv",header=TRUE)


##dat_res we create in 0_UNCERTAINTY...
##ress calls a cell from this data frame
## EIR_L we use to calibrate the model
## dide_no calls what is required from sites file

malsim_smc_cali_f = function(EIR_L,row_drawn,top_up){
  
  
  year <- 365
  month <- 30
  sim_length <- 9 * year
  ## This is spanning Jan 2016 - Dec 2024
  ## received NNP nets 2020
  ## 
  
  ## and finally we wish to see the difference for these places running forward from 2020-2023 (1 or 3 years)
  human_population <- 10000
  starting_EIR <- EIR_L
  
  
  simparams <- get_parameters(
    list(
      human_population = human_population,
      # irs_correlation = 
      
      prevalence_rendering_min_ages = c(0.5, 0, 0, 5,  0, 2) * 365, ## Prev in 6 months to 14 years measured
      prevalence_rendering_max_ages = c(14,  5,10,10,100,10) * 365,
      
      clinical_incidence_rendering_min_ages = c(0.5) * 365, ## 6-months to 10 years clin_inc
      clinical_incidence_rendering_max_ages = c(10) * 365,
      
      model_seasonality = TRUE, ## Seasonality to match study site inputs [sites_13]
      g0 = sites$seasonal_a0[1],
      g = c(sites$seasonal_a1[1], sites$seasonal_a2[1], sites$seasonal_a3[1]),
      h = c(sites$seasonal_b1[1], sites$seasonal_b2[1], sites$seasonal_b3[1]),
      
      individual_mosquitoes = FALSE
      # bednets = TRUE,
      # drugs = TRUE,
      # smc = TRUE
    )
  )
  
  # set species
  fun_mint_params <- fun_params # fun
  gamb_mint_params <- gamb_params # gamb sl
  
  
  fun_mint_params['species'] <- "fun_misungwe"
  fun_mint_params['blood_meal_rates'] <- 1/3 # 1/duration of gonothropic cycle
  fun_mint_params['foraging_time'] <- 0.69 # time spent foraging
  fun_mint_params['Q0'] <- 0.94 # human blood index
  fun_mint_params['phi_bednets'] <- 0.9 # proportion biting in bed
  fun_mint_params['phi_indoors'] <- 0.98 # proportion biting indoors
  fun_mint_params['mum'] <- 0.112 # death rate or 1/life expectancy
  
  gamb_mint_params['species'] <- "gamb_misungwe"
  gamb_mint_params['blood_meal_rates'] <- 1/3 # 1/duration of gonothropic cycle
  gamb_mint_params['foraging_time'] <- 0.69   # time spent foraging
  gamb_mint_params['Q0'] <- 0.92              # human blood index
  gamb_mint_params['phi_bednets'] <- 0.81     # proportion biting in bed
  gamb_mint_params['phi_indoors'] <- 0.89     # proportion biting indoors
  gamb_mint_params['mum'] <- 0.132            # death rate or 1/life expectancy
  
  simparams <- set_species(simparams,
                           list(fun_mint_params,gamb_mint_params),
                           c(data1$prop_fun[1],1-data1$prop_fun[1]))                         
  
  # set treatment
  simparams <- set_drugs(simparams, list(AL_params, SP_AQ_params, DHA_PQP_params))
  simparams <- set_clinical_treatment(simparams, 
                                      drug=1,
                                      time=c(100),
                                      coverage=c(0.258))
  simparams <- set_clinical_treatment(simparams, 
                                      drug=2,
                                      time=c(100),
                                      coverage=c(0.233))
  simparams <- set_clinical_treatment(simparams, 
                                      drug=3,
                                      time=c(100),
                                      coverage=c(0))
  
  ## estimate equilibrium now
  simparams <- set_equilibrium(simparams, starting_EIR)
  
  # set mosquito nets
  bednetparams <- simparams
  
  ## as done
  bednet_events = data.frame(
    timestep = round(c(0,0.5,
                       1,1.5,2,2.5,3,3.5, ## baseline with top up of nets each 6 months (always pyrethroid only nets)
                       4,4.5,5,5.5,6,6.5,       ## trial nets introduced year 4, jan 26-28th, with top up each 6 months of pyr-only
                       7,7.5,8,8.5,9,9.5) * year + ## post trial years 
                       c(0,0,27,rep(0,5),      
                         27,rep(0,5),
                         27,rep(0,5)),0),
    name=c("2015","2015.5","2016","2016.5","2017","2017.5","2018","2018.5",
           "2019","2019.5","2020","2020.5","2021","2021.5",
           "2022","2022.5","2023","2023.5","2024","2024.5")
  )
  
  # each net will be changing the row from 1 to 3
  bednetparams_1 <- set_bednets(
    bednetparams,
    
    timesteps = bednet_events$timestep,
    
    coverages = c(0.12,0.12,data1$historic_ITN_use[row_drawn],rep(0.12,5),
                  data1$itn_use[row_drawn],rep(top_up,5), ## From Misungwe trial
                  
                  data1$itn_use[row_drawn],rep(top_up,5)),   ## planned for 2020 - ** Assuming the distribution coverage matched the RCT estimate
    
    retention = data1$itn_leave_dur[row_drawn] * year, ## from trial (try matching this for all as well...)
    
    ## each row needs to show the efficacy parameter across years (and cols are diff mosquito)
    dn0 = matrix(c(rep(data1$itn_kill_pyrethroid_ITN[row_drawn],8),
                   data1$itn_kill_TRIAL_NET[row_drawn],
                   rep(data1$itn_kill_pyrethroid_ITN[row_drawn],5),
                   data1$itn_kill_TRIAL_NET[row_drawn],
                   rep(data1$itn_kill_pyrethroid_ITN[row_drawn],5),
                   rep(data1$itn_kill_pyrethroid_ITN[row_drawn],8),
                   data1$itn_kill_TRIAL_NET[row_drawn],
                   rep(data1$itn_kill_pyrethroid_ITN[row_drawn],5),
                   data1$itn_kill_TRIAL_NET[row_drawn],
                   rep(data1$itn_kill_pyrethroid_ITN[row_drawn],5)),
                 nrow=20, ncol=2),
    rn =  matrix(c(rep(data1$itn_repel_pyrethroid_ITN[row_drawn],8),
                   data1$itn_repel_TRIAL_NET[row_drawn],
                   rep(data1$itn_repel_pyrethroid_ITN[row_drawn],5),
                   data1$itn_repel_TRIAL_NET[row_drawn],
                   rep(data1$itn_repel_pyrethroid_ITN[row_drawn],5),
                   rep(data1$itn_repel_pyrethroid_ITN[row_drawn],8),
                   data1$itn_repel_TRIAL_NET[row_drawn],
                   rep(data1$itn_repel_pyrethroid_ITN[row_drawn],5),
                   data1$itn_repel_TRIAL_NET[row_drawn],
                   rep(data1$itn_repel_pyrethroid_ITN[row_drawn],5)),
                 nrow=20, ncol=2),
    rnm = matrix(.24, nrow=20, ncol=2),
    gamman = c(rep(data1$itn_half_life_pyrethroid_ITN[row_drawn] * 365, 8),
               data1$itn_half_life_TRIAL_NET[row_drawn] * 365,
               rep(data1$itn_half_life_pyrethroid_ITN[row_drawn] * 365, 5),
               data1$itn_half_life_TRIAL_NET[row_drawn] * 365,
               rep(data1$itn_half_life_pyrethroid_ITN[row_drawn] * 365, 5))
  )
  
  correlationsb1 <- get_correlation_parameters(bednetparams_1)
  correlationsb1$inter_round_rho('bednets', 0)
  
  ## Run the simulations
  output1 <- run_simulation(sim_length, bednetparams_1,correlationsb1)
  
  output1$pv_182.5_5110 = output1$n_detect_182.5_5110/output1$n_182.5_5110
  output1$pv_0_1825 = output1$n_detect_0_1825/output1$n_0_1825
  output1$pv_1825_3650 = output1$n_detect_1825_3650/output1$n_1825_3650
  output1$pv_0_3650 = output1$n_detect_0_3650/output1$n_0_3650
  output1$pv_730_3650 = output1$n_detect_730_3650/output1$n_730_3650
  output1$pv_0_36500 = output1$n_detect_0_36500/output1$n_0_36500
  
  output1$clin_inc_182.5_3650 = output1$n_inc_clinical_182.5_3650/output1$n_182.5_3650
  
  # Define target, here two prevalence measures:
  target <- data1$prev_baseline[row_drawn]
  # Time points at which to match target
  target_tt <- c(3*365+31+28+31+30+31+30+31+31+30+15)
  
  summary_pfprev_6m_14y = function(x){
    prev_6m_14y <- x$n_detect_182.5_5110/x$n_182.5_5110
    return(prev_6m_14y)
  }
  
  set.seed(123)
  out_net1 <- calibrate(parameters = bednetparams_1,
                        target = target,
                        target_tt = target_tt,
                        summary_function = summary_pfprev_6m_14y,
                        tolerance = 0.02, 
                        interval = c(1, 250))##upper bound needs to be high enough so negative differences are not returned in uniroot
  
  return(out_net1$root)
  
  # return(data.frame(timestep = output1$timestep,
  #                   prev_05_14 = output1$pv_182.5_5110,
  #                   prev_0_5 = output1$pv_0_1825,
  #                   prev_5_10 =  output1$pv_1825_3650,
  #                   prev_0_10 =  output1$pv_0_3650,
  #                   prev_all_age = output1$pv_0_36500,
  #                   prev_2_10 =  output1$pv_730_3650,
  #                   clin_inc_6m_10y = output1$clin_inc_182.5_3650))
}

pyr_nets_cali = malsim_smc_cali_f(EIR_L = 150,row_drawn =1,top_up = 0.1577571)
pbo_nets_cali = malsim_smc_cali_f(EIR_L = 90, row_drawn =2,top_up = 0.1334457)
pp_nets_cali =  malsim_smc_cali_f(EIR_L = 130,row_drawn =3,top_up = 0.1451762)



malsim_smc_f = function(EIR_L,row_drawn,top_up){
  
  
  year <- 365
  month <- 30
  sim_length <- 9 * year
  ## This is spanning Jan 2016 - Dec 2024
  ## received NNP nets 2020
  ## 
  
  ## and finally we wish to see the difference for these places running forward from 2020-2023 (1 or 3 years)
  human_population <- 10000
  starting_EIR <- EIR_L
  
  
  simparams <- get_parameters(
    list(
      human_population = human_population,
      # irs_correlation = 
      
      prevalence_rendering_min_ages = c(0.5, 0, 0, 5,  0, 2) * 365, ## Prev in 6 months to 14 years measured
      prevalence_rendering_max_ages = c(14,  5,10,10,100,10) * 365,
      
      clinical_incidence_rendering_min_ages = c(0.5) * 365, ## 6-months to 10 years clin_inc
      clinical_incidence_rendering_max_ages = c(10) * 365,
      
      model_seasonality = TRUE, ## Seasonality to match study site inputs [sites_13]
      g0 = sites$seasonal_a0[1],
      g = c(sites$seasonal_a1[1], sites$seasonal_a2[1], sites$seasonal_a3[1]),
      h = c(sites$seasonal_b1[1], sites$seasonal_b2[1], sites$seasonal_b3[1]),
      
      individual_mosquitoes = FALSE
      # bednets = TRUE,
      # drugs = TRUE,
      # smc = TRUE
    )
  )
  
  # set species
  fun_mint_params <- fun_params # fun
  gamb_mint_params <- gamb_params # gamb sl
  
  
  fun_mint_params['species'] <- "fun_misungwe"
  fun_mint_params['blood_meal_rates'] <- 1/3 # 1/duration of gonothropic cycle
  fun_mint_params['foraging_time'] <- 0.69 # time spent foraging
  fun_mint_params['Q0'] <- 0.94 # human blood index
  fun_mint_params['phi_bednets'] <- 0.9 # proportion biting in bed
  fun_mint_params['phi_indoors'] <- 0.98 # proportion biting indoors
  fun_mint_params['mum'] <- 0.112 # death rate or 1/life expectancy
  
  gamb_mint_params['species'] <- "gamb_misungwe"
  gamb_mint_params['blood_meal_rates'] <- 1/3 # 1/duration of gonothropic cycle
  gamb_mint_params['foraging_time'] <- 0.69   # time spent foraging
  gamb_mint_params['Q0'] <- 0.92              # human blood index
  gamb_mint_params['phi_bednets'] <- 0.81     # proportion biting in bed
  gamb_mint_params['phi_indoors'] <- 0.89     # proportion biting indoors
  gamb_mint_params['mum'] <- 0.132            # death rate or 1/life expectancy
  
  simparams <- set_species(simparams,
                           list(fun_mint_params,gamb_mint_params),
                           c(data1$prop_fun[1],1-data1$prop_fun[1]))                         
  
  # set treatment
  simparams <- set_drugs(simparams, list(AL_params, SP_AQ_params, DHA_PQP_params))
  simparams <- set_clinical_treatment(simparams, 
                                      drug=1,
                                      time=c(100),
                                      coverage=c(0.258))
  simparams <- set_clinical_treatment(simparams, 
                                      drug=2,
                                      time=c(100),
                                      coverage=c(0.233))
  simparams <- set_clinical_treatment(simparams, 
                                      drug=3,
                                      time=c(100),
                                      coverage=c(0))
  
  ## estimate equilibrium now
  simparams <- set_equilibrium(simparams, starting_EIR)
  
  # set mosquito nets
  bednetparams <- simparams
  
  ## as done
  bednet_events = data.frame(
    timestep = round(c(0,0.5,
                       1,1.5,2,2.5,3,3.5, ## baseline with top up of nets each 6 months (always pyrethroid only nets)
                       4,4.5,5,5.5,6,6.5,       ## trial nets introduced year 4, jan 26-28th, with top up each 6 months of pyr-only
                       7,7.5,8,8.5,9,9.5) * year + ## post trial years 
                       c(0,0,27,rep(0,5),      
                         27,rep(0,5),
                         27,rep(0,5)),0),
    name=c("2015","2015.5","2016","2016.5","2017","2017.5","2018","2018.5",
           "2019","2019.5","2020","2020.5","2021","2021.5",
           "2022","2022.5","2023","2023.5","2024","2024.5")
  )
  
  # each net will be changing the row from 1 to 3
  bednetparams_1 <- set_bednets(
    bednetparams,
    
    timesteps = bednet_events$timestep,
    
    coverages = c(0.12,0.12,data1$historic_ITN_use[row_drawn],rep(0.12,5),
                  data1$itn_use[row_drawn],rep(top_up,5), ## From Misungwe trial
                  
                  data1$itn_use[row_drawn],rep(top_up,5)),   ## planned for 2020 - ** Assuming the distribution coverage matched the RCT estimate
    
    retention = data1$itn_leave_dur[row_drawn] * year, ## from trial (try matching this for all as well...)
    
    ## each row needs to show the efficacy parameter across years (and cols are diff mosquito)
    dn0 = matrix(c(rep(data1$itn_kill_pyrethroid_ITN[row_drawn],8),
                   data1$itn_kill_TRIAL_NET[row_drawn],
                   rep(data1$itn_kill_pyrethroid_ITN[row_drawn],5),
                   data1$itn_kill_TRIAL_NET[row_drawn],
                   rep(data1$itn_kill_pyrethroid_ITN[row_drawn],5),
                   rep(data1$itn_kill_pyrethroid_ITN[row_drawn],8),
                   data1$itn_kill_TRIAL_NET[row_drawn],
                   rep(data1$itn_kill_pyrethroid_ITN[row_drawn],5),
                   data1$itn_kill_TRIAL_NET[row_drawn],
                   rep(data1$itn_kill_pyrethroid_ITN[row_drawn],5)),
                 nrow=20, ncol=2),
    rn =  matrix(c(rep(data1$itn_repel_pyrethroid_ITN[row_drawn],8),
                   data1$itn_repel_TRIAL_NET[row_drawn],
                   rep(data1$itn_repel_pyrethroid_ITN[row_drawn],5),
                   data1$itn_repel_TRIAL_NET[row_drawn],
                   rep(data1$itn_repel_pyrethroid_ITN[row_drawn],5),
                   rep(data1$itn_repel_pyrethroid_ITN[row_drawn],8),
                   data1$itn_repel_TRIAL_NET[row_drawn],
                   rep(data1$itn_repel_pyrethroid_ITN[row_drawn],5),
                   data1$itn_repel_TRIAL_NET[row_drawn],
                   rep(data1$itn_repel_pyrethroid_ITN[row_drawn],5)),
                 nrow=20, ncol=2),
    rnm = matrix(.24, nrow=20, ncol=2),
    gamman = c(rep(data1$itn_half_life_pyrethroid_ITN[row_drawn] * 365, 8),
               data1$itn_half_life_TRIAL_NET[row_drawn] * 365,
               rep(data1$itn_half_life_pyrethroid_ITN[row_drawn] * 365, 5),
               data1$itn_half_life_TRIAL_NET[row_drawn] * 365,
               rep(data1$itn_half_life_pyrethroid_ITN[row_drawn] * 365, 5))
  )
  
  correlationsb1 <- get_correlation_parameters(bednetparams_1)
  correlationsb1$inter_round_rho('bednets', 0)
  
  ## Run the simulations
  output1 <- run_simulation(sim_length, bednetparams_1,correlationsb1)
  
  output1$pv_182.5_5110 = output1$n_detect_182.5_5110/output1$n_182.5_5110
  output1$pv_0_1825 = output1$n_detect_0_1825/output1$n_0_1825
  output1$pv_1825_3650 = output1$n_detect_1825_3650/output1$n_1825_3650
  output1$pv_0_3650 = output1$n_detect_0_3650/output1$n_0_3650
  output1$pv_730_3650 = output1$n_detect_730_3650/output1$n_730_3650
  output1$pv_0_36500 = output1$n_detect_0_36500/output1$n_0_36500
  
  output1$clin_inc_182.5_3650 = output1$n_inc_clinical_182.5_3650/output1$n_182.5_3650
  
  return(data.frame(timestep = output1$timestep,
                    prev_05_14 = output1$pv_182.5_5110,
                    prev_0_5 = output1$pv_0_1825,
                    prev_5_10 =  output1$pv_1825_3650,
                    prev_0_10 =  output1$pv_0_3650,
                    prev_all_age = output1$pv_0_36500,
                    prev_2_10 =  output1$pv_730_3650,
                    clin_inc_6m_10y = output1$clin_inc_182.5_3650))
}

## Things to try - using the range of uncertainty from the resistance estimates
## Setting net retention to be equivalent for the 3 nets
## Updating the drug treatment assumptions if trial team has more info
## simulating no top up of nets
## simulating top up with trial specific nets in each arm
## outputing net use?

# pyr_nets_cali = 133.6573
# pbo_nets_cali = 70.83473
# pp_nets_cali = 96.35578
  
pyr_nets = malsim_smc_f(EIR_L = pyr_nets_cali,row_drawn =1,top_up = 0.1577571)
pbo_nets = malsim_smc_f(EIR_L = pbo_nets_cali,row_drawn =2,top_up = 0.1334457)
pp_nets =  malsim_smc_f(EIR_L = pp_nets_cali,row_drawn =3,top_up = 0.1451762)

pyr_nets_min = malsim_smc_f(EIR_L = pyr_nets_cali,row_drawn =4,top_up = 0.1577571)
pbo_nets_min = malsim_smc_f(EIR_L = pbo_nets_cali,row_drawn =5,top_up = 0.1334457)
pp_nets_min =  malsim_smc_f(EIR_L = pp_nets_cali,row_drawn =6,top_up = 0.1451762)

pyr_nets_max = malsim_smc_f(EIR_L = pyr_nets_cali,row_drawn =7,top_up = 0.1577571)
pbo_nets_max = malsim_smc_f(EIR_L = pbo_nets_cali,row_drawn =8,top_up = 0.1334457)
pp_nets_max =  malsim_smc_f(EIR_L = pp_nets_cali,row_drawn =9,top_up = 0.1451762)

par(mfrow=c(2,2))
par(mar=c(4,4,3,2))

plot(pyr_nets$prev_05_14 ~ pyr_nets$timestep,
     main = "",ylim=c(0,0.8),
     ylab = "Prevalence in children 0.5 to 14-yrs (%)",
     cex.lab=1.2,cex.axis=1.2,
     col="white",yaxt="n",
     xlab = "",type="l",xlim=c(3,7)*365,xaxt="n")
axis(1,at=c(365*c(2,3,4,5,6,7,8,9)),
     labels = c("Jan 2017","Jan 2018","Jan 2019","Jan 2020","Jan 2021","Jan 2022","Jan 2023","Jan 2024"))
axis(2,las=2,at=c(0,0.2,0.4,0.6,0.8),labels=seq(0,80,20),cex=1.2,cex.axis=1.2)

polygon(c(pyr_nets$timestep,rev(pyr_nets$timestep)),
        c(pyr_nets_min$prev_05_14,rev(pyr_nets_max$prev_05_14)),
        border=NA,col=adegenet::transp("darkgrey",0.3))
lines(pyr_nets$prev_05_14 ~ pyr_nets$timestep,col="black",lty=1,lwd=1) ## pyr-only


polygon(c(pbo_nets$timestep,rev(pbo_nets$timestep)),
        c(pbo_nets_min$prev_05_14,rev(pbo_nets_max$prev_05_14)),
        border=NA,col=adegenet::transp("purple",0.3))
lines(pbo_nets$prev_05_14 ~ pbo_nets$timestep,col="purple",lty=1,lwd=1) ## pyr-PBO

polygon(c(pp_nets$timestep,rev(pp_nets$timestep)),
        c(pp_nets_min$prev_05_14,rev(pp_nets_max$prev_05_14)),
        border=NA,col=adegenet::transp("aquamarine3",0.4))
lines(pp_nets$prev_05_14 ~ pp_nets$timestep,col="aquamarine4",lty=1,lwd=1) ## IG2

### Check and use what is in Mosha et al 2022
# Post 12 months (Jan 2019) median 17.6 (min 0 and max 73.3)
# Post 18 months (Aug 2019) median 48.5 (min 3.0 and max 79.2)
# Post 24 months (Jan 2020) median 35.5 (min 2.4 and max 76.7)
points(c(0.459,0.42,0.427)~rep(c(3*365+31+28+31+30+31+30+31+31+30+15),3),
       col=c("black","purple","aquamarine3"),pch=19,cex=2)

year = 365
timestep = round(c(0, 4, 7) * year + c(27, ## from NNP Sept 2017
                                       27,##assumed target for NNP is 1 month after baseline survey
                                       27),0)
abline(v=timestep,col="grey",lty=2)

time_obs_epi = c(5*year+1,5*year+31+28+31+30+31+30+31+1,6*year+1)
points(c(0.31,0.51,0.44) ~ time_obs_epi,cex=1.3,pch=19) ## pyrethroid only
points(c(0.20,0.41,0.40) ~ time_obs_epi,pch=19,cex=1.3,col="purple") ## pyrethroid PBO
points(c(0.154,0.39,0.23) ~ time_obs_epi,pch=19,cex=1.3,col="aquamarine4") ## pyrethroid pyrrole


##############
##
## lm of arm preds
##
obs_val = c(0.31,0.51,0.44,0.20,0.41,0.40,0.154,0.39,0.23)

pred_val_med = c(pyr_nets$prev_05_14[c(time_obs_epi)],
                 pbo_nets$prev_05_14[c(time_obs_epi)],
                 pp_nets$prev_05_14[c(time_obs_epi)])

pred_val_min = c(pyr_nets_min$prev_05_14[c(time_obs_epi)],
                 pbo_nets_min$prev_05_14[c(time_obs_epi)],
                 pp_nets_min$prev_05_14[c(time_obs_epi)])

pred_val_max = c(pyr_nets_max$prev_05_14[c(time_obs_epi)],
                 pbo_nets_max$prev_05_14[c(time_obs_epi)],
                 pp_nets_max$prev_05_14[c(time_obs_epi)])

m1 = lm(obs_val ~ pred_val_med + 0)
summary.lm(m1)

plot(obs_val ~ pred_val_med,
     ylim = c(0,0.6),xlim=c(0,0.6),
     cex.lab=1.2,cex.axis=1.2,
     ylab = "Empirical data prevalence (%)",
     xlab = "Model simulated prevalence (%)",
     pch=19,col=rep(c("black","purple","aquamarine4"),each=3))
abline(a=0,b=1,lty=2,col="grey")
x=seq(0,1,0.01)
y = summary.lm(m1)$coef[1,1] * x + 0
y_l = (summary.lm(m1)$coef[1,1]+summary.lm(m1)$coef[1,2]) * x + 0
y_u = (summary.lm(m1)$coef[1,1]-summary.lm(m1)$coef[1,2]) * x + 0
lines(y[1:60] ~ x[1:60])
polygon(c(x[1:60],rev(x[1:60])),
        c(y_l[1:60],rev(y_u[1:60])),col=adegenet::transp("grey",0.5),
        border=NA)
points(obs_val ~ pred_val_med,pch=19,col=rep(c("black","purple","aquamarine4"),each=3))
segments(x0=pred_val_min,x1 = pred_val_max,
         y0=obs_val,y1=obs_val,
         col=rep(c("black","purple","aquamarine4"),each=3))
#################
##
## Clinical incidence 

annual_obs_incy1_pyr = sum(pyr_nets$clin_inc_6m_10y[1487:1851])
annual_obs_incy1_pbo = sum(pbo_nets$clin_inc_6m_10y[1487:1851])
annual_obs_incy1_pyrrole = sum(pp_nets$clin_inc_6m_10y[1487:1851])

annual_obs_incy2_pyr = sum(pyr_nets$clin_inc_6m_10y[1851:2216])
annual_obs_incy2_pbo = sum(pbo_nets$clin_inc_6m_10y[1851:2216])
annual_obs_incy2_pyrrole = sum(pp_nets$clin_inc_6m_10y[1851:2216])

minannual_obs_incy1_pyr = sum(pyr_nets_min$clin_inc_6m_10y[1487:1851])
minannual_obs_incy1_pbo = sum(pbo_nets_min$clin_inc_6m_10y[1487:1851])
minannual_obs_incy1_pyrrole = sum(pp_nets_min$clin_inc_6m_10y[1487:1851])

minannual_obs_incy2_pyr = sum(pyr_nets_min$clin_inc_6m_10y[1851:2216])
minannual_obs_incy2_pbo = sum(pbo_nets_min$clin_inc_6m_10y[1851:2216])
minannual_obs_incy2_pyrrole = sum(pp_nets_min$clin_inc_6m_10y[1851:2216])

maxannual_obs_incy1_pyr = sum(pyr_nets_max$clin_inc_6m_10y[1487:1851])
maxannual_obs_incy1_pbo = sum(pbo_nets_max$clin_inc_6m_10y[1487:1851])
maxannual_obs_incy1_pyrrole = sum(pp_nets_max$clin_inc_6m_10y[1487:1851])

maxannual_obs_incy2_pyr = sum(pyr_nets_max$clin_inc_6m_10y[1851:2216])
maxannual_obs_incy2_pbo = sum(pbo_nets_max$clin_inc_6m_10y[1851:2216])
maxannual_obs_incy2_pyrrole = sum(pp_nets_max$clin_inc_6m_10y[1851:2216])



barplot(c(annual_obs_incy1_pyr,annual_obs_incy2_pyr,NA,
          annual_obs_incy1_pbo,annual_obs_incy2_pbo,NA,
          annual_obs_incy1_pyrrole,annual_obs_incy2_pyrrole),
        # col=adegenet::transp(c("darkgrey","darkgrey",NA,
        #                      "purple","purple",NA,
        #                      "aquamarine3","aquamarine3"),0.6),
        col=rep("white",8), bty="n",
        yaxt="n",ylim = c(0,2.5),border=NA,
        cex.lab=1.2,cex.axis=1.2,
        ylab = "Clinical incidence per child-year (0.5 to 10-yrs)")
axis(2,las=2,seq(0,1.6,0.2),cex=1.2,cex.axis=1.2)

segments(x0=seq(0.8,9,length=8),x1=seq(0.8,9,length=8),
         y0=c(minannual_obs_incy1_pyr,minannual_obs_incy2_pyr,NA,
              minannual_obs_incy1_pbo,minannual_obs_incy2_pbo,NA,
              minannual_obs_incy1_pyrrole,minannual_obs_incy2_pyrrole),
         y1=c(maxannual_obs_incy1_pyr,maxannual_obs_incy2_pyr,NA,
              maxannual_obs_incy1_pbo,maxannual_obs_incy2_pbo,NA,
              maxannual_obs_incy1_pyrrole,maxannual_obs_incy2_pyrrole),
         col=c("black","black",NA,
               "purple","purple",NA,
               "aquamarine4","aquamarine4"))
points(c(annual_obs_incy1_pyr,annual_obs_incy2_pyr,NA,
         annual_obs_incy1_pbo,annual_obs_incy2_pbo,NA,
         annual_obs_incy1_pyrrole,annual_obs_incy2_pyrrole)~seq(0.8,9,length=8),
       col=c("black","black",NA,
             "purple","purple",NA,
             "aquamarine4","aquamarine4"),
       pch=15)
axis(1,at=seq(0.8,9,length=8),
     labels=c("Yr 1","Yr 2","","Yr 1","Yr 2","","Yr 1","Yr 2"),cex=1.2,cex.axis=1.2)

points(c(0.32,0.57,NA,0.13,0.48,NA,0.13,0.31)~seq(0.9,9.1,length=8),
       col=c("black","black",NA,
             "purple","purple",NA,
             "aquamarine4","aquamarine4"),pch=19)

legend("topright",title = "Net type",
       legend = c("Pyrethroid-only",
                  "Pyrethroid-PBO",
                  "Pyrethroid-pyrrole"),
       col=c("black","purple","aquamarine3"),
       cex=1,pch=15,bty="n")


## relative difference in clinical incidence from pyrethroid-only nets
pbo_rel_inc_y1 = (annual_obs_incy1_pyr-annual_obs_incy1_pbo)/annual_obs_incy1_pyr
pbo_rel_inc_y2 = (annual_obs_incy2_pyr-annual_obs_incy2_pbo)/annual_obs_incy2_pyr

pp_rel_inc_y1 = (annual_obs_incy1_pyr-annual_obs_incy1_pyrrole)/annual_obs_incy1_pyr
pp_rel_inc_y2 = (annual_obs_incy2_pyr-annual_obs_incy2_pyrrole)/annual_obs_incy2_pyr

maxpbo_rel_inc_y1 = (annual_obs_incy1_pyr-maxannual_obs_incy1_pbo)/annual_obs_incy1_pyr
maxpbo_rel_inc_y2 = (annual_obs_incy2_pyr-maxannual_obs_incy2_pbo)/annual_obs_incy2_pyr

maxpp_rel_inc_y1 = (annual_obs_incy1_pyr-maxannual_obs_incy1_pyrrole)/annual_obs_incy1_pyr
maxpp_rel_inc_y2 = (annual_obs_incy2_pyr-maxannual_obs_incy2_pyrrole)/annual_obs_incy2_pyr

minpbo_rel_inc_y1 = (annual_obs_incy1_pyr-minannual_obs_incy1_pbo)/annual_obs_incy1_pyr
minpbo_rel_inc_y2 = (annual_obs_incy2_pyr-minannual_obs_incy2_pbo)/annual_obs_incy2_pyr

minpp_rel_inc_y1 = (annual_obs_incy1_pyr-minannual_obs_incy1_pyrrole)/annual_obs_incy1_pyr
minpp_rel_inc_y2 = (annual_obs_incy2_pyr-minannual_obs_incy2_pyrrole)/annual_obs_incy2_pyr

c(0.32,0.57,NA,0.13,0.48,NA,0.13,0.31)
pbo_trial_relative_inc_y1 = (0.32 - 0.13)/0.32
pbo_trial_relative_inc_y2 = (0.32 - 0.48)/0.32
pyrrole_trial_relative_inc_y1 = (0.57 - 0.13)/0.57
pyrrole_trial_relative_inc_y2 = (0.57 - 0.31)/0.57


barplot(c(pbo_rel_inc_y1,pbo_rel_inc_y2,NA,
          pp_rel_inc_y1,pp_rel_inc_y2),
        col=rep("white",5), bty="n",
        yaxt="n",ylim = c(-1,1),border=NA,
        cex.lab=1.2,cex.axis=1.2,
        ylab = "Relative reduction in incidence")
axis(2,las=2,seq(0,1,0.2),labels=seq(0,100,20),cex=1.2,cex.axis=1.2)

segments(x0=seq(0.8,5.7,length=5),x1=seq(0.8,5.7,length=5),
         y0=c(minpbo_rel_inc_y1,minpbo_rel_inc_y2,NA,
              minpp_rel_inc_y1,minpp_rel_inc_y2),
         y1=c(maxpbo_rel_inc_y1,maxpbo_rel_inc_y2,NA,
              maxpp_rel_inc_y1,maxpp_rel_inc_y2),
         col=c("purple","purple",NA,
               "aquamarine4","aquamarine4"))
points(c(pbo_rel_inc_y1,pbo_rel_inc_y2,NA,
         pp_rel_inc_y1,pp_rel_inc_y2)~seq(0.8,5.7,length=5),
       col=c("purple","purple",NA,
             "aquamarine4","aquamarine4"),
       pch=15)
axis(1,at=seq(0.8,5.7,length=5),
     labels=c("Yr 1","Yr 2","","Yr 1","Yr 2"),cex=1.2,cex.axis=1.2)

points(c(pbo_trial_relative_inc_y1,pbo_trial_relative_inc_y2,NA,
         pyrrole_trial_relative_inc_y1,pyrrole_trial_relative_inc_y2)~seq(0.9,5.8,length=5),
       col=c("purple","purple",NA,
             "aquamarine4","aquamarine4"),pch=19)

par(xpd=NA,cex = 1)

text(x = -8.6, y = 4,"(A)")
text(x = -1, y = 4,"(B)")
text(x = -8.6, y = 1.2,"(C)")
text(x = -1, y = 1.2,"(D)")

#dim = 1000, 820