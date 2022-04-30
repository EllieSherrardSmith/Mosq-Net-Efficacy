library(malariasimulation)
library(malariaEquilibrium)
library(reshape2)
library(ggplot2)




set.seed(20130)
###########################################################################
##
## RCT with 4 arms
##  * Standard pyrethroid nets
##  * Pyrethroid PBO nets
##  * Standard pyrethroid nets + Actellic IRS
##  * Pyrethroid PBO nets + Actellic IRS

## Protopopoff 2018 doi:10.1016/S0140-6736(18)30427-6


## sites
sites_13 = read.csv("validation/data/site_file_13.csv",header=TRUE)
## specifics
params_13 = read.csv("validation/data/input_params_13_update_net_parameters2_all_hutdataB.csv",header=TRUE) 

# params_13$ITN_1 = 0
# params_13$IRS_1
## Site data,  
# sites_13$itn_cov4 = sites_13$itn_cov5 = sites_13$itn_cov6 = 
#   sites_13$itn_cov7 = sites_13$itn_cov8 = sites_13$itn_cov9 = sites_13$itn_cov10 = 
#   sites_13$itn_cov11=sites_13$itn_cov12 = sites_13$itn_cov13 = sites_13$itn_cov14 = 
#   sites_13$itn_cov15 = sites_13$itn_cov16 = sites_13$itn_cov17 = sites_13$itn_cov18 = 
#   sites_13$itn_cov19 = sites_13$itn_cov20 = 0


##
##
## MODEL SET UP - for arm 1 and 3 and distinct for 2 and 4 (0.68, 0.67 versus 0.61 and 0.64)
year <- 365
month <- 30
sim_length <- 10 * year
human_population <- 10000
starting_EIR_U <- 240 ## upper estimate
starting_EIR_L <- 180 ## lower estimate

simparams <- get_parameters(
  list(
    human_population = human_population,
    # irs_correlation = 
    
    prevalence_rendering_min_ages = 0.5 * 365, ## Prev in 6 months to 14 years measured
    prevalence_rendering_max_ages = 14 * 365,
    
    model_seasonality = TRUE, ## Seasonality to match study site inputs [sites_13]
    g0 = sites_13$seasonal_a0[1],
    g = c(sites_13$seasonal_a1[1], sites_13$seasonal_a2[1], sites_13$seasonal_a3[1]),
    h = c(sites_13$seasonal_b1[1], sites_13$seasonal_b2[1], sites_13$seasonal_b3[1]),
    
    Q0 = sites_13$gamb_ss_Q0[1],
    individual_mosquitoes = FALSE ## Update next
  )
)

simparamsU <- set_equilibrium(simparams, starting_EIR_U)
simparamsL <- set_equilibrium(simparams, starting_EIR_L)

# set species
simparamsU <- set_species(simparamsU,
                          species=list(gamb_params, arab_params, fun_params),
                          proportions=c(sites_13$prop_gamb_ss[1],sites_13$prop_arab[1],sites_13$prop_fun[1]))
# set treatment
simparamsU <- set_drugs(simparamsU, list(AL_params, SP_AQ_params, DHA_PQP_params))
simparamsU <- set_clinical_treatment(simparamsU, 
                                     drug=1,
                                     time=c(100),
                                     coverage=c(sites_13$drug_cov_0_0[1]))
simparamsU <- set_clinical_treatment(simparamsU, 
                                     drug=2,
                                     time=c(100),
                                     coverage=c(sites_13$drug_cov_1_0[1]))


simparamsL <- set_species(simparamsL,
                          species=list(gamb_params, arab_params, fun_params),
                          proportions=c(sites_13$prop_gamb_ss[1001],sites_13$prop_arab[1001],sites_13$prop_fun[1001]))
# set treatment
simparamsL <- set_drugs(simparamsL, list(AL_params, SP_AQ_params, DHA_PQP_params))
simparamsL <- set_clinical_treatment(simparamsL, 
                                     drug=1,
                                     time=c(100),
                                     coverage=c(sites_13$drug_cov_0_0[1001]))
simparamsL <- set_clinical_treatment(simparamsL, 
                                     drug=2,
                                     time=c(100),
                                     coverage=c(sites_13$drug_cov_1_0[1001]))



## Set up the independent arms of the trial

bednetparamsU <- simparamsU
bednetparamsL <- simparamsL

bednet_events = data.frame(
  timestep = c(1, 4, 7) * year,
  name=c("Bednets start1", "redistribute","redistribute")
)


## Pyrethroid only nets
bednetparams_1 <- set_bednets(
  bednetparamsU,
  
  timesteps = bednet_events$timestep,
  
  coverages = c(sites_13$itn_cov20[1], ## historic 
                sites_13$itn_cov20[1], ## historic  
                params_13$ITN_1[1]),   ## TRIAL NETS
  
  retention = params_13$itn_leave_dur[1] * year,
  
  ## April 2022 update - WHO Recommended Nets alone: RCT Validation
  ## 83% pyrethroid resistance
  ##use parameters/pyrethroid_only_nets.csv
  dn0 = matrix(c(0.245997991, 0.245997991, 0.245997991), nrow=3, ncol=3),
  rn = matrix(c(0.706149167,0.706149167,0.706149167), nrow=3, ncol=3),
  rnm = matrix(c(.24, .24, .24), nrow=3, ncol=3),
  gamman = c(2.288886539, 2.288886539, 2.288886539) * 365
)


## Pyrethroid PBO nets
bednetparams_2 <- set_bednets(
  bednetparamsL,
  
  timesteps = bednet_events$timestep,
  
  coverages = c(sites_13$itn_cov20[1], ## historic 
                sites_13$itn_cov20[1], ## historic  
                params_13$ITN_1[1001]),   ## TRIAL NETS
  
  retention = params_13$itn_leave_dur[1001] * year,
  
  ## April 2022 update - WHO Recommended Nets alone: RCT Validation
  ## 83& pyrethroid resistance
  ## Pyrethroid-pbo nets
  dn0 = matrix(c(0.389611527, 
                 0.389611527, 
                 0.389611527), nrow=3, ncol=3),
  rn = matrix(c(0.596244584,
                0.596244584,
                0.596244584), nrow=3, ncol=3),
  rnm = matrix(c(.24, .24, .24), nrow=3, ncol=3),
  gamman = c(1.90343459,1.90343459,1.90343459) * 365 ## TRIAL NETS
)

##
##
##
## For arms 3 and 4 there is also nets

sprayingparams1 <- bednetparams_1

peak <- peak_season_offset(sprayingparams1)

spraying_events = data.frame(
  timestep = c(7) * year +30, ## spraying done end of Jan
  name=c("Spraying (1 year only)")
)

sprayingparams_1 <- set_spraying(
  sprayingparams1, ## details same as pyrethroid only nets
  
  timesteps = spraying_events$timestep,
  coverages = rep(params_13$IRS_1[2001], 1),

  ls_theta = matrix(c(params_13$irs_decay_mort1[2001]), nrow = 1, ncol = 3),
  ls_gamma = matrix(c(params_13$irs_decay_mort2[2001]), nrow = 1, ncol = 3),
  ks_theta = matrix(c(params_13$irs_decay_succ1[2001]), nrow = 1, ncol = 3),
  ks_gamma = matrix(c(params_13$irs_decay_succ2[2001]), nrow = 1, ncol = 3),
  ms_theta = matrix(c(params_13$irs_decay_det1[2001]), nrow = 1, ncol = 3),
  ms_gamma = matrix(c(params_13$irs_decay_mort2[2001]), nrow = 1, ncol = 3) ## same as mortality as unable to measure deterrence
)



sprayingparams2 <- bednetparams_2

peak <- peak_season_offset(sprayingparams2)

sprayingparams_2 <- set_spraying(
  sprayingparams2,
  
  timesteps = spraying_events$timestep,
  coverages = rep(params_13$IRS_1[3001], 1),
  
  ls_theta = matrix(c(params_13$irs_decay_mort1[3001]), nrow = 1, ncol = 3),
  ls_gamma = matrix(c(params_13$irs_decay_mort2[3001]), nrow = 1, ncol = 3),
  ks_theta = matrix(c(params_13$irs_decay_succ1[3001]), nrow = 1, ncol = 3),
  ks_gamma = matrix(c(params_13$irs_decay_succ2[3001]), nrow = 1, ncol = 3),
  ms_theta = matrix(c(params_13$irs_decay_det1[3001]), nrow = 1, ncol = 3),
  ms_gamma = matrix(c(params_13$irs_decay_mort2[3001]), nrow = 1, ncol = 3) ## same as mortality as unable to measure deterrence
)


correlationsb1 <- get_correlation_parameters(bednetparams_1)
correlationsb2 <- get_correlation_parameters(bednetparams_2)
correlationsb1s <- get_correlation_parameters(sprayingparams_1)
correlationsb2s <- get_correlation_parameters(sprayingparams_2)

correlationsb1$inter_round_rho('bednets', 1)
correlationsb2$inter_round_rho('bednets', 1)
correlationsb1s$inter_intervention_rho('bednets', 'spraying', 1)
correlationsb2s$inter_intervention_rho('bednets', 'spraying', 1)
# 
# #...
# run_simulation(timesteps, parameters, correlations)

output1 <- run_simulation(sim_length, bednetparams_1,correlationsb1)
output2 <- run_simulation(sim_length, bednetparams_2,correlationsb2)
output3 <- run_simulation(sim_length, sprayingparams_1,correlationsb1s)
output4 <- run_simulation(sim_length, sprayingparams_2,correlationsb2s)

output1$pv_182.5_5110 = output1$n_detect_182.5_5110/output1$n_182.5_5110
output2$pv_182.5_5110 = output2$n_detect_182.5_5110/output2$n_182.5_5110
output3$pv_182.5_5110 = output3$n_detect_182.5_5110/output3$n_182.5_5110
output4$pv_182.5_5110 = output4$n_detect_182.5_5110/output4$n_182.5_5110

par(mfrow=c(1,1))
plot(output1$pv_182.5_5110 ~ output1$timestep,
     xlim=c(6.8,10)*365,
     ylim=c(0,1),pch="",
     col="darkred",ylab = "Prevalence 6 months to 14 years")
lines(output1$pv_182.5_5110 ~ output1$timestep,col="darkred")
lines(output2$pv_182.5_5110 ~ output1$timestep,col="aquamarine3")
lines(output3$pv_182.5_5110 ~ output1$timestep,col="orange")
lines(output4$pv_182.5_5110 ~ output1$timestep,col="darkblue")

loop_func_arm_1 = function(row_val){
  ## Pyrethroid only nets
  bednetparams_1 <- set_bednets(
    bednetparamsU,
    
    timesteps = bednet_events$timestep,
    
    coverages = c(sites_13$itn_cov20[row_val], ## historic 
                  sites_13$itn_cov20[row_val], ## historic  
                  params_13$ITN_1[row_val]),   ## TRIAL NETS
    
    retention = params_13$itn_leave_dur[row_val] * year,
    
    dn0 = matrix(c(0.245997991, 0.245997991, 0.245997991), nrow=3, ncol=3),
    rn = matrix(c(0.706149167,0.706149167,0.706149167), nrow=3, ncol=3),
    rnm = matrix(c(.24, .24, .24), nrow=3, ncol=3),
    gamman = c(2.288886539, 2.288886539, 2.288886539) * 365
  )
  correlationsb1 <- get_correlation_parameters(bednetparams_1)
  correlationsb1$inter_round_rho('bednets', 1)
  
  output1 <- run_simulation(sim_length, bednetparams_1,correlationsb1)
  
  store1 = output1$pv_182.5_5110
  
  return(store1)
}

loop_func_arm_2 = function(row_val){
  ## Pyrethroid PBO nets
  bednetparams_2 <- set_bednets(
    bednetparamsL,
    
    timesteps = bednet_events$timestep,
    
    coverages = c(sites_13$itn_cov20[row_val], ## historic 
                  sites_13$itn_cov20[row_val], ## historic  
                  params_13$ITN_1[1000+row_val]),   ## TRIAL NETS
    
    retention = params_13$itn_leave_dur[1000+row_val] * year,
    
    ## Pyrethroid-pbo nets
    dn0 = matrix(c(0.389611527, 
                   0.389611527, 
                   0.389611527), nrow=3, ncol=3),
    rn = matrix(c(0.596244584,
                  0.596244584,
                  0.596244584), nrow=3, ncol=3),
    rnm = matrix(c(.24, .24, .24), nrow=3, ncol=3),
    gamman = c(1.90343459,1.90343459,1.90343459) * 365 ## TRIAL NETS
  )
  
  correlationsb2 <- get_correlation_parameters(bednetparams_2)
  correlationsb2$inter_round_rho('bednets', 1)
  
  output2 <- run_simulation(sim_length, bednetparams_2,correlationsb2)
  
  store2 = output2$pv_182.5_5110
  
  return(store2)
}



loop_func_arm_3 = function(row_val){
  ## Pyrethroid only nets
  bednetparams_1 <- set_bednets(
    bednetparamsU,
    
    timesteps = bednet_events$timestep,
    
    coverages = c(sites_13$itn_cov20[row_val], ## historic 
                  sites_13$itn_cov20[row_val], ## historic  
                  params_13$ITN_1[2000+row_val]),   ## TRIAL NETS
    
    retention = params_13$itn_leave_dur[2000+row_val] * year,
    
    dn0 = matrix(c(0.245997991, 0.245997991, 0.245997991), nrow=3, ncol=3),
    rn = matrix(c(0.706149167,0.706149167,0.706149167), nrow=3, ncol=3),
    rnm = matrix(c(.24, .24, .24), nrow=3, ncol=3),
    gamman = c(2.288886539, 2.288886539, 2.288886539) * 365
  )
  
  ## For arms 3 and 4 there is also nets
  
  sprayingparams1 <- bednetparams_1
  
  peak <- peak_season_offset(sprayingparams1)
  
  spraying_events = data.frame(
    timestep = c(7) * year +30, ## spraying done end of Jan
    name=c("Spraying (1 year only)")
  )
  
  sprayingparams_1 <- set_spraying(
    sprayingparams1, ## details same as pyrethroid only nets
    
    timesteps = spraying_events$timestep,
    coverages = rep(params_13$IRS_1[2000+row_val], 1),
    
    ls_theta = matrix(c(params_13$irs_decay_mort1[2000+row_val]), nrow = 1, ncol = 3),
    ls_gamma = matrix(c(params_13$irs_decay_mort2[2000+row_val]), nrow = 1, ncol = 3),
    ks_theta = matrix(c(params_13$irs_decay_succ1[2000+row_val]), nrow = 1, ncol = 3),
    ks_gamma = matrix(c(params_13$irs_decay_succ2[2000+row_val]), nrow = 1, ncol = 3),
    ms_theta = matrix(c(params_13$irs_decay_det1[2000+row_val]), nrow = 1, ncol = 3),
    ms_gamma = matrix(c(params_13$irs_decay_mort2[2000+row_val]), nrow = 1, ncol = 3) ## same as mortality as unable to measure deterrence
  )
  correlationsb1s <- get_correlation_parameters(sprayingparams_1)
  correlationsb1s$inter_intervention_rho('bednets', 'spraying', 1)
  
  
  output3 <- run_simulation(sim_length, sprayingparams_1,correlationsb1s)
  
  store3 = output3$pv_182.5_5110
  
  return(store3)
}


loop_func_arm_4 = function(row_val){
  ## Pyrethroid PBO nets
  bednetparams_2 <- set_bednets(
    bednetparamsL,
    
    timesteps = bednet_events$timestep,
    
    coverages = c(sites_13$itn_cov20[row_val], ## historic 
                  sites_13$itn_cov20[row_val], ## historic  
                  params_13$ITN_1[1000+row_val]),   ## TRIAL NETS
    
    retention = params_13$itn_leave_dur[1000+row_val] * year,
    
    ## Pyrethroid-pbo nets
    dn0 = matrix(c(0.389611527, 
                   0.389611527, 
                   0.389611527), nrow=3, ncol=3),
    rn = matrix(c(0.596244584,
                  0.596244584,
                  0.596244584), nrow=3, ncol=3),
    rnm = matrix(c(.24, .24, .24), nrow=3, ncol=3),
    gamman = c(1.90343459,1.90343459,1.90343459) * 365 ## TRIAL NETS
  )
  
  sprayingparams2 <- bednetparams_2
  
  peak <- peak_season_offset(sprayingparams2)
  
  sprayingparams_2 <- set_spraying(
    sprayingparams2,
    
    timesteps = spraying_events$timestep,
    coverages = rep(params_13$IRS_1[3000+row_val], 1),
    
    ls_theta = matrix(c(params_13$irs_decay_mort1[3000+row_val]), nrow = 1, ncol = 3),
    ls_gamma = matrix(c(params_13$irs_decay_mort2[3000+row_val]), nrow = 1, ncol = 3),
    ks_theta = matrix(c(params_13$irs_decay_succ1[3000+row_val]), nrow = 1, ncol = 3),
    ks_gamma = matrix(c(params_13$irs_decay_succ2[3000+row_val]), nrow = 1, ncol = 3),
    ms_theta = matrix(c(params_13$irs_decay_det1[3000+row_val]), nrow = 1, ncol = 3),
    ms_gamma = matrix(c(params_13$irs_decay_mort2[3000+row_val]), nrow = 1, ncol = 3) ## same as mortality as unable to measure deterrence
  )
  correlationsb2s <- get_correlation_parameters(sprayingparams_2)
  correlationsb2s$inter_intervention_rho('bednets', 'spraying', 1)
  
  output4 <- run_simulation(sim_length, sprayingparams_2,correlationsb2s)
  
  store4 = output4$pv_182.5_5110
  
  return(store4)
}

rand_num = sample(1:1000,11,replace=FALSE)
outputs_arm1 = outputs_arm2 = outputs_arm3 = outputs_arm4 = array(dim=c(nrow=length(output1$timestep),ncol=11))
for(i in 1:2){
  
  outputs_arm1[,i+1] = loop_func_arm_1(rand_num[i])
  outputs_arm2[,i+1] = loop_func_arm_2(rand_num[i])
  outputs_arm3[,i+1] = loop_func_arm_3(rand_num[i])
  outputs_arm4[,i+1] = loop_func_arm_4(rand_num[i])
  
}



plot(output1$pv_182.5_5110 ~ output1$timestep,
     xlim=c(6.8,10)*365,
     ylim=c(0,1),pch="",
     col="darkred",ylab = "Prevalence 6 months to 14 years")
for(i in 2:3){
  lines(outputs_arm1[,i] ~ output1$timestep,col=adegenet::transp("darkred",0.3))  
  lines(outputs_arm2[,i] ~ output1$timestep,col=adegenet::transp("aquamarine3",0.3))  
  lines(outputs_arm3[,i] ~ output1$timestep,col=adegenet::transp("gold2",0.3))  
  lines(outputs_arm4[,i] ~ output1$timestep,col=adegenet::transp("darkblue",0.3))  
}


abline(v=7*365,lty=2,col="grey")
abline(v=7*365+1*365,lty=2,col="grey")
abline(v=7*365+2*365,lty=2,col="grey")
abline(v=7*365+3*365,lty=2,col="grey")

lines(rowMeans(outputs_arm1[,2:11]) ~ output1$timestep,col="darkred",lwd=2)
lines(rowMeans(outputs_arm2[,2:11]) ~ output1$timestep,col="aquamarine3",lwd=2)
lines(rowMeans(outputs_arm3[,2:11]) ~ output1$timestep,col="gold2",lwd=2)
lines(rowMeans(outputs_arm4[,2:11]) ~ output1$timestep,col="darkblue",lwd=2)

# lines(output1$pv_182.5_5110 ~ output1$timestep,col="darkred",lwd=2)
# lines(output2$pv_182.5_5110 ~ output2$timestep,col="aquamarine3",lwd=2)
# lines(output3$pv_182.5_5110 ~ output3$timestep,col="gold2",lwd=2)
# lines(output4$pv_182.5_5110 ~ output4$timestep,col="darkblue",lwd=2)

### Add trial points

Tlo_arms = array(dim=c(7,4))
Tlo_arms[,1] = c(0.68, 0.555, 0.553, 0.530, 0.68,  0.827, 0.648)
Tlo_arms[,2] = c(0.61, 0.458, 0.311, 0.348, 0.459, 0.6835,0.490)
Tlo_arms[,3] = c(0.67, 0.386, 0.264, 0.399, 0.551, NA,    NA)
Tlo_arms[,4] = c(0.64, 0.375, 0.286, 0.291, 0.437, NA,    NA)
time_match_L1 = c(-0.1666667, 0.33,  0.75,  1.33,  1.75,  2.33,2.75) * 365 + (7*365) 


## add baselines
cols_13 = c("darkred","aquamarine3","gold2","darkblue")
for(i in 1:4){
  points(Tlo_arms[1:5,i] ~ time_match_L1[1:5], col = cols_13[i], pch=19,cex=1.5)
}










