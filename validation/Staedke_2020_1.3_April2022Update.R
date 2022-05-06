library(malariasimulation)
library(malariaEquilibrium)
library(reshape2)
library(ggplot2)




set.seed(20130)


#############################
## RCT with 2 arms
##
## Pyrethroid nets
## Pyrethroid-PBO nets

## Staedke 2020 The Lancet


## sites
sites_14 = read.csv("validation/data/site_file_14.csv",header=TRUE)
## specifics
params_14 = read.csv("validation/data/input_params_14_net_parameters2_all_hutdataB.csv",header=TRUE) 

dat_res_pyr = read.csv("parameters/pyrethroid_only_nets.csv",header=TRUE) 
dat_res_pbo = read.csv("parameters/pyrethroid_pbo_nets.csv",header=TRUE) 
dat_res_pp = read.csv("parameters/pyrethroid_pyrrole_nets.csv",header=TRUE) 



malsim_actual_f = function(dat_row){
  
  year <- 365
  month <- 30
  sim_length <- 11 * year ## initially just til 2020 as need extra resistance data!
  ## This is spanning Jan 2014 - Dec 2025
  ## Assume that all places have received nets at start of 2014
  ## then relative to jan 2017 (see uganda_what_if...)
  ## then as noted for 2020
  
  
  ## and finally we wish to see the difference for these places running forward from 2020-2023 (1 or 3 years)
  human_population <- 10000
  starting_EIR <- 12
    
  simparams <- get_parameters(
    list(
      human_population = human_population,
      # irs_correlation = 
      
      prevalence_rendering_min_ages = 2 * 365, ## Prev in 6 months to 14 years measured
      prevalence_rendering_max_ages = 10 * 365,
      
      clinical_incidence_rendering_min_ages = 0 * 365, ## All age clin_inc
      clinical_incidence_rendering_max_ages = 100 * 365,
      
      model_seasonality = TRUE, ## Seasonality to match study site inputs [sites_13]
      g0 = sites_14$seasonal_a0[dat_row],
      g = c(sites_14$seasonal_a1[dat_row], sites_14$seasonal_a2[dat_row], sites_14$seasonal_a3[dat_row]),
      h = c(sites_14$seasonal_b1[dat_row], sites_14$seasonal_b2[dat_row], sites_14$seasonal_b3[dat_row]),
      
      Q0 = sites_14$gamb_ss_Q0[dat_row],
      individual_mosquitoes = FALSE ## Update next
    )
  )
  
  simparams <- set_equilibrium(simparams, starting_EIR)
  
  # set species
  simparams <- set_species(simparams,
                           species=list(gamb_params, arab_params, fun_params),
                           proportions=c(sites_14$prop_gamb_ss[dat_row],
                                         sites_14$prop_arab[dat_row],
                                         sites_14$prop_fun[dat_row]))
  # set treatment
  simparams <- set_drugs(simparams, list(AL_params, SP_AQ_params, DHA_PQP_params))
  simparams <- set_clinical_treatment(simparams, 
                                      drug=1,
                                      time=c(100),
                                      coverage=c(sites_14$drug_cov_0_0[dat_row]))
  simparams <- set_clinical_treatment(simparams, 
                                      drug=2,
                                      time=c(100),
                                      coverage=c(sites_14$drug_cov_1_0[dat_row]))
  
  
  ## Set up bed nets
  
  bednetparams <- simparams
  
  ## as done
  bednet_events = data.frame(
    timestep = c(0, 3, 6) * year + c(196, 
                                     196,##july as an average 2017
                                     196),
    name=c("background", 
           "trial_nets",
           "2020_nets")
  )
  
  
  # pyrethroid nets (rows 1:1000)
  bednetparams_1 <- set_bednets(
    bednetparams,
    
    timesteps = bednet_events$timestep,
    
    coverages = c(sites_14$itn_cov20[dat_row]*2, ## historic prior to RCT - it is 40% the year prior... so assuming 3 year decay in adherrence
                  params_14$ITN_1[dat_row], ## historic during RCT (on deployment) 
                  params_14$ITN_1[dat_row]),   ## planned for 2020 - ** Assuming the distribution coverage matched the RCT estimate
    
    retention = params_14$itn_leave_dur[dat_row] * year, ## Keeping this as it was observed during RCT
    ## 70% mortality 2017
    ## 48% mortality 2020 onwards
    dn0 = matrix(c(dat_res_pyr[which(dat_res_pyr$bioassay_mortality == 0.7),3],
                   dat_res_pyr[which(dat_res_pyr$bioassay_mortality == 0.48),3],
                   dat_res_pyr[which(dat_res_pyr$bioassay_mortality == 0.48),3]), nrow=3, ncol=3),
    rn = matrix(c(dat_res_pyr[which(dat_res_pyr$bioassay_mortality == 0.7),6],
                  dat_res_pyr[which(dat_res_pyr$bioassay_mortality == 0.48),6],
                  dat_res_pyr[which(dat_res_pyr$bioassay_mortality == 0.48),6]), nrow=3, ncol=3),
    rnm = matrix(c(.24, .24, .24), nrow=3, ncol=3),
    gamman = c(dat_res_pyr[which(dat_res_pyr$resistance == 0.7),9],
               dat_res_pyr[which(dat_res_pyr$resistance == 0.48),9],
               dat_res_pyr[which(dat_res_pyr$resistance == 0.48),9]) * 365
  )
  
  # pyrethroid PBO nets (rows 1:1000)
  bednetparams_2 <- set_bednets(
    bednetparams,
    
    timesteps = bednet_events$timestep,
    
    coverages = c(sites_14$itn_cov20[dat_row+1000]*2, ## historic prior to RCT - it is 40% the year prior... so assuming 3 year decay in adherrence
                  params_14$ITN_1[dat_row+1000], ## historic during RCT (on deployment) 
                  params_14$ITN_1[dat_row+1000]),   ## planned for 2020 - ** Assuming the distribution coverage matched the RCT estimate
    
    retention = params_14$itn_leave_dur[dat_row+1000] * year, ## Keeping this as it was observed during RCT
    
    ## 70% mortality 2017
    ## 48% mortality 2020 onwards
    dn0 = matrix(c(dat_res_pyr[which(dat_res_pyr$bioassay_mortality == 0.7),3],
                   dat_res_pbo[which(dat_res_pbo$bioassay_mortality == 0.48),3],
                   dat_res_pbo[which(dat_res_pbo$bioassay_mortality == 0.48),3]), nrow=3, ncol=3),
    rn = matrix(c(dat_res_pyr[which(dat_res_pyr$bioassay_mortality == 0.7),6],
                  dat_res_pbo[which(dat_res_pbo$bioassay_mortality == 0.48),6],
                  dat_res_pbo[which(dat_res_pbo$bioassay_mortality == 0.48),6]), nrow=3, ncol=3),
    rnm = matrix(c(.24, .24, .24), nrow=3, ncol=3),
    gamman = c(dat_res_pyr[which(dat_res_pyr$resistance == 0.7),9],
               dat_res_pyr[which(dat_res_pbo$resistance == 0.48),9],
               dat_res_pyr[which(dat_res_pbo$resistance == 0.48),9]) * 365
  )
  
  
  ## assume the same people are getting nets each round
  correlationsb1 <- get_correlation_parameters(bednetparams_1)
  correlationsb2 <- get_correlation_parameters(bednetparams_2)
  correlationsb1$inter_round_rho('bednets', 1)
  correlationsb2$inter_round_rho('bednets', 1)
  
  
  
  ## Run the simulations
  output1 <- run_simulation(sim_length, bednetparams_1,correlationsb1)
  output2 <- run_simulation(sim_length, bednetparams_2,correlationsb2)
  
  output1$pv_730_3650 = output1$n_detect_730_3650/output1$n_730_3650
  output2$pv_730_3650 = output2$n_detect_730_3650/output2$n_730_3650
  
  
  return(data.frame(timestep = output1$timestep,
                    
                    pyr = output1$pv_730_3650,
                    pbo = output2$pv_730_3650
  ))
}

rand_val = sample(1:1000,10,replace=FALSE)

dat = list()

for(i in 1:3){
  dat[[i]] = malsim_actual_f(rand_val[i])
}

par(mfrow=c(1,2))
plot(dat[[1]][,2] ~ dat[[1]][,1],
     xlim=c(3,5.5)*365,
     ylim=c(0,0.5),type="l",
     ylab = "Prevalence 6 months to 14 years",
     col="white",
     yaxt="n",xaxt="n",
     xlab = "Time in months")
axis(1, at =c(3,4,5)*365,labels=c("Jan 2017","Jan 2018","Jan 2019"))
axis(2, las=2,at=seq(0,1,0.2),labels=seq(0,100,20))

for(i in 1:3){
  lines(dat[[i]][,2] ~ dat[[i]][,1],col=adegenet::transp("black",0.4))  
  lines(dat[[i]][,3] ~ dat[[i]][,1],col=adegenet::transp("purple",0.7))  
}

Tqo_arm1 = c(0.193,0.145,0.13,0.14)
Tqo_arm2 = c(0.192,0.107,0.106,0.118)
prev_measured_time = c(3*year+157, c(4 * year + c(15, ## 6 months later (jan 2018)
                                 196,## 12 months later (july 2018)
                                 380)))## 18 months later (jan 2019)

abline(v=3*year+196,lty=2,col="darkgrey")

Tqo_arm_min = c(0.099,0.035)
Tqo_arm_max = c(0.416,0.401)
segments(x0=c(prev_measured_time[1],prev_measured_time[1]+14),
         x1=c(prev_measured_time[1],prev_measured_time[1]+14),
         y0 = Tqo_arm_min,
         y1 = Tqo_arm_max,col=c("darkgrey","purple"),lwd=2)

## Baseline all 4 arms
## 23.7% (12.8 ? 42)	25.1% (7.2 ? 41.7)	15.9% (1.5 ? 37.5)	8.2% (2.8 ? 34.1)

# abline(v=0.677,lty=2)

points(Tqo_arm1 ~ prev_measured_time,pch=19,cex=1,col="darkgrey")
points(Tqo_arm2 ~ c(prev_measured_time+14),pch=19,cex=1,col="purple")

obs_uganda = c(Tqo_arm1[2:4],Tqo_arm2[2:4])

time_uganda = rep(prev_measured_time[2:4],2)

preds_uganda = c(dat[[1]][c(prev_measured_time[2:4]),2],
                 dat[[1]][c(prev_measured_time[2:4]),3])

plot(obs_uganda~preds_uganda,ylim=c(0,0.4),xlim=c(0,0.4),
     ylab = "Empirical data prevalence (%)",
     xlab = "Model simulated prevalence (%)",
     col = rep(c("grey","purple"),each=3),
     pch=19)
abline(a=0,b=1,lty=2,col="grey")

legend("topleft",title = "Intervention",
       legend = c("Pyrethroid-only nets",
                  "Pyrethroid-PBO nets"),
       col=c("grey","purple"),pch=19,bty="n")
## dim 1000,550
