library(malariasimulation)
library(malariaEquilibrium)
library(reshape2)
library(ggplot2)




set.seed(20130)
###########################################################################
##
## Comparing
##  * Pyrethroid only nets
##  * Pyrethroid PBO nets
##  * Pyrethroid-pyrrole nets

## Piggy backing on set up for a previous RCT
## sites
sites_13 = read.csv("validation/data/site_file_13.csv",header=TRUE)
## specifics
params_13 = read.csv("validation/data/input_params_13_update_net_parameters2_all_hutdataB.csv",header=TRUE) 

##
##
## MODEL SET UP - for arm 1 and 3 and distinct for 2 and 4 (0.68, 0.67 versus 0.61 and 0.64)
year <- 365
month <- 30
sim_length <- 10 * year
human_population <- 10000
starting_EIR_U <- 40 

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




## Set up the independent arms of the trial

bednetparamsU <- simparamsU

bednet_events = data.frame(
  timestep = c(1, 4, 7) * year,
  name=c("Bednets start1", "redistribute","redistribute")
)


## Pyrethroid only nets
bednetparams_1_0 <- set_bednets(
  bednetparamsU,
  
  timesteps = bednet_events$timestep,
  
  coverages = c(sites_13$itn_cov20[1], ## historic 
                sites_13$itn_cov20[1], ## historic  
                params_13$ITN_1[1]),   ## TRIAL NETS
  
  retention = params_13$itn_leave_dur[1] * year,
  
  ## April 2022 update - WHO Recommended Nets alone: Comparison
  ## 0 pyrethroid resistance
  ##use parameters/pyrethroid_only_nets.csv
  dn0 = matrix(0.313552075, nrow=3, ncol=3),
  rn = matrix(0.658532942, nrow=3, ncol=3),
  rnm = matrix(c(.24, .24,.24), nrow=3, ncol=3),
  gamman = c(2.64,2.64,2.64) * 365
)

bednetparams_1_0low <- set_bednets(
  bednetparamsU,
  
  timesteps = bednet_events$timestep,
  
  coverages = c(sites_13$itn_cov20[1], ## historic 
                sites_13$itn_cov20[1], ## historic  
                params_13$ITN_1[1]),   ## TRIAL NETS
  
  retention = params_13$itn_leave_dur[1] * year,
  
  ## April 2022 update - WHO Recommended Nets alone: Comparison
  ## 0 pyrethroid resistance
  dn0 = matrix(0.28143939, nrow=3, ncol=3),
  rn = matrix(0.629340769, nrow=3, ncol=3),
  rnm = matrix(c(.24, .24,.24), nrow=3, ncol=3),
  gamman = c(2.64,2.64,2.64) * 365
)

bednetparams_1_0upp <- set_bednets(
  bednetparamsU,
  
  timesteps = bednet_events$timestep,
  
  coverages = c(sites_13$itn_cov20[1], ## historic 
                sites_13$itn_cov20[1], ## historic  
                params_13$ITN_1[1]),   ## TRIAL NETS
  
  retention = params_13$itn_leave_dur[1] * year,

  ## April 2022 update - WHO Recommended Nets alone: Comparison
  ## 0 pyrethroid resistance
  dn0 = matrix(0.350317658, nrow=3, ncol=3),
  rn = matrix(0.68227717, nrow=3, ncol=3),
  rnm = matrix(c(.24, .24,.24), nrow=3, ncol=3),
  gamman = c(2.64,2.64,2.64) * 365
)

bednetparams_1_40 <- set_bednets(
  bednetparamsU,
  
  timesteps = bednet_events$timestep,
  
  coverages = c(sites_13$itn_cov20[1], ## historic 
                sites_13$itn_cov20[1], ## historic  
                params_13$ITN_1[1]),   ## TRIAL NETS
  
  retention = params_13$itn_leave_dur[1] * year,
  
  ## April 2022 update - WHO Recommended Nets alone: Comparison
  ## 0 pyrethroid resistance
  ##use parameters/pyrethroid_only_nets.csv
  dn0 = matrix(0.294764746, nrow=3, ncol=3),
  rn = matrix(0.672645681, nrow=3, ncol=3),
  rnm = matrix(c(.24, .24,.24), nrow=3, ncol=3),
  gamman = c(2.51923897,2.51923897,2.51923897) * 365
)


bednetparams_1_80 <- set_bednets(
  bednetparamsU,
  
  timesteps = bednet_events$timestep,
  
  coverages = c(sites_13$itn_cov20[1], ## historic 
                sites_13$itn_cov20[1], ## historic  
                params_13$ITN_1[1]),   ## TRIAL NETS
  
  retention = params_13$itn_leave_dur[1] * year,
  
  ## April 2022 update - WHO Recommended Nets alone: Comparison
  ## 0 pyrethroid resistance
  ##use parameters/pyrethroid_only_nets.csv
  dn0 = matrix(0.25177413, nrow=3, ncol=3),
  rn = matrix(0.702446273, nrow=3, ncol=3),
  rnm = matrix(c(.24, .24,.24), nrow=3, ncol=3),
  gamman = c(2.313621542,2.313621542,2.313621542) * 365
)

bednetparams_1_100 <- set_bednets(
  bednetparamsU,
  
  timesteps = bednet_events$timestep,
  
  coverages = c(sites_13$itn_cov20[1], ## historic 
                sites_13$itn_cov20[1], ## historic  
                params_13$ITN_1[1]),   ## TRIAL NETS
  
  retention = params_13$itn_leave_dur[1] * year,
  
  ## April 2022 update - WHO Recommended Nets alone: Comparison
  ## 0 pyrethroid resistance
  ##use parameters/pyrethroid_only_nets.csv
  dn0 = matrix(0, nrow=3, ncol=3),
  rn = matrix(0.745631384, nrow=3, ncol=3),
  rnm = matrix(c(.24, .24,.24), nrow=3, ncol=3),
  gamman = c(1.827237287,1.827237287,1.827237287) * 365
)


## Pyrethroid PBO nets
bednetparams_2_0 <- set_bednets(
  bednetparamsU,
  
  timesteps = bednet_events$timestep,
  
  coverages = c(sites_13$itn_cov20[1], ## historic 
                sites_13$itn_cov20[1], ## historic  
                params_13$ITN_1[1001]),   ## TRIAL NETS
  
  retention = params_13$itn_leave_dur[1001] * year,

  ## April 2022 update - WHO Recommended Nets alone: Comparison
  dn0 = matrix(0.470220097, nrow=3, ncol=3),
  rn = matrix(0.523949565, nrow=3, ncol=3),
  rnm = matrix(c(.24, .24, .24), nrow=3, ncol=3),
  gamman = c(2.64,2.64,2.64) * 365 ## TRIAL NETS
)

bednetparams_2_40 <- set_bednets(
  bednetparamsU,
  
  timesteps = bednet_events$timestep,
  
  coverages = c(sites_13$itn_cov20[1], ## historic 
                sites_13$itn_cov20[1], ## historic  
                params_13$ITN_1[1001]),   ## TRIAL NETS
  
  retention = params_13$itn_leave_dur[1001] * year,
  
  ## April 2022 update - WHO Recommended Nets alone: Comparison
  dn0 = matrix(0.449544064, nrow=3, ncol=3),
  rn = matrix(0.542941503, nrow=3, ncol=3),
  rnm = matrix(c(.24, .24, .24), nrow=3, ncol=3),
  gamman = c(2.404391848,2.404391848,2.404391848) * 365 ## TRIAL NETS
)

bednetparams_2_80 <- set_bednets(
  bednetparamsU,
  
  timesteps = bednet_events$timestep,
  
  coverages = c(sites_13$itn_cov20[1], ## historic 
                sites_13$itn_cov20[1], ## historic  
                params_13$ITN_1[1001]),   ## TRIAL NETS
  
  retention = params_13$itn_leave_dur[1001] * year,
  
  ## April 2022 update - WHO Recommended Nets alone: Comparison
  dn0 = matrix(0.398090444, nrow=3, ncol=3),
  rn = matrix(0.588888977, nrow=3, ncol=3),
  rnm = matrix(c(.24, .24, .24), nrow=3, ncol=3),
  gamman = c(1.958622999,1.958622999,1.958622999) * 365 ## TRIAL NETS
)

bednetparams_2_100 <- set_bednets(
  bednetparamsU,
  
  timesteps = bednet_events$timestep,
  
  coverages = c(sites_13$itn_cov20[1], ## historic 
                sites_13$itn_cov20[1], ## historic  
                params_13$ITN_1[1001]),   ## TRIAL NETS
  
  retention = params_13$itn_leave_dur[1001] * year,
  
  ## April 2022 update - WHO Recommended Nets alone: Comparison
  dn0 = matrix(0.082935739, nrow=3, ncol=3),
  rn = matrix(0.766553383, nrow=3, ncol=3),
  rnm = matrix(c(.24, .24, .24), nrow=3, ncol=3),
  gamman = c(1.102581665,1.102581665,1.102581665) * 365 ## TRIAL NETS
)


## Pyrethroid pyrrole nets
bednetparams_3_0 <- set_bednets(
  bednetparamsU,
  
  timesteps = bednet_events$timestep,
  
  coverages = c(sites_13$itn_cov20[1], ## historic 
                sites_13$itn_cov20[1], ## historic  
                params_13$ITN_1[1001]),   ## TRIAL NETS
  
  retention = params_13$itn_leave_dur[1001] * year,
  
  ## April 2022 update - WHO Recommended Nets alone: Comparison
  dn0 = matrix(0.528938884, nrow=3, ncol=3),
  rn = matrix(0.468795373, nrow=3, ncol=3),
  rnm = matrix(c(.24, .24, .24), nrow=3, ncol=3),
  gamman = c(2.64,2.64,2.64) * 365 ## TRIAL NETS
)

bednetparams_3_40 <- set_bednets(
  bednetparamsU,
  
  timesteps = bednet_events$timestep,
  
  coverages = c(sites_13$itn_cov20[1], ## historic 
                sites_13$itn_cov20[1], ## historic  
                params_13$ITN_1[1001]),   ## TRIAL NETS
  
  retention = params_13$itn_leave_dur[1001] * year,
  
  ## April 2022 update - WHO Recommended Nets alone: Comparison
  dn0 = matrix(0.518145305, nrow=3, ncol=3),
  rn = matrix(0.479048086, nrow=3, ncol=3),
  rnm = matrix(c(.24, .24, .24), nrow=3, ncol=3),
  gamman = c(2.490349351,2.490349351,2.490349351) * 365 ## TRIAL NETS
)

bednetparams_3_80 <- set_bednets(
  bednetparamsU,
  
  timesteps = bednet_events$timestep,
  
  coverages = c(sites_13$itn_cov20[1], ## historic 
                sites_13$itn_cov20[1], ## historic  
                params_13$ITN_1[1001]),   ## TRIAL NETS
  
  retention = params_13$itn_leave_dur[1001] * year,
  
  ## April 2022 update - WHO Recommended Nets alone: Comparison
  dn0 = matrix(0.488890709, nrow=3, ncol=3),
  rn = matrix(0.506589695, nrow=3, ncol=3),
  rnm = matrix(c(.24, .24, .24), nrow=3, ncol=3),
  gamman = c(2.143357997,2.143357997,2.143357997) * 365 ## TRIAL NETS
)

bednetparams_3_100 <- set_bednets(
  bednetparamsU,
  
  timesteps = bednet_events$timestep,
  
  coverages = c(sites_13$itn_cov20[1], ## historic 
                sites_13$itn_cov20[1], ## historic  
                params_13$ITN_1[1001]),   ## TRIAL NETS
  
  retention = params_13$itn_leave_dur[1001] * year,

  ## April 2022 update - WHO Recommended Nets alone: Comparison
  dn0 = matrix(0.190141348, nrow=3, ncol=3),
  rn = matrix(0.737480169, nrow=3, ncol=3),
  rnm = matrix(c(.24, .24, .24), nrow=3, ncol=3),
  gamman = c(0.907427092,0.907427092,0.907427092) * 365 ## TRIAL NETS
)

correlationsb1_0 <- get_correlation_parameters(bednetparams_1_0)
correlationsb1_0l <- get_correlation_parameters(bednetparams_1_0low)
correlationsb1_0u <- get_correlation_parameters(bednetparams_1_0upp)
correlationsb1_40 <- get_correlation_parameters(bednetparams_1_40)
correlationsb1_80 <- get_correlation_parameters(bednetparams_1_80)
correlationsb1_100 <- get_correlation_parameters(bednetparams_1_100)

correlationsb2_0 <- get_correlation_parameters(bednetparams_2_0)
correlationsb2_40 <- get_correlation_parameters(bednetparams_2_40)
correlationsb2_80 <- get_correlation_parameters(bednetparams_2_80)
correlationsb2_100 <- get_correlation_parameters(bednetparams_2_100)

correlationsb3_0 <- get_correlation_parameters(bednetparams_3_0)
correlationsb3_40 <- get_correlation_parameters(bednetparams_3_40)
correlationsb3_80 <- get_correlation_parameters(bednetparams_3_80)
correlationsb3_100 <- get_correlation_parameters(bednetparams_3_100)

output1_0 <- run_simulation(sim_length, bednetparams_1_0,correlationsb1_0)
output1_0l <- run_simulation(sim_length, bednetparams_1_0low,correlationsb1_0l)
output1_0u <- run_simulation(sim_length, bednetparams_1_0upp,correlationsb1_0u)
output1_40 <- run_simulation(sim_length, bednetparams_1_40,correlationsb1_40)
output1_80 <- run_simulation(sim_length, bednetparams_1_80,correlationsb1_80)
output1_100 <- run_simulation(sim_length, bednetparams_1_100,correlationsb1_100)

output2_0 <- run_simulation(sim_length, bednetparams_2_0,correlationsb2_0)
output2_40 <- run_simulation(sim_length, bednetparams_2_40,correlationsb2_40)
output2_80 <- run_simulation(sim_length, bednetparams_2_80,correlationsb2_80)
output2_100 <- run_simulation(sim_length, bednetparams_2_100,correlationsb2_100)

output3_0 <- run_simulation(sim_length, bednetparams_3_0,correlationsb3_0)
output3_40 <- run_simulation(sim_length, bednetparams_3_40,correlationsb3_40)
output3_80 <- run_simulation(sim_length, bednetparams_3_80,correlationsb3_80)
output3_100 <- run_simulation(sim_length, bednetparams_3_100,correlationsb3_100)

output1_0$pv_182.5_5110 = output1_0$n_detect_182.5_5110/output1_0$n_182.5_5110
output1_0l$pv_182.5_5110 = output1_0l$n_detect_182.5_5110/output1_0l$n_182.5_5110
output1_0u$pv_182.5_5110 = output1_0u$n_detect_182.5_5110/output1_0u$n_182.5_5110
output1_40$pv_182.5_5110 = output1_40$n_detect_182.5_5110/output1_40$n_182.5_5110
output1_80$pv_182.5_5110 = output1_80$n_detect_182.5_5110/output1_80$n_182.5_5110
output1_100$pv_182.5_5110 = output1_100$n_detect_182.5_5110/output1_100$n_182.5_5110

output2_0$pv_182.5_5110 = output2_0$n_detect_182.5_5110/output2_0$n_182.5_5110
output2_40$pv_182.5_5110 = output2_40$n_detect_182.5_5110/output2_40$n_182.5_5110
output2_80$pv_182.5_5110 = output2_80$n_detect_182.5_5110/output2_80$n_182.5_5110
output2_100$pv_182.5_5110 = output2_100$n_detect_182.5_5110/output2_100$n_182.5_5110

output3_0$pv_182.5_5110 = output3_0$n_detect_182.5_5110/output3_0$n_182.5_5110
output3_40$pv_182.5_5110 = output3_40$n_detect_182.5_5110/output3_40$n_182.5_5110
output3_80$pv_182.5_5110 = output3_80$n_detect_182.5_5110/output3_80$n_182.5_5110
output3_100$pv_182.5_5110 = output3_100$n_detect_182.5_5110/output3_100$n_182.5_5110



par(mfrow=c(2,2))
plot(output1_0$pv_182.5_5110 ~ output1_0$timestep,
     xlim=c(6.8,10)*365,
     ylim=c(0,1),pch="",
     col="darkred",ylab = "Prevalence 6 months to 14 years")
polygon(c(output1_0$timestep,rev(output1_0$timestep)),
        c(output1_0l$pv_182.5_5110,rev(output1_0u$pv_182.5_5110)),border=NA,col = adegenet::transp("darkred",0.4))
lines(output1_40$pv_182.5_5110 ~ output1_40$timestep,col="red",lty=2)
lines(output1_80$pv_182.5_5110 ~ output1_80$timestep,col="darkorange",lty=3)
lines(output1_100$pv_182.5_5110 ~ output1_100$timestep,col="orange",lty=4)

plot(output2_0$pv_182.5_5110 ~ output2_0$timestep,
     xlim=c(6.8,10)*365,
     ylim=c(0,1),pch="",
     col="darkblue",ylab = "Prevalence 6 months to 14 years")
lines(output2_40$pv_182.5_5110 ~ output2_40$timestep,col="blue",lty=1)
lines(output2_80$pv_182.5_5110 ~ output2_80$timestep,col="aquamarine4",lty=1)
lines(output2_100$pv_182.5_5110 ~ output2_100$timestep,col="aquamarine",lty=1)


plot(output2_0$pv_182.5_5110 ~ output3_0$timestep,
     xlim=c(6.8,10)*365,
     ylim=c(0,1),pch="",
     col="darkgreen",ylab = "Prevalence 6 months to 14 years")
lines(output3_40$pv_182.5_5110 ~ output3_40$timestep,col="green",lty=1)
lines(output3_80$pv_182.5_5110 ~ output3_80$timestep,col="lightgreen",lty=1)
lines(output3_100$pv_182.5_5110 ~ output3_100$timestep,col="yellow",lty=1)


plot(output1_0$pv_182.5_5110 ~ output1_0$timestep,
     xlim=c(6.8,10)*365,
     ylim=c(0,1),pch="",
     col="darkred",ylab = "Prevalence 6 months to 14 years")
lines(output1_0$pv_182.5_5110 ~ output1_0$timestep,col="darkred",lty=1)
lines(output2_0$pv_182.5_5110 ~ output2_0$timestep,col="darkblue",lty=1)
lines(output3_0$pv_182.5_5110 ~ output3_0$timestep,col="darkgreen",lty=1)


lines(output1_100$pv_182.5_5110 ~ output1_100$timestep,col="darkred",lwd=3)
lines(output2_100$pv_182.5_5110 ~ output2_100$timestep,col="darkblue",lwd=3)
lines(output3_100$pv_182.5_5110 ~ output3_100$timestep,col="darkgreen",lwd=3)
