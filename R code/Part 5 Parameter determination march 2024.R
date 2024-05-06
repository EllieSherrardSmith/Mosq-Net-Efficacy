#############################################################
##
## Final script to determine parameters for net efficacy
##
##############################################################

###########################
## Limitations
## DECISION: to assume all nets are working equivalently
##           when it comes to expt hut associations
##           Use median estimates for these from all data
##           Include uncertainty for the mortality associations
##           That is bioassay and hut
##           And added mortality from new nets
##
##           In the absence of durability data assume 
##           half life equivalent to pyrethroid-only nets
#####################################

####################################################
##

## Critical is to keep, for each fit
## the same row of parameters for that simulation from stan models
setwd("C:/Users/esherrar/Documents/Rprojects/Mosq-Net-Efficacy")

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## keep rows associated from the bayes posterior draws
data_picker = sample(1:4000,size = 1000,replace=TRUE)

## draw from the posterior distribution with the respective inputs

## Part 1 Susceptibility assay as proxy for resistance

# Global set

ll_1 = readRDS("stan model outputs/log_logistic_fit.rds")
# ll_1 = readRDS("stan model outputs/log_logistic_fit_origcheck.RDS")
LL_fit <- rstan::extract(ll_1, permuted = TRUE)

#function shape:
x = seq(0,1,length=101) ## this is mosquito survival 

f_LOG_logistic <- function(x, b, a){
  surv = 1 - (1/(1+((1-x)/b)^(-a)))
  mort = (1-surv)
  
  return(mort)
}
## returns hut mortality

## keep uncertainty
LL_b_full <- LL_fit$b[data_picker]
LL_a_full <- LL_fit$a[data_picker]

## 
## Part 2 Benefit of pyrethroid-PBO and pyrethroid-pyrrole nets

simulated_individuals = 1e+06

benefitpbo <- readRDS("stan model outputs/ento_pbo_beta_binomial_benefit_inv_link_LancetGH2024.RDS")
pbo_bene <- rstan::extract(benefitpbo, permuted = TRUE)

## Narrow uncertainty (binomial)
fitbene_1 <- pbo_bene$alpha_1[data_picker]
fitbene_2 <- pbo_bene$alpha_2[data_picker]

mt_samples2 <- (pbo_bene$d_m * 1.0) / (simulated_individuals * 1.0)
mt_sorted2 <- apply(mt_samples2, 2, sort)
N_samples2 <- dim(mt_sorted2)[1]

testpbo = expand.grid(id = 1:101)
rr = seq(0.05, 0.95, length = 1000)
for(i in 1:1000){
  LB_ID <- round(N_samples2*rr[i])
  testpbo[,i] <- mt_sorted2[round(LB_ID,0),]
  
}

a1_test = a2_test = numeric(4000)
rand_draws_betabin = sample(1:1000,size=4000,replace=TRUE)
for(i in 1:4000){
  # hut_mort_LLPBOnew <- mt_samples2[rand_draws_betabin[i],]
  lo <- loess(testpbo[,rand_draws_betabin[i]]~seq(0,1,length=101))
  hut_mort_LLPBOnew <- predict(lo,se = TRUE)
  new_count = round(100*hut_mort_LLPBOnew$fit,0)
  new_count_inverse = 100 - new_count
  
  new_count_r = ifelse(new_count > 100,100,
                       ifelse(new_count < 0, 0,new_count))
  new_count_inverse_r = 100 - new_count_r
  
  estPBO <- cbind(new_count_r,
                  new_count_inverse_r) ~ seq(0,1,length=101) 
  glm_2 <- glm(estPBO, family = binomial())
  role.fitted2 <- predict(glm_2, se.fit = TRUE, type = "response")
  a1_test[i] = summary(glm_2)$coeff[1,1]
  a2_test[i] = summary(glm_2)$coeff[2,1]
} 

## full uncertainty (beta-binomial)
fitbene_1betabin <- c(a1_test)
fitbene_2betabin <- c(a2_test)

## Benefit of pyrethroid-pyrrole nets
## using all data and 72 hour mortality
benefit_pp <- readRDS("stan model outputs/ento_g2_beta_binomial_benefit_inv_link_LancetGH2024.RDS") ## ALL

pp_bene <- rstan::extract(benefit_pp, permuted = TRUE)

## Narrow uncertainty (binomial)
fitbene_1a <- pp_bene$alpha_1[data_picker]
fitbene_2a <- pp_bene$alpha_2[data_picker]

mt_samples <- (pp_bene$d_m * 1.0) / (simulated_individuals * 1.0)
mt_sorted <- apply(mt_samples, 2, sort)
N_samples <- dim(mt_sorted)[1]

testpp = expand.grid(id = 1:101)
rr = seq(0.05, 0.95, length = 1000)
for(i in 1:1000){
  LB_ID <- round(N_samples*rr[i])
  testpp[,i] <- mt_sorted[round(LB_ID,0),]
  
}

a1_test = a2_test = numeric(4000)
rand_draws_betabin = sample(1:1000,size=4000,replace=TRUE)
for(i in 1:4000){
  # hut_mort_LLPPnew <- mt_samples[rand_draws_betabin[i],]
  lo <- loess(testpp[,rand_draws_betabin[i]]~seq(0,1,length=101))
  hut_mort_LLPPnew <- predict(lo,se = TRUE)
  new_count = round(100*hut_mort_LLPPnew$fit,0)
  new_count_inverse = 100 - new_count
  
  new_count_r = ifelse(new_count > 100,100,
                       ifelse(new_count < 0, 0,new_count))
  new_count_inverse_r = 100 - new_count_r
  
  estPP <- cbind(new_count_r,
                 new_count_inverse_r) ~ seq(0,1,length=101) 
  glm_2 <- glm(estPP, family = binomial())
  role.fitted2 <- predict(glm_2, se.fit = TRUE, type = "response")
  a1_test[i] = summary(glm_2)$coeff[1,1]
  a2_test[i] = summary(glm_2)$coeff[2,1]
}

fitbene_1a_betabin <- c(a1_test)
fitbene_2a_betabin <- c(a2_test)

## Part 3 MORTALITY TO DETERRENCE
# load fit - just using All nets
# fit1_a <- readRDS("stan model outputs/April_2022_ento_deterrence_AllRec.rds")
# fit1_a_fit <- rstan::extract(fit1_a, permuted = TRUE)
# 
# dt_det = data.frame(c = fit1_a_fit$c,d = fit1_a_fit$d,e = fit1_a_fit$e)
# 
# saveRDS(dt_det,"stan model outputs/April_2022_ento_deterrence_AllRec_extract.rds")
fit1_a_fit = readRDS("stan model outputs/April_2022_ento_deterrence_AllRec_extract.rds")

## All nets (WHO recommended as per Okumu & Finda 2021) 
fit1_a_c <- fit1_a_fit$c[c(data_picker)]
fit1_a_d <- fit1_a_fit$d[c(data_picker)]
fit1_a_e <- fit1_a_fit$e[c(data_picker)]

fit1_a_med_c = median(fit1_a_fit$c)
fit1_a_med_d = median(fit1_a_fit$d)
fit1_a_med_e = median(fit1_a_fit$e)


##  Part 3 MORTALITY TO feed
# fit3_a <- readRDS("stan model outputs/April_2022_ento_feeding_0washes_AllRec.rds")
# fit3_a_fit <- rstan::extract(fit3_a, permuted = TRUE)
# 
# dt_fed = data.frame(a = fit3_a_fit$a, b = fit3_a_fit$b)
# saveRDS(dt_fed,"stan model outputs/April_2022_ento_feeding_0washes_AllRec_extract.rds")
fit3_a_fit = readRDS("stan model outputs/April_2022_ento_feeding_0washes_AllRec_extract.rds")

fit3_a_f <- median(fit3_a_fit$a)
fit3_a_g <- median(fit3_a_fit$b)

## assuming a log-logistic association between susc bioassay and hut surv
## beta binomial for mortality
## median estimates for feeding and deterrence

resistance_ITN_default_params_2_f = function(product, 
                                             data_picker_rand){ ## random draw from posterior pred  
  
  ## PARAMETERS
  
  # give the uncertainty for the log-logistic function
  #Assay to hut mortality conversion - median estimates	
  param_b = LL_b_full[data_picker_rand] 
  param_a = LL_a_full[data_picker_rand] 
  
  # x = res
  
  f_LOG_logistic <- function(x, b, a){
    surv = 1 - (1/(1+((1-x)/b)^(-a)))
    mort = (1-surv)
    
    return(mort)
  }
  
  # mort = 1 - res ## added benefit is relative to mortality in pyr-only nets
  
  f_LOG_logistic_with_benefit <- function(mort, alp1, alp2){
    mort_newnet = 1 / (1 + exp(-alp1 - alp2 * mort))
    # surv = (1-mort)
    return(mort_newnet)
  }
  ## returns hut mortality
  
  hut_mort_LL <- f_LOG_logistic(x = seq(0,1,length=101),#** 
                                param_b, param_a)
  
  ## using all data for pyrethroid-pbo
  hut_mort_LLPBOtemp <- f_LOG_logistic_with_benefit(mort = hut_mort_LL,
                                                    alp1 = fitbene_1betabin[data_picker_rand],
                                                    alp2 = fitbene_2betabin[data_picker_rand])
  
  hut_mort_LLPBOtempNARROW <- f_LOG_logistic_with_benefit(mort = hut_mort_LL,
                                                alp1 = fitbene_1[data_picker_rand],
                                                alp2 = fitbene_2[data_picker_rand])

  ## using all data for pyrethroid-pyrrole
  hut_mort_LLpptemp <- f_LOG_logistic_with_benefit(mort = hut_mort_LL,#seq(0,1,length=101),
                                                   alp1 = fitbene_1a_betabin[data_picker_rand],
                                                   alp2 = fitbene_2a_betabin[data_picker_rand])
  hut_mort_LLpptempNARROW <- f_LOG_logistic_with_benefit(mort = hut_mort_LL,#seq(0,1,length=101),
                                                alp1 = fitbene_1a[data_picker_rand],
                                                alp2 = fitbene_2a[data_picker_rand])

  
  #specify whichever net is used in the RCT
  #this will determine what the mortality is in the hut trial
  mort_huta = if(product==0) hut_mort_LL else if(product==1) hut_mort_LLPBOtemp else if(product==2) hut_mort_LLpptemp 
  mort_hut = mort_huta
  
  
  hut_surv = 1 - mort_hut
  # ff, ff1, and ff2 all match up for the associations for net det and fed
  # when using the generic parameters as here
  
  ## The maximum successful feeding probability per feeding attempt 
  ## (feeding and not dying) in the absence of interventions 
  kp0=0.699 ## derived from Lines et al 1987 and Curtis et al 1990 
  
  ## Now we work through the probability steps to determine the key input parameters for the model
  ## These probability relationships are determined by Rebecca Nash, Ben and Tom see email notes above
  
  fit1_a_med_c = median(fit1_a_fit$c)
  fit1_a_med_d = median(fit1_a_fit$d)
  fit1_a_med_e = median(fit1_a_fit$e)
  
  ## This is association with hut survival
  det_hut = fit1_a_med_e * (exp(fit1_a_med_d * (1 - exp(fit1_a_med_c * hut_surv)) / fit1_a_med_c))
  
  fit3_a_f <- median(fit3_a_fit$a)
  fit3_a_g <- median(fit3_a_fit$b)
  
  ## This is association with hut survival
  suc_hut = (1 - (exp(fit3_a_g * (1 - exp(fit3_a_f * hut_surv))/fit3_a_f)))
  
  ## This is association with hut survival
  rep_hut = (1 - suc_hut - mort_hut)
  
  xx = data.frame(hut_surv,mort_hut,suc_hut,rep_hut,det_hut)
  ## Combine to estimate the 3 key probable outcomes of feeding attempts
  ## Here we adjust for those mosquitoes not entering treated huts (determined by deterrence)
  n1n0 = 1-xx$det_hut
  kp1  = n1n0*xx$suc_hut
  lp1  = n1n0*xx$mort_hut
  jp1  = n1n0*xx$rep_hut+(1-n1n0)
  
  kp1 = ifelse(kp1 > kp0,kp0,kp1) ## Capping impact so max feeding is no bigger than assumed
  # # max feeding for no interventions (kp0 = 0.699, Griffin et al 2010)
  # ## (time = 0 time steps after net implementation)
  # 
  # kp0 = 1
  # 
  r_ITN0  = (1-kp1/kp0)*(jp1/(lp1+jp1))	#probability of repeating behaviour
  d_ITN0  = (1-kp1/kp0)*(lp1/(lp1+jp1))	#probability of dying with an encounter with ITN
  s_ITN0  = 1-d_ITN0-r_ITN0             #probability of successfully feeding (surviving and feeding)
  
  # plot(r_ITN0 ~ c(1-mort),ylim=c(0,1),xlim=c(0,1),xlab = "Susc bioassay survial",type="l",col="orange") 
  
  ## The half-life of the net relative to it's capacity to kill mosquitoes
  ## with the insecticide active ingredient (a pyrethroid) when there is 
  ## no resistance in mosquitoes. 
  
  ## Repeat these to determin the maximum effect which combine to help determine ITN half life
  ## We will stick to pyr-params for half life and update in 2023 once new data are available
  mort_maxA   =  if(product == 0) f_LOG_logistic(x = 0,#** this is surv i.e. mort max when surv =0
                                                 param_b, param_a) else if(product == 1) f_LOG_logistic_with_benefit(mort = hut_mort_LL, 
                                                                                                                     alp1 = fitbene_1[data_picker_rand],
                                                                                                                     alp2 = fitbene_2[data_picker_rand]) else if(product == 2) f_LOG_logistic_with_benefit(mort = hut_mort_LL,#seq(0,1,length=101), 
                                                                                                                                                                                                           alp1 = fitbene_1a[data_picker_rand],
                                                                                                                                                                                                           alp2 = fitbene_2a[data_picker_rand])
  
  
  
  mort_max = max(mort_maxA)
  
  #Decay in insecticide pyrethroid-only net		
  # p = -2.36 rhop = -3.05 ## with binomial fits... 
  mup =	-2.36  #2.66157
  rhop = -3.05 #-4.05578 #NEW -2.581591 ###hlf_rho# #rho_p ##  sample(test$b,size=1)	##array(c(rep(-3.007,3),rep(-3.74,3),rep(-2.295,3)),c(3,3)) ## ... gam.medians[10]
  net_halflife=2.64
  
  #{halflife}
  my_max_washes   = log(2)/(1/(1+exp(-(mup +rhop*(mort_max)))))
  
  ## Uncertainty
  net_half_life_min = 2
  net_half_life_max = 3
  
  ## FOR NON PYR-ONLY NETS THIS WILL RETURN HIGH HALF LIFE
  ## WE RECOMMEND ONLY USING PY-ONLY HALF LIFE UNTIL WE CAN
  ## VALIDATE OTHER NETS 
  wash_decay_rate_a = mup +rhop*(mort_hut)
  # wash_decay_rate   = log(2)/(exp(wash_decay_rate_a)/(1+exp(wash_decay_rate_a)))
  wash_decay_rate   = log(2)/(1/(1+exp(-wash_decay_rate_a)))
  itn_half_life     = wash_decay_rate/my_max_washes*net_halflife
  itn_half_life_max     = wash_decay_rate/my_max_washes*net_half_life_max  
  itn_half_life_min     = wash_decay_rate/my_max_washes*net_half_life_min
  ##No need to re-adjust these anymore
  ##Final Parameter estimates for the transmission model
  ERG_d_ITN0 <- d_ITN0
  ERG_s_ITN0 <- s_ITN0
  ERG_r_ITN0 <- 1-ERG_d_ITN0-ERG_s_ITN0
  
  #########################
  ##
  ## now repeating with narrow binomial uncertainty for additional benefit from ngen ITN
  ##
  ##
  #specify whichever net is used in the RCT
  #this will determine what the mortality is in the hut trial
  mort_huta2 = if(product==0) hut_mort_LL else if(product==1) hut_mort_LLPBOtempNARROW else if(product==2) hut_mort_LLpptempNARROW 
  mort_hut2 = mort_huta2
  
  
  hut_surv2 = 1 - mort_hut2
  # ff, ff1, and ff2 all match up for the associations for net det and fed
  # when using the generic parameters as here
  
  ## The maximum successful feeding probability per feeding attempt 
  ## (feeding and not dying) in the absence of interventions 
  kp0=0.699 ## derived from Lines et al 1987 and Curtis et al 1990 
  
  ## Now we work through the probability steps to determine the key input parameters for the model
  ## These probability relationships are determined by Rebecca Nash, Ben and Tom see email notes above
  
  fit1_a_med_c = median(fit1_a_fit$c)
  fit1_a_med_d = median(fit1_a_fit$d)
  fit1_a_med_e = median(fit1_a_fit$e)
  
  ## This is association with hut survival
  det_hut2 = fit1_a_med_e * (exp(fit1_a_med_d * (1 - exp(fit1_a_med_c * hut_surv2)) / fit1_a_med_c))
  
  fit3_a_f <- median(fit3_a_fit$a)
  fit3_a_g <- median(fit3_a_fit$b)
  
  ## This is association with hut survival
  suc_hut2 = (1 - (exp(fit3_a_g * (1 - exp(fit3_a_f * hut_surv2))/fit3_a_f)))
  
  ## This is association with hut survival
  rep_hut2 = (1 - suc_hut2 - mort_hut2)
  
  xx2 = data.frame(hut_surv2,mort_hut2,suc_hut2,rep_hut2,det_hut2)
  ## Combine to estimate the 3 key probable outcomes of feeding attempts
  ## Here we adjust for those mosquitoes not entering treated huts (determined by deterrence)
  n1n0_2 = 1-xx2$det_hut2
  kp2  = n1n0_2*xx2$suc_hut2
  lp2  = n1n0_2*xx2$mort_hut2
  jp2  = n1n0_2*xx2$rep_hut2+(1-n1n0_2)
  
  kp2 = ifelse(kp2 > kp0,kp0,kp2) ## Capping impact so max feeding is no bigger than assumed
  # # max feeding for no interventions (kp0 = 0.699, Griffin et al 2010)
  # ## (time = 0 time steps after net implementation)
  # 
  # kp0 = 1
  # 
  r2_ITN0  = (1-kp2/kp0)*(jp2/(lp2+jp2))	#probability of repeating behaviour
  d2_ITN0  = (1-kp2/kp0)*(lp2/(lp2+jp2))	#probability of dying with an encounter with ITN
  s2_ITN0  = 1-d2_ITN0-r2_ITN0             #probability of successfully feeding (surviving and feeding)
  
  # plot(r_ITN0 ~ c(1-mort),ylim=c(0,1),xlim=c(0,1),xlab = "Susc bioassay survial",type="l",col="orange") 
  
  ## The half-life of the net relative to it's capacity to kill mosquitoes
  ## with the insecticide active ingredient (a pyrethroid) when there is 
  ## no resistance in mosquitoes. 
  
  ## Repeat these to determin the maximum effect which combine to help determine ITN half life
  ## We will stick to pyr-params for half life and update in 2023 once new data are available
  ## WE RECOMMEND ONLY USING PY-ONLY HALF LIFE UNTIL WE CAN
  ## VALIDATE OTHER NETS 
  wash_decay_rate_a2 = mup +rhop*(mort_hut2)
  # wash_decay_rate   = log(2)/(exp(wash_decay_rate_a)/(1+exp(wash_decay_rate_a)))
  wash_decay_rate2   = log(2)/(1/(1+exp(-wash_decay_rate_a2)))
  itn_half_life2     = wash_decay_rate2/my_max_washes*net_halflife
  itn_half_life_max2     = wash_decay_rate2/my_max_washes*net_half_life_max  
  itn_half_life_min2     = wash_decay_rate2/my_max_washes*net_half_life_min
  ##No need to re-adjust these anymore
  ##Final Parameter estimates for the transmission model
  ERG_d2_ITN0 <- d2_ITN0
  ERG_s2_ITN0 <- s2_ITN0
  ERG_r2_ITN0 <- 1-ERG_d2_ITN0-ERG_s2_ITN0
 
  ## Print out these estimates to a data.frame as the function output
  uncertainty_resistance_params_nets = data.frame(ERG_d_ITN0,ERG_r_ITN0,itn_half_life,
                                                  itn_half_life_min,itn_half_life_max,
                                                  ERG_d2_ITN0,ERG_r2_ITN0,itn_half_life2,
                                                  itn_half_life_min2,itn_half_life_max2,
                                                  wash_decay_rate,
                                                  det_hut=xx$det_hut,suc_hut=xx$suc_hut,mort_hut=xx$mort_hut,rep_hut=xx$rep_hut,
                                                  det_hut2=xx2$det_hut2,suc_hut2=xx2$suc_hut2,mort_hut2=xx2$mort_hut2,rep_hut2=xx2$rep_hut2,
                                                  n1n0,kp1,lp1,jp1,
                                                  n1n0_2,kp2,lp2,jp2,
                                                  bioassay_surv = seq(0,1,length=101))
  
  return(uncertainty_resistance_params_nets)
}

#################################
##
## Global estimates for all of us to use
## Now we have a process to draw any uncertainty estimate and return 
## parameters for ITNs

## We want to fairly represent the median uncertainty

vec = 1:1000
matrix_dn0 = matrix_rn0 = matrix_halflife = array(dim=c(nrow(test),length(vec)))
matrix_dn0NARROW = matrix_rn0NARROW = matrix_halflife_NARROW = array(dim=c(nrow(test),length(vec)))


for(i in 1:1000){
  test = resistance_ITN_default_params_2_f(product = 0, ## PYRETHROID ONLY LLIN 
                                           data_picker_rand = i) ## any number from 1 to 1000
  matrix_dn0[,i] = test$ERG_d_ITN0
  matrix_rn0[,i] = test$ERG_r_ITN0
  matrix_halflife[,i] = test$itn_half_life
  

  matrix_dn0NARROW[,i] = test$ERG_d2_ITN0
  matrix_rn0NARROW[,i] = test$ERG_r2_ITN0
  matrix_halflife_NARROW[,i] = test$itn_half_life2
}

data_base_1 = data.frame(id = 1:101000,
                         resistance = rep(seq(0,1,length=101),1000),
                         dn0 = round(c(matrix_dn0),4),
                         rn_pyr =  c(matrix_rn0),
                         mean_duration = c(matrix_halflife))

saveRDS(data_base_1, "parameters/pyrethroid_uncertainty_LancetGH2024")

## pyrethroid PBO

vec = 1:1000
matrix_dn0_pbo = matrix_rn0_pbo = matrix_halflife_pbo = array(dim=c(nrow(test),length(vec)))
for(i in 1:1000){
  test = resistance_ITN_default_params_2_f(product = 1, ## PYRETHROID ONLY LLIN 
                                           data_picker_rand = i) ## any number from 1 to 1000
  matrix_dn0[,i] = test$ERG_d_ITN0
  matrix_rn0[,i] = test$ERG_r_ITN0
  matrix_halflife[,i] = test$itn_half_life
  
  
  matrix_dn0NARROW[,i] = test$ERG_d2_ITN0
  matrix_rn0NARROW[,i] = test$ERG_r2_ITN0
  matrix_halflife_NARROW[,i] = test$itn_half_life2
}

## For gamman, we assume that the mortality given standard ITNs is linked to a given gamman
## in the absence of empirical data to see how mortality wanes for newer net types, 
## we make the assumption that a given mortality level is linked to the gamman
## so we need to align mortalities between standard ITNs and new to select the gamman
## for any level of resistance (as newer nets are shifted, relative to standard, toward higher mortality on this scale)
## the process below enables this:
data_base_2 = data.frame(id = 1:101000,
                         resistance = rep(seq(0,1,length=101),1000),
                         dn0 = round(c(matrix_dn0),4),
                         rn_pbo =  c(matrix_rn0),
                         mean_duration = c(matrix_halflife))

library(dplyr)
merge <- left_join(data_base_2, data_base_1, by = "dn0", relationship = "many-to-many")
DD2 = merge[order(merge$id.x), ]
DD2b = DD2[order(-DD2$dn0),]
DD2b$mn_dur = numeric(length=nrow(DD2b))
for(i in 1:nrow(DD2b)){
  DD2b$mn_dur[i] = max(c(DD2b$mean_duration.x[i],DD2b$mean_duration.y[i]),na.rm = TRUE)
  }
DD2b = distinct(DD2b, id.x, .keep_all = TRUE)
DD2b = DD2b[order(DD2b$id.x),]

saveRDS(DD2b, "parameters/pbo_uncertainty_using_pyrethroid_dn0_for_mn_durability_LancetGH2024")



## pyrethroid-pyrrole

vec = 1:1000
matrix_dn0_pyrrole = matrix_rn0_pyrrole = matrix_halflife_pyrrole = array(dim=c(nrow(test),length(vec)))
for(i in 1:1000){
  test = resistance_ITN_default_params_2_f(product = 2, ## PYRETHROID ONLY LLIN 
                                           data_picker_rand = i) ## any number from 1 to 1000
  matrix_dn0[,i] = test$ERG_d_ITN0
  matrix_rn0[,i] = test$ERG_r_ITN0
  matrix_halflife[,i] = test$itn_half_life
  
  
  matrix_dn0NARROW[,i] = test$ERG_d2_ITN0
  matrix_rn0NARROW[,i] = test$ERG_r2_ITN0
  matrix_halflife_NARROW[,i] = test$itn_half_life2
}


data_base_3 = data.frame(id = 1:101000,
                         resistance = rep(seq(0,1,length=101),1000),
                         dn0 = round(c(matrix_dn0),4),
                         rn_pbo =  c(matrix_rn0),
                         mean_duration = c(matrix_halflife))


merge <- left_join(data_base_3, data_base_1, by = "dn0", relationship = "many-to-many")
DD2 = merge[order(merge$id.x), ]
DD2b = DD2[order(-DD2$dn0),]
DD2b$mn_dur = numeric(length=nrow(DD2b))
for(i in 1:nrow(DD2b)){
  DD2b$mn_dur[i] = max(c(DD2b$mean_duration.x[i],DD2b$mean_duration.y[i]),na.rm = TRUE)
}
DD2b = distinct(DD2b, id.x, .keep_all = TRUE)
DD2b = DD2b[order(DD2b$id.x),]

saveRDS(DD2b, "parameters/pyrrole_uncertainty_using_pyrethroid_dn0_for_mn_durability_LancetGH2024")


