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

## Just confirm nothing spurious here
hut_mort_LL <- f_LOG_logistic(x = seq(0,1,length=101), LL_b_full[5], LL_a_full[5])
plot(hut_mort_LL ~ seq(0,1,length=101),ylim=c(0,1),
     xlab="bioassay survival"
)

## When bioassay mortality is 0, hut mortality is 0.57 ish
## When bioassay mortality is 1. hut mortlity is 1

## 
## Part 2 Benefit of pyrethroid-PBO and pyrethroid-pyrrole nets

simulated_individuals = 1e+06
benefitpbo <- readRDS("stan model outputs/ento_pbo_beta_binomial_benefit_Jan2024data.RDS")
pbo_bene <- rstan::extract(benefitpbo, permuted = TRUE)
mt_samples2 <- pbo_bene$rho / simulated_individuals

fitbene_1 <- pbo_bene$a1[data_picker]
fitbene_2 <- pbo_bene$a2[data_picker]

f_LOG_logistic_with_benefit <- function(mort, alp1, alp2){
  mort_newnet = 1 / (1 + exp(-alp1 - alp2 * mort))
  # surv = (1-mort)
  return(mort_newnet)
}

hut_mort_LLPBO <- f_LOG_logistic_with_benefit(mort = hut_mort_LL, 
                                              alp1 = fitbene_1[5],
                                              alp2 = fitbene_2[5])
## for a given resistance level (ncol),
## select a value for the mortality induced by new net
# hut_mort_LLPBO = mt_samples2[7,]
# lo <- loess(mt_samples2[7,]~seq(0,1,length=101))

plot(hut_mort_LLPBO ~ seq(0,1,length=101),
     xlab = "hut_mort_LL",
     ylab = "New net hut mortality",
     ylim=c(0,1),xlim=c(0,1))

a1_test = a2_test = 
  a1_testLO = a2_testLO = 
  a1_testUP = a2_testUP = numeric(1334)
rand_draws_betabin = sample(1:4000,size=1334,replace=FALSE)
for(i in 1:1334){
  hut_mort_LLPBOnew <- mt_samples2[rand_draws_betabin[i],]
  lo <- loess(hut_mort_LLPBOnew~seq(0,1,length=101))
  hut_mort_LLPBOnew <- predict(lo,se = TRUE)
  new_count = round(100*hut_mort_LLPBOnew$fit,0)
  new_count_inverse = 100 - new_count
  
  new_count_LO = round(100*(hut_mort_LLPBOnew$fit-1.96*hut_mort_LLPBOnew$se.fit),0)
  new_count_inverse_LO = 100 - new_count_LO

  new_count_UP = round(100*(hut_mort_LLPBOnew$fit+1.96*hut_mort_LLPBOnew$se.fit),0)
  new_count_inverse_UP = 100 - new_count_UP
  
  new_count_r = ifelse(new_count > 100,100,
                       ifelse(new_count < 0, 0,new_count))
  new_count_inverse_r = 100 - new_count_r
  
  estPBO <- cbind(new_count_r,
                  new_count_inverse_r) ~ seq(0,1,length=101) 
  glm_2 <- glm(estPBO, family = binomial())
  role.fitted2 <- predict(glm_2, se.fit = TRUE, type = "response")
  a1_test[i] = summary(glm_2)$coeff[1,1]
  a2_test[i] = summary(glm_2)$coeff[2,1]

  new_count_rLO = ifelse(new_count_LO > 100,100,
                       ifelse(new_count_LO < 0, 0,new_count_LO))
  new_count_inverse_rLO = 100 - new_count_rLO

  new_count_rUP = ifelse(new_count_UP > 100,100,
                         ifelse(new_count_UP < 0, 0,new_count_UP))
  new_count_inverse_rUP = 100 - new_count_rUP
  
  estPBOL <- cbind(new_count_rLO,
                  new_count_inverse_rLO) ~ seq(0,1,length=101) 
  glm_2L <- glm(estPBOL, family = binomial())
  role.fitted2L <- predict(glm_2L, se.fit = TRUE, type = "response")
  a1_testLO[i] = summary(glm_2L)$coeff[1,1]
  a2_testLO[i] = summary(glm_2L)$coeff[2,1]
  
  estPBOU <- cbind(new_count_rUP,
                   new_count_inverse_rUP) ~ seq(0,1,length=101) 
  glm_2U <- glm(estPBOU, family = binomial())
  role.fitted2U <- predict(glm_2U, se.fit = TRUE, type = "response")
  a1_testUP[i] = summary(glm_2U)$coeff[1,1]
  a2_testUP[i] = summary(glm_2U)$coeff[2,1]
  
}
fitbene_1betabin <- c(a1_test,a1_testLO,a1_testUP)
fitbene_2betabin <- c(a2_test,a2_testLO,a2_testUP)

length(fitbene_1betabin)
length(fitbene_2betabin)

hut_mort_LLPBOnew = hut_mort_LLPBOold = expand.grid(mort = seq(0,1,length=101))
for(i in sample(1:4000,40,replace = FALSE)){
  hut_mort_LLPBOold <- f_LOG_logistic_with_benefit(mort = hut_mort_LL,
                                                   alp1 = fitbene_1[i],
                                                   alp2 = fitbene_2[i])
  
  hut_mort_LLPBOnew <- f_LOG_logistic_with_benefit(mort = hut_mort_LL,
                                                        alp1 = fitbene_1betabin[i],
                                                        alp2 = fitbene_2betabin[i])
  lines(hut_mort_LLPBOnew ~ seq(0,1,length=101),col = adegenet::transp("blue",0.4))
  lines(hut_mort_LLPBOold ~ seq(0,1,length=101),col = adegenet::transp("purple",0.4))
}

## Benefit of pyrethroid-pyrrole nets
## using all data and 72 hour mortality
benefit_pp <- readRDS("stan model outputs/ento_g2_beta_binomial_benefit_Jan2024data.RDS") ## ALL

pp_bene <- rstan::extract(benefit_pp, permuted = TRUE)
mt_samples <- pp_bene$rho / simulated_individuals

fitbene_1a <- pp_bene$a1[data_picker]
fitbene_2a <- pp_bene$a2[data_picker]

hut_mort_LLPP <- f_LOG_logistic_with_benefit(mort = hut_mort_LL,
                                             alp1 = fitbene_1a[5],
                                             alp2 = fitbene_2a[5])

## Benefit of new nets against the hut mortality
lines(hut_mort_LLPP ~ seq(0,1,length=101),col="aquamarine3", lwd=2)

a1_test = a2_test = 
  a1_testLO = a2_testLO = 
  a1_testUP = a2_testUP = numeric(1334)
rand_draws_betabin = sample(1:4000,size=1334,replace=FALSE)
for(i in 1:1334){
  hut_mort_LLPPnew <- mt_samples[rand_draws_betabin[i],]
  lo <- loess(hut_mort_LLPPnew~seq(0,1,length=101))
  hut_mort_LLPPnew <- predict(lo,se = TRUE)
  new_count = round(100*hut_mort_LLPPnew$fit,0)
  new_count_inverse = 100 - new_count
  
  new_count_LO = round(100*(hut_mort_LLPPnew$fit-1.96*hut_mort_LLPPnew$se.fit),0)
  new_count_inverse_LO = 100 - new_count_LO
  
  new_count_UP = round(100*(hut_mort_LLPPnew$fit+1.96*hut_mort_LLPPnew$se.fit),0)
  new_count_inverse_UP = 100 - new_count_UP
  
  new_count_r = ifelse(new_count > 100,100,
                       ifelse(new_count < 0, 0,new_count))
  new_count_inverse_r = 100 - new_count_r
  
  estPP <- cbind(new_count_r,
                  new_count_inverse_r) ~ seq(0,1,length=101) 
  glm_2 <- glm(estPP, family = binomial())
  role.fitted2 <- predict(glm_2, se.fit = TRUE, type = "response")
  a1_test[i] = summary(glm_2)$coeff[1,1]
  a2_test[i] = summary(glm_2)$coeff[2,1]
  
  new_count_rLO = ifelse(new_count_LO > 100,100,
                         ifelse(new_count_LO < 0, 0,new_count_LO))
  new_count_inverse_rLO = 100 - new_count_rLO
  
  new_count_rUP = ifelse(new_count_UP > 100,100,
                         ifelse(new_count_UP < 0, 0,new_count_UP))
  new_count_inverse_rUP = 100 - new_count_rUP
  
  estPPL <- cbind(new_count_rLO,
                   new_count_inverse_rLO) ~ seq(0,1,length=101) 
  glm_2L <- glm(estPPL, family = binomial())
  role.fitted2L <- predict(glm_2L, se.fit = TRUE, type = "response")
  a1_testLO[i] = summary(glm_2L)$coeff[1,1]
  a2_testLO[i] = summary(glm_2L)$coeff[2,1]
  
  estPPU <- cbind(new_count_rUP,
                   new_count_inverse_rUP) ~ seq(0,1,length=101) 
  glm_2U <- glm(estPPU, family = binomial())
  role.fitted2U <- predict(glm_2U, se.fit = TRUE, type = "response")
  a1_testUP[i] = summary(glm_2U)$coeff[1,1]
  a2_testUP[i] = summary(glm_2U)$coeff[2,1]
  
}
fitbene_1a_betabin <- c(a1_test,a1_testLO,a1_testUP)
fitbene_2a_betabin <- c(a2_test,a2_testLO,a2_testUP)

hut_mort_LLPPnew = hut_mort_LLPPold = expand.grid(mort = seq(0,1,length=101))
for(i in sample(1:4000,40,replace = FALSE)){
  hut_mort_LLPPold <- f_LOG_logistic_with_benefit(mort = hut_mort_LL,
                                                   alp1 = fitbene_1a[i],
                                                   alp2 = fitbene_2a[i])
  
  hut_mort_LLPPnew <- f_LOG_logistic_with_benefit(mort = hut_mort_LL,
                                                   alp1 = fitbene_1a_betabin[i],
                                                   alp2 = fitbene_2a_betabin[i])
  lines(hut_mort_LLPPnew ~ seq(0,1,length=101),col = adegenet::transp("darkgreen",0.4))
  lines(hut_mort_LLPPold ~ seq(0,1,length=101),col = adegenet::transp("aquamarine",0.4))
}

## ONLY USING MEDIAN ESTIMATES FROM FULL FITS
## i.e. all data from WHO recommended nets
## and assuming these associations hold for all nets


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

quantile(fit1_a_fit$c,c(0.1,0.9))
which(fit1_a_fit$c > 2.08 & fit1_a_fit$c < 2.087)

mort = seq(0,1,length=101)
deterrence = fit1_a_med_e * (exp(fit1_a_med_d * (1 - exp(fit1_a_med_c * mort)) / fit1_a_med_c))

plot(deterrence ~ mort,xaxt="n",yaxt="n",type="l",
     xlab = "Hut mortality of any WHO Recommended net (%)",
     ylab = "Deterrence due to presence of that net (%)",
     ylim=c(0,1),xlim=c(0,1))
axis(1,at=seq(0,1,0.2),labels=seq(0,100,20))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))

## 
##  Part 3 MORTALITY TO feed
# fit3_a <- readRDS("stan model outputs/April_2022_ento_feeding_0washes_AllRec.rds")
# fit3_a_fit <- rstan::extract(fit3_a, permuted = TRUE)
# 
# dt_fed = data.frame(a = fit3_a_fit$a, b = fit3_a_fit$b)
# saveRDS(dt_fed,"stan model outputs/April_2022_ento_feeding_0washes_AllRec_extract.rds")
fit3_a_fit = readRDS("stan model outputs/April_2022_ento_feeding_0washes_AllRec_extract.rds")

fit3_a_f <- median(fit3_a_fit$a)
fit3_a_g <- median(fit3_a_fit$b)

feeding = (1 - (exp(fit3_a_g * (1 - exp(fit3_a_f * mort))/fit3_a_f)))

plot(feeding ~ mort,xaxt="n",yaxt="n",type="l",
     xlab = "Hut mortality of any WHO Recommended net (%)",
     ylab = "Success feeding in presence of that net (%)",
     ylim=c(0,1),xlim=c(0,1))
axis(1,at=seq(0,1,0.2),labels=seq(0,100,20))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))

## Update in 2023 once data are available
## Half life we could stick to the original - but this is then
## assuming a logistic association between susc bioassay and hut surv

resistance_ITN_default_params_2_f = function(product, 
                                             data_picker_rand){ ## random draw from posterior pred  
  
  ## The parameters included here are from the statistical analysis 
  ## following on from Nash et al 2021 
  
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
  
  # hut_mort_LLPBOtemp <- f_LOG_logistic_with_benefit(mort = hut_mort_LL,
  #                                               alp1 = fitbene_1[data_picker_rand],
  #                                               alp2 = fitbene_2[data_picker_rand])
  # 
  ## using all data for pyrethroid-pyrrole
  hut_mort_LLpptemp <- f_LOG_logistic_with_benefit(mort = hut_mort_LL,#seq(0,1,length=101),
                                                   alp1 = fitbene_1a_betabin[data_picker_rand],
                                                   alp2 = fitbene_2a_betabin[data_picker_rand])
  # hut_mort_LLpptemp <- f_LOG_logistic_with_benefit(mort = hut_mort_LL,#seq(0,1,length=101),
  #                                               alp1 = fitbene_1a[data_picker_rand],
  #                                               alp2 = fitbene_2a[data_picker_rand])
  # 
  
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
  mup =	2.66157
  rhop =-4.05578 #NEW -2.581591 ###hlf_rho# #rho_p ##  sample(test$b,size=1)	##array(c(rep(-3.007,3),rep(-3.74,3),rep(-2.295,3)),c(3,3)) ## ... gam.medians[10]
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
  
  ## Print out these estimates to a data.frame as the function output
  uncertainty_resistance_params_nets = data.frame(ERG_d_ITN0,ERG_r_ITN0,itn_half_life,
                                                  itn_half_life_min,itn_half_life_max,
                                                  det_hut=xx$det_hut,suc_hut=xx$suc_hut,mort_hut=xx$mort_hut,rep_hut=xx$rep_hut,
                                                  n1n0,kp1,lp1,jp1,
                                                  bioassay_surv = seq(0,1,length=101))
  
  return(uncertainty_resistance_params_nets)
}

test = resistance_ITN_default_params_2_f(product = 0, ## pyrethroid only nets 
                                         data_picker_rand = 45) ## any number from 1 to 1000
head(test)
tail(test)

test = resistance_ITN_default_params_2_f(product = 1, ## pyrethroid PBO nets 
                                         data_picker_rand = 45) ## any number from 1 to 1000
head(test)
tail(test)


test = resistance_ITN_default_params_2_f(product = 2, ## PYRETHROID pyrrole LLIN 
                                         data_picker_rand = 45) ## any number from 1 to 1000
head(test)
tail(test)
########################
## To draw the comparison plots


rand_draws = function(product,##0 for pyrethroid only, 1 for PBO, 2 for G2, 3 for G2 without new data
                      rand_draw, ## any value betweem 1 and 1000
                      type_of_net){ ## "name of net type"
  
  test = resistance_ITN_default_params_2_f(product = product, ## PYRETHROID ONLY LLIN 
                                           data_picker_rand = rand_draw) ## any number from 1 to 1000
  
  ## Standard nets
  printed = type_of_net#"Pyrethroid-only nets"
  ## inputs
  
  my.cols<-c(RColorBrewer::brewer.pal(12, "Paired"),"white",
             RColorBrewer::brewer.pal(8, "Dark2"))
  
  n1n0 = 1-test$det_hut
  kp1  = test$n1n0*test$suc_hut
  lp1  = test$n1n0*test$mort_hut
  jp1  = test$n1n0*test$rep_hut+(1-test$n1n0)
  all_mosq = test$n1n0 + test$kp1 + test$jp1 + test$lp1
  
  my.success.kill = #
    100 - (lp1*100)
  # my.success.kill = # this is everything but the mortality
  
  my.kill.det = #
    100 - (lp1*100) - (test$det_hut*100)
  # my.kill.det = # this is all but kill and deter
  
  my.success.a = # this is just fed
    100*kp1
  # my.success.a = # this is just fed
  
  # par(mfrow=c(2,2))
  my.x100 = seq(100,0,length=101) ## just 1 big square
  plot(my.x100,my.x100,type="n",main=printed,cex.lab=1,cex.axis=1,cex.main=1,
       xlim=c(0,100),ylim=c(0,100),las=1,
       xlab="Bioassay (% survival at susceptibility test)",
       ylab="Probable outcome of feeding attempt (%)")
  polygon(c(my.x100,rev(my.x100)),c(rep(0,length(my.x100)),rep(100,length(my.x100))),col=my.cols[1],border=my.cols[1])
  polygon(c(my.x100,rev(my.x100)),c(rep(0,length(my.x100)),my.success.kill),col=my.cols[3],border=my.cols[3])
  polygon(c(my.x100,rev(my.x100)),c(rep(0,length(my.x100)),my.kill.det),col=my.cols[7],border=my.cols[7])
  polygon(c(my.x100,rev(my.x100)),c(rep(0,length(my.x100)),my.success.a),col=my.cols[5],border=my.cols[5])
  
  
  
}

par(mfrow=c(2,4))
for(i in sample(1:1000,size=8,replace=FALSE)){
  rand_draws(product = 0,##0 for pyrethroid only, 1 for PBO, 2 for G2, 3 for G2 without new data
             rand_draw = i, ## any value betweeN 1 and 1000
             type_of_net = "Pyrethroid-only nets") ## "name of net type"
  print(i)  
}
  

for(i in sample(1:1000,size=8,replace=FALSE)){
  rand_draws(product = 1,##0 for pyrethroid only, 1 for PBO, 2 for G2, 3 for G2 without new data
             rand_draw = i, ## any value betweeN 1 and 1000
             type_of_net = "Pyrethroid PBO nets") ## "name of net type"
  print(i)
}

for(i in sample(1:1000,size=8,replace=FALSE)){
  rand_draws(product = 2,##0 for pyrethroid only, 1 for PBO, 2 for G2, 3 for G2 without new data
             rand_draw = i, ## any value betweeN 1 and 1000
             type_of_net = "Pyrethroid-pyrrole nets") ## "name of net type"
  print(i)
}

#################################
##
## Global estimates for all of us to use

vec = 1:1000
matrix_dn0 = matrix_rn0 = matrix_halflife = array(dim=c(nrow(test),length(vec)))
for(i in 1:1000){
  test = resistance_ITN_default_params_2_f(product = 0, ## PYRETHROID ONLY LLIN 
                                           data_picker_rand = i) ## any number from 1 to 1000
  matrix_dn0[,i] = test$ERG_d_ITN0
  matrix_rn0[,i] = test$ERG_r_ITN0
  matrix_halflife[,i] = test$itn_half_life
  
}

dn0MEAN = rowMeans(matrix_dn0)
rn0MEAN = rowMeans(matrix_rn0)
halflifeMEAN = rowMeans(matrix_halflife)

dn0 = rn0 = gamman = array(dim=c(nrow(test),3))
for(j in 1:nrow(test)){
  dn0[j,] = c(as.numeric(quantile(matrix_dn0[j,],c(0.05,0.5,0.95))))
  rn0[j,] = c(as.numeric(quantile(matrix_rn0[j,],c(0.05,0.5,0.95))))
  gamman[j,] = c(as.numeric(quantile(matrix_halflife[j,],c(0.05,0.5,0.95))))
}

pyrethroidOnlyNets = data.frame(dn0_lo10 = dn0[,1],dn0_med = dn0[,2],dn0_up90 = dn0[,3],
                                rn0_lo10 = rn0[,1],rn0_med = rn0[,2],rn0_up90 = rn0[,3],
                                gamman_lo10 = gamman[,1],gamman_med = gamman[,2],gamman_up90 = gamman[,3])
pyrethroidOnlyNets$bioassay_mortality = seq(1,0,length=101)
head(pyrethroidOnlyNets)

## pyrethroid PBO

vec = 1:1000
matrix_dn0_pbo = matrix_rn0_pbo = matrix_halflife_pbo = array(dim=c(nrow(test),length(vec)))
for(i in 1:1000){
  test = resistance_ITN_default_params_2_f(product = 1, ## PYRETHROID ONLY LLIN 
                                           data_picker_rand = i) ## any number from 1 to 1000
  matrix_dn0_pbo[,i] = test$ERG_d_ITN0
  matrix_rn0_pbo[,i] = test$ERG_r_ITN0
  matrix_halflife_pbo[,i] = test$itn_half_life
  
}

dn0MEAN_pbo = rowMeans(matrix_dn0_pbo)
rn0MEAN_pbo = rowMeans(matrix_rn0_pbo)
halflifeMEAN_pbo = rowMeans(matrix_halflife_pbo)

dn0 = rn0 = gamman = array(dim=c(nrow(test),3))
for(j in 1:nrow(test)){
  dn0[j,] = c(as.numeric(quantile(matrix_dn0_pbo[j,],c(0.05,0.5,0.95))))
  rn0[j,] = c(as.numeric(quantile(matrix_rn0_pbo[j,],c(0.05,0.5,0.95))))
  gamman[j,] = c(as.numeric(quantile(matrix_halflife_pbo[j,],c(0.05,0.5,0.95))))
}

pyrethroidPBONets = data.frame(dn0_lo10 = dn0[,1],dn0_med = dn0[,2],dn0_up90 = dn0[,3],
                                rn0_lo10 = rn0[,1],rn0_med = rn0[,2],rn0_up90 = rn0[,3],
                                gamman_lo10 = gamman[,1],gamman_med = gamman[,2],gamman_up90 = gamman[,3])
pyrethroidPBONets$bioassay_mortality = seq(1,0,length=101)
head(pyrethroidPBONets)


## pyrethroid-pyrrole

vec = 1:1000
matrix_dn0_pyrrole = matrix_rn0_pyrrole = matrix_halflife_pyrrole = array(dim=c(nrow(test),length(vec)))
for(i in 1:1000){
  test = resistance_ITN_default_params_2_f(product = 2, ## PYRETHROID ONLY LLIN 
                                           data_picker_rand = i) ## any number from 1 to 1000
  matrix_dn0_pyrrole[,i] = test$ERG_d_ITN0
  matrix_rn0_pyrrole[,i] = test$ERG_r_ITN0
  # matrix_halflife[,i] = test$itn_half_life
  
}

dn0MEAN = rowMeans(matrix_dn0_pyrrole)
rn0MEAN = rowMeans(matrix_rn0_pyrrole)
halflifeMEAN = rowMeans(matrix_halflife)

dn0 = rn0 = gamman = array(dim=c(nrow(test),3))
for(j in 1:nrow(test)){
  dn0[j,] = c(as.numeric(quantile(matrix_dn0_pyrrole[j,],c(0.05,0.5,0.95))))
  rn0[j,] = c(as.numeric(quantile(matrix_rn0_pyrrole[j,],c(0.05,0.5,0.95))))
  gamman[j,] = c(as.numeric(quantile(matrix_halflife[j,],c(0.05,0.5,0.95))))
}

ppNets = data.frame(dn0_lo10 = dn0[,1],dn0_med = dn0[,2],dn0_up90 = dn0[,3],
                               rn0_lo10 = rn0[,1],rn0_med = rn0[,2],rn0_up90 = rn0[,3],
                               gamman_lo10 = gamman[,1],gamman_med = gamman[,2],gamman_up90 = gamman[,3])
ppNets$bioassay_mortality = seq(1,0,length=101)
head(ppNets)
tail(ppNets)

write.csv(pyrethroidOnlyNets,"parameters/loglogistic_beta_binomial_pyrethroid_only_nets.csv")
write.csv(pyrethroidPBONets,"parameters/loglogistic_beta_binomial_pyrethroid_pbo_nets.csv")
write.csv(ppNets,"parameters/loglogistic_beta_binomial_pyrethroid_pyrrole_nets.csv")
