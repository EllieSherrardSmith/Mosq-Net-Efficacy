#############################
##
##  4 fitting half life following Churcher et al 2016
##
#############################
library(rstan)
options(mc.cores=4)
rstan_options(auto_write = TRUE)
library(dplyr)
library(tidyr)

## Estimating half life

df <- read.csv("data/All Hut Studies April 2022.csv") # mac

# Remove Ifakara and those without totals or number dead


df2 <- df %>%
  filter(N_dead_comb!="NA" & Totals_known!="0" & N_total!="NA" & Intervention!="Untreated" & Hut!="Ifakara")
# all the data

# # select the nets with 0 washes
selection0 <- c(0)
dfAll <- df2 %>% filter(No_washes%in%selection0)


selectiona <- c("Olyset Net", "Interceptor", 
                "PermaNet 2.0","PermaNet2", 
                "MiraNet",
                "DuraNet","Duranet",
                "Yahe",
                "MAGNet", "DawaPlus 2.0", "Yorkool",
                "Panda Net 2.0","PandaNet 2.0",
                "Royal Sentry")#,
## We do not know if it is representative to use other nets     
                # "Olyset Plus","OlysetPlus",
                # "PermaNet 3.0","PermaNet3", 
                # "Veeralin", 
                # "DawaPlus 3.0","DawaPlus 4.0","Dawaplus 4.0", 
                # 
                # "Interceptor G2","InterceptorG2",
                # "Phantom SC + Fendona SC",
                # "PermaNetP191")
dfa <- dfAll %>% filter(Intervention%in%selectiona)

unique(df2$No_washes)
selection20 <- c(20)
dfa20 <- df2 %>% filter(No_washes%in%selection20)
dfa20 <- dfa20 %>% filter(Intervention%in%selectiona)


dfa_match = merge(dfa,dfa20,by="Orig_ID")


dfa_match = dfa_match %>% drop_na(N_total.x)
dfa_match = dfa_match %>% drop_na(N_dead_comb.x)
dfa_match = dfa_match %>% drop_na(N_total.y)
dfa_match = dfa_match %>% drop_na(N_dead_comb.y)

dfa_match = subset(dfa_match,dfa_match$N_total.x != 0)
dim(dfa_match)

## Simple model first
setup_inputs = function(df2){
  
  # Prep data
  S <- as.numeric(length(df2$N_total.x))
  prop_dead <- c(df2$N_dead_comb.x)/df2$N_total.x
  X_halflife <- c(df2$N_total.y - df2$N_dead_comb.y)
  X_caught_halflife <- df2$N_total.y
  
  
  # Need to pass it the data:
  data_stan <- list(S=S, 
                    prop_dead=prop_dead,
                    X_halflife=X_halflife,
                    N_caught_halflife=X_caught_halflife)#,
                    # nsite = nsite,
                    # site = site)
  
  
  return(data_stan)
  
}

# Compile model

full_model<- stan_model("R code/stan models/Model_halflife_0.stan") # flex params


stan_base <- rstan::stan(file="R code/stan models/Model_halflife_0.stan", 
                         data=setup_inputs(dfa_match), 
                         warmup=1000,
                         control = list(adapt_delta = 0.8,
                                        max_treedepth = 20),
                         iter=2000, chains=4)
base <- rstan::extract(stan_base)
median(base$a)
median(base$b)

run_f = function(df2,dataset){
  data_stan = setup_inputs(df2)
  
  # Run model
  fit_full <- sampling(full_model, data=data_stan, iter=2000, chains=4)
  # params_a = extract(fit_full)
  
  # save fit
  saveRDS(fit_full, paste0("stan model outputs/ento_fit1_half_life_",dataset,".rds"))
}
run_f(dfa_match,"pyrethroid_nets")#all data


# load fit
fit_a <- readRDS("stan model outputs/ento_fit1_half_life_pyrethroid_nets.rds")


#########################################
##
## Plotting the data against the half-life predictions
##
##
###########################################

## proportion of mosquitoes killed for 0 washes

prop_killed_0_wash = setup_inputs(dfa_match)$prop_dead
prop_killed20_wash = setup_inputs(dfa_match)$X_halflife/setup_inputs(dfa_match)$N_caught_halflife

test = base
names(test)
sim_X = seq(0,1,0.01)
P_hl1 <- 1/(1 + exp(-(median(test$a) + median(test$b) * sim_X)))
P_hl1_low <- 1/(1 + exp(-(quantile(test$a,0.025) + quantile(test$b,0.025) * sim_X)))
P_hl1_upp <- 1/(1 + exp(-(quantile(test$a,0.975) + quantile(test$b,0.975) * sim_X)))
  
par(mfrow=c(1,1))
  
lp_tau = seq(0,1,0.01)-0.5
mort1 = seq(0,1,length=length(lp_tau))
mu1 = -2.36
rho1 = -3.05
original_eLife = 1/(1 + exp(-(-mu1 - rho1 * lp_tau)))
  
plot(original_eLife ~ c(1-mort1),ylim=c(0,1),xlim=c(0,1),
       ylab = "EHT Mortality after 20 washes (%)",
       xlab = "EHT Mortality after 0 washes (%)",
       xaxt="n",yaxt="n")
axis(1,at=seq(0,1,0.2),labels=seq(0,100,20))
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20))
  
lines(P_hl1 ~ sim_X,col="red")
polygon(c(sim_X,rev(sim_X)),
          c(P_hl1_low,rev(P_hl1_upp)),border=NA,col=adegenet::transp("red",0.4))
points(prop_killed_0_wash, prop_killed20_wash,
         ylim=c(0,1),xlim=c(0,1),
         col=adegenet::transp("darkred",0.6),pch=19)
  
Hw = log(2)/(1/(1+exp(-(median(test$a) + median(test$b) * sim_X))))

##############################
##
## Check with model code
##
##############################

#{halflife}
#Assay to hut mortality conversion - median estimates	
param_a = 0.89 ## All 0.89 ## West 0.72 ## East 1.05
param_b = 0.47 ## All 0.47 ## West 0.88 ## East 0.36

resistance = seq(0,1,0.01)
mort_assay = 1 - resistance
mort_huta   = 1/(1+(mort_assay/param_a)^-param_b)
mort_maxA   = 1/(1+((1)/param_a)^-param_b)

## ORIGINAL_ELIFE WITH OLDER FUNCTION FOR Bioasay and EHT mortality
alpha1 = 0.63
beta1 = 4.00
mort_hutaORIG   = (1/(1+exp(-(alpha1 +beta1*(mort_assay-0.5)	))))
mort_maxAORIG   = (1/(1+exp(-(alpha1 +beta1*(1-0.5)	))))

# my_max_washes_a = mu1 +rho1*(mort_maxA-0.5)		
# my_max_washes   = log(2)/(exp(my_max_washes_a)/(1+exp(my_max_washes_a)))
# 
## OR
# my_max_washes = log(2) /(1/(1+exp(-(mu1 +rho1*(mort_maxA-0.5)	))))
my_max_washes = log(2) /(1/(1+exp(-(mu1 +rho1*(mort_maxA)	))))
my_max_washesORIG = log(2) /(1/(1+exp(-(mu1 +rho1*(mort_maxAORIG-0.5)	))))

# wash_decay_rate_a = mu1 +rho1*(mort_huta-0.5)
# wash_decay_rate   = log(2)/(exp(wash_decay_rate_a)/(1+exp(wash_decay_rate_a)))

## OR
# wash_decay_rate = log(2) /(1/(1+exp(-(mu1 +rho1*(mort_huta-0.5)	))))
wash_decay_rate = log(2) /(1/(1+exp(-(mu1 +rho1*(mort_huta)	))))
wash_decay_rateORIG = log(2) /(1/(1+exp(-(mu1 +rho1*(mort_hutaORIG-0.5)	))))
net_halflife = 2.64

itn_half_life     = wash_decay_rate/my_max_washes*net_halflife
itn_half_lifeORIG = wash_decay_rateORIG/my_max_washesORIG*net_halflife


plot(itn_half_life ~ resistance,ylim=c(0,2.7),type="l",
     yaxt = "n",ylab = "Insecticidal impact half life (years)",
     xaxt = "n",xlab = "Survival in bioassay (%)")
axis(1,at=seq(0,1,0.2),labels=seq(0,100,20))
axis(2,las=2,at=seq(0,2.5,0.5))


P_hl <- 1/(1 + exp(-(mean(test$a) + mean(test$b) * mort_huta)))
P_hlMAX <- 1/(1 + exp(-(mean(test$a) + mean(test$b) * mort_maxA )))

P_hl_low <- 1/(1 + exp(-(mean(test$a) + mean(test$b) * mort_huta)))
P_hl_upp <- 1/(1 + exp(-(mean(test$a) + mean(test$b) * mort_huta)))
P_hlMAX_low <- 1/(1 + exp(-(quantile(test$a,0.025) + quantile(test$b,0.025) * mort_maxA)))
P_hlMAX_upp <- 1/(1 + exp(-(quantile(test$a,0.975) + quantile(test$b,0.975) * mort_maxA)))

P_hlORIG <- 1/(1 + exp(-(mean(test$a) + mean(test$b) * mort_hutaORIG)))
P_hlMAXORIG <- 1/(1 + exp(-(mean(test$a) + mean(test$b) * mort_maxAORIG )))

P_hl_lowORIG <- 1/(1 + exp(-(mean(test$a) + mean(test$b) * mort_hutaORIG)))
P_hl_uppORIG <- 1/(1 + exp(-(mean(test$a) + mean(test$b) * mort_hutaORIG)))
P_hlMAX_lowORIG <- 1/(1 + exp(-(quantile(test$a,0.025) + quantile(test$b,0.025) * mort_maxAORIG)))
P_hlMAX_uppORIG <- 1/(1 + exp(-(quantile(test$a,0.975) + quantile(test$b,0.975) * mort_maxAORIG)))

wash_decay_rateNEW = log(2) /P_hl
wash_decay_rateNEW_low = log(2) /P_hl_low
wash_decay_rateNEW_upp = log(2) /P_hl_upp

my_max_washesNEW = log(2)/P_hlMAX
my_max_washesNEW_low = log(2)/P_hlMAX_low
my_max_washesNEW_upp = log(2)/P_hlMAX_upp
itn_half_lifeNEW     = wash_decay_rateNEW/my_max_washesNEW*net_halflife

itn_half_lifeNEW_low     = wash_decay_rateNEW_low/my_max_washesNEW_low*net_halflife
itn_half_lifeNEW_upp     = wash_decay_rateNEW_upp/my_max_washesNEW_upp*net_halflife

lines(itn_half_lifeNEW ~ resistance,col="red")
lines(itn_half_lifeORIG ~ resistance,col="grey25",lty=2)

polygon(c(resistance,rev(resistance)),
        c(itn_half_lifeNEW_low,rev(itn_half_lifeNEW_upp)),border=NA,col=adegenet::transp("red",0.4))

wash_decay_rateNEWORIG = log(2) /P_hlORIG
wash_decay_rateNEW_lowORIG = log(2) /P_hl_lowORIG
wash_decay_rateNEW_uppORIG = log(2) /P_hl_uppORIG

my_max_washesNEWORIG = log(2)/P_hlMAXORIG
my_max_washesNEW_lowORIG = log(2)/P_hlMAX_lowORIG
my_max_washesNEW_uppORIG = log(2)/P_hlMAX_uppORIG
itn_half_lifeNEWORIG     = wash_decay_rateNEWORIG/my_max_washesNEWORIG*net_halflife

itn_half_lifeNEW_lowORIG     = wash_decay_rateNEW_lowORIG/my_max_washesNEW_lowORIG*net_halflife
itn_half_lifeNEW_uppORIG     = wash_decay_rateNEW_uppORIG/my_max_washesNEW_uppORIG*net_halflife

lines(itn_half_lifeNEWORIG ~ resistance,col="red",lty=2)

polygon(c(resistance,rev(resistance)),
        c(itn_half_lifeNEW_lowORIG,rev(itn_half_lifeNEW_uppORIG)),border=NA,col=adegenet::transp("red",0.4))









