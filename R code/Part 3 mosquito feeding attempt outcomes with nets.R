
################################
##
## Part 3: Probable outcomes of feeding attempts in huts
##
#################################

##################################################################
##
## Read in an updated file for ALL net studies 

## All data available up to April 2022
df <- read.csv("data/All Hut Studies April 2022.csv",header = TRUE)


# Remove Ifakara, UNTREATED and those without totals or number dead
# given Nash et al 2021 = Ifakara huts produced distinct data patterns
# also the best parameter set for validations is using West/East hut data (Sherrard-Smith et al. 2022 Nat Communs)
library(dplyr)
df2 <- df %>%
  filter(N_dead_comb!="NA" & Totals_known!="0" & N_total!="NA" & Intervention!="Untreated" & Hut!="Ifakara")
# all the data

# # select the nets with 0 washes
selection0 <- c(0)
dfAll <- df2 %>% filter(No_washes%in%selection0)
# dfAll we are interested in the associations before washing/aging nets here

## Pyrethroid only nets included are;
## (all WHO Recommended according to Okumu & Finda 2021 book chapter)
# Olyset Net
# Interceptor
# Royal Sentry     
# Royal Sentry 2.0 ## No data
# PermaNet 2.0
# DuraNet LLIN
# MiraNet      
# MAGNet
# YaheLN       
# DawaPlus 2.0 LLIN
# SafeNet      # No data
# Yorkool LN
# Panda Net2.0 LLIN      


## Pyrethroid-PBO nets included are;
## (all WHO Recommended according to Okumu & Finda 2021 book chapter)

# Olyset Plus
# PermaNet 3.0
# Veeralin
# Dawa Plus 3.0 (ALSO 4.0)
# Tsara Boost # No data


## Pyrethroid-pyrrole nets 
# Interceptor G2
# Phantom SC + Fendona SC
# PermaNetP191


# All WHO recommended nets together
selectiona <- c("Olyset Net", "Interceptor", 
                "PermaNet 2.0","PermaNet2", 
                "MiraNet",
                "DuraNet","Duranet",
                "Yahe",
                "MAGNet", "DawaPlus 2.0", "Yorkool",
                "Panda Net 2.0","PandaNet 2.0",
                "Royal Sentry",
                
                "Olyset Plus","OlysetPlus",
                "PermaNet 3.0","PermaNet3", 
                "Veeralin", 
                "DawaPlus 3.0","DawaPlus 4.0","Dawaplus 4.0", 
                
                "Interceptor G2","InterceptorG2",
                "Phantom SC + Fendona SC",
                "PermaNetP191")
dfAllRec <- dfAll %>% filter(Intervention%in%selectiona)


add_pyr_baselines_f = function(df2){
  
  nsite <- max(df2$Site_number)
  site <- df2$Site_number
  
  for(i in 1:length(df2$R_ID)){
    df2$pyr_base_count[i] = round(mean(df2$N_total[df2$LLIN_Gen == 1 & df2$R_ID == df2$R_ID[i]],na.rm = TRUE),0)
    df2$pyr_base_killed[i] = round(mean(df2$N_dead_comb[df2$LLIN_Gen == 1 & df2$R_ID == df2$R_ID[i]],na.rm = TRUE),0)
    
  }
  
  df2$pyr_base_countN = ifelse(df2$LLIN_Gen == 1,df2$N_total,df2$pyr_base_count)  
  df2$pyr_base_killN = ifelse(df2$LLIN_Gen == 1,df2$N_dead_comb,df2$pyr_base_killed)
  
  df3 <- df2[complete.cases(df2$pyr_base_countN), ] 
  data2 <- df3[complete.cases(df3$pyr_base_killN), ] 
  
  return(data2)
}

dfAlltest = add_pyr_baselines_f(dfAll)


setup_inputs_ALL = function(df2){
  
  # Prep data
  S <- as.numeric(length(df2$N_total))
  X_survive <- df2$N_total-df2$N_dead_comb
  N_caught <- df2$N_total
  nsite <- max(df2$Site_number)
  site <- df2$Site_number
  X_control <- df2$N_Control
  S_sim <- 300
  theta_sim <- seq(0, 1, length.out = S_sim)
  
  # Need to pass it the data:
  data_stan <- list(S=S, 
                    X_survive=X_survive,
                    N_caught=N_caught,
                    nsite=nsite,
                    site=site,
                    X_control=X_control,
                    S_sim=S_sim,
                    theta_sim=theta_sim)

  
  df2 <- df2 %>%
    mutate(X_survive=N_total-N_dead_comb,
           frac_surv=X_survive / N_caught,
           Perc_surv=frac_surv*100,
           change=(X_control-N_caught) / X_control,
           Perc_int_vs_control=change*100)
  
  return(data_stan)
  
}


setup_inputs_pyr = function(df2){

  # Prep data
  S <- as.numeric(length(df2$pyr_base_countN))
  X_survive <- df2$pyr_base_countN-df2$pyr_base_killN
  N_caught <- df2$pyr_base_countN
  nsite <- max(df2$Site_number)
  site <- df2$Site_number
  X_control <- df2$N_Control
  S_sim <- 300
  theta_sim <- seq(0, 1, length.out = S_sim)
  
  # Need to pass it the data:
  data_stan <- list(S=S, 
                    X_survive=X_survive,
                    N_caught=N_caught,
                    nsite=nsite,
                    site=site,
                    X_control=X_control,
                    S_sim=S_sim,
                    theta_sim=theta_sim)
  
  df2 <- df2 %>%
    mutate(X_survive=N_total-N_dead_comb,
           frac_surv=X_survive / N_caught,
           Perc_surv=frac_surv*100,
           change=(X_control-N_caught) / X_control,
           Perc_int_vs_control=change*100)
  
  return(data_stan)
  
}


# Compile model
full_model<- stan_model("R code/stan models/Deter_model_NoRE_DETER_opt1.stan") # flex params

run1_f = function(df2,dataset){
  data_stan = setup_inputs_ALL(df2)
  
  # Run model
  fit_full <- sampling(full_model, data=data_stan, iter=2000, chains=4)
  
  # save fit
  saveRDS(fit_full, paste0("stan model outputs/April_2022_ento_deterrence_",dataset,".rds"))
}

run1_f(dfAllRec,"AllRec")

fit1_all <- readRDS("stan model outputs/April_2022_ento_deterrence_AllRec.rds")
print(fit1_all, pars=c("c","d","e", "kappa","P_control")) 

# Extract ratio
plot_ratio_prob <- function(P_param, fit){
  P <- rstan::extract(fit, P_param)[[1]]
  P_med <- apply(P, 2, median)
  return(P_med*100)
}

P_ratio_all <- plot_ratio_prob("P_ratio", fit1_all)

create_output_f = function(fit_full){
  # Extract upper credible intervals
  plot_upper_funct <- function(P_param, fit){
    P <- rstan::extract(fit, P_param)[[1]]
    P_med <- apply(P, 2, median)
    P_upper <- apply(P, 2, function(x) quantile(x,0.025))
    return(P_upper)
  }
  
  P_ratio_upper <- plot_upper_funct("P_ratio", fit_full)
  P_det_upper <- plot_upper_funct("P_deter_sim", fit_full)
  P_inside_upper <- plot_upper_funct("P_inside_sim", fit_full)
  
  
  plot_med_funct <- function(P_param, fit){
    P <- rstan::extract(fit, P_param)[[1]]
    P_med <- apply(P, 2, median)
    return(P_med)
  }
  
  P_ratio_med <- plot_med_funct("P_ratio", fit_full)
  P_det_med <- plot_med_funct("P_deter_sim", fit_full)
  P_inside_med <- plot_med_funct("P_inside_sim", fit_full)
  
  
  # Extract lower credible intervals
  plot_lower_funct <- function(P_param, fit){
    P <- rstan::extract(fit, P_param)[[1]]
    P_med <- apply(P, 2, median)
    P_lower <- apply(P, 2, function(x) quantile(x,0.975))
    return(P_lower) 
  }
  
  P_ratio_lower <- plot_lower_funct("P_ratio", fit_full)
  P_det_lower <- plot_lower_funct("P_deter_sim", fit_full)
  P_inside_lower <- plot_lower_funct("P_inside_sim", fit_full)
  
  
  a_theta <- seq(0,1,length=300)
  my_surv <- a_theta*100
  
  cred_int_data <- data.frame(P_ratio_med, P_ratio_upper, P_ratio_lower,my_surv)
  det_cred <- data.frame(P_det_med, P_det_upper, P_det_lower,a_theta)
  inside_cred <- data.frame(P_inside_med, P_inside_upper, P_inside_lower,a_theta)
  
  return(list(cred_int_data,
              det_cred,
              inside_cred))  
}

mod_data_allRec = create_output_f(fit1_all)

S_sim <- 300
theta_sim <- seq(0, 1, length.out = S_sim)

par(mfrow=c(2,2))
plot(mod_data_allRec[[2]][,1] ~ theta_sim,type = "l",
     ylim=c(0,1),xlim=c(0,1),yaxt="n",
     col="grey45",
     main="Using matched net data for x-axis",
     ylab="Probability of being deterred",
     xlab="Probability of hut survival")
axis(2,las=2,at=seq(0,1,0.20),labels=seq(0,1,0.2))
##    a. All data (all nets as per Nash et al)
polygon(c(theta_sim,rev(theta_sim)),c(mod_data_allRec[[2]][,2],rev(mod_data_allRec[[2]][,3])),border = NA,col=adegenet::transp("grey75",0.4))
lines(mod_data_allRec[[2]][,2] ~ theta_sim,lty=2,col="grey75")
lines(mod_data_allRec[[2]][,3] ~ theta_sim,lty=2,col="grey75")


# extract P_control
extract_control <- function(P_param, fit_full){
  P <- rstan::extract(fit_full, P_param)[[1]]
  return(median(P))
}

P_control_all <- extract_control("P_control", fit1_all)


plot(mod_data_allRec[[3]][,1] ~ theta_sim,type = "l",
     ylim=c(0,1),xlim=c(0,1),yaxt="n",
     col="grey55",
     ylab="Probability of being caught in control hut",
     xlab="Probability of hut survival")
axis(2,las=2,at=seq(0,1,0.20),labels=seq(0,1,0.2))
##    a. All data (all nets as per Nash et al)
polygon(c(theta_sim,rev(theta_sim)),c(mod_data_allRec[[3]][,2],rev(mod_data_allRec[[3]][,3])),border = NA,col=adegenet::transp("grey75",0.4))
lines(mod_data_allRec[[3]][,2] ~ theta_sim,lty=2,col="grey75")
lines(mod_data_allRec[[3]][,3] ~ theta_sim,lty=2,col="grey75")
abline(h=P_control_all,lwd=2,col="grey75")


########################
##
## Mortality and successful feeding

dfAllRec = dfAllRec[!is.na(dfAllRec$N_fed),]
dfAllRec = dfAllRec[!is.na(dfAllRec$N_dead_comb),]
dfAllRec = dfAllRec[!is.na(dfAllRec$N_total),]
dfAllRec = subset(dfAllRec,dfAllRec$N_total != 0)

setup_inputs_fed = function(df2){

    # Data for stan
  S <- as.numeric(length(df2$N_total))
  X_survive <- df2$N_total-df2$N_dead_comb
  N_caught <- df2$N_total
  X_fed <- df2$N_fed
  X_survfed <- as.integer(X_fed * (1-(df2$N_dead_comb/df2$N_total))) # No. fed * proportion survived
  Perc_survfed <- X_survfed/df2$N_total*100
  nsite <- max(df2$Site_number)
  site <- df2$Site_number
  S_sim <- 300
  theta_sim <- seq(0, 1, length.out = S_sim)
  
  # Need to pass it the data:
  data_stan <- list(S=S, 
                    X_survive=X_survive,
                    N_caught=N_caught,
                    X_succfed=X_survfed,
                    nsite=nsite,
                    site=site,
                    S_sim = S_sim,
                    theta_sim = theta_sim)
  
  return(data_stan)
}


full_model<- stan_model("R code/stan models/Model_succfed_WithREs.stan")

run_f2 = function(df2,dataset){
  data_stan = setup_inputs_fed(df2)
  
  # Fit model
  fit_hut_ALL <- sampling(full_model, data=data_stan, iter=2000, chains=4)
  
  # save fit
  saveRDS(fit_hut_ALL, paste0("stan model outputs/April_2022_ento_feeding_0washes_",dataset,".rds"))
}



run_f2(dfAllRec,"AllRec")

fit3_a <- readRDS("stan model outputs/April_2022_ento_feeding_0washes_AllRec.rds")

print(fit3_a, pars=c("a","b")) 

create_output2_f = function(fit_hut_ALL){
  # Extract prob of succ fed
  plot_fed_prob <- function(P_param, fit){
    P <- rstan::extract(fit, P_param)[[1]]
    P_med <- apply(P, 2, median)
    return(P_med*100)
  }
  
  P_succfed <- plot_fed_prob("P_fed_sim", fit_hut_ALL)
  
  
  # Extract upper and lower cred intervals
  plot_upper_funct <- function(P_param, fit){
    P <- rstan::extract(fit, P_param)[[1]]
    P_upper <- apply(P, 2, function(x) quantile(x,0.025))
    return(P_upper)
  }
  
  plot_lower_funct <- function(P_param, fit){
    P <- rstan::extract(fit, P_param)[[1]]
    P_lower <- apply(P, 2, function(x) quantile(x,0.975))
    return(P_lower) 
  }
  
  P_fed_upper <- plot_upper_funct("P_fed_sim", fit_hut_ALL)*100
  P_fed_lower <- plot_lower_funct("P_fed_sim", fit_hut_ALL)*100
  
  
  a_theta <- seq(0,1,length=300)
  my_surv <- a_theta*100
  
  cred_int_data <- data.frame(P_succfed, P_fed_upper, P_fed_lower,my_surv)
  
  return(cred_int_data)  
}
cred_int_data_a = create_output2_f(fit3_a)

dfAllRec$pc_surv <- (dfAllRec$N_total-dfAllRec$N_dead_comb)/dfAllRec$N_total*100
dfAllRec$prop_surv_fed <- dfAllRec$pc_fed/100*dfAllRec$pc_surv/100
dfAllRec$suc.u <- dfAllRec$prop_surv_fed*dfAllRec$N_total
dfAllRec$pc_surv_fed <- 100* dfAllRec$suc.u/dfAllRec$N_total
dfAllRec$prop_dead <- dfAllRec$N_dead_comb/dfAllRec$N_total
dfAllRec$repelled<-(1-dfAllRec$prop_surv_fed-dfAllRec$prop_dead)*100  # proportion repelled * 100 = percentage repelled

dfAllRec$pt_cex = c(0.5 + dfAllRec$N_total/1000)

outline.hut <- ifelse(dfAllRec$Hut=="West","#00008b",
                      ifelse(dfAllRec$Hut=="East","#8b0000",
                             ifelse(dfAllRec$Hut=="Ifakara","#006400",NA)))

fed.dataa <- dfAllRec %>% filter(Intervention%in%selectiona)

plot(dfAllRec$pc_surv_fed ~ dfAllRec$pc_surv, cex=dfAllRec$pt_cex,
     col=adegenet::transp("white",0.4),
     pch=19,
     ylab = "Successfully blood fed (%)",
     xlab = "Mosquito survival in EHT (%)",
     yaxt="n",ylim=c(0,100))
axis(2,las=2,at=c(0,20,40,60,80,100))
points(fed.dataa$pc_surv_fed ~ fed.dataa$pc_surv,
       cex=fed.dataa$pt_cex,
       col=adegenet::transp("grey",0.4),
       pch=19,ylim=c(0,100))

lines(cred_int_data_a$P_succfed ~ cred_int_data_a$my_surv, col = "grey55")
polygon(c(cred_int_data_a$my_surv,rev(cred_int_data_a$my_surv)),
        c(cred_int_data_a$P_fed_upper,rev(cred_int_data_a$P_fed_lower)),
        border=NA,col=adegenet::transp("grey55",0.4))



