## Part 1

## fitting the stan model for mosquito survival in susc bioassay 
## to predict the survival in experimental huts

library(rstan)

# For execution on a local, multicore CPU with excess RAM we recommend calling
options(mc.cores = parallel::detectCores())
# To avoid recompilation of unchanged Stan programs, we recommend calling
rstan_options(auto_write = TRUE)

## data
mt = read.csv("data/mortality_data_from_Nash2021.csv",header=TRUE)

## remove Ifakara huts
mt = subset(mt, mt$Hut != "Ifakara")
library(dplyr)
selectiona <- c("Olyset Net", "Interceptor", 
                "PermaNet 2.0","PermaNet", 
                "MiraNet",
                "DuraNet","Duranet",
                "Yahe",
                "MAGNet", "DawaPlus 2.0", "Yorkool",
                "Panda Net 2.0","PandaNet 2.0",
                "Royal Sentry")
mt <- mt %>% filter(Net%in%selectiona)

data_list_MT = list(S = nrow(mt),
                    N_b = mt$Bioassay_N_tested,
                    N_h = mt$N_total,
                    X_b = mt$Bioassay_N_tested - mt$Bioassay_N_dead,
                    X_h = mt$N_total - mt$N_dead,
                    nsite = length(unique(mt$Site)),
                    site = as.numeric(factor(mt$Site)),
                    S_test = 101,
                    theta_b_test = seq(0,1,length=101))

stan_base <- rstan::stan(file="R code/stan models/Bioassay/A1 logistic Model.stan", 
                         data=data_list_MT, 
                         warmup=1000,
                         control = list(adapt_delta = 0.8,
                                        max_treedepth = 20),
                         iter=2000, chains=4)
base <- rstan::extract(stan_base)

median(base$c)
median(base$a)

f_logistic <- function(x, a, c){
        mort = 1/(1+exp(-(x-c)*a))
        surv = (1-mort)
        return(mort)
}

x = seq(0,1,length=101)
surv = f_logistic(x=x, a=median(base$a), c=median(base$c))

# surv_low = f_logistic(x=x, a=as.numeric(quantile(base$a,0.05)), c=as.numeric(quantile(base$c,0.95)))
# surv_upp = f_logistic(x=x, a=as.numeric(quantile(base$a,0.95)), c=as.numeric(quantile(base$c,0.05)))
plot(surv ~ x,type="l",
     xlab = "Survival at discriminatory dose bioassay (%)",
     ylab = "Survival in experimental hut trial (%)",
     xaxt="n",yaxt="n",ylim=c(0,1),xlim=c(0,1))
axis(1,at=seq(0,1,by=0.2),labels=seq(0,100,by=20))
axis(2,las=2,at=seq(0,1,by=0.2),labels=seq(0,100,by=20))

# aaL = which(base$a < as.numeric(quantile(base$a,0.025)) & base$a > as.numeric(quantile(base$a,0.025))-0.02)[1]
# aaU = which(base$a > as.numeric(quantile(base$a,0.975)) & base$a < as.numeric(quantile(base$a,0.975))+0.02)[1]
# 
# ccL = which(base$c < as.numeric(quantile(base$c,0.025)) & base$c > as.numeric(quantile(base$c,0.025))-0.02)[1]
# ccU = which(base$c > as.numeric(quantile(base$c,0.975)) & base$c < as.numeric(quantile(base$c,0.975))+0.02)[1]

for(i in c(20,3000)){
  surv = (1/(1+exp(-(x-base$c[i])*base$a[i])))
  lines(surv ~ x,col=adegenet::transp("grey",0.8))
}

surv_low = f_logistic(x=x, a=base$a[20], c=base$c[20])
surv_upp = f_logistic(x=x, a=base$a[3000], c=base$c[3000])
polygon(c(x,rev(x)),c(surv_low,rev(surv_upp)),col=adegenet::transp("green",0.4),border=NA)




points(c((mt$Bioassay_N_tested - mt$Bioassay_N_dead)/mt$Bioassay_N_tested),
       c((mt$N_total - mt$N_dead)/mt$N_total),col=adegenet::transp("orange",0.7),pch=19,
       cex = c(log(50*(mt$N_total/max(mt$N_total))+1)))

## saveRDS(stan_base,"stan model outputs/log_logistic_fit.RDS")
