##############################################################################
##
## Additional benefit of next-generation nets
## April 2022
##
###############################################################################


###################################################
##
## G2: New nets analysis, ESSENTIALS showing the relationships for the
## new nets are consistent with the standard and PBO-nets defined in 
## Churcher et al 2016 / Nash et al 2021 / Sherrard-Smith et al 2021
##
#########################################################################


library(rstan)

# For execution on a local, multicore CPU with excess RAM we recommend calling
options(mc.cores = parallel::detectCores())
# To avoid recompilation of unchanged Stan programs, we recommend calling
rstan_options(auto_write = TRUE)

g2_dat = read.csv("data/pyrethroid-pyrrole nets April 2022.csv",header=TRUE)

tail(g2_dat)
library(tidyr)
g2_dat = g2_dat %>% drop_na(n_dead_72h)
g2_dat = g2_dat %>% drop_na(n_total)
g2_dat = g2_dat %>% drop_na(n_dead_72h_PYRNET)
g2_dat = g2_dat %>% drop_na(n_total_PYRNET)

###########################################
##
## Part 1:
##
## Added benefit of the nets

##Relationship 1 pyrethroid-pyrrole benefit on top of pyrethroid-only nets 
##LOGISTIC BENEFIT
d_t = c(g2_dat$n_dead_72h) # number of mosquitoes dying IRS HUTS
n_t = c(g2_dat$n_total) # Number of mosquitoes entering IRS huts
x1 = c(g2_dat$n_dead_72h_PYRNET)/
  c(g2_dat$n_total_PYRNET)



data_G2 = data.frame(d_t,n_t,x1)
data_G2$Proportion_dead = data_G2$d_t/data_G2$n_t
x = seq(1,0,length=100)
fmG2 <- cbind(data_G2$d_t,data_G2$n_t-data_G2$d_t) ~ data_G2$x1 
glm_2 <- glm(fmG2, family = binomial())
role.fitted2 <- predict(glm_2, se.fit = TRUE, type = "response")
summary(glm_2)$coeff[2,1]
summary(glm_2)$coeff[1,1]

data_list_g2 = list(N = nrow(data_G2),
                     n_t = data_G2$n_t,
                     d_t = data_G2$d_t,
                     x = data_G2$x1)

stan_base <- stan(file="R code/stan models/binomial_fit_nets.stan", 
                  data=data_list_g2, 
                  warmup=1000,
                  control = list(adapt_delta = 0.8,
                                 max_treedepth = 20),
                  iter=2000, chains=4)
base <- rstan::extract(stan_base)

# saveRDS(stan_base,"stan model outputs/ento_g2_benefit_Apr2022data.RDS")

## diagnostics
# traceplot(stan_base)
# stan_diag(stan_base,
#           information = c("sample","stepsize", "treedepth","divergence"),
#           chain = 0)
# stan_hist(stan_base, include = TRUE, unconstrain = FALSE,
          # inc_warmup = FALSE)
# stan_rhat(stan_base,bins = 10)
mean(base$alpha2)
mean(base$alpha1)


preds_g2 = array(dim=c(100,4000))
for(i in 1:4000){
  preds_g2[,i] = 1 / (1 + exp(-base$alpha1[i] - base$alpha2[i]*x))
}

median_G2_pred = 1 / (1 + exp(-quantile(base$alpha1,0.5) - quantile(base$alpha2,0.5)*x))
upper_G2_pred = 1 / (1 + exp(-quantile(base$alpha1,0.975) - quantile(base$alpha2,0.975)*x))
lower_G2_pred = 1 / (1 + exp(-quantile(base$alpha1,0.025) - quantile(base$alpha2,0.025)*x))

##
## Drawing the figure
##
par(mfrow=c(1,1))
par(mar=c(5,6,3,3))

plot(median_G2_pred ~ x, ylab="Mosquito mortality next generation ITN (%)",
     ylim=c(0,1),xlim=c(0,1),pch="",yaxt="n",xaxt="n",
     xlab = "Mosquito mortality pyrethroid-only LLIN (%)",cex.lab=1.2,cex.axis=1.2)
axis(1,at=seq(0,1,0.2),labels=seq(0,100,20),cex.axis=1.2)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20),cex.axis=1.2)

g2_dat$x = g2_dat$n_dead_72h_PYRNET/g2_dat$n_total_PYRNET
g2_dat$prop_killed_g2 = g2_dat$n_dead_72h/g2_dat$n_total
points(g2_dat$prop_killed_g2[g2_dat$data_visible == 1] ~ 
         g2_dat$x[g2_dat$data_visible == 1],
       cex=1.6,pch=19)#,col="aquamarine3")
polygon(c(x,rev(x)),c(upper_G2_pred,rev(lower_G2_pred)),col=adegenet::transp("darkgreen",0.4),border=NA)


lines(median_G2_pred ~ x,col="darkgreen",lwd=2,lty=1)


abline(0,1,lty=2,col="blue")


############################
## 
## Add PBO line

##Relationship 2 PBO benefit on top of Normal LN 
##LOGISTIC BENEFIT
pbo_dat = read.csv("data/pyrethroid-pbo nets April 2022.csv")
names(pbo_dat)

pbo_dat = pbo_dat %>% drop_na(n_dead_24h)
pbo_dat = pbo_dat %>% drop_na(n_total)
pbo_dat = pbo_dat %>% drop_na(n_dead_24h_PYRNET)
pbo_dat = pbo_dat %>% drop_na(n_total_PYRNET)


d_t = c(pbo_dat$n_dead_24h) # number of mosquitoes dying IRS HUTS
n_t = c(pbo_dat$n_total) # Number of mosquitoes entering IRS huts
x1 = c(pbo_dat$n_dead_24h_PYRNET)/
  c(pbo_dat$n_total_PYRNET)

data_PBO = data.frame(d_t,n_t,x1)
data_PBO = data_PBO[complete.cases(data_PBO),]
data_PBO$Proportion_dead = data_PBO$d_t/data_PBO$n_t
x = seq(1,0,length=100)
fmPBO <- cbind(data_PBO$d_t,data_PBO$n_t-data_PBO$d_t) ~ data_PBO$x1 
glm_2 <- glm(fmPBO, family = binomial())
role.fitted2 <- predict(glm_2, se.fit = TRUE, type = "response")
summary(glm_2)$coeff[2,1]
summary(glm_2)$coeff[1,1]

data_list_PBO = list(N = nrow(data_PBO),
                     n_t = data_PBO$n_t,
                     d_t = data_PBO$d_t,
                     x = data_PBO$x1)

stan_base <- stan(file="R code/stan models/binomial_fit_nets.stan", 
                  data=data_list_PBO, 
                  warmup=1000,
                  control = list(adapt_delta = 0.8,
                                 max_treedepth = 20),
                  iter=2000, chains=4)
base <- rstan::extract(stan_base)

mean(base$alpha2)
mean(base$alpha1)

# saveRDS(stan_base,"stan model outputs/ento_pbo_benefit_Apr2022data.RDS")

preds_pbo = array(dim=c(100,4000))
for(i in 1:4000){
  preds_pbo[,i] = 1 / (1 + exp(-base$alpha1[i] - base$alpha2[i]*x))
}
median_pbo_pred = 1 / (1 + exp(-quantile(base$alpha1,0.5) - quantile(base$alpha2,0.5)*x))
upper_pbo_pred = 1 / (1 + exp(-quantile(base$alpha1,0.975) - quantile(base$alpha2,0.975)*x))
lower_pbo_pred = 1 / (1 + exp(-quantile(base$alpha1,0.025) - quantile(base$alpha2,0.025)*x))


pbo_option2 = 1 / (1 + exp(-summary(glm_2)$coeff[2,1]*x - summary(glm_2)$coeff[1,1]))
pbo_option2low = 1 / (1 + exp(-(summary(glm_2)$coeff[2,1]-1.96*summary(glm_2)$coeff[2,2])*x - (summary(glm_2)$coeff[1,1]-1.96*summary(glm_2)$coeff[1,2])))
pbo_option2upp = 1 / (1 + exp(-(summary(glm_2)$coeff[2,1]+1.96*summary(glm_2)$coeff[2,2])*x - (summary(glm_2)$coeff[1,1]+1.96*summary(glm_2)$coeff[1,2])))

pbo_dat$prop_killed_pbo = data_list_PBO$d_t/data_list_PBO$n_t
pbo_dat$x = data_list_PBO$x

points(pbo_dat$prop_killed_pbo[pbo_dat$data_visible == 1] ~ 
         pbo_dat$x[pbo_dat$data_visible == 1],
       cex=1.6,pch=19,col="purple")

polygon(c(x,rev(x)),c(upper_pbo_pred,rev(lower_pbo_pred)),col=adegenet::transp("purple",0.4),border=NA)


lines(median_pbo_pred ~ x,col="purple",lwd=2)

legend("bottomright",legend=c("Pyr-pyrrole (72 hr)",
                              "Pyr-PBO (24 hr)"),ncol=2,
       pch=c(19,19),col=c("aquamarine3","purple"),lty=1,bty="n")


##############################################
##
## Adding uncertainty by fitting the beta_binomial
library(shinystan)

#Logical indicators
diagnostic_plots <- FALSE
save_plot <- FALSE

#Create fine scaled data frame
fine_df <- data.frame("resistance" = seq(0.0,
                                   1.0,
                                   0.01))

#rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
n_chains <- 4
tot_iter <- 3000
burnin <- 2000
betbin_options <- list(adapt_delta = 0.999,
                       stepsize = 0.01,
                       max_treedepth = 20)
simulated_individuals <- 1e6

#Extract data for Stan
data_list_g2_beta = list(N = nrow(data_G2),
                         n_t = data_G2$n_t,
                         d_t = data_G2$d_t,
                         x = data_G2$x1,
                         
                         R = dim(fine_df)[1],
                         tau = fine_df$resistance,
                         N_M = simulated_individuals)


betbin_fit <- stan(file = 'R code/stan models/beta_binomial_fit_nets.stan',
                   data = data_list_g2_beta,
                   chains = n_chains,
                   iter = tot_iter,
                   warmup = burnin,
                   init_r = 1e-2,
                   control = betbin_options)

# saveRDS(betbin_fit,"stan model outputs/ento_g2_beta_binomial_benefit_Jan2024data.RDS")

#Investigate diagnostic plots
if (diagnostic_plots) {
  launch_shinystan(betbin_fit)
}



#Extract samples
betbin_samples <- rstan::extract(betbin_fit)

## overlay the original
preds2_g2 = array(dim=c(100,4000))
for(i in 1:4000){
  preds2_g2[,i] = 1 / (1 + exp(-betbin_samples$alpha1[i] - betbin_samples$alpha2[i]*x))
}

median_G2_pred2 = 1 / (1 + exp(-quantile(betbin_samples$alpha1,0.5) - quantile(betbin_samples$alpha2,0.5)*x))
upper_G2_pred2 = 1 / (1 + exp(-quantile(betbin_samples$alpha1,0.975) - quantile(betbin_samples$alpha2,0.975)*x))
lower_G2_pred2 = 1 / (1 + exp(-quantile(betbin_samples$alpha1,0.025) - quantile(betbin_samples$alpha2,0.025)*x))

lines(median_G2_pred2 ~ x)
polygon(c(x,rev(x)),
        c(lower_G2_pred2,rev(upper_G2_pred2)),col=adegenet::transp("aquamarine",0.4),border=NA)

mean(betbin_samples$a1)
mean(betbin_samples$a2)

mean(betbin_samples$alpha1)
mean(betbin_samples$alpha2)

mt_samples <- betbin_samples$rho / simulated_individuals
mt_sorted <- apply(mt_samples, 2, sort)
N_samples <- dim(mt_sorted)[1]
LB_ID <- round(N_samples*0.025)
UB_ID <- round(N_samples*0.975)
fine_df$mt_LB <- mt_sorted[LB_ID,]
fine_df$mt_UB <- mt_sorted[UB_ID,]
fine_df$mt_median <- apply(mt_samples, 2, median)

for(i in 1:nrow(mt_samples)){
  points(mt_samples[i,] ~ fine_df$resistance,col = adegenet::transp("aquamarine",0.3))
}

# polygon(c(fine_df$resistance,rev(fine_df$resistance)),
#         c(fine_df$mt_UB,rev(fine_df$mt_LB)),col=adegenet::transp("aquamarine",0.2),border=NA)
lines(fine_df$mt_UB ~ fine_df$resistance,lwd=2,col="aquamarine3",lty=2)
lines(fine_df$mt_LB ~ fine_df$resistance,lwd=2,col="aquamarine3",lty=2)

lines(colMeans(betbin_samples$p) ~ fine_df$resistance)
# lines(fine_df$mt_median ~ fine_df$resistance)



data_list_PBO = list(N = nrow(data_PBO),
                     n_t = data_PBO$n_t,
                     d_t = data_PBO$d_t,
                     x = data_PBO$x1,
                     R = dim(fine_df)[1],
                     tau = fine_df$resistance,
                     N_M = simulated_individuals)


betbin_fit2 <- stan(file = 'R code/stan models/beta_binomial_fit_nets.stan',
                   data = data_list_PBO,
                   chains = n_chains,
                   iter = tot_iter,
                   warmup = burnin,
                   init_r = 1e-2,
                   control = betbin_options)

# saveRDS(betbin_fit2,"stan model outputs/ento_pbo_beta_binomial_benefit_Jan2024data.RDS")

#Extract samples
betbin_samples2 <- rstan::extract(betbin_fit2)

## overlay the original
preds2_pbo = array(dim=c(100,4000))
for(i in 1:4000){
  preds2_pbo[,i] = 1 / (1 + exp(-betbin_samples2$alpha1[i] - betbin_samples2$alpha2[i]*x))
}

median_pbo_pred2 = 1 / (1 + exp(-quantile(betbin_samples2$alpha1,0.5) - quantile(betbin_samples2$alpha2,0.5)*x))
upper_pbo_pred2 = 1 / (1 + exp(-quantile(betbin_samples2$alpha1,0.975) - quantile(betbin_samples2$alpha2,0.975)*x))
lower_pbo_pred2 = 1 / (1 + exp(-quantile(betbin_samples2$alpha1,0.025) - quantile(betbin_samples2$alpha2,0.025)*x))

lines(median_pbo_pred2 ~ x)
polygon(c(x,rev(x)),
        c(lower_pbo_pred2,rev(upper_pbo_pred2)),col=adegenet::transp("purple",0.4),border=NA)

mean(betbin_samples2$a1)
mean(betbin_samples2$a2)

mean(betbin_samples2$alpha1)
mean(betbin_samples2$alpha2)

mt_samples2 <- betbin_samples2$rho / simulated_individuals
mt_sorted2 <- apply(mt_samples2, 2, sort)
N_samples2 <- dim(mt_sorted2)[1]
LB_ID <- round(N_samples2*0.025)
UB_ID <- round(N_samples2*0.975)
fine_df$mt_LB2 <- mt_sorted2[LB_ID,]
fine_df$mt_UB2 <- mt_sorted2[UB_ID,]
fine_df$mt_median2 <- apply(mt_samples2, 2, median)

# polygon(c(fine_df$resistance,rev(fine_df$resistance)),
#         c(fine_df$mt_UB,rev(fine_df$mt_LB)),col=adegenet::transp("purple",0.4),border=NA)
lines(fine_df$mt_UB2 ~ fine_df$resistance,lwd=2,col="purple",lty=2)
lines(fine_df$mt_LB2 ~ fine_df$resistance,lwd=2,col="purple",lty=2)

lines(colMeans(betbin_samples$p) ~ fine_df$resistance)
