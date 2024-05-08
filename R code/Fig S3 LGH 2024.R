##LOGISTIC BENEFIT
pbo_dat = read.csv("data/pyrethroid-pbo nets April 2024.csv")
names(pbo_dat)

pbo_dat = subset(pbo_dat,pbo_dat$washed == 0)

d_t = c(pbo_dat$N_dead_comb_B) # number of mosquitoes dying IRS HUTS
n_t = c(pbo_dat$N_total_B) # Number of mosquitoes entering IRS huts
x1 = c(pbo_dat$N_dead_comb_pyrethroid)/
  c(pbo_dat$N_total_pyrethroid)

data_PBO = data.frame(d_t,n_t,x1)
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

mean(base$alpha2);quantile(base$alpha2,c(0.05,0.5,0.95))
mean(base$alpha1);quantile(base$alpha1,c(0.05,0.5,0.95))

# saveRDS(stan_base,"stan model outputs/ento_pbo_benefit_LancetGH2024.RDS")

preds_pbo = array(dim=c(100,4000))
for(i in 1:4000){
  preds_pbo[,i] = 1 / (1 + exp(-base$alpha1[i] - base$alpha2[i]*x))
}
median_pbo_pred = 1 / (1 + exp(-quantile(base$alpha1,0.5) - quantile(base$alpha2,0.5)*x))
upper_pbo_pred = 1 / (1 + exp(-quantile(base$alpha1,0.975) - quantile(base$alpha2,0.975)*x))
lower_pbo_pred = 1 / (1 + exp(-quantile(base$alpha1,0.025) - quantile(base$alpha2,0.025)*x))
mean_pbo_pred = 1 / (1 + exp(-mean(base$alpha1) - mean(base$alpha2)*x))


pbo_dat$prop_killed_pbo = pbo_dat$N_dead_comb_B/pbo_dat$N_total_B
pbo_dat$x = pbo_dat$N_dead_comb_pyrethroid/pbo_dat$N_total_pyrethroid
pbo_dat$size = sqrt(pbo_dat$N_total_B)/6

par(mfrow=c(1,2))
par(mar=c(5,6,3,3))


plot(median_pbo_pred ~ x, ylab="Mosquito mortality pyrethroid-PBO ITN (%)",
     ylim=c(0,1),xlim=c(0,1),pch="",yaxt="n",xaxt="n",
     xlab = "Mosquito mortality pyrethroid-only LLIN (%)",cex.lab=1.2,cex.axis=1.2)
axis(1,at=seq(0,1,0.2),labels=seq(0,100,20),cex.axis=1.2)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20),cex.axis=1.2)

points(pbo_dat$prop_killed_pbo ~ 
         pbo_dat$x,
       cex=pbo_dat$size,pch=19,col=adegenet::transp("purple",0.5))

# polygon(c(x,rev(x)),c(upper_pbo_pred,rev(lower_pbo_pred)),col=adegenet::transp("purple",0.4),border=NA)
# 
lines(median_pbo_pred ~ x,col="purple",lwd=2)
lines(mean_pbo_pred ~ x,col="purple",lwd=2)

# legend("bottomright",legend = c("West hut design",
#                                 "East hut design"),
#        title = "Hut type",
#        pch = c(19,17),
#        col = "purple",
#        lty=1,bty = "n")

## Fitting the beat binomial
data_list_PBO = list(N = nrow(data_PBO),
                     n = data_PBO$n_t,
                     d = data_PBO$d_t,
                     x = data_PBO$x1,
                     
                     M = dim(fine_df)[1],
                     x_m = fine_df$resistance,
                     m = simulated_individuals)


betbin_fit2 <- stan(file = 'R code/stan models/beta_binomial_fit_nets_inv_logit_link.stan',
                    data = data_list_PBO,
                    chains = n_chains,
                    iter = tot_iter,
                    warmup = burnin,
                    init_r = 1e-2,
                    control = betbin_options)

# saveRDS(betbin_fit2,"stan model outputs/ento_pbo_beta_binomial_benefit_inv_link_LancetGH2024.RDS")

#Extract samples
betbin_samples2 <- rstan::extract(betbin_fit2)
# saveRDS(betbin_samples2,"stan model outputs/ento_pbo_beta_binomial_benefit_inv_link_LancetGH2024_extracted_samples.RDS")
mean(betbin_samples2$alpha_1);quantile(betbin_samples2$alpha_1,c(0.05,0.5,0.95))
mean(betbin_samples2$alpha_2);quantile(betbin_samples2$alpha_2,c(0.05,0.5,0.95))
mean(betbin_samples2$k);quantile(betbin_samples2$k,c(0.05,0.5,0.95))

mt_samples2 <- (betbin_samples2$d_m * 1.0) / (simulated_individuals * 1.0)
mt_sorted2 <- apply(mt_samples2, 2, sort)
N_samples2 <- dim(mt_sorted2)[1]
LB_ID2 <- round(N_samples2*0.05)
UB_ID2 <- round(N_samples2*0.95)
fine_df$mt_LB2 <- mt_sorted2[LB_ID2,]
fine_df$mt_UB2 <- mt_sorted2[UB_ID2,]
fine_df$mt_median2 <- apply(mt_samples2, 2, median)

lines(fine_df$mt_UB2 ~ fine_df$resistance,lwd=2,col="purple",lty=2)
lines(fine_df$mt_LB2 ~ fine_df$resistance,lwd=2,col="purple",lty=2)

abline(0,1,lty=2,col="blue")

# lines(colMeans(betbin_samples2$p) ~ fine_df$resistance)

####################################################
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

g2_dat = read.csv("data/pyrethroid-pyrrole nets April 2024.csv",header=TRUE)

g2_dat = subset(g2_dat,g2_dat$is_washed == 0)
g2_dat = subset(g2_dat,g2_dat$Notes == 0)

###########################################
##
## Part 1:
##
## Added benefit of the nets

##Relationship 1 pyrethroid-pyrrole benefit on top of pyrethroid-only nets 
##LOGISTIC BENEFIT
d_t = c(g2_dat$N_dead_comb_B) # number of mosquitoes dying IRS HUTS
n_t = c(g2_dat$N_total_B) # Number of mosquitoes entering IRS huts
x1 = c(g2_dat$N_dead_comb_pyrethroid)/
  c(g2_dat$N_total_pyrethroid)


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

## Fitting a binomial regression to the data
stan_base <- stan(file="R code/stan models/binomial_fit_nets.stan", 
                  data=data_list_g2, 
                  warmup=1000,
                  control = list(adapt_delta = 0.8,
                                 max_treedepth = 20),
                  iter=2000, chains=4)
base <- rstan::extract(stan_base)

# saveRDS(stan_base,"stan model outputs/ento_g2_benefit_LancetGH2024.RDS")

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
mean_G2_pred = 1 / (1 + exp(-mean(base$alpha1) - mean(base$alpha2)*x))

##
## Drawing the figure
##


plot(median_G2_pred ~ x, ylab="Mosquito mortality next generation ITN (%)",
     ylim=c(0,1),xlim=c(0,1),pch="",yaxt="n",xaxt="n",
     xlab = "Mosquito mortality pyrethroid-only LLIN (%)",cex.lab=1.2,cex.axis=1.2)
axis(1,at=seq(0,1,0.2),labels=seq(0,100,20),cex.axis=1.2)
axis(2,las=2,at=seq(0,1,0.2),labels=seq(0,100,20),cex.axis=1.2)

g2_dat$x2 = g2_dat$N_dead_comb_pyrethroid/g2_dat$N_total_pyrethroid
g2_dat$prop_killed_g2 = g2_dat$N_dead_comb_B/g2_dat$N_total_B
g2_dat$size = sqrt(g2_dat$N_total_B)/4

# points(g2_dat$xx_B ~ g2_dat$x)
points(g2_dat$prop_killed_g2 ~ 
         g2_dat$x2,
       cex=g2_dat$size,pch=19,col=adegenet::transp("aquamarine2",0.5))

# polygon(c(x,rev(x)),c(upper_G2_pred,rev(lower_G2_pred)),col=adegenet::transp("darkgreen",0.4),border=NA)


lines(median_G2_pred ~ x,col="darkgreen",lwd=2,lty=1)
lines(mean_G2_pred ~ x,col="darkgreen",lwd=2,lty=1)
lines(fine_df$mt_UB2 ~ fine_df$resistance,lwd=2,col="purple",lty=2)
lines(fine_df$mt_LB2 ~ fine_df$resistance,lwd=2,col="purple",lty=2)
lines(median_pbo_pred ~ x,col=adegenet::transp("purple",0.5),lwd=2)


abline(0,1,lty=2,col="blue")

## Fitting a beta-binomial 

##############################################
##
## Adding uncertainty by fitting the beta_binomial
library(shinystan)

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
                         n = data_G2$n_t,
                         d = data_G2$d_t,
                         x = data_G2$x1,
                         
                         M = dim(fine_df)[1],
                         x_m = fine_df$resistance,
                         m = simulated_individuals)


betbin_fit <- stan(file = 'R code/stan models/beta_binomial_fit_nets_inv_logit_link.stan',
                   data = data_list_g2_beta,
                   chains = n_chains,
                   iter = tot_iter,
                   warmup = burnin,
                   init_r = 1e-2,
                   control = betbin_options)

# saveRDS(betbin_fit,"stan model outputs/ento_g2_beta_binomial_benefit_inv_link_LancetGH2024.RDS")


#Investigate diagnostic plots
diagnostic_plots <- FALSE
if (diagnostic_plots) {
  shinystan::launch_shinystan(betbin_fit)
}

#Extract samples
betbin_samples <- rstan::extract(betbin_fit)

## overlay the original
mean(betbin_samples$alpha_1)
mean(betbin_samples$alpha_2)

mt_samples <- (betbin_samples$d_m * 1.0) / (simulated_individuals * 1.0)
mt_sorted <- apply(mt_samples, 2, sort)
N_samples <- dim(mt_sorted)[1]
LB_ID <- round(N_samples*0.05)
UB_ID <- round(N_samples*0.95)
fine_df$mt_LB <- mt_sorted[LB_ID,]
fine_df$mt_UB <- mt_sorted[UB_ID,]
fine_df$mt_median <- apply(mt_samples, 2, median)

lines(fine_df$mt_UB ~ fine_df$resistance,lwd=2,col="aquamarine3",lty=2)
lines(fine_df$mt_LB ~ fine_df$resistance,lwd=2,col="aquamarine3",lty=2)

points(pbo_dat$prop_killed_pbo ~ 
         pbo_dat$x,
       cex=pbo_dat$size,pch=1,col=adegenet::transp("purple",0.5))


# Useful if plotting on the same figure
legend("bottomright",legend=c("Pyr-pyrrole (72 hr)",
                              "Pyr-PBO (24 hr)"),ncol=1,
       pch=c(19,1),col=c("aquamarine3","purple"),lty=1,bty="n")



par(xpd=NA,cex = 1.11)

text(x = -1.85, y = 1.15,"A")
text(x = -.32, y = 1.15,"B")

