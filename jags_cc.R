#setwd('/work/MiaJacobsen#3200/exam')

install.packages("pacman")
pacman::p_load(R2jags, parallel, polspline, tidyverse, glue, patchwork, bayestestR)


########### INFO ABOUT THE GAME ##########
groupSize <- 4
ntrials <- 10
pi <- 1.6 # multiplication factor in game
ntokens <- 20
vals <- seq(0,ntokens,1) #possible values to contribute - from 0 to 20 tokens


########## LOADING DATA AND MERGING ###########
### trust measures
# note: we have summed 'somewhat agree' and 'completely agree'
V23 = c(24.6, 39.1, 51.2, 30.0, 24.5, 4.8, 49.3, 28.0, 33.8, 45.6)
V126 = c(66.4, 77.9, 85.4, 77.7, 68, 73.1, 83.7, 72.2, 74.3, 80.8)
V128 = c(13.5, 39.4, 50.7, 45.6, 18.3, 14.5, 9.9, 14.9, 25.2, 48.2)
# confidence in government
V138 = c(42.8, 36.8, 65.1, 32.4, 23.6, 61.7, 87.6, 45.6, 22.7, 38.9)

source('pgg_funcs.R')
df <- PGG_data_setup('/work/216377/Module5/data/HerrmannThoeniGaechterDATA.csv',
                     V23, V126, V128, V138)

########## DATA PREPARATION ########## 
ngroups <- length(unique(df$groupid))

## group average contribution on each trial 
Gga_np <- df %>% 
  group_by(groupid, period) %>% 
  summarise(mean_send = mean(senderscontribution)) 

## group average contribution on each trial for each subject
# this fixes problem with some groups not having 4 subject on each trial 
Ggas_np <- Gga_np %>% 
  slice(rep(row_number(), 4)) %>% 
  arrange(groupid, period)

## Nation for every group
Nation_df <- df %>% 
  group_by(groupid) %>% 
  summarise(nation = mean(nation))

## calculating winnings 
source('pgg_funcs.R')
winnings <- calculate_winnings_df(df, pi, Nation_df)


### making arrays for JAGS
Ggas <- array(data = Ggas_np$mean_send, 
              dim = c(groupSize, ntrials, ngroups))
Nation <- array(data = unlist(Nation_df$nation))

########## CC MODEL ONE VARIBALE AT A TIME ###########
nnations <- length(unique(df$nation))
ngroups <- length(unique(df$groupid))

Ga <- Ggas

df <- arrange(df, groupid, period)
c <- array(data = df$senderscontribution, dim = c(groupSize, ntrials, ngroups))

data <- list("groupSize", "ngroups", "ntrials", "nnations","c","Ga","X","Nation","invSigma") 
params <- c("beta0_alpha","betaX_alpha","beta0_rho","betaX_rho","beta0_omega","betaX_omega") 

# V23
X <- V23
X <- (X-mean(X))/sd(X)
invSigma <- solve(t(X)%*%X) # required for JZS priors

CC.samples.23 <- jags.parallel(data, inits=NULL, params,
                               model.file ="CC_corr.txt",
                               n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster=4)
# V126
X <- V126
X <- (X-mean(X))/sd(X)
invSigma <- solve(t(X)%*%X) # required for JZS priors

CC.samples.126 <- jags.parallel(data, inits=NULL, params,
                                model.file ="CC_corr.txt",
                                n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster=4)

# V128
X <- V128
X <- (X-mean(X))/sd(X)
invSigma <- solve(t(X)%*%X) # required for JZS priors

CC.samples.128 <- jags.parallel(data, inits=NULL, params,
                                model.file ="CC_corr.txt",
                                n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster=4)

# V138
X <- V138
X <- (X-mean(X))/sd(X)
invSigma <- solve(t(X)%*%X) # required for JZS priors

CC.samples.138 <- jags.parallel(data, inits=NULL, params,
                                model.file ="CC_corr.txt",
                                n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster=4)

####### plotting ########
# V23
cc_samples_23 <- data.frame(CC.samples.23$BUGSoutput$sims.list$betaX_alpha,
                            CC.samples.23$BUGSoutput$sims.list$betaX_omega,
                            CC.samples.23$BUGSoutput$sims.list$betaX_rho)
names(cc_samples_23) <- c("Initial belief", "Belief learning rate", 
                          "Cooperation preference")

a1 <-  cc_samples_23 %>% 
  ggplot()+
  geom_density(aes(x = `Initial belief`))+
  geom_point(aes(x = mean(`Initial belief`), y=0), 
             col ="black")+
  geom_segment(aes(x = ci(`Initial belief`, ci = 0.95)$CI_low, y = 0, 
                   xend = ci(`Initial belief`, ci = 0.95)$CI_high, yend = 0), 
               col = 'black')+
  geom_vline(xintercept = 0.0, linetype = 'dotted', color = 'black', size = .4)+
  theme_bw()+
  scale_x_continuous(limits = c(-1, 1))+
  scale_y_continuous(limits = c(0,11))+
  ggtitle('Most people can be trusted')


a2 <-  cc_samples_23 %>% 
  ggplot()+
  geom_density(aes(x = `Belief learning rate`))+
  geom_point(aes(x = mean(`Belief learning rate`), y=0), 
             col ="black")+
  geom_segment(aes(x = ci(`Belief learning rate`, ci = 0.95)$CI_low, y = 0, 
                   xend = ci(`Belief learning rate`, ci = 0.95)$CI_high, yend = 0), 
               col = 'black')+
  geom_vline(xintercept = 0.0, linetype = 'dotted', color = 'black', size = .4)+
  theme_bw()+
  scale_x_continuous(limits = c(-1, 1))+
  scale_y_continuous(limits = c(0,11))

a3 <-  cc_samples_23 %>% 
  ggplot()+
  geom_density(aes(x = `Cooperation preference`))+
  geom_point(aes(x = mean(`Cooperation preference`), y=0), 
             col ="black")+
  geom_segment(aes(x = ci(`Cooperation preference`, ci = 0.95)$CI_low, y = 0, 
                   xend = ci(`Cooperation preference`, ci = 0.95)$CI_high, yend = 0), 
               col = 'black')+
  geom_vline(xintercept = 0.0, linetype = 'dotted', color = 'black', size = .4)+
  theme_bw()+
  scale_x_continuous(limits = c(-1, 1))+
  scale_y_continuous(limits = c(0,11))

a1+a2+a3

#V126
cc_samples_126 <- data.frame(CC.samples.126$BUGSoutput$sims.list$betaX_alpha,
                             CC.samples.126$BUGSoutput$sims.list$betaX_omega,
                             CC.samples.126$BUGSoutput$sims.list$betaX_rho)
names(cc_samples_126) <- c("Initial belief", "Belief learning rate", 
                           "Cooperation preference")

b1 <-  cc_samples_126 %>% 
  ggplot()+
  geom_density(aes(x = `Initial belief`))+
  geom_point(aes(x = mean(`Initial belief`), y=0), 
             col ="black")+
  geom_segment(aes(x = ci(`Initial belief`, ci = 0.95)$CI_low, y = 0, 
                   xend = ci(`Initial belief`, ci = 0.95)$CI_high, yend = 0), 
               col = 'black')+
  geom_vline(xintercept = 0.0, linetype = 'dotted', color = 'black', size = .4)+
  theme_bw()+  
  scale_x_continuous(limits = c(-1, 1))+
  scale_y_continuous(limits = c(0,11))+
  ggtitle('Trust in your neighborhood')


b2 <-  cc_samples_126 %>% 
  ggplot()+
  geom_density(aes(x = `Belief learning rate`))+
  geom_point(aes(x = mean(`Belief learning rate`), y=0), 
             col ="black")+
  geom_segment(aes(x = ci(`Belief learning rate`, ci = 0.95)$CI_low, y = 0, 
                   xend = ci(`Belief learning rate`, ci = 0.95)$CI_high, yend = 0), 
               col = 'black')+
  geom_vline(xintercept = 0.0, linetype = 'dotted', color = 'black', size = .4)+
  theme_bw() +
  scale_x_continuous(limits = c(-1, 1))+
  scale_y_continuous(limits = c(0,11))

b3 <-  cc_samples_126 %>% 
  ggplot()+
  geom_density(aes(x = `Cooperation preference`))+
  geom_point(aes(x = mean(`Cooperation preference`), y=0), 
             col ="dodgerblue3")+
  geom_segment(aes(x = ci(`Cooperation preference`, ci = 0.95)$CI_low, y = 0, 
                   xend = ci(`Cooperation preference`, ci = 0.95)$CI_high, yend = 0), 
               col = 'dodgerblue3')+
  geom_vline(xintercept = 0.0, linetype = 'dotted', color = 'black', size = .4)+
  theme_bw()+
  scale_x_continuous(limits = c(-1, 1))+
  scale_y_continuous(limits = c(0,11))

b1+b2+b3

#V128
cc_samples_128 <- data.frame(CC.samples.128$BUGSoutput$sims.list$betaX_alpha,
                             CC.samples.128$BUGSoutput$sims.list$betaX_omega,
                             CC.samples.128$BUGSoutput$sims.list$betaX_rho)
names(cc_samples_128) <- c("Initial belief", "Belief learning rate", 
                           "Cooperation preference")

c1 <-  cc_samples_128 %>% 
  ggplot()+
  geom_density(aes(x = `Initial belief`))+
  geom_point(aes(x = mean(`Initial belief`), y=0), 
             col ="black")+
  geom_segment(aes(x = ci(`Initial belief`, ci = 0.95)$CI_low, y = 0, 
                   xend = ci(`Initial belief`, ci = 0.95)$CI_high, yend = 0), 
               col = 'black')+
  geom_vline(xintercept = 0.0, linetype = 'dotted', color = 'black', size = .4)+
  theme_bw()+
  scale_x_continuous(limits = c(-1, 1))+
  scale_y_continuous(limits = c(0,11))+
  ggtitle('Trust in people you meet for the 1st time')


c2 <-  cc_samples_128 %>% 
  ggplot()+
  geom_density(aes(x = `Belief learning rate`))+
  geom_point(aes(x = mean(`Belief learning rate`), y=0), 
             col ="black")+
  geom_segment(aes(x = ci(`Belief learning rate`, ci = 0.95)$CI_low, y = 0, 
                   xend = ci(`Belief learning rate`, ci = 0.95)$CI_high, yend = 0), 
               col = 'black')+
  geom_vline(xintercept = 0.0, linetype = 'dotted', color = 'black', size = .4)+
  theme_bw()+
  scale_x_continuous(limits = c(-1, 1))+
  scale_y_continuous(limits = c(0,11))

c3 <-  cc_samples_128 %>% 
  ggplot()+
  geom_density(aes(x = `Cooperation preference`))+
  geom_point(aes(x = mean(`Cooperation preference`), y=0), 
             col ="dodgerblue3")+
  geom_segment(aes(x = ci(`Cooperation preference`, ci = 0.95)$CI_low, y = 0, 
                   xend = ci(`Cooperation preference`, ci = 0.95)$CI_high, yend = 0), 
               col = 'dodgerblue3')+
  geom_vline(xintercept = 0.0, linetype = 'dotted', color = 'black', size = .4)+
  theme_bw()+
  scale_x_continuous(limits = c(-1, 1))+
  scale_y_continuous(limits = c(0,11))

c1+c2+c3

#V138
cc_samples_138 <- data.frame(CC.samples.138$BUGSoutput$sims.list$betaX_alpha,
                             CC.samples.138$BUGSoutput$sims.list$betaX_omega,
                             CC.samples.138$BUGSoutput$sims.list$betaX_rho)
names(cc_samples_138) <- c("Initial belief", "Belief learning rate", 
                           "Cooperation preference")

d1 <-  cc_samples_138 %>% 
  ggplot()+
  geom_density(aes(x = `Initial belief`))+
  geom_point(aes(x = mean(`Initial belief`), y=0), 
             col ="black")+
  geom_segment(aes(x = ci(`Initial belief`, ci = 0.95)$CI_low, y = 0, 
                   xend = ci(`Initial belief`, ci = 0.95)$CI_high, yend = 0), 
               col = 'black')+
  geom_vline(xintercept = 0.0, linetype = 'dotted', color = 'black', size = .4)+
  theme_bw()+
  scale_x_continuous(limits = c(-1, 1))+
  scale_y_continuous(limits = c(0,11))+
  ggtitle('Confidence in government')


d2 <-  cc_samples_138 %>% 
  ggplot()+
  geom_density(aes(x = `Belief learning rate`))+
  geom_point(aes(x = mean(`Belief learning rate`), y=0), 
             col ="black")+
  geom_segment(aes(x = ci(`Belief learning rate`, ci = 0.95)$CI_low, y = 0, 
                   xend = ci(`Belief learning rate`, ci = 0.95)$CI_high, yend = 0), 
               col = 'black')+
  geom_vline(xintercept = 0.0, linetype = 'dotted', color = 'black', size = .4)+
  theme_bw()+
  scale_x_continuous(limits = c(-1, 1))+
  scale_y_continuous(limits = c(0,11))

d3 <-  cc_samples_138 %>% 
  ggplot()+
  geom_density(aes(x = `Cooperation preference`))+
  geom_point(aes(x = mean(`Cooperation preference`), y=0), 
             col ="black")+
  geom_segment(aes(x = ci(`Cooperation preference`, ci = 0.95)$CI_low, y = 0, 
                   xend = ci(`Cooperation preference`, ci = 0.95)$CI_high, yend = 0), 
               col = 'black')+
  geom_vline(xintercept = 0.0, linetype = 'dotted', color = 'black', size = .4)+
  theme_bw()+
  scale_x_continuous(limits = c(-1, 1))+
  scale_y_continuous(limits = c(0,11))

d1+d2+d3

(a1+a2+a3)/(b1+b2+b3)/(c1+c2+c3)/(d1+d2+d3)

###### CC MODEL ALL VARIABLES ########
fullControl.CC <- function (X1,X2,X3,X4,groupSize,ngroups,ntrials,nnations,c,Ga,Nation) {
  
  # standardise covariates
  X1 <- (X1-mean(X1))/sd(X1)
  X2 <- (X2-mean(X2))/sd(X2)
  X3 <- (X3-mean(X3))/sd(X3)
  X4 <- (X4-mean(X4))/sd(X4)
  
  X <- cbind(X1,X2,X3,X4)
  V <- solve(t(X)%*%X)
  
  data <- list("groupSize", "ngroups", "ntrials", "nnations","c","Ga","X","V","Nation") #data inputted into jags
  params <- c("beta0_alpha","betaX_alpha","beta0_rho","betaX_rho","beta0_omega","betaX_omega") #parameters we'll track in jags
  
  # - run jags code
  CC.samples <- jags(data, inits=NULL, params,
                     model.file ="CC_corr_partial.txt",
                     n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1)
  
  return(list(CC.samples))
  
}

Full_control.CC <- fullControl.CC(V23,V126,V128,V138,groupSize,ngroups,ntrials,nnations,c,Ga,Nation)

########### plotting #######
full_cc_df <- data.frame(Full_control.CC[[1]]$BUGSoutput$summary)
full_cc_df <- cbind(parameter = rownames(full_cc_df), full_cc_df)
rownames(full_cc_df) <- 1:nrow(full_cc_df)

labels = c("Confidence in government",
           "Trust in people you meet for the first time",
           "Trust in neighborhood",
           "Most people can be trusted")

cc1 <- full_cc_df %>% 
  filter(grepl("betaX_alpha",parameter)) %>% 
  ggplot(aes(x = mean, y = parameter, colour = X2.5.> 0 | X97.5. < 0))+
  geom_pointrange(aes(xmin = X2.5., xmax = X97.5.), size = 0.8)+
  geom_vline(xintercept = 0.0, linetype = 'dotted', color = 'black', size = .4)+
  scale_colour_manual(values = c("black", "dodgerblue3"))+
  scale_y_discrete(labels = labels, position = "right")+
  labs(y = '',
       x = 'Standardized Estimates',
       title = 'Initial Belief')+
  theme_bw()+
  scale_x_continuous(limits = c(-6, 6))+
  theme(legend.position = "none")

# COMMENT ON THE FACT THAT SCALES ARE DIFFERENT
cc2 <- full_cc_df %>% 
  filter(grepl("betaX_omega",parameter)) %>% 
  ggplot(aes(x = mean, y = parameter, colour = X2.5.> 0 | X97.5. < 0))+
  geom_pointrange(aes(xmin = X2.5., xmax = X97.5.), size = 0.8)+
  geom_vline(xintercept = 0.0, linetype = 'dotted', color = 'black', size = .4)+
  scale_colour_manual(values = c("black", "dodgerblue3"))+
  scale_y_discrete(labels = labels, position = "right")+
  labs(y = '',
       x = 'Standardized Estimates',
       title = 'Belief Learning Rate')+
  theme_bw()+
  scale_x_continuous(limits = c(-1, 1))+
  theme(legend.position = "none")

cc3 <- full_cc_df %>% 
  filter(grepl("betaX_rho",parameter)) %>% 
  ggplot(aes(x = mean, y = parameter, colour = X2.5.> 0 | X97.5. < 0))+
  geom_pointrange(aes(xmin = X2.5., xmax = X97.5.), size = 0.8)+
  geom_vline(xintercept = 0.0, linetype = 'dotted', color = 'black', size = .4)+
  scale_colour_manual(values = c("black", "dodgerblue3"))+
  scale_y_discrete(labels = labels, position = "right")+
  labs(y = '',
       x = 'Standardized Estimates',
       title = 'Cooperation Preference')+
  theme_bw()+
  scale_x_continuous(limits = c(-1, 1))+
  theme(legend.position = "none")

cc1/cc2/cc3
