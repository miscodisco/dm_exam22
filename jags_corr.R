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
Nation <- array(data = unlist(Nation_df$nation))
# also standardizing winnings
Y = array((winnings$win - mean(winnings$win)) / sd(winnings$win))

## looking at problematic groups
df_gga <- merge(Gga_np, df, all = TRUE)
count_probs <- df_gga %>%
  group_by(groupid) %>%
  count(subjectid) %>%
  filter(n != 10) %>%
  mutate(diff = abs(n - 10))

problem_groups <- unique(count_probs$groupid)
length(problem_groups)
sum(count_probs$diff)


########### JAGS - WIN CORRELATION #######

ngroups = length(unique(df$groupid))
nnations = length(unique(df$nation))

data <- list("ngroups", "Y", "nnations","X","Nation","invSigma") 
params <- c("beta0","betaX") 

# run jags code
X = array((V23-mean(V23))/sd(V23))
invSigma <- solve(t(X)%*%X)
win.samples.23 <- jags.parallel(data, inits=NULL, params,
                                model.file ="/work/MiaJacobsen#3200/livecoding/win_corr.txt",
                                n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster=4)

X = array((V126-mean(V126))/sd(V126))
invSigma <- solve(t(X)%*%X)
win.samples.126 <- jags.parallel(data, inits=NULL, params,
                                 model.file ="/work/MiaJacobsen#3200/livecoding/win_corr.txt",
                                 n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster=4)

X = array((V128-mean(V23))/sd(V128))
invSigma <- solve(t(X)%*%X)
win.samples.128 <- jags.parallel(data, inits=NULL, params,
                                 model.file ="/work/MiaJacobsen#3200/livecoding/win_corr.txt",
                                 n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster=4)

X = array((V138-mean(V23))/sd(V138))
invSigma <- solve(t(X)%*%X)
win.samples.138 <- jags.parallel(data, inits=NULL, params,
                                 model.file ="/work/MiaJacobsen#3200/livecoding/win_corr.txt",
                                 n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster=4)

################# JAGS - PARTIAL CORRELATION FOR ALL TRUST VARIABLES ################# 
fullControl.win <- function (X1,X2,X3,X4,ngroups,nnations,winnings,Nation) {
  
  # standardise covariates
  X1 <- array((X1-mean(X1))/sd(X1))
  X2 <- array((X2-mean(X2))/sd(X2))
  X3 <- array((X3-mean(X3))/sd(X3))
  X4 <- array((X4-mean(X4))/sd(X4))
  
  X <- cbind(X1,X2,X3,X4)
  V <- solve(t(X)%*%X)
  
  Y = array((winnings$win - mean(winnings$win)) / sd(winnings$win))
  
  data <- list("ngroups", "Y", "nnations","X","V","Nation") #data inputted into jags
  params <- c("beta0","betaX","prior_T") #parameters we'll track in jags
  
  # - run jags code
  win.samples <- jags(data, inits=NULL, params,
                      model.file ="win_corr_partial.txt",
                      n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1)
  
  return(list(win.samples))
  
}

Full_control.win <- fullControl.win(V23, V126, V128, V138, ngroups,nnations,winnings,Nation)


############# PLOTTING ################
# win correlation for V23
samples_23 <- data.frame(win.samples.23$BUGSoutput$sims.list$betaX)
names(samples_23) <- "V1"
summary_23 <- data.frame(win.samples.23$BUGSoutput$summary)

w1 <- samples_23 %>% 
  ggplot(aes(x = V1))+
  geom_density()+
  geom_point(aes(x=summary_23[2,1], y=0), 
             col ="black")+
  geom_segment(aes(x = summary_23[2,3], y = 0, 
                   xend = summary_23[2,7], yend = 0), 
               col = 'black')+
  geom_vline(xintercept = 0.0, linetype = 'dotted', color = 'black', size = .4)+
  theme_bw()+
  scale_x_continuous(limits = c(-0.41, 0.41))+
  scale_y_continuous(limits = c(0,6.5))+
  labs(x = "'Most people can be trusted', standardized",
       title = 'A)')

# win correlation for V126
samples_126 <- data.frame(win.samples.126$BUGSoutput$sims.list$betaX)
names(samples_126) <- "V1"
summary_126 <- data.frame(win.samples.126$BUGSoutput$summary)

w2 <- samples_126 %>% 
  ggplot(aes(x = V1))+
  geom_density()+
  geom_point(aes(x=summary_126[2,1], y=0), 
             col ="black")+
  geom_segment(aes(x = summary_126[2,3], y = 0, 
                   xend = summary_126[2,7], yend = 0), 
               col = 'black')+
  geom_vline(xintercept = 0.0, linetype = 'dotted', color = 'black', size = .4)+
  theme_bw()+
  scale_x_continuous(limits = c(-0.41, 0.41))+
  scale_y_continuous(limits = c(0,6.5))+
  labs(x = "'Trust in neighborhood', standardized",
       title = 'B)')

# win correlation for V128
samples_128 <- data.frame(win.samples.128$BUGSoutput$sims.list$betaX)
names(samples_128) <- "V1"
summary_128 <- data.frame(win.samples.128$BUGSoutput$summary)

w3 <- samples_128 %>% 
  ggplot(aes(x = V1))+
  geom_density()+
  geom_point(aes(x=summary_128[2,1], y=0), 
             col ="black")+
  geom_segment(aes(x = summary_128[2,3], y = 0, 
                   xend = summary_128[2,7], yend = 0), 
               col = 'black')+
  geom_vline(xintercept = 0.0, linetype = 'dotted', color = 'black', size = .4)+
  theme_bw()+
  scale_x_continuous(limits = c(-0.41, 0.41))+
  scale_y_continuous(limits = c(0,6.5))+
  labs(x = "'Trust in people you meet for the first time', standardized",
       title = 'C)')

# win correlation for V138
samples_138 <- data.frame(win.samples.138$BUGSoutput$sims.list$betaX)
names(samples_138) <- "V1"
summary_138 <- data.frame(win.samples.138$BUGSoutput$summary)

w4 <- samples_138 %>% 
  ggplot(aes(x = V1))+
  geom_density()+
  geom_point(aes(x=summary_138[2,1], y=0), 
             col ="black")+
  geom_segment(aes(x = summary_138[2,3], y = 0, 
                   xend = summary_138[2,7], yend = 0), 
               col = 'black')+
  geom_vline(xintercept = 0.0, linetype = 'dotted', color = 'black', size = .4)+
  theme_bw()+
  scale_x_continuous(limits = c(-0.41, 0.41))+
  scale_y_continuous(limits = c(0,6.5))+
  labs(x = "'Confidence in government', standardized",
       title = 'D)')

##### partial correlation for all trust measures
# summarize posterior estimates for controlled winnings model
win.model.summary <- cbind(Full_control.win[[1]]$BUGSoutput$summary[2:6,3], #2.5%
                           Full_control.win[[1]]$BUGSoutput$summary[2:6,1], #mean
                           Full_control.win[[1]]$BUGSoutput$summary[2:6,7]) #97.5%


win_summary <- data.frame(Full_control.win[[1]]$BUGSoutput$summary)
win_summary <- cbind(parameter = rownames(win_summary), win_summary)
rownames(win_summary) <- 1:nrow(win_summary)

parameters = win_summary$parameter[2:5]

labels = c("Confidence in government",
           "Trust in people you meet for the first time",
           "Trust in neighborhood",
           "Most people can be trusted")

w5 <- win_summary %>% 
  filter(parameter %in% parameters) %>% 
  mutate(parameter = fct_relevel(parameter, 
                                 "betaX[4]", "betaX[3]", "betaX[2]", 
                                 "betaX[1]")) %>%
  ggplot(aes(x = mean, y = parameter, colour = X2.5.> 0 | X97.5. < 0))+
  geom_pointrange(aes(xmin = X2.5., xmax = X97.5.), size = 0.8)+
  geom_vline(xintercept = 0.0, linetype = 'dotted', color = 'black', size = .4)+
  scale_colour_manual(values = c("black", "dodgerblue3"))+
  scale_y_discrete(labels = labels, position = "right")+
  labs(y = '',
       x = 'Standardized Estimates')+
  theme_bw()+
  scale_x_continuous(limits = c(-0.9, 0.9))+
  theme(legend.position = "none")
w5
(w1+w2+w3+w4)

