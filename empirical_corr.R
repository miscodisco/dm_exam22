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
Nation_df <- df %>% 
  group_by(groupid) %>% 
  summarise(nation = mean(nation))

## calculating winnings 
source('pgg_funcs.R')
winnings <- calculate_winnings_df(df, pi, Nation_df)


########## EMPIRICAL CORRELATION ##########
# for plotting purposes
nation_df <- df %>% 
  distinct(nation, .keep_all = T)
group_df <- df %>% 
  distinct(groupid, .keep_all = T)

# mean winnings at nation-level with standard deviation
empirical_win <- winnings %>%
  group_by(nation) %>%
  summarise(mean_win = mean(win),
            sd_win = sd(win)) %>%
  mutate(mean_minus_sd = mean_win - sd_win,
         mean_plus_sd = mean_win + sd_win)

# for plotting nation-level winnings over each trust measure
empirical_win <- merge(empirical_win, nation_df)

# winnings over V23 
p1 <- empirical_win %>% 
  ggplot(aes(V23, mean_win))+
  geom_pointrange(aes(ymin = mean_minus_sd, ymax = mean_plus_sd))+
  theme_bw()+
  geom_smooth(method = 'lm', alpha = 0, col = 'red')+
  xlab("'Most people can be trusted'")+
  ylab('Winnings')+
  ggtitle('A)')

# winnings over V126
p2 <- empirical_win %>% 
  ggplot(aes(V126, mean_win))+
  geom_pointrange(aes(ymin = mean_minus_sd, ymax = mean_plus_sd))+
  geom_smooth(method = 'lm', alpha = 0, col = 'blue')+
  theme_bw()+
  xlab("'Trust in your neighborhood'")+
  ylab('Winnings')+
  ggtitle('B)')

# winnings over V128
p3 <- empirical_win %>% 
  ggplot(aes(V128, mean_win))+
  geom_pointrange(aes(ymin = mean_minus_sd, ymax = mean_plus_sd))+
  geom_smooth(method = 'lm', alpha = 0, col = 'dark green')+
  theme_bw()+
  xlab("'Trust in people you meet for the first time'")+
  ylab('Winnings')+
  ggtitle('C)')

# winnings over V138
p4 <- empirical_win %>% 
  ggplot(aes(V138, mean_win))+
  geom_pointrange(aes(ymin = mean_minus_sd, ymax = mean_plus_sd))+
  theme_bw()+
  geom_smooth(method = 'lm', alpha = 0, col = 'black')+
  xlab("'Confidence in government'")+
  ylab('Winnings')+
  ggtitle('D)')

# all measures 
colors <- c("A)" = "red", "B)" = "blue", "C)" = 'dark green', "D)" = 'black')
p5 <- empirical_win %>% 
  mutate(V23 = (V23-mean(V23))/sd(V23),
         V126 = (V126-mean(V126))/sd(V126),
         V128 = (V128-mean(V128))/sd(V128),
         V138 = (V138-mean(V138))/sd(V138),) %>% 
  ggplot()+
  geom_smooth(aes(V23, mean_win, col = "A)"), fill = 'red', method = 'lm', alpha = 0.1)+
  geom_smooth(aes(V126, mean_win, col = "B)"), fill = 'blue', method = 'lm', alpha = 0.1)+
  geom_smooth(aes(V128, mean_win, col = "C)"), fill = 'dark green', method = 'lm', alpha = 0.1)+
  geom_smooth(aes(V138, mean_win, col = "D)"), fill = 'grey', method = 'lm', alpha = 0.1)+
  theme_bw()+
  labs(x ='Standardized estimates',
       y = 'Winnings',
       color = "Legend",
       title = 'E)')+
  scale_color_manual(values = colors)

(p1+p2+p3+p4)/p5






