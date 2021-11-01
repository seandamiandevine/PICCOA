#SETUP ####
# function to load packages or install if not already installed
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("dplyr", "ggplot2", "ggpubr", "lme4", "reshape2", "broom.mixed", "rtdists", 'effsize', "tidyr", "purrr")
ipak(packages)

rm(list = ls())
options(pillar.sigfig = 4)

se <- function(x) sd(x, na.rm = T)/sqrt(length(x)) # custom se function 

# Load from directory based on system type 
if(as.list(Sys.info())$sysname == 'Linux') dir = '~/Desktop/PICCAgeing/Analysis&Data/' else dir = "C:/Users/Sean/Desktop/LDM Lab/PICCAgeing" 
setwd(dir)

# LOAD EXP 1 DOTS DATA ####
files = list.files(path='Analysis&Data/Exp1/Dots',pattern = ".csv")
d <- do.call(rbind,  lapply(paste0(dir, '/Analysis&Data/Exp1/Dots/', files), read.csv, header = T, stringsAsFactors = F))

# CLEAN  ####
d <- 
  d %>% 
  filter(block != 'Practice') %>% 
  mutate(trialintensity = trialintensity-154, 
         response = ifelse(response == 'a', 1, 0), 
         id = factor(id), 
         age_group = factor(ifelse(age < 60, 'Young', 'Old'), levels = c("Young", 'Old')), 
         condition = recode(factor(condition),  
                      '0' = 'Stable', 
                      '1' = 'Decreasing'), 
         block = factor(block, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 
                                          16, 17)),
         colour0 = trialintensity/100, 
         trial0 = trial/max(trial)
         )

# BIN DATA ####
d <-
  d %>% 
  mutate(strength_bin = ntile(trialintensity, 5), 
         timebin4 = ntile(trial, 4))

# summarise/visualise if you want 
binsum <- 
  d %>% 
  mutate(timebin4 = factor(timebin4)) %>%
  filter(timebin4 == 1 | timebin4 == 4) %>% 
  group_by(id, age_group, condition, strength_bin, timebin4) %>% 
  summarise(pblue = mean(response)) %>% 
  ggplot(aes(x = strength_bin, y = pblue, colour=timebin4)) + 
  stat_summary(fun.y = mean, geom='point', position=position_dodge(0.9)) + 
  stat_summary(fun.data = mean_se, geom='errorbar', position=position_dodge(0.9)) + 
  stat_summary(fun.y = mean, geom='line', position=position_dodge(0.9)) + 
  facet_grid(condition ~ age_group)
  

# DIFFUSION ANALYSIS W/ RTDISTS ####
# Following https://cran.r-project.org/web/packages/rtdists/vignettes/reanalysis_rr98.html
# Maximum likelihood estimation
d$strength_bin <- factor(d$strength_bin)
ndrift <- length(levels(d$strength_bin))

d$response_ddm <- ifelse(d$response == 1, 'upper', 'lower')



# Function that creates random start values 


# Function that tries different random start values until it works:


# Fit! 


# VISUALISE ####
base_pars_plot <-
  d.fit[c(1, 3:9)] %>% 
  melt(id.vars = c('id', 'age_group', 'condition')) %>% 
  ggplot(aes(y = value, x = age_group, fill = condition)) + 
  stat_summary(fun.y = mean, geom='bar', position = position_dodge(0.9)) + 
  stat_summary(fun.data = mean_se, geom='errorbar', position = position_dodge(0.9)) + 
  geom_point(position = position_dodge(0.9), alpha=0.2) + 
  facet_wrap(~variable, scales = 'free_y') + 
  labs(x = '', y='Parameter Value', fill = 'Condition') + 
  theme_classic()

drift_rate_plot <- 
  d.fit[c(1, 3, 4, 10:14)] %>%
  melt(id.vars = c('id', 'age_group', 'condition')) %>% 
  mutate(variable = as.numeric(gsub('v_', "", variable, fixed = T))) %>% 
  ggplot(aes(x = variable, y = value, color = condition)) + 
  stat_summary(fun.y = mean, geom='point', position = position_dodge(0.9)) + 
  stat_summary(fun.data = mean_se, geom='errorbar', position = position_dodge(0.9)) + 
  stat_summary(fun.y = mean, geom='line', position = position_dodge(0.9)) + 
  geom_point(position = position_dodge(0.9), alpha=0.2) + 
  facet_wrap(~age_group) + 
  labs(x = 'Stimulus Strength', y='Drift Rate', color = 'Condition') + 
  theme_classic()

ggarrange(base_pars_plot, drift_rate_plot, labels = LETTERS[1:2],
          common.legend = T, legend = "bottom")

# If you divide subjects by parameter difference, instead of age difference, is PICC recovered?
d.fit$a_split <- ntile(d.fit$a, 2) # median split on threshold
for(id in unique(d.fit$id)) { 
  
  d[d$id==id, 'a_split'] <- d.fit$a_split[d.fit$id==id]
  
}

a_picc_plot_all <- 
  d %>%
  mutate(a_split = factor(ifelse(a_split==1, 'Low a', 'High a'), levels = c('Low a', 'High a')), 
         timebin4 = factor(timebin4)) %>% 
  group_by(a_split, condition, trialintensity, timebin4) %>% 
  summarise(pBlue = length(response[response == 1])/length(response)) %>%
  rename(colour = trialintensity, trialbins = timebin4) %>% 
  filter(trialbins == 1 | trialbins == 4) %>% 
  ggplot(aes(x = colour, y = pBlue, color = trialbins)) + geom_point() + 
  geom_line(stat="smooth",method = "glm", method.args = list(
    family="binomial"), se=FALSE,size=1.8,alpha=.7) + 
  ylab('% Dots Judged as Blue') + xlab('') +
  scale_colour_manual(labels=c("Initial 200 trials", "Final 200 trials"), name=NULL,values=c("#0066CC", '#990000')) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = c(5, 96), 
                     labels = c('Very Purple', 'Very Blue')) + 
  theme_bw() + 
  ylim(c(0.00, 1.00)) + 
  facet_wrap(~a_split+condition, nrow = 2, ncol = 2) + 
  theme(axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), 
        legend.position = 'none', 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, colour = 'black'), 
        axis.title.y = element_text(size = 12), 
        strip.placement = "outside", 
        strip.text = element_text(size = 12))

a_picc_plot_YA <- 
  d %>%
  filter(age_group == 'Young') %>% 
  mutate(a_split = factor(ifelse(a_split==1, 'Low a', 'High a'), levels = c('Low a', 'High a')), 
         timebin4 = factor(timebin4)) %>% 
  group_by(a_split, condition, trialintensity, timebin4) %>% 
  summarise(pBlue = length(response[response == 1])/length(response)) %>%
  rename(colour = trialintensity, trialbins = timebin4) %>% 
  filter(trialbins == 1 | trialbins == 4) %>% 
  ggplot(aes(x = colour, y = pBlue, color = trialbins)) + geom_point() + 
  geom_line(stat="smooth",method = "glm", method.args = list(
    family="binomial"), se=FALSE,size=1.8,alpha=.7) + 
  ylab('% Dots Judged as Blue') + xlab('') +
  ggtitle("Young Adults Only") + 
  scale_colour_manual(labels=c("Initial 200 trials", "Final 200 trials"), name=NULL,values=c("#0066CC", '#990000')) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = c(5, 96), 
                     labels = c('Very Purple', 'Very Blue')) + 
  theme_bw() + 
  ylim(c(0.00, 1.00)) + 
  facet_wrap(~a_split+condition, nrow = 2, ncol = 2) + 
  theme(axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), 
        legend.position = 'none', 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, colour = 'black'), 
        axis.title.y = element_text(size = 12), 
        strip.placement = "outside", 
        strip.text = element_text(size = 12))


a_picc_plot_OA <- 
  d %>%
  filter(age_group == 'Old') %>% 
  mutate(a_split = factor(ifelse(a_split==1, 'Low a', 'High a'), levels = c('Low a', 'High a')), 
         timebin4 = factor(timebin4)) %>% 
  group_by(a_split, condition, trialintensity, timebin4) %>% 
  summarise(pBlue = length(response[response == 1])/length(response)) %>%
  rename(colour = trialintensity, trialbins = timebin4) %>% 
  filter(trialbins == 1 | trialbins == 4) %>% 
  ggplot(aes(x = colour, y = pBlue, color = trialbins)) + geom_point() + 
  geom_line(stat="smooth",method = "glm", method.args = list(
    family="binomial"), se=FALSE,size=1.8,alpha=.7) + 
  ylab('% Dots Judged as Blue') + xlab('') +
  ggtitle("Old Adults Only") + 
  scale_colour_manual(labels=c("Initial 200 trials", "Final 200 trials"), name=NULL,values=c("#0066CC", '#990000')) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = c(5, 96), 
                     labels = c('Very Purple', 'Very Blue')) + 
  theme_bw() + 
  ylim(c(0.00, 1.00)) + 
  facet_wrap(~a_split+condition, nrow = 2, ncol = 2) + 
  theme(axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), 
        legend.position = 'none', 
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, colour = 'black'), 
        axis.title.y = element_text(size = 12), 
        strip.placement = "outside", 
        strip.text = element_text(size = 12))

ggarrange(a_picc_plot_all, 
          ggarrange(a_picc_plot_YA, a_picc_plot_OA, labels = c("B", "C")), 
          nrow = 2, labels = 'A')