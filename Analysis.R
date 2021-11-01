# Setup -------------------------------------------------------------------
library(dplyr); library(reshape2) # for data manipulation
library(ggplot2);library(ggpubr)  # for plotting
library(lme4); library(emmeans)   # for MLM
library(rtdists)                  # for DDM
library(car)                      # for model comparison and misc.


if(as.list(Sys.info())$sysname == 'Linux')
  dir = '~/Desktop/Ongoing/PICCAgeing/Analysis&Data/' 
if(as.list(Sys.info())$sysname == 'Windows') 
  dir = "C://Users//LSDMlab_RA//Desktop//PICCAgeing//Analysis&Data/" 
setwd(dir)

# Dots ####
files = list.files(paste0('Exp1/Dots'), pattern = ".csv")
tables <- lapply(paste0('Exp1/Dots/', files), read.csv, header = T, stringsAsFactors = F)
d.dots <- do.call(rbind, tables)

backwardmap.dots <- c(68, 109) # subjects who used the wrong key mapping
d.dots <-
  d.dots %>% 
  filter(block != 'Practice') %>% 
  mutate(trialintensity = trialintensity-154, 
         response = ifelse(id %in% backwardmap.dots, 
                           ifelse(response=='a', 0, 1),
                           ifelse(response=='a', 1, 0)), 
         responsenum = as.numeric(response)-1, 
         id = as.factor(id), 
         age_group = factor(ifelse(age < 60, 'YA', 'OA'), levels = c('YA', 'OA')), 
         condition = factor(ifelse(condition == 0, 'Stable', 'Decreasing'), 
                            levels = c('Stable', 'Decreasing')), 
         block = factor(block, levels = c(1:16)), 
         colour0 = trialintensity/max(trialintensity), 
         trial0 = trial/max(trial),
         timebin4 = ntile(trial, 4))

# * Visualise ####
dotspiccplot <- 
  d.dots %>% 
  mutate(age_group = factor(ifelse(age_group=='YA', 'Young', 'Old'), 
                            levels=c('Young', 'Old')), 
         timebin4 = factor(timebin4)) %>% 
  group_by(age_group, condition, trialintensity, timebin4) %>% 
  summarise(pBlue = length(response[response == 1])/length(response)) %>%
  filter(timebin4 == 1 | timebin4 == 4) %>% 
  ggplot(aes(x = trialintensity, y = pBlue, color = timebin4)) + 
  geom_point() + 
  geom_line(stat="smooth",method = "glm", method.args = list(
    family="binomial"), se=FALSE,size=1.8,alpha=.7) + 
  labs(x= '', y = '% Dots Judged as Blue', color = '') +
  scale_colour_manual(labels=c("Initial 200 trials", "Final 200 trials"), 
                      name=NULL,values=c("#0066CC", '#990000')) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                     limits = c(0, 1)) +
  scale_x_continuous(breaks = c(5, 96), 
                     labels = c('Very\nPurple', 'Very\nBlue')) + 
  theme_bw() + 
  facet_wrap(~age_group+condition, nrow = 2, ncol = 2) + 
  theme(axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5), 
        legend.position = 'bottom')

dotspiccplot_zoom <- 
  d.dots %>%
  mutate(age_group = factor(ifelse(age_group=='YA', 'Young', 'Old'), 
                            levels=c('Young', 'Old')), 
         stimbin=ntile(trialintensity, 10)) %>% 
  group_by(id, age_group, condition, stimbin, timebin4) %>% 
  summarise(pBlue = mean(response), 
            n=n()) %>% 
  mutate(timebin4 = factor(timebin4)) %>% 
  filter(timebin4 == 1 | timebin4 == 4, 
         stimbin==10) %>% 
  group_by(id, age_group, condition, stimbin) %>% 
  summarise(change = pBlue[timebin4==4] - pBlue[timebin4==1]) %>% 
  ggplot(aes(x = age_group, y = change, fill=condition)) + 
  stat_summary(fun.y = mean, geom = 'bar', position = position_dodge(0.9), color='black') + 
  stat_summary(fun.data = mean_se, geom = 'errorbar', position = position_dodge(0.9), width=0) + 
  labs(x = '',
       y = '% Change in Colour Judgements\n(First 200 Trials \u2013 Last 200 Trials)', 
       fill='Condition', 
       title='Bluest Dot') + 
  scale_y_continuous(labels=scales::percent_format(accuracy = 1)) + 
  scale_fill_manual(values=c('#4A3E3D', '#C16527')) + 
  theme_bw() 

# * Analyse ####
m0.dots = glmer(response ~ 1 + (1|id) , family = binomial, control=glmerControl(optimizer="bobyqa"), data = d.dots)
m1.dots = glmer(response ~ trial0*colour0 + (trial0:colour0| id) , family = binomial, control=glmerControl(optimizer="bobyqa"), data = d.dots)
m2.dots = glmer(response ~ condition*trial0*colour0 + (trial0:colour0| id) , family = binomial, control=glmerControl(optimizer="bobyqa"), data = d.dots)
m3.dots =  glmer(response ~ age_group*condition*trial0*colour0 + (trial0:colour0| id) , family = binomial, control=glmerControl(optimizer="bobyqa"), data = d.dots)
anova(m3.dots, m2.dots, m1.dots, m0.dots)

Anova(m3.dots)
contrasts.dots = pairs(emmeans(m3.dots, 
                               ~age_group:condition:trial0:colour0, 
                               at=list(trial0 = 0.875, colour0=0.5),
                               type='response'))
contrasts.dots.ci = confint(contrasts.dots)

# Ethics ------------------------------------------------------------------
files = list.files('Exp1/Ethics/', pattern = ".csv")
tables <- lapply(paste0('Exp1/Ethics/', files), read.csv, header = T, stringsAsFactors = F)
d.ethics <- do.call(rbind, tables)

backwardmap.ethics <- c(2, 4, 5, 6, 12, 17, 29, 32, 39, 50, 110, 116, 132, 
                        140, 170, 177, 73)
d.ethics <- 
  d.ethics %>% 
  filter(trial != 'practice', 
         !is.na(id)) %>% 
  mutate(id = factor(id), 
         block = factor(block, levels=1:10), 
         response = ifelse(id %in% backwardmap.ethics, 
                           ifelse(response == 'a', 0, 1), 
                           ifelse(response == 'a', 1, 0)),
         responsenum = as.numeric(response)-1, 
         age_group = factor(ifelse(age < 60, 'YA', 'OA'), levels = c('YA', 'OA')), 
         condition = factor(ifelse(condition == 0, 'Stable', 'Decreasing'), 
                            levels = c('Stable', 'Decreasing')), 
         norm_mean0 = norm_mean/max(norm_mean), 
         norm_mean = as.numeric(-norm_mean), 
         trial = as.numeric(trial), 
         timebin = ntile(trial, 5), 
         normbin = ntile(norm_mean, 24),
         trial0 = trial/max(trial))

# * Visualise ####
ethicspiccplot <-
  d.ethics %>%
  mutate(age_group = factor(ifelse(age_group=='YA', 'Young', 'Old'), 
                            levels=c('Young', 'Old'))) %>% 
  group_by(age_group, condition, normbin, timebin) %>% 
  summarise(pBlue = length(response[response == 1])/length(response)) %>% 
  mutate(timebin = factor(timebin)) %>% 
  filter(timebin == 1 | timebin == 5) %>% 
  ggplot(aes(x = normbin, y = pBlue, color = timebin)) + 
  geom_point() + 
  geom_line(stat="smooth",method = "glm", method.args = list(
    family="binomial"), se=FALSE,size=1.8,alpha=.7) + 
  labs(y = '% Unethical Judgements', x='', color='') +
  scale_colour_manual(labels=c("Initial 48 trials", "Final 48 trials"), name=NULL,values=c("#0066CC", '#990000')) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                     limits = c(0, 1)) +
  scale_x_continuous(breaks = c(2, 22), 
                     labels = c('Very\nEthical', 'Very\nUnethical')) + 
  theme_bw() + 
  facet_wrap(~age_group+condition, nrow = 2, ncol = 2) + 
  theme(axis.ticks.x = element_blank(), 
        plot.title = element_text(hjust = 0.5), 
        legend.position = 'bottom')

ethicspiccplot_zoom <- 
  d.ethics %>%
  mutate(age_group = factor(ifelse(age_group=='YA', 'Young', 'Old'), 
                            levels=c('Young', 'Old'))) %>% 
  group_by(id, age_group, condition, normbin, timebin) %>% 
  summarise(pNorm = mean(response), 
            n=n()) %>% 
  mutate(timebin = factor(timebin)) %>% 
  filter(timebin == 1 | timebin == 5, 
         normbin==24) %>% 
  group_by(id, age_group, condition, normbin) %>% 
  summarise(change = pNorm[timebin==5] - pNorm[timebin==1]) %>% 
  ggplot(aes(x = age_group, y = change, fill=condition)) + 
  stat_summary(fun.y = mean, geom = 'bar', position = position_dodge(0.9), color='black') + 
  stat_summary(fun.data = mean_se, geom = 'errorbar', position = position_dodge(0.9), width=0) + 
  labs(x = '',
       y = '% Change in Unethical Judgements\n(First 48 Trials \u2013 Last 48 Trials)', 
       fill='Condition', 
       title='Most Unethical Proposal') + 
  scale_y_continuous(labels=scales::percent_format(accuracy = 1)) + 
  scale_fill_manual(values=c('#4A3E3D', '#C16527')) + 
  theme_bw()

# fig. 1
ggarrange(dotspiccplot, dotspiccplot_zoom,
          ethicspiccplot, ethicspiccplot_zoom, 
          labels = letters[1:4], legend='bottom')

# * Analyse #### 
m0.ethics = glmer(response ~ 1 + (1 | id), family = binomial, control=glmerControl(optimizer="bobyqa"), data = d.ethics)
m1.ethics = glmer(response ~ trial0:norm_mean0 + (trial0:norm_mean0 | id), family = binomial, control=glmerControl(optimizer="bobyqa"), data = d.ethics)
m2.ethics = glmer(response ~ condition*trial0*norm_mean0 + (trial0:norm_mean0 | id), family = binomial, control=glmerControl(optimizer="bobyqa"), data = d.ethics)
m3.ethics = glmer(response ~ age_group*condition*trial0*norm_mean0 + (trial0:norm_mean0 | id), family = binomial, control=glmerControl(optimizer="bobyqa"), data = d.ethics)

anova(m3.ethics, m2.ethics, m1.ethics, m0.ethics)
Anova(m3.ethics)

contrasts.ethics = pairs(emmeans(m3.ethics, 
                                 ~age_group:condition:trial0:norm_mean0,
                                 at = list(trial0 = 0.875, norm_mean0 = 0.2904039),
                                 type='response'))
contrasts.ethics.ci = confint(contrasts.ethics)

# Sequential Choice Model (Wilson, 2018) ####
source('Models/seq/fitpicc.R')
# * Fit ####
dots.fit = fitpicc(d.dots, 'dots', iterList = backwardmap.dots)
ethics.fit = fitpicc(d.ethics, 'ethics', iterList = backwardmap.ethics)

# skip fitting
# dots.fit <- read.csv('Models/seq/dots.fit.csv')
# ethics.fit <- read.csv('Models/seq/ethics.fit.csv')

# * Visualise ####
pars = c('B0', 'Bf', 'BF', 'Bc', 'lF', 'lc')
ranges = list(c(-5, 5), c(0, 20), c(-1, 1), c(-1, 1), c(-1, 1), c(-1, 1))
plots = list()
for(p in 1:length(pars)) { 
  
  if(pars[p]=='B0') {
    ylab = 'Parameter\nValue' 
    plotlab = 'a'
  } else {
    ylab = ''
    plotlab = ''
  }
  thisdotplot <- 
    dots.fit %>%
    select(id, age_group, condition, pars[p]) %>%
    melt(id.vars=c('id', 'age_group', 'condition')) %>%
    mutate(age_group = factor(ifelse(age_group == 'OA', 'Old', 'Young'), levels = c('Young', 'Old')),  
           condition = factor(condition, levels = c('Stable', 'Decreasing')), 
           variable = factor(variable, 
                             levels = c('B0', 'Bf', 'BF', 'Bc', 'lF', 'lc'), 
                             ordered = T, 
                             labels = c(expression('\u03B2'[0]), 
                                        expression('\u03B2'[f]), 
                                        expression('\u03B2'[F]), 
                                        expression('\u03B2'[c]), 
                                        expression('\u03BB'[F]), 
                                        expression('\u03BB'[c])))) %>%
    ggplot(aes(x = age_group, y = value, colour = condition)) +
    stat_summary(fun.y = mean, geom='point') +
    stat_summary(fun.data = mean_se, geom='errorbar', width=.2) +
    geom_hline(yintercept = 0, linetype = 'dotted', size= 1) + 
    facet_wrap(~variable, labeller = label_parsed) +
    scale_y_continuous(limits = ranges[[p]], 
                       breaks = seq(ranges[[p]][1], ranges[[p]][2], length.out = 3)) + 
    labs(x = '', y = ylab, colour = 'Condition', tag = plotlab) +
    scale_colour_manual(values=c('#4A3E3D', '#C16527')) + 
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5), 
          strip.text = element_text(size=14)) 
  
  n = length(plots)+1
  plots[[n]] <- thisdotplot
  
}

for(p in 1:length(pars)) { 
  
  if(pars[p]=='B0') {
    ylab = 'Parameter\nValue' 
    plotlab = 'b'
  } else {
    ylab = ''
    plotlab = ''
  }
  thisethicsplot <- 
    ethics.fit %>%
    select(id, age_group, condition, pars[p]) %>%
    melt(id.vars=c('id', 'age_group', 'condition')) %>%
    mutate(age_group = factor(ifelse(age_group == 'OA', 'Old', 'Young'), levels = c('Young', 'Old')),  
           condition = factor(condition, levels = c('Stable', 'Decreasing')), 
           variable = factor(variable, 
                             levels = c('B0', 'Bf', 'BF', 'Bc', 'lF', 'lc'), 
                             ordered = T, 
                             labels = c(expression('\u03B2'[0]), 
                                        expression('\u03B2'[f]), 
                                        expression('\u03B2'[F]), 
                                        expression('\u03B2'[c]), 
                                        expression('\u03BB'[F]), 
                                        expression('\u03BB'[c])))) %>%
    ggplot(aes(x = age_group, y = value, colour = condition)) +
    stat_summary(fun.y = mean, geom='point') +
    stat_summary(fun.data = mean_se, geom='errorbar', width=.2) +
    geom_hline(yintercept = 0, linetype = 'dotted', size= 1) + 
    facet_wrap(~variable, labeller = label_parsed) +
    scale_y_continuous(limits = ranges[[p]], 
                       breaks = seq(ranges[[p]][1], ranges[[p]][2], length.out = 3)) + 
    labs(x = '', y = ylab, colour = 'Condition', tag=plotlab) +
    scale_colour_manual(values=c('#4A3E3D', '#C16527')) + 
    theme_bw() +
    theme(plot.title = element_text(hjust=0.5), 
          strip.text = element_text(size=14)) 
  
  n = length(plots)+1
  plots[[n]] <- thisethicsplot
}

ggarrange(plotlist = plots, common.legend = T, legend = 'bottom', nrow = 2, ncol = 6)

# * Analyse ####
# multivariate multiple regression
dots.seq.mmr = lm(cbind(B0, Bf, BF, Bc, lF, lc) ~ age_group*condition, data=dots.fit)
ethics.seq.mmr = lm(cbind(B0, Bf, BF, Bc, lF, lc) ~ age_group*condition, data=ethics.fit)

# Age Related RT Differences ####
RTsum <- rbind(
  d.dots %>% group_by(id, age_group, condition) %>% summarise(meanRT = mean(RT), task = 'Dots'), 
  d.ethics %>% group_by(id, age_group, condition) %>% summarise(meanRT = mean(RT), task = 'Ethics')
)
RTsum$age_group <- factor(ifelse(RTsum$age_group == 'YA', 'Young', 'Old'), 
                          levels = c('Young', 'Old'))

# * Visualise ####
RTexp1plot <- 
  RTsum %>% 
  ggplot(aes(x = age_group, y = meanRT, fill = condition)) + 
  stat_summary(fun.y = mean, geom='bar', position = position_dodge(0.9), colour='black') +
  stat_summary(fun.data = mean_se, geom='errorbar', position = position_dodge(0.9), width = 0.2) +
  geom_point(alpha = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values=c('#4A3E3D', '#C16527')) + 
  facet_wrap(~task, scales = 'free') +
  labs(x = '', y = 'Mean RT (in ms.)', fill= 'Condition') +
  theme_bw() 

# * Analyse ####
mRT.dots <- lm(meanRT ~ age_group*condition, data=RTsum[RTsum$task=='Dots',]) 
mRT.dots.confint <- confint(mRT.dots)
mRT.ethics <- lm(meanRT ~ age_group*condition, data=RTsum[RTsum$task=='Ethics',]) 
mRT.ethics.confint <- confint(mRT.ethics)

# Drift Diffusion Model of Dots Data ---------------------------------------------------
library(rtdists) # needed to compute drift diffusion likelihood
source('Models/DDM/fitddm.R')

dots.ddm <- fitddm(d.dots, nFit = 50, savedir = 'Models/DDM/out/exp1/')

# skip fitting.
# dots.ddm <- read.csv('Models/DDM/out/ddmfit1.csv')

# * Visualise ####
base_pars_plot <-
  dots.ddm %>%
  select(!contains('v_'), -loglik, -convergence) %>%
  melt(id.vars = c('id', 'age_group', 'condition')) %>% 
  mutate(age_group = factor(ifelse(age_group == 'OA', 'Old', 'Young'),
                            levels = c('Young', 'Old')), 
         condition = factor(condition, levels = c('Stable', 'Decreasing'))) %>% 
  ggplot(aes(y = value, x = age_group, fill = condition)) + 
  stat_summary(fun.y = mean, geom='bar', position = position_dodge(0.9)) + 
  stat_summary(fun.data = mean_se, geom='errorbar', position = position_dodge(0.9)) + 
  geom_point(position = position_dodge(0.9), alpha=0.2) + 
  facet_wrap(~variable, scales = 'free_y') + 
  scale_fill_manual(values=c('#4A3E3D', '#C16527')) + 
  labs(x = '', y='Parameter\nValue', fill = 'Condition') + 
  theme_bw()

drift_rate_plot <- 
  dots.ddm %>%
  select(id, age_group, condition, contains('v_')) %>% 
  melt(id.vars = c('id', 'age_group', 'condition')) %>% 
  mutate(age_group = factor(ifelse(age_group == 'OA', 'Old', 'Young'), 
                            levels = c('Young', 'Old'))) %>% 
  mutate(variable = as.numeric(gsub('v_', "", variable, fixed = T))) %>% 
  ggplot(aes(x = variable, y = value, color = condition)) + 
  stat_summary(fun.y = mean, geom='point', position = position_dodge(0.9)) + 
  stat_summary(fun.data = mean_se, geom='errorbar', position = position_dodge(0.9)) + 
  stat_summary(fun.y = mean, geom='line', position = position_dodge(0.9)) + 
  geom_point(position = position_dodge(0.9), alpha=0.2) + 
  facet_wrap(~age_group) + 
  scale_colour_manual(values=c('#4A3E3D', '#C16527')) + 
  scale_x_continuous(breaks = c(1, 5), 
                     labels = c('Very\nPurple', 'Very\nBlue')) + 
  labs(x = '', y='Drift Rate', color = 'Condition') + 
  theme_bw()

ggarrange(base_pars_plot, drift_rate_plot, labels = letters[1:2],
          common.legend = T, legend = "bottom")

# * Analyse ####
dots.ddm.mmr = lm(cbind(a, t0, sv, sz, z, v_1, v_2, v_3, v_4, v_5) ~ age_group*condition, data=dots.ddm)

# * Parameters vs. Random Effects ####
raneff.ddm =
  raneffs %>% 
  mutate(id = as.numeric(as.character(id))) %>% 
  full_join(., dots.ddm) %>% 
  filter(loglik < 1000) # these are the same subjects with really high mean RT 

# * Visualise ####
raneff_a_plot <- 
  raneff.ddm %>% 
  ggplot(aes(x = a, y = `trial0:colour0`, colour=condition)) + 
  geom_point(alpha=0.2) + 
  geom_smooth(method = 'lm', se=F) + 
  labs(x = expression(paste('Boundary Seperation (',italic(a), ')')), y = expression(U[0][id]~'\u002B \u03B2'[trial0:colour0]),
       colour='Condition') + 
  scale_color_brewer(palette = 'Dark2') + 
  theme_bw()

# * Analyse ####
# Simple multiple regression 
raneffddm.lm = lm(`trial0:colour0` ~ a*condition, data = raneff.ddm)