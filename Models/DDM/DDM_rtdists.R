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
if(as.list(Sys.info())$sysname == 'Linux') dir = '~/Desktop/Ongoing/PICCAgeing/Analysis&Data/' else dir = "C:/Users/Sean/Desktop/LDM Lab/PICCAgeing" 
setwd(dir)

# LOAD EXP 1 DOTS DATA ####
files = list.files(path='Exp1/Dots',pattern = ".csv")
d <- do.call(rbind,  lapply(paste0(dir, 'Exp1/Dots/', files), read.csv, header = T, stringsAsFactors = F))

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

# # summarise/visualise if you want 
# binsum <- 
#   d %>% 
#   mutate(timebin4 = factor(timebin4)) %>%
#   filter(timebin4 == 1 | timebin4 == 4) %>% 
#   group_by(id, age_group, condition, strength_bin, timebin4) %>% 
#   summarise(pblue = mean(response)) %>% 
#   ggplot(aes(x = strength_bin, y = pblue, colour=timebin4)) + 
#   stat_summary(fun.y = mean, geom='point', position=position_dodge(0.9)) + 
#   stat_summary(fun.data = mean_se, geom='errorbar', position=position_dodge(0.9)) + 
#   stat_summary(fun.y = mean, geom='line', position=position_dodge(0.9)) + 
#   facet_grid(condition ~ age_group)
#   

# DIFFUSION ANALYSIS W/ RTDISTS ####
# Following https://cran.r-project.org/web/packages/rtdists/vignettes/reanalysis_rr98.html
# Maximum likelihood estimation
d$strength_bin <- factor(d$strength_bin)
ndrift <- length(levels(d$strength_bin))

d$response_ddm <- ifelse(d$response == 1, 'upper', 'lower')

objective_diffusion_separate <- function(pars, rt, response, drift, ...) {
  non_v_pars <- grep("^v", names(pars), invert = TRUE, value = TRUE)
  base_par <- length(non_v_pars)  # number of non-drift parameters
  densities <- vector("numeric", length(rt))
  for (i in 1:ndrift) {
    densities[drift == levels(drift)[i]] <- 
      ddiffusion(rt[drift == levels(drift)[i]], 
                 response = response[drift == levels(drift)[i]], 
                 a=pars["a"], t0=pars["t0"],  
                 sv=pars["sv"],
                 sz=if ("sz" %in% non_v_pars) pars["sz"] else 0.1,
                 z=if ("z" %in% non_v_pars) pars["z"]*pars["a"] else 0.5*pars["a"],
                 st0=if ("st0" %in% non_v_pars) pars["st0"] else 0, 
                 v=pars[base_par+i])
  }
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}

# Function that creates random start values 
get_start <- function(base_par, n_drift = ndrift) {
  start1 <- c(
    a = runif(1, 0.5, 3),
    a_1 = runif(1, 0.5, 3), 
    a_2 = runif(1, 0.5, 3),
    t0 = runif(1, 0, 0.5), 
    z = runif(1, 0.4, 0.6), 
    sz = runif(1, 0, 0.5),
    sv = runif(1, 0, 0.5),
    st0 = runif(1, 0, 0.5),
    d = rnorm(1, 0, 0.05)
  )
  start2 <- sort(rnorm(n_drift), decreasing = FALSE)
  names(start2) <- paste0("v_", seq_len(n_drift))
  c(start1[base_par], start2)
}

# Function that tries different random start values until it works:
ensure_fit <-  function(data, start_function, objective_function,
                        base_pars, n_drift = 5, n_fits = 1, 
                        lower = c(rep(0, length(base_pars)), -Inf,
                                  rep(-Inf,length(start_function(base_pars))-length(base_pars)))) {
  
    best_fit <- list(objective = 1e+06)
    for (i in seq_len(n_fits)) {
      start_ll <- 1e+06
      cat("\nfinding viable start point.\n") 
      while(start_ll == 1e+06) {
        start <- start_function(base_pars, n_drift=ndrift)
        start_ll <- objective_function(start, 
                                       rt = data$RT, response = data$response_ddm, 
                                       drift = factor(data$strength_bin, seq_len(ndrift)), 
                                       timebin4 = data$timebin4)
      }
      cat("\nstart fitting.\n") 
      
      fit <- nlminb(start, objective_function, 
                    rt = data$RT, response = data$response_ddm, 
                    drift = factor(data$strength_bin, seq_len(ndrift)), 
                    timebin4 = data$timebin4,
                    lower = lower)
      
      if (fit$objective < best_fit$objective) best_fit <- fit
    }
    out <- matrix(c(fit$par, fit$objective), 1, length(fit$par)+1, byrow = T)
    colnames(out) <- c(names(fit$par), 'objective')
    out
  }

# Fit! 
d.fit <- data.frame()
for(n in 1:length(unique(d$id))) { 
  
  cat(paste0("\n-------- FITTING SUBJECT ",n,"/", length(unique(d$id))," ------------\n")) # keep track 
  
  id = unique(d$id)[n]
  
  fit <- ensure_fit(data = d[d$id==id,], 
                    start_function = get_start, 
                    objective_function = objective_diffusion_separate, 
                    base_pars = c("a", "t0", "sv", "sz", "z"))
  
  fit <- cbind(d[d$id==id, c('id', 'age', 'condition', 'age_group')][1,], fit)
  d.fit <- rbind(d.fit, data.frame(fit))
  
}

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


# ANALYSE ####
# multiple multivariate regression
pars_mmr = lm(cbind(a, t0, sv, sz, z, v_1, v_2, v_3, v_4, v_5) ~ age_group*condition, data= d.fit)
cohend_a = cohen.d(d  = d.fit$a, f= d.fit$age_group)
cohend_t0 = cohen.d(d  = d.fit$t0, f= d.fit$age_group)

# EXTRA ####
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