
library(lme4) # for MLM

# Dots --------------------------------------------------------------------

files   = list.files('data/Dots/', pattern = ".csv")
tables  = lapply(paste0('data/Dots/', files), read.csv, header = T, stringsAsFactors = F)
d.dots  = do.call(rbind, tables)
bw.dots =  c(68, 109) # subjects who used the wrong key mapping

# Clean 

d.dots = d.dots[d.dots$block!='Practice', ]
d.dots$trialintensity = d.dots$trialintensity-154
d.dots$response  = ifelse(d.dots$id %in% bw.dots, 
                         ifelse(d.dots$response=='a', 0, 1),
                         ifelse(d.dots$response=='a', 1, 0))
d.dots$age_group = factor(ifelse(d.dots$age < 60, 'YA', 'OA'), levels = c('YA', 'OA'))
d.dots$condition = factor(ifelse(d.dots$condition == 0, 'Stable', 'Decreasing'), levels = c('Stable', 'Decreasing'))
d.dots$colour0   = d.dots$trialintensity/max(d.dots$trialintensity)
d.dots$trial0    = d.dots$trial/max(d.dots$trial)
d.dots$timebin   = cut(d.dots$trial, 4, labels = F)
d.dots$stimbin   = cut(d.dots$trialintensity, 20, labels=F)

# Visualize
# Main choice plot
pdf('plots/DotsChoice.pdf', width=10, height=6)

pBlue = tapply(d.dots$response, list(d.dots$stimbin, d.dots$timebin, d.dots$condition, d.dots$age_group), mean)

layout(matrix(1:4, 2,2,byrow=T))
par(mar=c(5.1 ,6.1, 4.1, 2.1))

# Young
plot(pBlue[,1,'Stable','YA'], type='b', col='red',  xaxt='n', ylab = 'p(Choose Blue)', 
     ylim=c(0,1), xlab='', main = 'Stable Prevalence')
lines(pBlue[,4,'Stable','YA'], type='b', col='blue')
axis(1, at=range(d.dots$stimbin), labels=c('Very Purple', 'Very Blue'))
legend('bottomright', bty='n', pch=1, lty=1,col=c('red','blue'), legend=c('First 200 Trials', 'Last 200 Trials'))
mtext('Young', side=2, line=5, font=2)

plot(pBlue[,1,'Decreasing','YA'], type='b', col='red',  xaxt='n',
     ylab = 'p(Choose Blue)', ylim=c(0,1), xlab='',  main = 'Decreasing Prevalence')
lines(pBlue[,4,'Decreasing','YA'], type='b', col='blue')
axis(1, at=range(d.dots$stimbin), labels=c('Very Purple', 'Very Blue'))

# Old
plot(pBlue[,1,'Stable','OA'], type='b', col='red',  xaxt='n', ylab = 'p(Choose Blue)', ylim=c(0,1), xlab='')
lines(pBlue[,4,'Stable','OA'], type='b', col='blue')
axis(1, at=range(d.dots$stimbin), labels=c('Very Purple', 'Very Blue'))
mtext('Old', side=2, line=5,font=2)

plot(pBlue[,1,'Decreasing','OA'], type='b', col='red',  xaxt='n', ylab = 'p(Choose Blue)', ylim=c(0,1), xlab='')
lines(pBlue[,4,'Decreasing','OA'], type='b', col='blue')
axis(1, at=range(d.dots$stimbin), labels=c('Very Purple', 'Very Blue'))

dev.off()

# Most ambiguous dot
pdf('plots/DotsAmbig.pdf', width=10, height=6)
par(mar=c(5.1 ,6.1, 4.1, 2.1))

tmp      = d.dots[d.dots$stimbin==median(d.dots$stimbin),]
diffs   = tapply(tmp$response, list(tmp$id, tmp$timebin, tmp$condition, tmp$age_group), mean)

diffYS  = diffs[,4,'Stable','YA'] - diffs[,1,'Stable','YA']
diffYD  = diffs[,4,'Decreasing','YA'] - diffs[,1,'Decreasing','YA']
diffOS  = diffs[,4,'Stable','OA'] - diffs[,1,'Stable','OA']
diffOD  = diffs[,4,'Decreasing','OA'] - diffs[,1,'Decreasing','OA']

mdiffs  = matrix(c(
  mean(diffYS, na.rm=T), 
  mean(diffOS, na.rm=T), 
  mean(diffYD, na.rm=T), 
  mean(diffOD, na.rm=T)), 
  2, 2, byrow=T,dimnames = list(c('Stable','Decreasing'), c('Young', 'Old')))*100

sediffs  = matrix(c(
  plotrix::std.error(diffYS, na.rm=T), 
  plotrix::std.error(diffOS, na.rm=T), 
  plotrix::std.error(diffYD, na.rm=T), 
  plotrix::std.error(diffOD, na.rm=T)), 
  2, 2, byrow=T,dimnames = list(c('Stable','Decreasing'), c('Young', 'Old')))*100

thisbar = barplot(mdiffs, beside=T, col=c('orange','brown'), ylim=range(pretty(c(mdiffs-sediffs,mdiffs+sediffs))),
                  ylab = '% Change in Colour Judgements\n(Last 200 Trials-First 200 Trials)',
                  main = 'Most Ambiguous Dots',
                  legend.text = T, args.legend = list(x='topright',bty='n', title='Prevalence Condition'))
arrows(thisbar, mdiffs-sediffs, thisbar, mdiffs+sediffs, length=0)

dev.off()

# MLM
d.dots$age_groupc = ifelse(d.dots$age_group =='YA', -1, 1)
d.dots$conditionc = ifelse(d.dots$condition =='Stable', -1, 1)
d.dots$colour0c   = d.dots$colour0 - 0.5

m0.dots = glmer(response ~ 1 + (1|id) , family = binomial, control=glmerControl(optimizer="bobyqa"), data = d.dots)
m1.dots = glmer(response ~ trial0*colour0c + (trial0| id) , family = binomial, control=glmerControl(optimizer="bobyqa"), data = d.dots)
m2.dots = glmer(response ~ conditionc*trial0*colour0c + (trial0| id) , family = binomial, control=glmerControl(optimizer="bobyqa"), data = d.dots)
m3.dots =  glmer(response ~ age_groupc*conditionc*trial0*colour0c + (trial0| id) , family = binomial, control=glmerControl(optimizer="bobyqa"), data = d.dots)
modcomp = anova(m3.dots, m2.dots, m1.dots, m0.dots)
CI      = confint(m3.dots, method='Wald')

out = list(m0=m0.dots, m1=m1.dots, m2=m2.dots,m3=m3.dots,modcomp=modcomp,CI=CI)
saveRDS(out, 'output/dotsmodels.rds')
capture.output(list(model=summary(m3.dots), CI=CI), file='output/dotsMLM.txt')

# Ethics ------------------------------------------------------------------

files     = list.files('data/Ethics/', pattern = ".csv")
tables    = lapply(paste0('data/Ethics/', files), read.csv, header = T, stringsAsFactors = F)
d.ethics  = do.call(rbind, tables)
bw.ethics = c(2, 4, 5, 6, 12, 17, 29, 32, 39, 50, 110, 116, 132, 140, 170, 177, 73)

# Clean
d.ethics = d.ethics[d.ethics$trial!='practice' & !is.na(d.ethics$id),]
d.ethics$response = ifelse(d.ethics$id %in% bw.ethics, 
                           ifelse(d.ethics$response == 'a', 0, 1), 
                           ifelse(d.ethics$response == 'a', 1, 0))
d.ethics$age_group  = factor(ifelse(d.ethics$age < 60, 'YA', 'OA'), levels = c('YA', 'OA'))
d.ethics$condition  = factor(ifelse(d.ethics$condition == 0, 'Stable', 'Decreasing'), levels = c('Stable', 'Decreasing'))
d.ethics$norm_mean0 = d.ethics$norm_mean/max(d.ethics$norm_mean)
d.ethics$norm_mean  = as.numeric(-d.ethics$norm_mean) # for plotting
d.ethics$trial      = as.numeric(d.ethics$trial)
d.ethics$trial0     = d.ethics$trial/max(d.ethics$trial)
d.ethics$timebin    = cut(d.ethics$trial, 5, labels=F)
d.ethics$stimbin    = cut(d.ethics$norm_mean, 7, labels=F)

# Visualize
pdf('plots/EthicsChoice.pdf', width=10, height=6)

pEth = tapply(d.ethics$response, list(d.ethics$stimbin, d.ethics$timebin, d.ethics$condition, d.ethics$age_group), mean)

layout(matrix(1:4, 2,2,byrow=T))
par(mar=c(5.1 ,6.1, 4.1, 2.1))

# Young
plot(pEth[,1,'Stable','YA'], type='b', col='red',  xaxt='n', ylab = 'p(Choose Blue)', 
     ylim=c(0,1), xlab='', main = 'Stable Prevalence')
lines(pEth[,5,'Stable','YA'], type='b', col='blue')
axis(1, at=range(d.ethics$stimbin), labels=c('Very Purple', 'Very Blue'))
legend('bottomright', bty='n', pch=1, lty=1,col=c('red','blue'), legend=c('First 48 Trials', 'Last 48 Trials'))
mtext('Young', side=2, line=5, font=2)

plot(pEth[,1,'Decreasing','YA'], type='b', col='red',  xaxt='n',
     ylab = 'p(Choose Blue)', ylim=c(0,1), xlab='',  main = 'Decreasing Prevalence')
lines(pEth[,5,'Decreasing','YA'], type='b', col='blue')
axis(1, at=range(d.ethics$stimbin), labels=c('Very Purple', 'Very Blue'))

# Old
plot(pEth[,1,'Stable','OA'], type='b', col='red',  xaxt='n', ylab = 'p(Choose Blue)', ylim=c(0,1), xlab='')
lines(pEth[,5,'Stable','OA'], type='b', col='blue')
axis(1, at=range(d.ethics$stimbin), labels=c('Very Purple', 'Very Blue'))
mtext('Old', side=2, line=5,font=2)

plot(pEth[,1,'Decreasing','OA'], type='b', col='red',  xaxt='n', ylab = 'p(Choose Blue)', ylim=c(0,1), xlab='')
lines(pEth[,5,'Decreasing','OA'], type='b', col='blue')
axis(1, at=range(d.ethics$stimbin), labels=c('Very Purple', 'Very Blue'))

dev.off()

# Most ambiguous dot
pdf('plots/EthicsAmbig.pdf', width=10, height=6)
par(mar=c(5.1 ,6.1, 4.1, 2.1))

tmp     = d.ethics[d.ethics$stimbin==median(unique(d.ethics$stimbin)),]
diffs   = tapply(tmp$response, list(tmp$id, tmp$timebin, tmp$condition, tmp$age_group), mean)


diffYS  = diffs[,5,'Stable','YA'] - diffs[,1,'Stable','YA']
diffYD  = diffs[,5,'Decreasing','YA'] - diffs[,1,'Decreasing','YA']
diffOS  = diffs[,5,'Stable','OA'] - diffs[,1,'Stable','OA']
diffOD  = diffs[,5,'Decreasing','OA'] - diffs[,1,'Decreasing','OA']

mdiffs  = matrix(c(
  mean(diffYS, na.rm=T), 
  mean(diffOS, na.rm=T), 
  mean(diffYD, na.rm=T), 
  mean(diffOD, na.rm=T)), 
  2, 2, byrow=T,dimnames = list(c('Stable','Decreasing'), c('Young', 'Old')))*100


sediffs  = matrix(c(
  plotrix::std.error(diffYS, na.rm=T), 
  plotrix::std.error(diffOS, na.rm=T), 
  plotrix::std.error(diffYD, na.rm=T), 
  plotrix::std.error(diffOD, na.rm=T)), 
  2, 2, byrow=T,dimnames = list(c('Stable','Decreasing'), c('Young', 'Old')))*10
sediffs[is.na(sediffs)]=0

thisbar = barplot(mdiffs, beside=T, xpd=F,col=c('orange','brown'),ylim=range(pretty(c(mdiffs-sediffs,mdiffs+sediffs+2))),
                  ylab = '% Change in Ethical Judgements\n(Last 48 Trials-First 48 Trials)',
                  main = 'Most Ambiguous Ethics Proposals',
                  legend.text = T, args.legend = list(x='topright',bty='n', title='Prevalence Condition'))
arrows(thisbar, mdiffs-sediffs, thisbar, mdiffs+sediffs, length=0)

dev.off()

# MLM
d.ethics$age_groupc   = ifelse(d.ethics$age_group =='YA', -1, 1)
d.ethics$conditionc   = ifelse(d.ethics$condition =='Stable', -1, 1)
d.ethics$norm_mean0c  = d.ethics$norm_mean0 - 0.5

m0.ethics = glmer(response ~ 1 + (1|id) , family = binomial, control=glmerControl(optimizer="bobyqa"), data = d.ethics)
m1.ethics = glmer(response ~ trial0*norm_mean0c + (trial0| id) , family = binomial, control=glmerControl(optimizer="bobyqa"), data = d.ethics)
m2.ethics = glmer(response ~ conditionc*trial0*norm_mean0c + (trial0| id) , family = binomial, control=glmerControl(optimizer="bobyqa"), data = d.ethics)
m3.ethics =  glmer(response ~ age_groupc*conditionc*trial0*norm_mean0c + (trial0| id) , family = binomial, control=glmerControl(optimizer="bobyqa"), data = d.ethics)
modcomp   = anova(m3.ethics, m2.ethics, m1.ethics, m0.ethics)
CI        = confint(m3.ethics, method='Wald')

out = list(m0=m0.ethics, m1=m1.ethics, m2=m2.ethics,m3=m3.ethics,modcomp=modcomp,CI=CI)
saveRDS(out, 'output/ethicsmodels.rds')
capture.output(list(model=summary(m3.ethics), CI=CI), file='output/ethicsMLM.txt')



# RT ----------------------------------------------------------------------

# Visualize

pdf('plots/meanRTs.pdf', width=10, height=6)
layout(matrix(1:2,ncol=2))

# Dots
mRT  = tapply(d.dots$RT, list(d.dots$condition, d.dots$age_group), mean)
colnames(mRT) = c('Young', 'Old')
seRT = tapply(d.dots$RT, list(d.dots$condition, d.dots$age_group), plotrix::std.error)
colnames(seRT) = c('Young', 'Old')

ylim = range(pretty(c(mRT-seRT, mRT+seRT)))
thisbar = barplot(mRT, beside=T, xpd=F,ylim=ylim, col=c('brown', 'orange'),
                  ylab = 'Avg. Response Time (s.)', main='Dots Task',
                  legend.text = T, args.legend = list(x='topleft', bty='n', title='Prevalence Condition'))
arrows(thisbar, mRT-seRT, thisbar, mRT+seRT, length=0)

# Ethics
mRT  = tapply(d.ethics$RT, list(d.ethics$condition, d.ethics$age_group), mean)
colnames(mRT) = c('Young', 'Old')
seRT = tapply(d.ethics$RT, list(d.ethics$condition, d.ethics$age_group), plotrix::std.error)
colnames(seRT) = c('Young', 'Old')

ylim = range(pretty(c(mRT-seRT, mRT+seRT)))
thisbar = barplot(mRT, beside=T, xpd=F,ylim=ylim, col=c('brown', 'orange'),
                  ylab = 'Avg. Response Time (s.)', main='Ethics Task',
                  legend.text = T, args.legend = list(x='topleft', bty='n', title='Prevalence Condition'))
arrows(thisbar, mRT-seRT, thisbar, mRT+seRT, length=0)

dev.off()

# Regression 
# dots
tmp = data.frame()
for(id in unique(d.dots$id)) {
  mrt  = mean(d.dots$RT[d.dots$id==id])
  cond = d.dots$conditionc[d.dots$id==id][1]
  ag   = d.dots$age_groupc[d.dots$id==id][1]
  tmp = rbind(tmp, data.frame(id=id,condition=cond,age_group=ag,mrt=mrt))
}

RTlm_dots = lm(mrt ~ age_group*condition, data=tmp)
summary(RTlm_dots)
CI = confint(RTlm_dots)
capture.output(list(model=summary(RTlm_dots), CI=CI), file='output/dotsRT.txt')

# ethics
tmp = data.frame()
for(id in unique(d.ethics$id)) {
  mrt  = mean(d.ethics$RT[d.ethics$id==id])
  cond = d.ethics$conditionc[d.ethics$id==id][1]
  ag   = d.ethics$age_groupc[d.ethics$id==id][1]
  tmp = rbind(tmp, data.frame(id=id,condition=cond,age_group=ag,mrt=mrt))
}

RTlm_ethics = lm(mrt ~ age_group*condition, data=tmp)
summary(RTlm_ethics)
CI = confint(RTlm_ethics)
capture.output(list(model=summary(RTlm_ethics), CI=CI), file='output/ethicsRT.txt')
