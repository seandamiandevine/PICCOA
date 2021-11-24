
source('models.R')
library(optimParallel)


# Dots --------------------------------------------------------------------


d = read.csv('clean_dots.csv')
models = data.frame()
Niter = 1
alreadyfit = list.files('fits/')
alreadyfit = sapply(alreadyfit, function(x) strsplit(x,split = '\\.')[[1]][1])

for(id in unique(d$id)){
  if(paste0('dots_',id) %in% alreadyfit) next
  cat('----------', match(id,unique(d$id)), '/', length(unique(d$id)), '----------\n')
  sub     = d[d$id==id, ]
  sub     = sub[sub$RT<30,] # remove outlandish RTs
  choices = sub$response
  RT      = sub$RT
  stim    = sub$trialintensity
  
  # # Visualize real subject curve
  # stimbin = dplyr::ntile(sub$trialintensity, 20)
  # mResp = tapply(choices, list(stimbin, sub$timebin4), mean)
  # plot(mResp[,1], type='b', col='red', xaxt='n', ylab='p(Blue)', xlab='', main=paste0('Empirical\nid=',id), ylim=c(0,1))
  # lines(mResp[,4], type='b', col='blue')
  # axis(1, at=c(1,20), labels=c('Very Purple', 'Very Blue'))
  # legend('topleft', bty='n', lty=1, pch=1, col=c('red', 'blue'), legend=c('First 200 Trials', 'Last 200 Trials'))
  
  bestLL  = 1e6
  bestFit = list()
  for(iter in 1:Niter){
    cat('*** iteration', iter, '/', Niter, '***\n')
    
    # specify cost fx 
    obfunc = function(params) {
      LL = RTFDDM(stim, RT, choices, params)
      if(any(abs(LL)==Inf)) return(1e6)
      -sum(LL)
    }
    # specify x0 that works
    lim = 1000
    for(x in 1:lim) {
      a0      = rgamma(1, 2, .5)
      ter0    = max(1e-6, rnorm(1, 0.2, .5))
      nk0     = sample(5:100, size=1)
      x0      = c(a0, ter0, nk0)
      if(obfunc(x0)<1e6) break
      if(x>lim/2) cat('finding x0 taking a long time....')
      if(x==lim) cat('failed to find suitable x0')
    }
    
    # fit
    cat('optimizing...\n')
    cl = makeCluster(detectCores(), type='FORK')
    opt = optimParallel(par=x0, fn=obfunc, lower = c(0.1,.05,5), upper = c(10,10,800), 
                        # control = list(trace=6), 
                        parallel = list(cl=cl, loginfo=T))
    closeAllConnections()
    
    if(opt$value == 1e6) opt = optim(par=x0, fn=obfunc)
    
    if(opt$value < bestLL) {
      bestLL = opt$value
      bestFit = opt
    }
  }
  
  # Grid-search for nk
  est_a   = bestFit$par[1]
  est_ter = bestFit$par[2] 
  probs = c()
  grid = 5:100
  for(k in grid){
    cat('nk=',k,'||')
    probs= c(probs, obfunc(c(est_a, est_ter, k)))
  }
  cat('\n')
  probs[is.na(probs)] = 1e6
  plot(grid, probs, xlab='nk', ylab='LL', type='b', main=id) # check for convexity
  est_nk = grid[match(min(probs), probs)]
  legend('topright', bty='n', paste0('est nk=', est_nk))
  
  LL = min(probs)
  # Compute BIC
  bic = function(LL, k, n) k*log(n) - 2*LL
  aic = function(LL,k) -2*LL + 2*k
  
  npars = length(opt$par)
  RFTDDMBIC = bic(-LL, npars, nrow(sub))
  RFTDDMAIC = aic(-LL, npars)
  
  # Save
  subfit = data.frame(
    id = rep(id, npars), 
    condition = rep(sub$condition[1], npars), 
    model = rep('RFTDDM', npars), 
    parsname = c('a', 't0', 'nk'), 
    pars = c(est_a, est_ter, est_nk), 
    LL = rep(LL, npars), 
    BIC = rep(RFTDDMBIC, npars), 
    AIC = rep(RFTDDMBIC, npars)
  )
  saveRDS(list(opt=bestFit, nk=est_nk, df=subfit), paste0('fits/dots_', id,'.rds'))
  models = rbind(models, subfit)
}
write.csv(models, 'dots_fit.csv')


# Analyze parameters 

models$age_group = factor(ifelse(models$id<100, 'Old', 'Young'), levels=c('Young', 'Old'))
models$condition = factor(models$condition)
contrasts(models$age_group) = contr.sum(2)
contrasts(models$condition) = contr.sum(2)

# a 
pdf('plots/a_ddm.pdf', width=6, height=6)

apars = models[models$parsname=='a',]

ma  = tapply(apars$pars, list(apars$condition, apars$age_group), mean)
sea = tapply(apars$pars, list(apars$condition, apars$age_group), plotrix::std.error)

ylim = range(pretty(c(ma-sea,ma+sea)))
plot(ma['Stable',], col='brown', ylim=ylim, xaxt='n', pch=16, cex=2, ylab='Decision Threshold', xlab='', main='Dots Task')
arrows(c(1,2), ma['Stable',]-sea['Stable',], c(1,2), ma['Stable',]+sea['Stable',], length=0, col='brown')
points(ma['Decreasing',], col='orange', pch=16, cex=2)
arrows(c(1,2), ma['Decreasing',]-sea['Decreasing',], c(1,2), ma['Decreasing',]+sea['Decreasing',], length=0, col='orange')
axis(1, at=c(1,2), labels=c('Young', 'Old'))
legend('topleft', bty='n', col=c('brown', 'orange'), lty=1, pch=16, legend=c('Stable', 'Decreasing'), title='Prevalence Condition')

dev.off()

alm = lm(pars~age_group*condition, data=apars)
summary(alm)
CI = confint(alm)
capture.output(list(summary(alm), CI=CI), file='output/alm.txt')


# ter
pdf('plots/t0_ddm.pdf', width=10, height=6)

t0pars = models[models$parsname=='t0',]

mt0  = tapply(t0pars$pars, list(t0pars$condition, t0pars$age_group), mean)
set0 = tapply(t0pars$pars, list(t0pars$condition, t0pars$age_group), plotrix::std.error)

ylim = range(pretty(c(mt0-set0,mt0+set0)))
plot(mt0['Stable',], col='brown', ylim=ylim, xaxt='n', pch=16, cex=2, ylab='Non-Decision Time', xlab='', main='Dots Task')
arrows(c(1,2), mt0['Stable',]-set0['Stable',], c(1,2), mt0['Stable',]+set0['Stable',], length=0, col='brown')
points(mt0['Decreasing',], col='orange', pch=16, cex=2)
arrows(c(1,2), mt0['Decreasing',]-set0['Decreasing',], c(1,2), mt0['Decreasing',]+set0['Decreasing',], length=0, col='orange')
axis(1, at=c(1,2), labels=c('Young', 'Old'))
legend('topleft', bty='n', col=c('brown', 'orange'), lty=1, pch=16, legend=c('Stable', 'Decreasing'))

dev.off()

t0lm = lm(pars~age_group*condition, data=t0pars)
summary(t0lm)
CI = confint(t0lm)
capture.output(list(summary(t0lm), CI=CI), file='output/t0lm.txt')

# nk
pdf('plots/t0_ddm.pdf', width=10, height=6)

nkpars = models[models$parsname=='nk',]

mnk  = tapply(nkpars$pars, list(nkpars$condition, nkpars$age_group), mean)
senk = tapply(nkpars$pars, list(nkpars$condition, nkpars$age_group), plotrix::std.error)

ylim = range(pretty(c(mnk-senk,mnk+senk)))
plot(mnk['Stable',], col='brown', ylim=ylim, xaxt='n', pch=16, cex=2, ylab='RFT Window', xlab='', main='Dots Task')
arrows(c(1,2), mnk['Stable',]-senk['Stable',], c(1,2), mnk['Stable',]+senk['Stable',], length=0, col='brown')
points(mnk['Decreasing',], col='orange', pch=16, cex=2)
arrows(c(1,2), mnk['Decreasing',]-senk['Decreasing',], c(1,2), mnk['Decreasing',]+senk['Decreasing',], length=0, col='orange')
axis(1, at=c(1,2), labels=c('Young', 'Old'))
legend('topleft', bty='n', col=c('brown', 'orange'), lty=1, pch=16, legend=c('Stable', 'Decreasing'))

dev.off()

nklm = lm(pars~age_group*condition, data=nkpars)
summary(nklm)
CI = confint(nklm)
capture.output(list(summary(nklm), CI=CI), file='output/nklm.txt')



# Ethics task -------------------------------------------------------------

d = read.csv('clean_ethics.csv')
models = data.frame()
Niter = 1
alreadyfit = list.files('fits/')
alreadyfit = sapply(alreadyfit, function(x) strsplit(x,split = '\\.')[[1]][1])

for(id in unique(d$id)){
  if(paste0('ethics_',id) %in% alreadyfit) next
  cat('----------', match(id,unique(d$id)), '/', length(unique(d$id)), '----------\n')
  sub     = d[d$id==id, ]
  sub     = sub[sub$RT < 20,]
  choices = sub$response
  RT      = sub$RT
  stim    = sub$norm_mean0*100
  failed  = F
  
  # # Visualize real subject curve
  # stimbin = dplyr::ntile(sub$trialintensity, 20)
  # mResp = tapply(choices, list(stimbin, sub$timebin4), mean)
  # plot(mResp[,1], type='b', col='red', xaxt='n', ylab='p(Blue)', xlab='', main=paste0('Empirical\nid=',id), ylim=c(0,1))
  # lines(mResp[,4], type='b', col='blue')
  # axis(1, at=c(1,20), labels=c('Very Purple', 'Very Blue'))
  # legend('topleft', bty='n', lty=1, pch=1, col=c('red', 'blue'), legend=c('First 200 Trials', 'Last 200 Trials'))
  
  bestLL  = 1e6
  bestFit = list()
  for(iter in 1:Niter){
    cat('*** iteration', iter, '/', Niter, '***\n')
    
    # specify cost fx 
    obfunc = function(params) {
      LL = RTFDDM(stim, RT, choices, params)
      if(any(abs(LL)==Inf)) return(1e6)
      -sum(LL)
    }
    # specify x0 that works
    lim = 1000
    for(x in 1:lim) {
      a0      = rgamma(1, 10, 1.5)
      ter0    = max(1, rnorm(1, 5, 3))
      nk0     = sample(5:100, size=1)
      x0      = c(a0, ter0, nk0)
      if(obfunc(x0)<1e6) break
      if(x==lim/2) cat('finding x0 taking a long time....\n')
      if(x==lim) cat('failed to find suitable x0\n')
    }
    
    # fit
    cat('optimizing...\n')
    cl = makeCluster(detectCores(), type='FORK')
    opt = optimParallel(par=x0, fn=obfunc, lower = c(0.1,.05,5), upper = c(50,50,800), 
                        # control = list(trace=6), 
                        parallel = list(cl=cl, loginfo=T))
    closeAllConnections()
    
    if(opt$value == 1e6) opt = optim(par=x0, fn=obfunc)
    
    if(opt$value <= bestLL) {
      bestLL = opt$value
      bestFit = opt
    }
  }
  
  if(bestLL==1e6) failed=T
  # Grid-search for nk
  est_a   = bestFit$par[1]
  est_ter = bestFit$par[2] 
  probs = c()
  grid = 5:100
  for(k in grid){
    cat('nk=',k,'||')
    probs= c(probs, obfunc(c(est_a, est_ter, k)))
  }
  cat('\n')
  probs[is.na(probs)] = 1e6
  plot(grid, probs, xlab='nk', ylab='LL', type='b', main=id) # check for convexity
  est_nk = grid[match(min(probs), probs)]
  legend('topright', bty='n', paste0('est nk=', est_nk))
  
  LL = min(probs)
  # Compute BIC
  bic = function(LL, k, n) k*log(n) - 2*LL
  aic = function(LL,k) -2*LL + 2*k
  
  npars = length(opt$par)
  RFTDDMBIC = bic(-LL, npars, nrow(sub))
  RFTDDMAIC = aic(-LL, npars)
  
  # Save
  subfit = data.frame(
    id = rep(id, npars), 
    condition = rep(sub$condition[1], npars), 
    model = rep('RFTDDM', npars), 
    parsname = c('a', 't0', 'nk'), 
    pars = c(est_a, est_ter, est_nk), 
    LL = rep(LL, npars), 
    BIC = rep(RFTDDMBIC, npars), 
    AIC = rep(RFTDDMBIC, npars)
  )
  if(failed) {
    subfit$pars = NA
    subfit$LL   = 1e6
    subfit$BIC  = NA
    subfit$AIC  = NA
  }
  saveRDS(list(opt=bestFit, nk=est_nk, df=subfit), paste0('fits/ethics_', id,'.rds'))
  models = rbind(models, subfit)
}
write.csv(models, 'ethics_fit.csv')


# Analyze parameters 

models$age_group = factor(ifelse(models$id<100, 'Old', 'Young'), levels=c('Young', 'Old'))
models$condition = factor(models$condition)
contrasts(models$age_group) = contr.sum(2)
contrasts(models$condition) = contr.sum(2)


cutoff = mean(apars$LL, na.rm=T) + 3*sd(apars$LL, na.rm=T)
apars  = apars[apars$LL < cutoff,]

# N excluded
tmp = models[!duplicated(models$id),]
nexc = tapply(tmp$BIC, list(tmp$condition, tmp$age_group), function(x) sum(is.na(x)))
sum(nexc)

# a 
apars = models[models$parsname=='a',]
ma  = tapply(apars$pars, list(apars$condition, apars$age_group), mean, na.rm=T)
sea = tapply(apars$pars, list(apars$condition, apars$age_group), plotrix::std.error, na.rm=T)

ylim = range(pretty(c(ma-sea,ma+sea)))
plot(ma['Stable',], col='brown', ylim=ylim, xaxt='n', pch=16, cex=2, ylab='Decision Threshold', xlab='', main='Ethics Task')
arrows(c(1,2), ma['Stable',]-sea['Stable',], c(1,2), ma['Stable',]+sea['Stable',], length=0, col='brown')
points(ma['Decreasing',], col='orange', pch=16, cex=2)
arrows(c(1,2), ma['Decreasing',]-sea['Decreasing',], c(1,2), ma['Decreasing',]+sea['Decreasing',], length=0, col='orange')
axis(1, at=c(1,2), labels=c('Young', 'Old'))
legend('topleft', bty='n', col=c('brown', 'orange'), lty=1, pch=16, legend=c('Stable', 'Decreasing'))

alm = lm(pars~age_group*condition, data=apars)
summary(alm)

# ter
t0pars = models[models$parsname=='t0',]

mt0  = tapply(t0pars$pars, list(t0pars$condition, t0pars$age_group), mean, na.rm=T)
set0 = tapply(t0pars$pars, list(t0pars$condition, t0pars$age_group), plotrix::std.error, na.rm=T)

ylim = range(pretty(c(mt0-set0,mt0+set0)))
plot(mt0['Stable',], col='brown', ylim=ylim, xaxt='n', pch=16, cex=2, ylab='Non-Decision Time', xlab='', main='Ethics Task')
arrows(c(1,2), mt0['Stable',]-set0['Stable',], c(1,2), mt0['Stable',]+set0['Stable',], length=0, col='brown')
points(mt0['Decreasing',], col='orange', pch=16, cex=2)
arrows(c(1,2), mt0['Decreasing',]-set0['Decreasing',], c(1,2), mt0['Decreasing',]+set0['Decreasing',], length=0, col='orange')
axis(1, at=c(1,2), labels=c('Young', 'Old'))
legend('topleft', bty='n', col=c('brown', 'orange'), lty=1, pch=16, legend=c('Stable', 'Decreasing'))

t0lm = lm(pars~age_group*condition, data=t0pars)
summary(t0lm)

# nk
nkpars = models[models$parsname=='nk',]

mnk  = tapply(nkpars$pars, list(nkpars$condition, nkpars$age_group), mean, na.rm=T)
senk = tapply(nkpars$pars, list(nkpars$condition, nkpars$age_group), plotrix::std.error, na.rm=T)

ylim = range(pretty(c(mnk-senk,mnk+senk)))
plot(mnk['Stable',], col='brown', ylim=ylim, xaxt='n', pch=16, cex=2, ylab='RFT Window', xlab='', main='Ethics Task')
arrows(c(1,2), mnk['Stable',]-senk['Stable',], c(1,2), mnk['Stable',]+senk['Stable',], length=0, col='brown')
points(mnk['Decreasing',], col='orange', pch=16, cex=2)
arrows(c(1,2), mnk['Decreasing',]-senk['Decreasing',], c(1,2), mnk['Decreasing',]+senk['Decreasing',], length=0, col='orange')
axis(1, at=c(1,2), labels=c('Young', 'Old'))
legend('topleft', bty='n', col=c('brown', 'orange'), lty=1, pch=16, legend=c('Stable', 'Decreasing'))

nklm = lm(pars~age_group*condition, data=nkpars)
summary(nklm)

