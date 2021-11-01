
set.seed(2021)
setwd('~/Documents/PICCOA/Models/RTF')
source('models.R')
d = read.csv('clean_dots.csv')

# Simulate ----------------------------------------------------------------

Nsub = 25
rand_decrease = as.character(sample(unique(d$id[d$condition=='Decreasing']), Nsub))
rand_stable = as.character(sample(unique(d$id[d$condition=='Stable']), Nsub))
dsim = d[d$id%in%c(rand_decrease, rand_stable), ]

Ctl = data.frame()
RF = data.frame()
for(id in unique(dsim$id)){
  cat('Simulating choices for id', match(id, unique(dsim$id)), '/', Nsub*2, '\n')
  stim = dsim[dsim$id==id, 'trialintensity']
  N = length(stim)
  spectrum = min(stim):max(stim)
  
  # Control model
  tau = runif(1, 0.001, 100)
  sigma = runif(1, 0.001, 50)
  resps = c()
  probs = c()
  for(t in 1:N){
    p = CDF(stim[t], spectrum, tau, sigma)
    probs[t] = p
    resps[t] = Choice(p)
  }
  thissim = data.frame(id = match(id, unique(dsim$id)),
                       condition = as.character(dsim[dsim$id==id, 'condition'][1]),
                       true_tau = tau, 
                       true_sigma = sigma, 
                       stim = stim,
                       choices = resps, 
                       probs = probs,
                       model = 'Ctl')
  Ctl = rbind(Ctl, thissim)
  
  # RTF Model
  tau = runif(1, 0.01, 100)
  sigma = runif(1, 0.01, 50)
  nk = sample(5:75, size=1)
  w = .5
  resps = c()
  probs = c()
  for(t in 1:N){
    if(t <= 5) {
      probs[t] = CDF(stim[t], spectrum, tau, sigma) 
      resps[t] = Choice(probs[t])
      next 
    } else if(t < nk){
      thisRange = 1:t
    } else {
      thisRange = (t-nk):t
    }
    maxX = max(stim[thisRange])
    minX = min(stim[thisRange])
    idx = match(stim[t], stim[thisRange])
    rank = rank(stim[thisRange])[idx]
    
    range = (stim[t]-minX)/(maxX-minX)
    if(is.na(range)) range = 0
    freq = (rank-1)/(nk-1)
    y = w*range + (1-w)*freq
    y = round(y*max(stim))
    y = ifelse(y<min(stim), min(stim), ifelse(y>max(stim), max(stim), y))
    
    p = CDF(y, spectrum, tau, sigma)
    probs[t] = p
    resps[t] = Choice(p)
  }
  thissim = data.frame(id = match(id, unique(dsim$id)),
                       condition = as.character(dsim[dsim$id==id, 'condition'][1]), 
                       true_tau0 = tau, 
                       true_sigma = sigma, 
                       true_nk = nk,
                       stim=stim,
                       choices = resps,
                       probs = probs,
                       model = 'RTF')
  RF = rbind(RF, thissim)
}

# Visualize simulations 
Ctl$trial = unlist(by(Ctl$id, Ctl$id, FUN = function(x) 1:length(x)))
Ctl$timebin = dplyr::ntile(Ctl$trial, 4)
Ctl$stimbin = dplyr::ntile(Ctl$stim, 20)
RF$trial = unlist(by(RF$id, RF$id, FUN = function(x) 1:length(x)))
RF$timebin = dplyr::ntile(RF$trial, 4)
RF$stimbin = dplyr::ntile(RF$stim, 20)

pdf('Simulations.pdf')
layout(matrix(1:4, 2,2, byrow=T))
# Ctl model
# Stable
mProb = tapply(Ctl$probs, list(Ctl$stimbin, Ctl$timebin, Ctl$condition), mean)
plot(mProb[,1,'Stable'], type='b', col='red', xaxt='n', ylab='p(Blue)', xlab='', main='Control\nStable', ylim=c(0,1))
lines(mProb[,4,'Stable'], type='b', col='blue')
axis(1, at=c(1,20), labels=c('Very Purple', 'Very Blue'))
legend('topleft', bty='n', lty=1, pch=1, col=c('red', 'blue'), legend=c('First 200 Trials', 'Last 200 Trials'))

# Decreasing
plot(mProb[,1,'Decreasing'], type='b', col='red', xaxt='n', ylab='p(Blue)', xlab='', main='Control\nDecreasing', ylim=c(0,1))
lines(mProb[,4,'Decreasing'], type='b', col='blue')
axis(1, at=c(1,20), labels=c('Very Purple', 'Very Blue'))
#legend('topleft', bty='n', lty=1, pch=1, col=c('red', 'blue'), legend=c('First 200 Trials', 'Last 200 Trials'))

# RF model
# Stable
mProb = tapply(RF$probs, list(RF$stimbin, RF$timebin, RF$condition), mean)
plot(mProb[,1,'Stable'], type='b', col='red', xaxt='n', ylab='p(Blue)', xlab='', main='RTF\nStable', ylim=c(0,1))
lines(mProb[,4,'Stable'], type='b', col='blue')
axis(1, at=c(1,20), labels=c('Very Purple', 'Very Blue'))
#legend('topleft', bty='n', lty=1, pch=1, col=c('red', 'blue'), legend=c('First 200 Trials', 'Last 200 Trials'))

# Decreasing
plot(mProb[,1,'Decreasing'], type='b', col='red', xaxt='n', ylab='p(Blue)', xlab='', main='RTF\nDecreasing', ylim=c(0,1))
lines(mProb[,4,'Decreasing'], type='b', col='blue')
axis(1, at=c(1,20), labels=c('Very Purple', 'Very Blue'))
#legend('topleft', bty='n', lty=1, pch=1, col=c('red', 'blue'), legend=c('First 200 Trials', 'Last 200 Trials'))

dev.off()

# Fit ---------------------------------------------------------------------

Niter = 1

# Control model
cat('------------Fitting Simulated Control Data------------\n')
RecoveryCtl = data.frame()
for(id in unique(Ctl$id)) {
  cat('----------', match(id,unique(Ctl$id)), '/', length(unique(Ctl$id)), '----------\n')
  sub = Ctl[Ctl$id==id, ]
  true_tau = sub$true_tau[1]
  true_sigma = sub$true_sigma[1]
  choices = sub$choices
  stim = sub$stim
  
  bestCtlLL = 1e6
  bestCtlOpt = list()
  for(iter in 1:Niter){
    cat('*** iteration', iter, '/', Niter, '***\n')
    obfunc = function(x) {
      LL = -sum(log(ControlMod(choices, stim, x)))
      min(LL, 1e6)
    }
    while(1){
      tau0 = runif(1, 0.001, 100)
      sigma0 = runif(1, 0.001, 100)
      x0 = c(tau0, sigma0)
      if(obfunc(x0)!=1e6 & !is.na(obfunc(x0))) break
    }
    optControl = optim(par=x0, fn=obfunc)
    if(optControl$value<bestCtlLL) {
      bestCtlLL = optControl$value
      bestCtlOpt = optControl
    }
  }
  if(length(bestCtlOpt)==0) {
    cat('FAILED TO FIT!\n')
    next  # failed to fit 
  }
  thisRecov = data.frame(id=id, true_tau=true_tau, true_sigma=true_sigma, 
                         est_tau=bestCtlOpt$par[1], est_sigma=bestCtlOpt$par[2]) 
  RecoveryCtl = rbind(RecoveryCtl, thisRecov)
}

pdf('CtlRecoveryPlots.pdf')
r = cor(RecoveryCtl$true_tau, RecoveryCtl$est_tau)
plot(RecoveryCtl$true_tau, RecoveryCtl$est_tau, 
     xlab=expression('True'~tau), ylab=expression('Estimated'~tau), 
     main='Regular CDF') 
legend('topleft', bty='n', legend=paste0('r=',round(r,4)))

r = cor(RecoveryCtl$true_sigma, RecoveryCtl$est_sigma)
plot(RecoveryCtl$true_sigma, RecoveryCtl$est_sigma, 
     xlab=expression('True'~sigma), ylab=expression('Estimated'~sigma), 
     main='Regular CDF') 
legend('topleft', bty='n', legend=paste0('r=',round(r,4)))

dev.off()

# RTF model 
cat('------------Fitting Simulated RTF Data------------\n')
RecoveryRTF = data.frame()
for(id in unique(RF$id)) {
  cat('----------', match(id,unique(RF$id)), '/', length(unique(RF$id)), '----------\n')
  sub = RF[RF$id==id, ]
  true_tau = sub$true_tau0[1]
  true_sigma = sub$true_sigma[1]
  true_nk = sub$true_nk[1]
  choices = sub$choices
  stim = sub$stim
  
  bestRTFLL = 1e6
  bestRTFOpt = list()
  for(iter in 1:Niter){
    cat('*** iteration', iter, '/', Niter, '***\n')
    obfunc = function(params) {
      LL = -sum(log(RTF(choices, stim, params)))
      min(LL, 1e6)
    }
    while(1){
      tau0 = runif(1, 0.001, 100)
      sigma0 = runif(1, 0.001, 100)
      nk0 = sample(5:75, size=1)
      x0 = c(tau0, sigma0, nk0)
      if(obfunc(x0)!=1e6 & !is.na(obfunc(x0))) break
    }
    optRTF = optim(par=x0, fn=obfunc)
    if(optRTF$value < bestRTFLL) {
      bestRTFLL = optRTF$value
      bestRTFOpt = optRTF
    }
  }
  if(length(bestRTFOpt)==0) {
    cat('FAILED TO FIT!\n')
    next  # failed to fit 
  }
  
  # Grid-search for nk
  est_tau = bestRTFOpt$par[1]
  est_sigma = bestRTFOpt$par[2]
  probs = c()
  grid = 5:75
  for(k in grid){
    probs= c(probs, obfunc(c(est_tau, est_sigma, k)))
  }
  probs[is.na(probs)] = 1e6
  plot(grid, probs, xlab='nk', ylab='LL', type='b') # check for convexity
  est_nk = grid[match(min(probs), probs)]
  legend('topright', bty='n',legend=paste0('true=',true_nk,'\nest=',est_nk))
  thisRecov = data.frame(id=id, true_tau=true_tau, true_sigma=true_sigma, true_nk=true_nk, 
                         est_tau=est_tau, est_sigma=est_sigma, est_nk=est_nk) 
  RecoveryRTF = rbind(RecoveryRTF, thisRecov)
}

pdf('RTFRecoveryPlots.pdf')

r = cor(RecoveryRTF$true_tau, RecoveryRTF$est_tau)
plot(RecoveryRTF$true_tau, RecoveryRTF$est_tau, 
     xlab=expression('True'~tau), ylab=expression('Estimated'~tau), 
     main='RTF') 
legend('topleft', bty='n', legend=paste0('r=',round(r,4)))

r = cor(RecoveryRTF$true_sigma, RecoveryRTF$est_sigma)
plot(RecoveryRTF$true_sigma, RecoveryRTF$est_sigma, 
     xlab=expression('True'~sigma), ylab=expression('Estimated'~sigma), 
     main='RTF') 
legend('topleft', bty='n', legend=paste0('r=',round(r,4)))

r = cor(RecoveryRTF$true_nk, round(RecoveryRTF$est_nk))
plot(RecoveryRTF$true_nk, round(RecoveryRTF$est_nk), 
     xlab=expression('True'~nk), ylab=expression('Estimated'~nk), 
     main='RTF') 
legend('topleft', bty='n', legend=paste0('r=',round(r,4)))

dev.off()
