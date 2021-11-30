library(rtdists)
set.seed(2021)
source('models.R')
d = read.csv('clean_dots.csv')

# Simulate ----------------------------------------------------------------

Nsub = 25
rand_decrease = as.character(sample(unique(d$id[d$condition=='Decreasing']), Nsub))
rand_stable = as.character(sample(unique(d$id[d$condition=='Stable']), Nsub))
dsim = d[d$id%in%c(rand_decrease, rand_stable), ]

RFDDM = data.frame()
for(id in unique(dsim$id)){
  cat('Simulating choices for id', match(id, unique(dsim$id)), '/', Nsub*2, '\n')
  stim = dsim[dsim$id==id, 'trialintensity']
  N = length(stim)
  
  # RTF DDM
  a = rgamma(1, 2, 1.5)
  ter = abs(rnorm(1, 0.2, .1))
  nk = sample(5:100, size=1)
  w = .5
  N = length(stim)
  resps = c()
  RTs   = c()
  minV  = -5
  maxV  = 5
  
  for(t in 1:N) {
    # range normalization
    if(t <= nk){
      y = stim[t]
    } else {
      thisRange = (t-nk):t
      maxX = max(stim[thisRange-1])
      minX = min(stim[thisRange-1])
      rank  = match(stim[t], sort(c(stim[t], stim[thisRange])))
      
      range = (stim[t]-minX)/(maxX-minX-1)
      if(is.na(range)) range = 0
      freq = (rank-1)/(nk-1)
      y = w*range + (1-w)*freq
      y = round(y*max(stim))
      y = ifelse(y<min(stim), min(stim), ifelse(y>max(stim), max(stim), y))
    }
    
    # Diffusion process
    v   = ((stim[t]-min(stim))/(max(stim)-min(stim))) * (maxV-minV) + minV
    pz  = (y-min(stim))/(max(stim)-min(stim))
    pz  = ifelse(pz<0.05, .05, ifelse(pz>.95, .95, pz))
    z   = a*pz
    while(1) {
      ddm = rdiffusion(1, a=a, v=v, t0=ter, z=z)
      if(ddm[['rt']] > ter) break
    }
    
    RTs   = c(RTs, ddm[['rt']])
    resps = c(resps, ifelse(ddm[['response']]=='upper', 1, 0))
  }
  
  thissim = data.frame(id = match(id, unique(dsim$id)),
                         condition = as.character(dsim[dsim$id==id, 'condition'][1]), 
                         true_a = a, 
                         true_ter = ter, 
                         true_nk = nk,
                         stim=stim,
                         choices = resps,
                         RT = RTs,
                         model = 'RTFDDM')
  RFDDM = rbind(RFDDM, thissim)
}

# write.csv(RFDDM, 'simulation_data.csv', row.names = F)

# Visualize simulations ---------------------------------------------------

# RFDDM = read.csv('simulation_data.csv')

pdf('plots/simulations1.pdf', width = 12, height=6 )

# RFDDM model
# overall
layout(matrix(1:2, 1, 2,byrow=T))

RFDDM$trial = unlist(by(RFDDM$id, RFDDM$id, FUN = function(x) 1:length(x)))
RFDDM$timebin = dplyr::ntile(RFDDM$trial, 4)
RFDDM$stimbin = dplyr::ntile(RFDDM$stim, 10)
RFDDM$abin = cut(RFDDM$true_a, 2, labels=c('Low', 'High'))

mResp = tapply(RFDDM$choices, list(RFDDM$stimbin, RFDDM$timebin, RFDDM$condition), mean)

# Stable
plot(rownames(mResp), mResp[,1,'Stable'], type='b', ylab = 'p(Blue)', xlab='', xaxt='n', col='red', 
     main = 'Stable Prevalence', ylim=c(0,1))
axis(1, at=c(1,max(RFDDM$stimbin)), labels = c('Very Purple', 'Very Blue'))
lines(rownames(mResp), mResp[,4,'Stable'], type='b', col='blue')
legend('topleft', bty='n', lty=1, pch=1, col=c('red', 'blue'), legend=c('First 200 trials', 'Last 200 trials'))

# Decreasing
plot(rownames(mResp), mResp[,1,'Decreasing'], type='b', ylab = 'p(Blue)', xlab='', xaxt='n', col='red', 
     main = 'Decreasing Prevalence', ylim=c(0,1))
axis(1, at=c(1,max(RFDDM$stimbin)), labels = c('Very Purple', 'Very Blue'))
lines(rownames(mResp), mResp[,4,'Decreasing'], type='b', col='blue')
#legend('topleft', bty='n', lty=1, pch=1, col=c('red', 'blue'), legend=c('First 200 trials', 'Last 200 trials'))

dev.off()

# Per abin
pdf('plots/simulations2.pdf', width = 8, height=6 )

layout(matrix(1:4, 2, 2,byrow=T))

mResp = tapply(RFDDM$choices, list(RFDDM$stimbin, RFDDM$timebin, RFDDM$condition, RFDDM$abin), mean)

# Low
# Stable
plot(rownames(mResp), mResp[,1,'Stable', 'Low'], type='b', ylab = 'p(Blue)', xlab='', xaxt='n', col='red', 
     main = 'Stable\n Low a', ylim=c(0,1))
axis(1, at=c(1,max(RFDDM$stimbin)), labels = c('Very Purple', 'Very Blue'))
lines(rownames(mResp), mResp[,4,'Stable', 'Low'], type='b', col='blue')
legend('topleft', bty='n', lty=1, pch=1, col=c('red', 'blue'), legend=c('First 200 trials', 'Last 200 trials'))

# Decreasing
plot(rownames(mResp), mResp[,1,'Decreasing', 'Low'], type='b', ylab = 'p(Blue)', xlab='', xaxt='n', col='red', 
     main = 'Decreasing \nLow a', ylim=c(0,1))
axis(1, at=c(1,max(RFDDM$stimbin)), labels = c('Very Purple', 'Very Blue'))
lines(rownames(mResp), mResp[,4,'Decreasing', 'Low'], type='b', col='blue')
#legend('topleft', bty='n', lty=1, pch=1, col=c('red', 'blue'), legend=c('First 200 trials', 'Last 200 trials'))

# High
# Stable
plot(rownames(mResp), mResp[,1,'Stable', 'High'], type='b', ylab = 'p(Blue)', xlab='', xaxt='n', col='red', 
     main = 'Stable\nHigh a', ylim=c(0,1))
axis(1, at=c(1,max(RFDDM$stimbin)), labels = c('Very Purple', 'Very Blue'))
lines(rownames(mResp), mResp[,4,'Stable', 'High'], type='b', col='blue')
#legend('topleft', bty='n', lty=1, pch=1, col=c('red', 'blue'), legend=c('First 200 trials', 'Last 200 trials'))

# Decreasing
plot(rownames(mResp), mResp[,1,'Decreasing', 'High'], type='b', ylab = 'p(Blue)', xlab='', xaxt='n', col='red', 
     main = 'Decreasing\nHigh a', ylim=c(0,1))
axis(1, at=c(1,max(RFDDM$stimbin)), labels = c('Very Purple', 'Very Blue'))
lines(rownames(mResp), mResp[,4,'Decreasing', 'High'], type='b', col='blue')
#legend('topleft', bty='n', lty=1, pch=1, col=c('red', 'blue'), legend=c('First 200 trials', 'Last 200 trials'))

dev.off()

# Fit ---------------------------------------------------------------------

library(optimParallel) # for fast optimization

Niter = 1

# RTF model 
cat('------------Fitting Simulated RTFDDM Data------------\n')
Recovery = data.frame()
for(id in unique(RFDDM$id)) {
  cat('\n----------', match(id,unique(RFDDM$id)), '/', length(unique(RFDDM$id)), '----------\n')
  sub         = RFDDM[RFDDM$id==id, ]
  true_a      = sub$true_a[1]
  true_ter    = sub$true_ter[1]
  true_nk     = sub$true_nk[1]
  choices     = sub$choices
  RT          = sub$RT
  stim        = sub$stim
  
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
    while(1) {
      a0      = rgamma(1, 2, 1.5)
      ter0    = abs(rnorm(1, 0.2, .1))
      nk0     = sample(5:100, size=1)
      x0      = c(a0, ter0, nk0)
      if(obfunc(x0)<1e6) break
    }
    
    # fit
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
  plot(grid, probs, xlab='nk', ylab='LL', type='b') # check for convexity
  est_nk = grid[match(min(probs), probs)]
  legend('topright', bty='n',legend=paste0('true=',true_nk,'\nest=',est_nk))
  
  thisRecov = data.frame(id=id,
                         true_a=true_a, true_ter=true_ter, true_nk=true_nk, 
                         est_a = est_a, est_ter = est_ter, est_nk = est_nk, 
                         LL = min(probs)
                         )
  Recovery = rbind(Recovery, thisRecov)
}

write.csv(Recovery, 'paramRecovery.csv', row.names = F)

# Visualize recovery ---------------------------------------------------------------------

pdf('plots/recovery.pdf')

Recovery$condition = ifelse(Recovery$id<=Nsub, 'Stable', 'Decreasing')

# a
plot(Recovery$true_a, Recovery$est_a, xlab=expression(a), ylab=expression(hat(a)))
r = cor(Recovery$true_a, Recovery$est_a)
legend('topleft', bty='n', legend=paste('r=',round(r,4)))

# ter
plot(Recovery$true_ter, Recovery$est_ter, xlab=expression(t0), ylab=expression(hat(t0)))
r = cor(Recovery$true_ter, Recovery$est_ter)
legend('topleft', bty='n', legend=paste('r=',round(r,4)))

# nk
layout(matrix(c(1,1,2,3), nrow=2,ncol=2,byrow = T))
# overall
plot(Recovery$true_nk, Recovery$est_nk, xlab=expression(nk), ylab=expression(hat(nk)), main='Overall')
r = cor(Recovery$true_nk, Recovery$est_nk)
legend('topleft', bty='n', legend=paste('r=',round(r,4)))

# by condition
recov1 = Recovery[Recovery$condition=='Stable', ]
plot(recov1$true_nk, recov1$est_nk, xlab=expression(nk), ylab=expression(hat(nk)), main='Stable')
r = cor(recov1$true_nk, recov1$est_nk)
legend('topleft', bty='n', legend=paste('r=',round(r,4)))

recov2 = Recovery[Recovery$condition=='Decreasing', ]
plot(recov2$true_nk, recov2$est_nk, xlab=expression(nk), ylab=expression(hat(nk)),main='Decreasing')
r = cor(recov2$true_nk, recov2$est_nk)
legend('topleft', bty='n', legend=paste('r=',round(r,4)))

dev.off()
