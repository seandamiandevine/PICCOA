
setwd('~/Documents/PICCOA/Models/RTF')
source('models.R')

# Dots --------------------------------------------------------------------

cat('********************** DOTS **********************\n')

d = read.csv('clean_dots.csv')
models = data.frame()
Niter = 1
alreadyfit = list.files('out/')
alreadyfit = sapply(alreadyfit, function(x) strsplit(x,split = '\\.')[[1]][1])

for(id in unique(d$id)){
  if(paste0('dots_',id) %in% alreadyfit) next
  cat('----------', match(id,unique(d$id)), '/', length(unique(d$id)), '----------\n')
  sub = d[d$id==id, ]
  choices = sub$response
  stim = sub$trialintensity
  
  # # Visualize real subject curve
  # stimbin = dplyr::ntile(sub$trialintensity, 20)
  # mResp = tapply(choices, list(stimbin, sub$timebin4), mean)
  # plot(mResp[,1], type='b', col='red', xaxt='n', ylab='p(Blue)', xlab='', main=paste0('Empirical\nid=',id), ylim=c(0,1))
  # lines(mResp[,4], type='b', col='blue')
  # axis(1, at=c(1,20), labels=c('Very Purple', 'Very Blue'))
  # legend('topleft', bty='n', lty=1, pch=1, col=c('red', 'blue'), legend=c('First 200 Trials', 'Last 200 Trials'))

  bestCtlLL = 1e6
  bestCtlOpt = list()
  bestRTFLL = 1e6
  bestRTFOpt = list()
  for(iter in 1:Niter) {
    cat('*** iteration', iter, '/', Niter, '***\n')
    # Control model 
    cat('fitting Ctl model || ')
    
    obfunc = function(params) {
      LL = -sum(log(ControlMod(choices, stim, params=params)))
      min(LL, 1e6)
    }
    while(1){
      tau0 = runif(1, 0.001, 100)
      sigma0 = runif(1, 0.001, 100)
      x0 = c(tau0, sigma0)
      if(!is.na(obfunc(x0)) & obfunc(x0)<1e6) break
    }
    optControl = optim(par=x0, fn=obfunc)
    if(optControl$value<bestCtlLL) {
      bestCtlLL = optControl$value
      bestCtlOpt = optControl
    }
    
    # RTF model 
    cat('fitting RTF model \n')
    obfunc = function(params) {
      LL = -sum(log(RTF(choices, stim, params)))
      min(LL, 1e6)
    }
    while(1){
      tau0 = runif(1, 0.001, 100)
      sigma0 = runif(1, 0.001, 100)
      nk0 = runif(1, 1, 100)
      x0 = c(tau0, sigma0, nk0)
      if(!is.na(obfunc(x0)) & obfunc(x0)<1e6) break
    }
    optRTF = optim(par=x0, fn=obfunc)
    if(optRTF$value<bestRTFLL) {
      bestRTFLL = optRTF$value
      bestRTFOpt = optRTF
    }
    # Grid-search for nk
    est_tau = bestRTFOpt$par[1]
    est_sigma = bestRTFOpt$par[2]
    probs = c()
    for(k in 1:100){
      probs[k] = obfunc(c(est_tau, est_sigma, k))
    }
    plot(1:100, probs, xlab='nk', ylab='LL', type='b') # check for convexity
    probs[is.na(probs)] = 1e6
    est_nk = match(min(probs), probs)
    LLRTF = probs[est_nk]
  }
  
  # Compute BIC
  bic = function(LL, k, n) k*log(n) - 2*LL
  aic = function(LL,k) -2*LL + 2*k
  
  CtlBIC = bic(-bestCtlLL, length(bestCtlOpt$par), nrow(sub))
  RTFBIC = bic(-LLRTF, length(bestRTFOpt$par), nrow(sub))
  CtlAIC = aic(-bestCtlLL, length(bestCtlOpt$par))
  RTFAIC = aic(-LLRTF, length(bestRTFOpt$par))
  RelativeLik = (CtlAIC-RTFAIC)/2
  
  # Save
  models_sub = data.frame(
    id = rep(id, 5), 
    condition = rep(sub$condition[1], 5), 
    model = c(rep('Ctl', 2), rep('RTF', 3)), 
    parsname = c('tau', 'sigma_ctl', 'tau0', 'sigma_rtf', 'nk'), 
    pars = c(bestCtlOpt$par, c(est_tau, est_sigma, est_nk)), 
    LL = -c(rep(bestCtlOpt$value, 2), rep(LLRTF, 3)), 
    convergence = c(rep(bestCtlOpt$convergence, 2), rep(bestRTFOpt$convergence, 3)), 
    BIC = c(rep(CtlBIC, 2), rep(RTFBIC, 3)), 
    AIC = c(rep(CtlAIC, 2), rep(RTFAIC, 3)), 
    RelLik = rep(RelativeLik, 5)
  )
  saveRDS(list('Ctl'=bestCtlOpt, 'RTF'=bestRTFOpt), paste0('out/dots_', id,'.rds'))
  models = rbind(models, models_sub)
}
write.csv(models, 'dots_fit.csv')


# Ethics ------------------------------------------------------------------
