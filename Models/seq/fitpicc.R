fitpicc = function(data, stim, checkConverge = T, quickfit = F, verbose = T, maxIter = 10000, 
                   savedir = 'Models/seq/out/', startfrom = 1, iterList=NULL) { 
  
  source('Models/seq/loglik.R')
  source('Models/seq/getx0.R')
  
  if(quickfit) checkConverge = F # can't be too picky
  nParams = 6
  # if(quickfit) {
  #   paramSweep = t(rep(1, nParams))
  # } else {
  #   paramGrid = seq(-5, 5, length.out = nParams) 
  #   paramSweep = gtools::permutations(v=paramGrid, 
  #                                     r=nParams, 
  #                                     n = length(paramGrid),
  #                                     repeats.allowed = F)
  #   }
  # 
  
  nFit = ifelse(quickfit, 10, 100)
  out <- matrix(NA, length(iterList), 12)
  if(is.null(iterList)) iterList = unique(data$id)
  
  for(id in startfrom:length(iterList)) {  # compute parameters for each subject
    if(verbose==F) cat('computing parameters for subject ', id, '\n') 
    
    d.sub <- data[data$id == iterList[id], ]
    
    switch (stim,
      'dots' = {
        red = as.numeric(gsub(",", "", substr(d.sub$trialstim, 2, 3)))
        # green <- as.numeric(gsub(",", "", substr(d.sub$trialstim, 5, 6)))
        blue = d.sub$trialintensity + 154 
        f = 1 - 2*((red - min(red))/(max(red) - min(red))) 
        c = as.numeric(ifelse(d.sub$response == 1, 1, -1)) # 1 is unethical; -1 is ethical
      }, 
      'ethics' = {
        f = 2* ((d.sub$norm_mean - min(d.sub$norm_mean))/(max(d.sub$norm_mean) - min(d.sub$norm_mean))) - 1
        c = as.numeric(ifelse(d.sub$response == 1, 1, -1))
      }, 
      errorCondition('Type of stimulus not specified (dots or ethics).')
    )
  
    bestVal = Inf 
    bestConv = NA
    bestParams = c()
    bestBIC = NA
    LB = c(-100, -100, -100, -100, 0, 0) # lower bound
    UB = c(100, 100, 100, 100, 1, 1)     # upper bound
    
    obFunc <- function(x) LL(f, c, x) 
    bic = function(ll, n, k) -2 * ll + log(n) * k # calculates BIC
    
    for(paramIter in 1:nFit) { 
      x0 = getx0(f, c, nParams, LB, UB, obFunc, maxIter)
      Xfit = nlminb(start = x0, objective = obFunc, lower = LB, upper = UB) 
      if(checkConverge) {
        if(Xfit$convergence==0 & Xfit$objective < bestVal) {
          bestVal = Xfit$objective
          bestParams = Xfit$par
          bestConv = Xfit$convergence # should always be 0
          bestBIC = bic(Xfit$objective, length(f), nParams)
        }
      } else {
        if(Xfit$objective < bestVal) {
          bestVal = Xfit$objective
          bestParams = Xfit$par
          bestConv = Xfit$convergence # may not be 0
          bestBIC = bic(-Xfit$objective, length(f), nParams)
        }
      }
      foundParam = ifelse(Xfit$objective < 1e10, 'yes.', 'no.')
      cat('fit params for subject', id, '/', length(iterList), '|',
            'x0 combo', paramIter, '/', nFit, '|' ,
           'convergence:', Xfit$convergence, '\n')
      if(verbose) cat('--- x0:', paste(x0, collapse = ', '), '---\n\n')
    }
    
    out[id, ] = c(d.sub$id[1], as.character(d.sub$age_group[1]), as.character(d.sub$condition[1]), 
                   bestParams, bestVal, bestBIC, bestConv)
    
    if (savedir!=FALSE) {
      if(quickfit) {qf = '_quickfit'} else {qf = ''}
      thisid = d.sub$id[1]
      save(list=c("thisid", "bestParams", "bestVal",  "bestBIC","bestConv"), 
           file=paste0(savedir,'/',stim,'/', 'seq_fit_', 
                       as.character(d.sub$age_group[1]), '_',
                       as.character(d.sub$condition[1]), '_',
                       as.character(d.sub$id[1]), 
                       qf,".RData"))
    }
    
    on.exit(return(out)) # get something if terminated
    
  }
  out = data.frame(out, stringsAsFactors = F)
  colnames(out) =  c('id', 'age_group', 'condition', 
                     'B0', 'Bf', 'BF', 'Bc', 'lF', 'lc', 'loglik', 'BIC', 'convergence')
  out[, 4:ncol(out)] = sapply(out[,4:ncol(out)], as.numeric)
  
  return(out)
}
