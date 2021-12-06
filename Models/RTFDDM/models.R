RTFDDM = function(stim, RT, choices, params) {
  library(rtdists)
  
  a      = params[1]
  ter    = params[2]
  nk     = round(params[3])
  m      = params[4]
  w      = .5
  N      = length(stim)
  LL     = c()
  prec   = max(0, nchar(strsplit(as.character(stim), "\\.")[[1]][2]), na.rm=T)
  ys = c()
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
      y = round(y*max(stim), prec)
      y = ifelse(y<min(stim), min(stim), ifelse(y>max(stim), max(stim), y))
      ys[t] = y
    }
    v = ifelse(stim[t] < 50, stim[t]*-m, stim[t]*m)  
    pz  = (y-min(stim))/(max(stim)-min(stim))
    pz  = ifelse(pz<0.05, .05, ifelse(pz>.95, .95, pz))
    z   = a*pz
    
    # Compute likelihood
    resp  = ifelse(choices[t]==1, 'upper', 'lower')
    lik   = ddiffusion(rt = RT[t], response=resp, a=a, v=v, t0=ter, z=z)
    LL[t] = log(lik)
  }
  
  return(LL)
}

CtlDDM = function(stim, RT, choices, params) {
  library(rtdists)
  
  a      = params[1]
  ter    = params[2]
  m      = params[3]
  pz     = 0.5
  z      = pz*a
  N      = length(stim)
  LL     = c()
  for(t in 1:N) {
    v = ifelse(stim[t] < 50, stim[t]*-m, stim[t]*m)  
    
    # Compute likelihood
    resp  = ifelse(choices[t]==1, 'upper', 'lower')
    lik   = ddiffusion(rt = RT[t], response=resp, a=a, v=v, t0=ter, z=z)
    LL[t] = log(lik)
    
  }
  
  return(LL)
}
