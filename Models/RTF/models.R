CDF = function(x, tau, sigma) {
  pnorm((x-tau)/sigma)
}

Choice = function(p) {
  sample(c(0,1), prob = c(1-p,p), size=1)
}

ControlMod = function(resp, x, params) {
  tau=params[1]
  sigma=params[2]
  N = length(x)
  probs = c()
  for(t in 1:N){
    p = CDF(x[t], tau, sigma)
    if(resp[t]==0) p=1-p
    probs[t] = p
  }
  probs
}

MovingWindow = function(resp, x, params) {
  tau0 = 0.5 # params[1]
  sigma = params[1]
  alpha =  params[2]
  
  N = length(x)
  probs=c()
  tau = c(tau0)
  for(t in 1:N){
    tau[t+1] = tau[t] + alpha*(x[t]-tau[t])
    p = CDF(x[t], tau[t], sigma) # fat
    if(resp[t]==0) p = 1-p       # not fat
    probs[t] = p
  }
  probs
}

RTF = function(resp, x, params) {
  tau=params[1]
  sigma=params[2]
  nk=round(params[3])
  w=.5 # fixed for now, like David
  N = length(x)
  probs = c()
  for(t in 1:N){
    if(t <= nk) {
      probs[t] = CDF(x[t], tau, sigma) 
      next 
    }
    thisRange = (t-nk):t
    maxX = max(x[thisRange])
    minX = min(x[thisRange])
    idx = match(x[t], x[thisRange])
    rank = rank(x[thisRange])[idx]
    
    range = (x[t]-minX)/(maxX-minX)
    if(is.na(range)) range=0
    freq = (rank-1)/(length(thisRange)-1)
    y = w*range + (1-w)*freq
    
    p = CDF(y*100, tau, sigma)
    if(resp[t]==0) p=1-p
    
    probs[t] = p
  }
  probs 
}