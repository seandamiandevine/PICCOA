loglik = function(f, c, params, return='loglik') { # calculates log likelihood of choice in model-based analysis below
  
  Fbar = c(0) # list of exponentially weighted past stimuli 
  Cbar = c(0) # list of exponentially weighted past responses 
  CP = c()    # choice probabilities 
  
  for(t in 1:length(f)) {
    
    # predict choice probability from regression parameters
    dV = params[1] + params[2]*f[t] + params[3]*Fbar[t] + params[4]*Cbar[t] 
    P = 1/ (1 + exp(dV)) 
    
    CP[t] <- ifelse(c[t] == -1, P, 1-P)
    
    Fbar[t+1] <- params[5]*Fbar[t] + f[t] # update list weighing it by decay parameter lambda
    Cbar[t+1] <- params[6]*Cbar[t] + c[t]
  }
  logl = -sum(log(CP)) # use negative log to find minima later
  if(logl == Inf) logl = 1e10
  
  switch (return,
    'loglik' = {
      return(logl)
    }, 
    'cp' = {
      return(CP)
    }, 
    'all' = {
      return(list(cp=CP, loglik = loglik))
    }
  )
}

LL = function(f, c, params) tryCatch(loglik(f, c, params), error=function(e) cat('error computing model with params:', params))
