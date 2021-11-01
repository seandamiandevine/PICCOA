fitddm <- function(data, quickfit = F, verbose = T, checkConverge = T, 
                   savedir = 'Models/DDM/out/', nFit = 100, startfrom = 1,
                   iterList = NULL) {
  
  source('Models/DDM/ddm_loglik.R')
  source('Models/DDM/ensure_fit.R')
  source('Models/DDM/getx0.R')
  
  if(is.null(iterList)) iterList = unique(data$id)
  nsub = length(iterList)
  data$strength_bin = as.numeric(cut(data$trialintensity, 
                                breaks = quantile(data$trialintensity, 
                                                  probs = seq(0, 1, length = 6), 
                                                  na.rm = TRUE,
                                                  type = 2),
                                include.lowest = TRUE))
  data$strength_bin <- factor(data$strength_bin)
  ndrift <- length(levels(data$strength_bin))
  data$response_ddm <- ifelse(data$response == 1, 'upper', 'lower')
  
  d.fit <- data.frame()
  for(n in startfrom:nsub) { 
    
    cat(paste0("\n-------- FITTING SUBJECT ",n,"/", nsub," ------------\n")) # keep track 
    
    id = iterList[n]
    
    fit <- ensure_fit(data = data[data$id==id,], 
                      start_function = getx0, 
                      objective_function = ddm_loglik, 
                      base_pars = c("a", "t0", "sv", "sz", "z"), 
                      nFit=nFit,
                      quickfit = quickfit, verbose = verbose, checkConverge = checkConverge, 
                      savedir = savedir)
    
    fit <- cbind(data[data$id==id, c('id', 'age', 'condition', 'age_group')][1,], fit)
    d.fit <- rbind(d.fit, data.frame(fit, stringsAsFactors = F))
    
  }
  
  return(d.fit)
  
}