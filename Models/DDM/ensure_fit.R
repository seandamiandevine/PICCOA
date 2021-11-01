ensure_fit <-  function(data, start_function, objective_function, nFit, 
                        base_pars, n_drift = 5, quickfit = F, 
                        verbose = T, checkConverge = T, 
                        savedir = 'Models/DDM/out/',
                        lower = c(rep(0, length(base_pars)), -Inf,
                                  rep(-Inf,length(start_function(base_pars, 5))-length(base_pars)))) {
  
  n_fits = ifelse(quickfit, 1, nFit)
  if(quickfit) checkConverge = F # can't be too picky 
  best_fit <- list(objective = 1e+06)
  for (i in 1:n_fits) {
    start_ll <- 1e+06
    cat("\nfinding viable start point for new fit.\n") 
    while(start_ll == 1e+06) {
      start <- start_function(base_pars, n_drift=n_drift)
      start_ll <- objective_function(start,
                                     rt = data$RT, response = data$response_ddm,
                                     drift = factor(data$strength_bin, seq_len(n_drift)))
    }
    cat("start fitting: x0 =",start, ".\n") 
    
    fit <- try(nlminb(start, objective_function, 
                    rt = data$RT, response = data$response_ddm, 
                    drift = factor(data$strength_bin, seq_len(n_drift)), 
                    lower = lower), silent = T)
    if(inherits(fit, 'try-error')) {
      cat('bad starting valued. moving on.\n') 
      next
    }

    if(checkConverge) {
      if (fit$objective < best_fit$objective & fit$convergence==0) best_fit <- fit
      if(i == n_fits & best_fit$objective==1e+06) {
        best_fit <- fit
        cat('COULD NOT CONVERGE ON PARAMETERS FOR SUBJECT', data$id[1], '\n')
      }
    } else {
      if (fit$objective < best_fit$objective) best_fit <- fit
    }
    
    if(verbose) {
      converged = ifelse(fit$convergence==0, 'yes.', 'no.')
      newmin = ifelse(fit$objective < best_fit$objective, 'yes.', 'no.')
      cat('fit model', i, '/', n_fits, '|', 'converged?', converged, '|', 'new min?', newmin, '\n')
    }
  }
  out <- matrix(c(best_fit$par, best_fit$objective, best_fit$convergence), 1, length(best_fit$par)+2, byrow = T)
  colnames(out) <- c(names(fit$par), 'objective', 'convergence')
  
  if(savedir!=F) {
    best_fit$id = as.character(data$id[1])
    best_fit$age_group = as.character(data$age_group[1])
    best_fit$condition = as.character(data$condition[1])
    save(best_fit, file=paste0(savedir, 'ddmfit_', best_fit$age_group, '_',
                                  best_fit$condition, '_', best_fit$id, '.RData'))
  }
  
  return(out)
}