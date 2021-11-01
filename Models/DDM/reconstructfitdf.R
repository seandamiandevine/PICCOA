reconstructfitdf <- function(inpath, outpath=NULL, outfile = 'ddmfit.csv') { 
  
  files = list.files(inpath)
  fit = data.frame(stringsAsFactors = F)
  bic = function(ll, n, k) -2 * ll + log(n) * k # calculates BIC
  
  for(f in 1:length(files)) {
    
    thisfit = mget(load(paste0(inpath, files[f])))
    id = as.numeric(gsub('.RData', '' , strsplit(files[f], '_')[[1]][4]))
    age_group =  strsplit(files[f], '_')[[1]][2]
    condition =  strsplit(files[f], '_')[[1]][3]
    #bestbic = bic(-thisfit$best_fit$objective, 800, length(thisfit$best_fit$par))
    thisfitdf = data.frame(id, 
                           age_group, 
                           condition, 
                           thisfit$best_fit$par[1],
                           thisfit$best_fit$par[2],
                           thisfit$best_fit$par[3],
                           thisfit$best_fit$par[4],
                           thisfit$best_fit$par[5],
                           thisfit$best_fit$par[6],
                           thisfit$best_fit$par[7],
                           thisfit$best_fit$par[8],
                           thisfit$best_fit$par[9],
                           thisfit$best_fit$par[10],
                           thisfit$best_fit$objective, 
                           thisfit$best_fit$convergence, 
                           stringsAsFactors = F)
    
    fit = rbind(fit, thisfitdf)
    
  }
  colnames(fit) = c('id', 'age_group', 'condition', names(thisfit$best_fit$par), 'loglik', 'convergence')
  rownames(fit) <- NULL
  
  if(is.null(outpath)==F) write.csv(fit, paste0(outpath, outfile))

  return(fit)
}

dots.ddm  <- reconstructfitdf(inpath='Models/DDM/out/exp1/', outpath='Models/DDM/out/', outfile = 'ddmfit1.csv')
#dots3.ddm <- reconstructfitdf(inpath='Models/DDM/out/exp3/', outpath='Models/DDM/out/', outfile = 'ddmfit3.csv')
