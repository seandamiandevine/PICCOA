dir = 'Models/seq/out/ethics/'
files = list.files(dir)
fit = data.frame(stringsAsFactors = F)
bic = function(ll, n, k) -2 * ll + log(n) * k # calculates BIC

for(f in 1:length(files)) {
  
  thisfit = mget(load(paste0(dir, files[f])))
  id = as.numeric(gsub('.RData', '' , strsplit(files[f], '_')[[1]][5]))
  age_group =  strsplit(files[f], '_')[[1]][3]
  condition =  strsplit(files[f], '_')[[1]][4]
  bestbic = bic(-thisfit$bestVal, 800, 6)
  thisfitdf = data.frame(id, 
                         age_group, 
                         condition, 
                         thisfit$bestParams[1],
                         thisfit$bestParams[2],
                         thisfit$bestParams[3],
                         thisfit$bestParams[4],
                         thisfit$bestParams[5],
                         thisfit$bestParams[6],
                         thisfit$bestVal, 
                         bestbic, 
                         thisfit$bestConv, stringsAsFactors = F)
  
  fit = rbind(fit, thisfitdf)
  
}
colnames(fit) = c('id', 'age_group', 'condition', 'B0', 'Bf', 'BF', 'Bc', 'lF', 'lc', 'loglik', 'BIC', 'convergence')

#dots.fit <- fit
ethics.fit <- fit

# write.csv(dots.fit, 'Models/seq/dots.fit.csv')
# write.csv(ethics.fit, 'Models/seq/ethics.fit.csv')
