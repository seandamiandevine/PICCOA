
# Dots --------------------------------------------------------------------

#---- Model comparison
dotmod = read.csv('dots_fit.csv')
modcomp = by(dotmod, list(dotmod$id, dotmod$condition, dotmod$model), function(x)
             c(id=as.character(x$id[1]), condition=as.character(x$condition[1]), 
               model=as.character(x$model[1]), LL = as.numeric(x$LL[1]), bic=as.numeric(x$BIC[1])))
modcomp = data.frame(do.call(rbind, modcomp), stringsAsFactors = F)

modcomp$LL = as.numeric(modcomp$LL)
modcomp$bic = as.numeric(modcomp$bic)
modcomp$condition = factor(modcomp$condition, levels = c('Stable', 'Decreasing'))
modcomp$age_group = factor(ifelse(as.numeric(modcomp$id)<100, 'Old', 'Young'), levels = c('Young', 'Old'))

# LL
pdf('plots/dots_LLComp.pdf', width = 10, height = 6)
layout(matrix(1:2, ncol=2))
mLL = tapply(modcomp$LL, list(modcomp$model, modcomp$condition, modcomp$age_group), mean)
seLL = tapply(modcomp$LL, list(modcomp$model, modcomp$condition, modcomp$age_group), plotrix::std.error)

# YA
thisbar = barplot(mLL[,,1], beside=T, ylim=range(pretty(c(mLL[,,1]-seLL[,,1], mLL[,,1]+seLL[,,1]))), xpd=F, 
                  ylab = 'Avg. Log-Likelihood', xlab='Condition', legend.text = T, 
                  args.legend = list(x='bottomleft', title='Model', bty='n'), 
                  main='Young Adults')
axis(1, at=c(mean(thisbar[,1]),mean(thisbar[,2])), labels = c('',''))
arrows(thisbar, mLL[,,1]-seLL[,,1], thisbar, mLL[,,1]+seLL[,,1], length = 0)

# OA
thisbar = barplot(mLL[,,2], beside=T, ylim=range(pretty(c(mLL[,,2]-seLL[,,2], mLL[,,2]+seLL[,,2]))), xpd=F, 
                  ylab = 'Avg. Log-Likelihood', xlab='Condition', legend.text = T, 
                  args.legend = list(x='bottomleft', title='Model', bty='n'), 
                  main='Old Adults')
axis(1, at=c(mean(thisbar[,1]),mean(thisbar[,2])), labels = c('',''))
arrows(thisbar, mLL[,,2]-seLL[,,2], thisbar, mLL[,,2]+seLL[,,2], length = 0)

dev.off()

# Sum BIC
pdf('plots/dots_BICcomp.pdf', width=10, height=6)
layout(matrix(1:2, ncol=2))
sBIC =tapply(modcomp$bic, list(modcomp$model, modcomp$condition, modcomp$age_group), sum)

# YA
thisbar = barplot(sBIC[,,1], beside=T, ylim=range(pretty(sBIC[,,1])), xpd=F, 
                  ylab = 'Sum BIC', xlab='Condition', legend.text = T, 
                  args.legend = list(x='topleft', title='Model', bty='n'), 
                  yaxt='n', main='Young Adults')
axis(1, at=c(mean(thisbar[,1]),mean(thisbar[,2])), labels = c('',''))
axis(2, at=range(sBIC[,,1]), label=c('Better Fit','Worse Fit'))

# OA
thisbar = barplot(sBIC[,,2], beside=T, ylim=range(pretty(sBIC[,,2])), xpd=F, 
                  ylab = 'Sum BIC', xlab='Condition', legend.text = T, 
                  args.legend = list(x='topleft', title='Model', bty='n'), 
                  yaxt='n', main='Old Adults')
axis(1, at=c(mean(thisbar[,1]),mean(thisbar[,2])), labels = c('',''))
axis(2, at=range(sBIC[,,2]), label=c('Better Fit','Worse Fit'))

dev.off()

# Within-subject comparison
diff = data.frame()
for(id in unique(modcomp$id)){
  tmp = modcomp[modcomp$id==id,]
  moddiff = tmp$bic[tmp$model=='Ctl'] - tmp$bic[tmp$model=='RTF']
  diff = rbind(diff, data.frame(id=id, condition=tmp$condition[1], 
                                age_group=tmp$age_group[1], diff=moddiff))
}
diff=diff[order(diff$diff), ]

pdf('plots/dots_diff.pdf', width=10, height=6)

# YA
tmp = diff[diff$age_group=='Young',]
barplot(tmp$diff[tmp$condition=='Decreasing'], col='blue', ylab = 'BIC Difference', yaxt='n', 
        main='Young Adults')
barplot(tmp$diff[tmp$condition=='Stable'], col='red', add=T, yaxt='n')
axis(2, at=range(tmp$diff),labels=c('Control Better', 'RTF Better'))
legend('topleft', fill=c('red', 'blue'), bty='n', legend = c('Stable', 'Decreasing'),title='Model')
    
# OA
tmp = diff[diff$age_group=='Old',]
barplot(tmp$diff[tmp$condition=='Decreasing'], col='blue', ylab = 'BIC Difference', yaxt='n', 
        main='Old Adults')
barplot(tmp$diff[tmp$condition=='Stable'], col='red', add=T, yaxt='n')
axis(2, at=range(tmp$diff),labels=c('Control Better', 'RTF Better'))
legend('topleft', fill=c('red', 'blue'), bty='n', legend = c('Stable', 'Decreasing'),title='Model')  

dev.off()

#---- Parameter comparison
pdf('plots/ParamComp.pdf', width=10,height=6)

dotmod$age_group = factor(ifelse(as.numeric(dotmod$id)<100, 'Old', 'Young'), levels = c('Young', 'Old'))
contrasts(dotmod$age_group) = contr.sum(2)
contrasts(dotmod$condition) = contr.sum(2)

# tau 
tau = dotmod[dotmod$parsname=='tau0',]
mtau = tapply(tau$pars, list(tau$age_group, tau$condition), mean)
setau = tapply(tau$pars,list(tau$age_group, tau$condition), plotrix::std.error)
thisbar = barplot(mtau, beside=T, ylim=range(pretty(c(mtau-setau, mtau+setau+5))), xpd=F,
                  ylab = expression(tau), legend.text = T, 
                  args.legend = list(title='Age Groups', bty='n'))
arrows(thisbar, mtau-setau, thisbar, mtau+setau, length=0)
taulm = summary.lm(lm(pars~age_group*condition, data=tau))
b = round(coef(taulm)['age_group1','Estimate'], 2)
p = round(coef(taulm)['age_group1','Pr(>|t|)'], 4)
legend('topleft', bty='n', legend=paste0('b=',b,', p=', p))

# sigma 
sigma = dotmod[dotmod$parsname=='sigma_rtf',]
msigma = tapply(sigma$pars, list(sigma$age_group, sigma$condition), mean)
sesigma = tapply(sigma$pars,list(sigma$age_group, sigma$condition), plotrix::std.error)
thisbar = barplot(msigma, beside=T, ylim=range(pretty(c(msigma-sesigma, msigma+sesigma+5))), xpd=F,
                  ylab = expression(sigma), legend.text = T, 
                  args.legend = list(title='Age Groups', bty='n'))
arrows(thisbar, msigma-sesigma, thisbar, msigma+sesigma, length=0)
sigmalm = summary.lm(lm(pars~age_group*condition, data=sigma))
b = round(coef(sigmalm)['age_group1','Estimate'], 2)
p = round(coef(sigmalm)['age_group1','Pr(>|t|)'], 4)
legend('topleft', bty='n', legend=paste0('b=',b,', p=', p))

# nk
tau = dotmod[dotmod$parsname=='nk',]
mtau = tapply(tau$pars, list(tau$age_group, tau$condition), mean)
setau = tapply(tau$pars,list(tau$age_group, tau$condition), plotrix::std.error)
thisbar = barplot(mtau, beside=T, ylim=range(pretty(c(mtau-setau, mtau+setau+5))), xpd=F,
                  ylab = expression(nk), legend.text = T, 
                  args.legend = list(title='Age Groups', bty='n'))
arrows(thisbar, mtau-setau, thisbar, mtau+setau, length=0)
taulm = summary.lm(lm(pars~age_group*condition, data=tau))
b = round(coef(taulm)['age_group1','Estimate'], 2)
p = round(coef(taulm)['age_group1','Pr(>|t|)'], 4)
legend('topleft', bty='n', legend=paste0('b=',b,', p=', p))

dev.off()