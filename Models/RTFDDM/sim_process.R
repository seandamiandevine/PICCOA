# simulate one drift process with RTF model


# high/low a diffusion process ----------------------------------------------

pdf('plots/simprocess.pdf', width=14, height=4)

layout(matrix(1:2, ncol=2))

obj     = 40
sub     = 70
minstim = 1
maxstim = 100
minV    = -5
maxV    = 5
z       = 0.5*a
ter     = .2
s       = 1
t       = 1
td      = .001

v  = (obj-minstim)/(maxstim-minstim) * (maxV-minV) + minV
pz = (sub-minstim)/(maxstim-minstim)

# low a
set.seed(2022)

a  = 1
z  = pz*a
x = c(z)

tp = 1
while(1) {
  # accumulate evidence
  if (x[tp] > a || x[tp] < 0) break
  x[tp+1] = x[tp] + v*td + s*rnorm(1)*sqrt(td)
  
  # update 
  tp=tp+1
  }


x[x>a] = a; x[x<0]=0
plot(1:tp, x, type='l', xlab='Time (ms.)', ylab='', ylim=c(0,a), main='Low Decision Thresholds',
     xaxt='n', yaxt='n')
abline(h=z, lty='dashed')
text(tp*.95, z+.03, 'Biased start')
abline(h=c(0,a), col=c('purple','blue'), lwd=3)
abline(h=0.5*a, col='darkgrey', lty='dashed')
text(tp*.94, 0.5*a-.04, 'Unbiased start')
arrows(0, z-.1, 20, .2, col='red')
text(3.5, .4, 'Drift rate', col='red')
arrows(0, -.15, tp, -.15, xpd=T)

# high a
set.seed(2022)

a  = 3
z  = pz*a
x = c(z)

tp = 1
while(1) {
  # accumulate evidence
  if (x[tp] > a || x[tp] < 0) break
  x[tp+1] = x[tp] + v*td + s*rnorm(1)*sqrt(td)
  
  # update 
  tp=tp+1
}


x[x>a] = a; x[x<0]=0
plot(1:tp, x, type='l', xlab='Time (ms.)', ylab='', ylim=c(0,a), main='Higher Decision Thresholds', 
     xaxt='n', yaxt='n')
abline(h=c(0,a), col=c('purple','blue'), lwd=3)
abline(h=z, lty='dashed')
text(tp*.9, z+.1, 'Biased start')
abline(h = c(1.5,2.5), col=c(scales::alpha('purple', 0.25), scales::alpha('blue', 0.25)), lwd=3)
arrows(0, -.45, tp, -.45, xpd=T)

dev.off()



# pBlue over a ------------------------------------------------------------

as = seq(.5, 3, by=.5)
pblue = c()
seblue = c()
mrt = c()
sert = c()
for(a in as) {
  rts = c()
  resps = c()
  for(i in 1:100) {
    sub= 70 
    v  = (obj-minstim)/(maxstim-minstim) * (maxV-minV) + minV
    pz = (sub-minstim)/(maxstim-minstim)
    z  = pz*a
    x = c(z)
    tp = 1
    while(1) {
      # accumulate evidence
      if (x[tp] > a |x[tp] < 0) break
      x[tp+1] = x[tp] + v*td + s*rnorm(1)*sqrt(td)
      
      # update 
      tp=tp+1
      if(tp>10000) errorCondition('diffusion process ran too long!!')
    }

    rts[i] = tp/(1/td)+ter
    resps[i] = ifelse(x[tp]<=0, 0, 1)
  }
  pblue = c(pblue, mean(resps))
  seblue = c(seblue, plotrix::std.error(resps))
  
  mrt = c(mrt, mean(rts))
  sert = c(sert, plotrix::std.error(rts))
}

pdf('plots/simSlope.pdf', width=10, height=6)
layout(matrix(1:2, ncol=2))

plot(as, pblue, type='b', xlab = 'Decision Thresholds', ylab = 'p(Judge Dot as Blue)', pch=16, ylim=c(0,1), 
     main='Higher Thresholds = Fewer Errors')
arrows(as, pblue-seblue, as, pblue+seblue, length=0)

plot(as, mrt, type='b', xlab = 'Decision Thresholds', ylab = 'Avg. RT (s.)', pch=16, ylim=range(pretty(c(mrt-sert, mrt+sert))), 
     main = 'Higher Thresholds = Longer RT')
arrows(as, mrt-sert, as, mrt+sert, length=0)

dev.off()

