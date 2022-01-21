source("~/Documents/GitHub/favar.segment/common_seg.R")
source("~/Documents/GitHub/favar.segment/idio_seg.R")
source("~/Documents/GitHub/favar.segment/misc.R")

n.seq <- 250 * c(4, 8, 12, 2)
p.seq <- 25 * c(1, 2, 4, 8)
G.seq <- c(1/8, 1/5, 1/4)
q <- 3
d <- 1

sim <- 200

for(nn in 1){
  n <- n.seq[nn]
  
  for(pp in 1:length(p.seq)){
    p <- p.seq[pp]
    
    common.out <- idio.out <- 
      array(0, dim = c(sim, length(G.seq), 2, 3)) # (max scaled, max unscaled) x (x1, x2, xi)
    
    for(ii in 1:sim){
      for(jj in 1:3){
        
        if(jj == 1){
          ss <- sim.data(n, p, q, cp.common = c(), den.common = 0, 
                         type.common = c('ma', 'ar')[1], ma.order = 2,
                         cp.idio = c(), size.idio = 0, seed = ii)
          qq <- q
        }  
        if(jj == 2){
          ss <- sim.data(n, p, q, cp.common = c(), den.common = 0, 
                         type.common = c('ma', 'ar')[2], ma.order = 0,
                         cp.idio = c(), size.idio = 0, seed = ii)
          qq <- q
        }
        if(jj == 3){
          ss <- sim.data(n, p, 0, cp.common = c(), den.common = 0, 
                         type.common = c('ma', 'ar')[2], ma.order = 0,
                         cp.idio = c(), size.idio = 0, seed = ii)
          qq <- 0
        }
        xx <- ss$x - apply(ss$x, 1, mean)
        
        for(gg in 1:3){
          G <- ceiling(n * G.seq[gg])
          
          ll <- max(1, floor(G^(1/3)))
          common.thr <- p * max(sqrt(ll * log(n)/G), 1/ll, 1/p)
          cts <- common.two.step(xx, G, Inf, ll, 2, 'avg')
          common.out[ii, gg, 1, jj] <- max(cts$stat)/common.thr
          common.out[ii, gg, 2, jj] <- max(cts$stat)
          
          ##
          
          stat <- rep(0, n)
          
          idio.thr <- max(sqrt(ll * log(n * p) / G), 1/ll, 1/sqrt(p))
          
          pcfa <- post.cp.fa(xx, c(), qq, 5, max(1, 4 * floor((n/log(n))^(1/3))))
          Gamma_c <- pcfa$Gamma_c[,, 1:(ll + 1),, drop = FALSE]
          idx <- rep(1, n)
          
          int <- 1:G
          icv <- idio.cv(xx = xx[, int, drop = FALSE], Gamma_c = Gamma_c, idx = idx[int], var.order = d, 
                         path.length = 10, n.folds = 1, do.plot = !TRUE)  
          acv <- acv.x(xx[, int, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1), drop = FALSE]
          mg <- make.gg(acv - Gamma_c[,,, 1], d)
          beta <- idio.beta(mg$GG, mg$gg, icv$lambda)$beta
          null.norm <- max(abs(mg$GG %*% beta - mg$gg))
          
          diff.Gamma_x <- acv - 
            acv.x(xx[, int + G, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1), drop = FALSE]
          mgd <- make.gg(diff.Gamma_x, d)
          stat[G] <- max(abs(mgd$GG %*% beta - mgd$gg)) / null.norm 
          
          tt <- G + 1
          while(tt <= n - G){
            for(h in 0:ll){
              diff.Gamma_x[,, h + 1] <- diff.Gamma_x[,, h + 1] -
                xx[, tt - G, drop = FALSE] %*% t(xx[, tt - G + h, drop = FALSE]) / G +
                xx[, tt - h, drop = FALSE] %*% t(xx[, tt, drop = FALSE]) / G +
                xx[, tt, drop = FALSE] %*% t(xx[, tt + h, drop = FALSE]) / G -
                xx[, tt + G - h, drop = FALSE] %*% t(xx[, tt + G, drop = FALSE]) / G
            }
            mgd <- make.gg(diff.Gamma_x, d)
            stat[tt] <- max(abs(mgd$GG %*% beta - mgd$gg)) / null.norm
            tt <- tt + 1
          }  
          
          idio.out[ii, gg, 1, jj] <- max(stat)/idio.thr
          idio.out[ii, gg, 2, jj] <- max(stat)
        }
      }
    }
    
    save(common.out, file = paste('common_n', n, 'p', p, '.RData', sep = ''))
    save(idio.out, file = paste('idio_n', n, 'p', p, '.RData', sep = ''))
    
  }
}


n.seq <- 250 * c(1, 2, 4, 8, 12)
p.seq <- 25 * c(1, 2, 4, 8)
G.seq <- c(1/8, 1/5, 1/4)

n <- n.seq[3]
p <- p.seq[3]

load(file = paste('~/downloads/sim/sim_lanc/common_n', n, 'p', p, '.RData', sep = ''))
load(file = paste('~/downloads/sim/sim_lanc/idio_n', n, 'p', p, '.RData', sep = ''))

dimnames(common.out)[[4]] <- dimnames(idio.out)[[4]] <- c('ma', 'ar', 'none')
dimnames(common.out)[[3]] <- dimnames(idio.out)[[3]] <- c('scale', 'no')
dimnames(common.out)[[2]] <- dimnames(idio.out)[[2]] <- round(n * G.seq)

kk <- 2

signif(apply(common.out[,, kk,], c(2, 3), quantile, .9), 4)
signif(apply(idio.out[,, kk,], c(2, 3), quantile, .9), 4)

jj <- 2 # ma, ar, none
kk <- 2

par(mfcol = c(2, 3))
for(jj in 1:3){
  boxplot(common.out[,, kk, jj], main = paste('common_n', n, 'p', p, sep = ''), ylim = c(1, 2))
  boxplot(idio.out[,, kk, jj], main = paste('idio_n', n, 'p', p, sep = ''), ylim = c(1, 5))
}
