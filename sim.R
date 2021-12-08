source("~/Dropbox/v/half baked/var/cp/code/haeran/favar.segment-main/common_seg.R")
source("~/Dropbox/v/half baked/var/cp/code/haeran/favar.segment-main/idio_seg.R")
source("~/Dropbox/v/half baked/var/cp/code/haeran/favar.segment-main/misc.R")

n.seq <- 250 * c(1, 2, 3, 5, 8, 13)
p.seq <- 100 * c(1, 2, 3, 5, 8)
G.seq <- c(1/8, 1/5, 1/4)
q <- 3
d <- 1

sim <- 200

common.out <- idio.out <- 
  array(0, dim = c(sim, length(n.seq), length(p.seq), length(G.seq), 2, 2)) # (max scaled, max unscaled) x (x, xi)

for(nn in 1:length(n.seq)){
  n <- n.seq[nn]
  
  for(pp in 1:length(p.seq)){
    p <- p.seq[pp]
    
    for(ii in 1:sim){
      ss <- sim.data(n, p, q,  
                     cp.common = c(), den.common = 0, type.common = c('ma', 'ar')[1], ma.order = 0,
                     cp.idio = c(), size.idio = 0, burnin = 100, seed = ii)
      for(jj in 1:2){
        
        if(jj == 1) xx <- ss$x - apply(ss$x, 1, mean)
        if(jj == 2) xx <- ss$xi - apply(ss$xi, 1, mean)
        
        for(gg in 1:3){
          G <- round(n.seq[nn] * G.seq[gg])
          
          ll <- max(1, floor(G^(1/3)))
          common.thr <- p * max(sqrt(ll * log(n)/G), 1/ll, 1/p)
          cts <- common.two.step(xx, G, 0, ll, 1, 'm', 'avg')
          common.out[ii, nn, pp, gg, 1, jj] <- max(cts$stat)/common.thr
          common.out[ii, nn, pp, gg, 2, jj] <- max(cts$stat)
          
          ##
          
          stat <- rep(0, n)
          
          idio.thr <- max(sqrt(ll * log(n * p) / G), 1/ll, 1/sqrt(p))
          
          pcfa <- post.cp.factor.analysis(xx, c(), q, 5, ll)
          Gamma_c <- pcfa$Gamma_c[,, 1:(ll + 1), 1]
          idx <- rep(1, n)
          
          int <- 1:G
          icv <- idio.cv(xx = xx[, int, drop = FALSE], Gamma_c = Gamma_c, idx = idx, var.order = d, 
                         path.length = 10, n.folds = 1)  
          acv <- acv.x(xx[, int, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1), drop = FALSE]
          mg <- make.gg(acv - Gamma_c, d)
          beta <- idio.beta(mg$GG, mg$gg, icv$lambda)$beta
          
          diff.Gamma_x <- acv - 
            acv.x(xx[, int + G, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1), drop = FALSE]
          mg <- make.gg(diff.Gamma_x, d)
          null.norm <- max(abs(mg$GG %*% beta - mg$gg))
          stat[G] <- 1
          tt <- G + 1
          while(t <= n - G){
            for(h in 0:ll){
              diff.Gamma_x[,, h + 1] <- diff.Gamma_x[,, h + 1] -
                xx[, tt - G, drop = FALSE] %*% t(xx[, tt - G + h, drop = FALSE]) / G +
                xx[, tt - h, drop = FALSE] %*% t(xx[, tt, drop = FALSE]) / G +
                xx[, tt, drop = FALSE] %*% t(xx[, tt + h, drop = FALSE]) / G -
                xx[, tt + G - h, drop = FALSE] %*% t(xx[, tt + G, drop = FALSE]) / G
            }
            mg <- make.gg(diff.Gamma_x, d)
            stat[tt] <- max(abs(mg$GG %*% beta - mg$gg)) / null.norm
            tt <- tt + 1
          }  
          
          idio.out[ii, nn, pp, gg, 1, jj] <- max(stat)/idio.thr
          idio.out[ii, nn, pp, gg, 2, jj] <- max(stat)
        }
      }
    }
  }
}
