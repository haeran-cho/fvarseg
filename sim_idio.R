source("~/Documents/GitHub/favar.segment/common_seg.R")
source("~/Documents/GitHub/favar.segment/idio_seg.R")
source("~/Documents/GitHub/favar.segment/misc.R")

n.seq <- 250 * c(4, 8)
p.seq <- 25 * c(2, 4, 6)
G.seq <- c(1/8, 1/5, 1/4)
q <- 2

d <- 1
KK <- 1

sim <- 200

for(nn in 1:length(n.seq)){
  n <- n.seq[nn]
  
  for(pp in 1:length(p.seq)){
    p <- p.seq[pp]
    
    idio.out <- 
      array(NA, dim = c(sim, length(G.seq), 3, 2, 3)) # (acv0, spec0, lambda) x (1, beta1) x (x1, x2, xi)
    
    for(ii in 1:sim){
      for(jj in 1:3){
        
        if(jj == 1){
          ss <- sim.data2(n, p, q,
                          cp.common = c(), den.common = 1, type.common = c('ma', 'ar')[1],
                          cp.idio = c(), size.idio = 1, d = d, do.scale = !FALSE, seed = ii)
          qq <- q
          xx <- ss$x * KK
        }  
        if(jj == 2){
          ss <- sim.data2(n, p, q,
                          cp.common = c(), den.common = 1, type.common = c('ma', 'ar')[2],
                          cp.idio = c(), size.idio = 1, d = d, do.scale = !FALSE, seed = ii)
          qq <- q
          xx <- ss$x * KK
        }
        if(jj == 3){
          ss <- sim.data2(n, p, 0,
                          cp.common = cp.common, den.common = 1, type.common = c('ma', 'ar')[2],
                          cp.idio = c(), size.idio = 1, d = d, do.scale = !FALSE, seed = ii)
          qq <- 0
          xx <- ss$x * KK
        }
        xx <- xx - apply(xx, 1, mean)
        
        for(gg in 1:3){
          G <- ceiling(n * G.seq[gg])
          
          ll <- max(1, floor(G^(1/3)))
          
          stat <- rep(0, n)
          
          pcfa <- post.cp.fa(xx, c(), qq, 5, max(1, 4 * floor((G/log(G))^(1/3))))
          Gamma_c <- pcfa$Gamma_c[,, 1:(ll + 1),, drop = FALSE]
          idx <- rep(1, n)
          
          int <- 1:G
          
          dpca <- fnets:::dyn.pca(xx[, int], q = qq, ic.op = 5)
          # spec0 <- Re(dpca$spec$Sigma_i[,, 1])
          # acv0 <- dpca$acv$Gamma_i[,, 1]
          acv <- dpca$acv$Gamma_x[,, 1:(ll + 1), drop = FALSE]
          dpca.l <- fnets:::dyn.pca(xx[, int[1:round(G/2)]], q = qq, ic.op = 5)
          dpca.r <- fnets:::dyn.pca(xx[, int[-(1:round(G/2))]], q = qq, ic.op = 5)
          
          spec0 <- dpca.l$spec$Sigma_i[,, 1:(ll + 1)] - dpca.r$spec$Sigma_i[,, 1:(ll + 1)]
          acv0 <- dpca.l$acv$Gamma_i[,, 1:(ll + 1)] - dpca.r$acv$Gamma_i[,, 1:(ll + 1)]
          
          ycv <- fnets:::yw.cv(xx[, int], method = 'ds', var.order = d, q = dpca$q, do.plot = !TRUE)
          lambda.path <- ycv$lambda.path
          mg <- fnets:::make.gg(dpca$acv$Gamma_i, d)
          ive <- fnets:::var.dantzig(mg$GG, mg$gg, ycv$lambda)
          beta <- ive$beta
          
          # lambda.max <- max(abs(acv0))
          # lambda.path <- round(exp(seq(log(lambda.max), log(lambda.max * .05), length.out = 10)), digits = 10)
          # icv <- idio.cv(xx = xx[, int, drop = FALSE], lambda.max = NULL, Gamma_c = Gamma_c, idx = idx[int], var.order = d, path.length = 10, n.folds = 1, do.plot = TRUE)  
          # acv <- acv.x(xx[, int, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1), drop = FALSE]
          # mg <- make.gg(acv - Gamma_c[,,, 1], d)
          # beta <- idio.beta(mg$GG, mg$gg, min(icv$lambda, lambda.path[5]))$beta
          
          # mv <- max(abs(ss$A.list[[1]])); par(mfrow = c(1, 2)); fields::imagePlot(beta, nlevel = 12, breaks = seq(-mv, mv, length.out = 13)); fields::imagePlot(t(ss$A.list[[1]]), nlevel = 12, breaks = seq(-mv, mv, length.out = 13))
          
          diff.Gamma_x <- acv - 
            acv.x(xx[, int + G, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1), drop = FALSE]
          mgd <- make.gg(diff.Gamma_x, d)
          stat[G] <- max(abs(mgd$GG %*% beta - mgd$gg)) 
          
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
            stat[tt] <- max(abs(mgd$GG %*% beta - mgd$gg))
            tt <- tt + 1
          }
          
          # ts.plot(stat/max(abs(spec0)))
          # ts.plot(stat/max(abs(acv0)))
          
          idio.out[ii, gg, 1, 1, jj] <- max(stat)/max(abs(spec0))
          idio.out[ii, gg, 2, 1, jj] <- max(stat)/max(abs(acv0))
          idio.out[ii, gg, 3, 1, jj] <- max(stat)/lambda.path[4]
          if(norm(beta, "1") > 0) idio.out[ii, gg, , 2, jj] <- idio.out[ii, gg, , 1, jj]/norm(beta, "1")
        }
      }
    }
    
    save(idio.out, file = paste('idio2_n', n, 'p', p, 'd', d, 'K', KK, '.RData', sep = ''))
    
  }
}


n.seq <- 250 * c(4, 8)
p.seq <- 25 * c(2, 4, 6)
G.seq <- c(1/8, 1/5, 1/4)

n <- n.seq[2]
p <- p.seq[1]

ll <- 1
d <- c(1, 2)[ll]
KK <- c(1, 100)[ll]

load(file = paste('~/downloads/sim/sim_lanc/idio_new_n', n, 'p', p, 'd', d, 'K', KK, '.RData', sep = ''))

dimnames(idio.out)[[5]] <- c('ma', 'ar', 'none')
dimnames(idio.out)[[4]] <- c('1', 'beta')
dimnames(idio.out)[[3]] <- c('spec', 'acv', 'none')
dimnames(idio.out)[[2]] <- round(n * G.seq)

out0 <- idio.out
out1 <- idio.out

jj <- 1 # c('spec', 'acv', 'lambda')
kk <- 1 # c('1', 'beta')

qu <- .95

# out0/out1
apply(out0[,, jj, kk,], c(2, 3), quantile, qu, TRUE) /
apply(out1[,, jj, kk,], c(2, 3), quantile, qu, TRUE)

# ma/ar
apply(out0[,, jj, kk, 1], 2, quantile, qu, TRUE)/apply(out0[,, jj, kk, 2], 2, quantile, qu, TRUE)
apply(out1[,, jj, kk, 1], 2, quantile, qu, TRUE)/apply(out1[,, jj, kk, 2], 2, quantile, qu, TRUE)

apply(out0[,, jj, kk,], c(2, 3), quantile, qu, TRUE)
apply(out1[,, jj, kk,], c(2, 3), quantile, qu, TRUE)

par(mfcol = c(2, 3))
for(ll in 1:3){
  boxplot(out0[,, jj, kk, ll], main = paste('idio0_n', n, 'p', p, sep = ''), ylim = c(1, 5))
  boxplot(out1[,, jj, kk, ll], main = paste('idio1_n', n, 'p', p, sep = ''), ylim = c(1, 5))
}

