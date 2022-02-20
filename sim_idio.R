source("~/Documents/GitHub/fnets.segment/common_seg.R")
source("~/Documents/GitHub/fnets.segment/idio_seg.R")
source("~/Documents/GitHub/fnets.segment/misc.R")

n.seq <- 250 * c(4, 8, 16)
p.seq <- 25 * c(2, 4, 6)
q <- 2

d <- 2
KK <- 100

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
        
        if(jj <= 2) G.seq <- round(seq(2.5 * p, n / min(4, n / (3 * p)), length.out = 4))
        if(jj == 3) G.seq <- round(seq(2 * p, n / min(5, n / (2.5 * p)), length.out = 4))
        
        for(gg in 1:4){
          G <- G.seq[gg]
          
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
          
          lambda.max <- max(abs(xx[, int] %*% t(xx[, int])/G))
          lambda.path <- round(exp(seq(log(lambda.max), log(lambda.max * .0001), length.out = 10)), digits = 10)
          # ycv <- fnets:::yw.cv(xx[, int], method = 'ds', var.order = d, q = dpca$q, do.plot = !TRUE)
          # lambda.path <- ycv$lambda.path
          mg <- fnets:::make.gg(dpca$acv$Gamma_i, d)
          ive <- fnets:::var.dantzig(mg$GG, mg$gg, lambda.path[4])
          # ive <- fnets:::var.dantzig(mg$GG, mg$gg, ycv$lambda)
          beta <- ive$beta
          
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


n.seq <- 250 * c(4, 8, 16)
p.seq <- 25 * c(2, 4, 6)

# model

n <- n.seq[2]
p <- p.seq[2]
G.seq <- 1:4
for(ll in 1:2){
  d <- c(1, 2)[ll]
  KK <- c(1, 100)[ll]
  load(file = paste('~/downloads/sim/sim_lanc/idio2_n', n, 'p', p, 'd', d, 'K', KK, '.RData', sep = ''))

  dimnames(idio.out)[[5]] <- c('ma', 'ar', 'none')
  dimnames(idio.out)[[4]] <- c('1', 'beta')
  dimnames(idio.out)[[3]] <- c('spec', 'acv', 'none')
  dimnames(idio.out)[[2]] <- G.seq

  if(ll == 1) out0 <- idio.out
  if(ll == 2) out1 <- idio.out
}

# sample size

ll <- 1; d <- c(1, 2)[ll]; KK <- c(1, 100)[ll]
p <- p.seq[2]

for(nn in c(1, 3)){
  n <- n.seq[nn]
  load(file = paste('~/downloads/sim/sim_lanc/idio2_n', n, 'p', p, 'd', d, 'K', KK, '.RData', sep = ''))
  
  dimnames(idio.out)[[5]] <- c('ma', 'ar', 'none')
  dimnames(idio.out)[[4]] <- c('1', 'beta')
  dimnames(idio.out)[[3]] <- c('spec', 'acv', 'none')
  dimnames(idio.out)[[2]] <- 1:4
  
  if(nn == 1) out0 <- idio.out
  if(nn == 3) out1 <- idio.out
}

# dimensionality

ll <- 1; d <- c(1, 2)[ll]; KK <- c(1, 100)[ll]
n <- n.seq[2]

for(pp in c(1, 3)){
  p <- p.seq[pp]
  load(file = paste('~/downloads/sim/sim_lanc/idio2_n', n, 'p', p, 'd', d, 'K', KK, '.RData', sep = ''))
  
  dimnames(idio.out)[[5]] <- c('ma', 'ar', 'none')
  dimnames(idio.out)[[4]] <- c('1', 'beta')
  dimnames(idio.out)[[3]] <- c('spec', 'acv', 'none')
  dimnames(idio.out)[[2]] <- 1:4
  
  if(pp == 1) out0 <- idio.out
  if(pp == 3) out1 <- idio.out
}

######################

jj <- 2 # c('spec', 'acv', 'lambda')
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
  boxplot(out0[,, jj, kk, ll], main = paste('idio0_n', n, 'p', p, sep = ''), ylim = c(0, 3))
  boxplot(out1[,, jj, kk, ll], main = paste('idio1_n', n, 'p', p, sep = ''), ylim = c(0, 3))
}

load('idio_fit.RData')

n <- 2000
p <- 100
G <- round(seq(2.5 * p, n / min(4, n / (3 * p)), length = 4))[3]
exp(predict(idio.fit.list[[3]], list(n = n, p = p, G = G)))
exp(predict(fit, list(n = n, p = p, G = G)))

######################

n.seq <- 250 * c(4, 8, 16)
p.seq <- 25 * c(2, 4, 6)

kk <- 1

qu <- c(.8, .9, .95)
yy <- xx <- c()
for(nn in 1:3){
  n <- n.seq[nn]
  for(pp in 1:3){
    p <- p.seq[pp]
    if(kk == 1) G.seq <- sort(round(seq(2 * p, n / min(4, n/(p * 2.5)), length.out = 4)))
    if(kk == 2) G.seq <- sort(round(seq(2.5 * p, n / min(5, n/(p * 3)), length.out = 4)))
    
    for(gg in 1:4){
      G <- G.seq[gg]
      xx <- rbind(xx, c(n, p, G, d))
    }
    
    out <- c()
    for(ll in 1:2){
      d <- c(1, 2)[ll]
      KK <- c(1, 100)[ll]  
      # load(file = paste('~/downloads/sim/sim_lanc/idio2_n', n, 'p', p, 'd', d, 'K', KK, '.RData', sep = ''))
      load(file = paste('archive/idio2_n', n, 'p', p, 'd', d, 'K', KK, '.RData', sep = ''))
      if(kk == 1) out <- rbind(out, rbind(idio.out[,, 2, 1, 1], idio.out[,, 2, 1, 2]))
      if(kk == 2) out <- rbind(out, idio.out[,, 2, 1, 3])
    }
                   
    tmp <- c()
    for(jj in 1:3) tmp <- cbind(tmp, apply(out, 2, quantile, qu[jj], TRUE))
    yy <- rbind(yy, tmp)
  }
}

df <- data.frame(n = xx[, 1], p = xx[, 2], G = xx[, 3], y80 = yy[, 1], y90 = yy[, 2], y95 = yy[, 3])

y <- df$y90
fit <- lm(log(y) ~ 0 + I(log(log(n))) + I(log(log(p))) + I(log(G)), data = df)
fit <- lm(y ~ I((log(n)/G)^(1/3)) + I(1/p), data = df)

summary(fit)
plot(exp(fitted(fit)), y); abline(a = 0, b = 1, col = 1)

idio.fit.list <- list()
for(jj in 1:3){
  if(jj == 1) y <- df$y80
  if(jj == 2) y <- df$y90
  if(jj == 3) y <- df$y95
  fit <- lm(log(y) ~ 0 + I(log(log(n))) + I(log(log(p))) + I(log(G)), data = df)  
  idio.fit.list[[jj]] <- fit
  
  # if(jj == 1){
  #   plot(exp(fitted(fit)), y); abline(a = 0, b = 1, col = 1)
  # } else points(exp(fitted(fit)), y, col = jj)
}

save(idio.fit.list, file = 'idio_fit.RData')

idio.fit.list0 <- idio.fit.list
save(idio.fit.list0, file = 'idio_fit0.RData')

n <- 1000
G <- 50
p <- 150
exp(predict(fit, list(n = n, p = p, G = G)))


######################

library(quantreg)

n.seq <- 250 * c(4, 8, 16)
p.seq <- 25 * c(2, 4, 6)

kk <- 1

yx <- c()
for(nn in 1:3){
  n <- n.seq[nn]
  for(pp in 1:3){
    p <- p.seq[pp]
    if(kk == 1) G.seq <- sort(round(seq(2 * p, n / min(4, n/(p * 2.5)), length.out = 4)))
    if(kk == 2) G.seq <- sort(round(seq(2.5 * p, n / min(5, n/(p * 3)), length.out = 4)))
    
    for(ll in 1:2){
      d <- c(1, 2)[ll]
      KK <- c(1, 100)[ll]  
      load(file = paste('archive/idio2_n', n, 'p', p, 'd', d, 'K', KK, '.RData', sep = ''))
      for(gg in 1:4){
        G <- G.seq[gg]
        if(kk == 1){
          tmp <- idio.out[, gg, 2, 1, 1]
          tmp <- c(tmp, idio.out[, gg, 2, 1, 2])
        } else tmp <- idio.out[, gg, 2, 1, 3]
        tmp <- cbind(tmp, n)
        tmp <- cbind(tmp, p)
        tmp <- cbind(tmp, G)
        
        yx <- rbind(yx, tmp)
      }
    }
  }
}

df <- data.frame(y = yx[, 1], n = yx[, 2], p = yx[, 3], G = yx[, 4])

fit <- rq(log(y) ~ 0 + I(log(log(n))) + I(log(log(p))) + I(log(G)), tau = .95, data = df)
summary(fit)

plot(exp(fitted(fit)), df$y); abline(a = 0, b = 1, col = 2)

rho <- function(u, tau = .5) u*(tau - (u < 0))
1 - sum(rho(fit$resid, fit$tau))/ sum(rho(log(df$y), fit$tau))

idio.fit.list <- list()
qu <- c(.8, .9, .95)
for(jj in 1:3){
  fit <- rq(log(y) ~ 0 + I(log(log(n))) + I(log(log(p))) + I(log(G)), tau = qu[jj], data = df)
  idio.fit.list[[jj]] <- fit
}

save(idio.fit.list, file = 'new_idio_fit.RData')

idio.fit.list0 <- idio.fit.list
save(idio.fit.list0, file = 'new_idio_fit0.RData')

n <- 2000
p <- 100
G.seq <- sort(round(seq(2.5 * p, n / min(4, n/(p * 3)), length.out = 4)))
G.seq <- sort(round(seq(2 * p, n / min(5, n/(p * 2.5)), length.out = 4)))

for(gg in 1:4) print(c(exp(predict(fit, list(n = n, p = p, G = G.seq[gg]))),
                       exp(predict(idio.fit.list[[3]], list(n = n, p = p, G = G.seq[gg])))))


load(file = paste('archive/idio2_n', n, 'p', p, 'd', 1, 'K', 1, '.RData', sep = ''))
tmp <- rbind(idio.out[,, 2, 1, 1], idio.out[,, 2, 1, 2])
load(file = paste('archive/idio2_n', n, 'p', p, 'd', 2, 'K', 100, '.RData', sep = ''))
tmp <- rbind(tmp, rbind(idio.out[,, 2, 1, 1], idio.out[,, 2, 1, 2]))
apply(tmp, 2, quantile, .95)

load(file = paste('archive/idio2_n', n, 'p', p, 'd', 1, 'K', 1, '.RData', sep = ''))
tmp <-idio.out[,, 2, 1, 3]
load(file = paste('archive/idio2_n', n, 'p', p, 'd', 2, 'K', 100, '.RData', sep = ''))
tmp <- rbind(tmp, idio.out[,, 2, 1, 3])
apply(tmp, 2, quantile, .95)


