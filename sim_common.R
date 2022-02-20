source("~/Documents/GitHub/fnets.segment/common_seg.R")
source("~/Documents/GitHub/fnets.segment/idio_seg.R")
source("~/Documents/GitHub/fnets.segment/misc.R")

n.seq <- 250 * c(4, 8, 16)
p.seq <- 25 * c(2, 4, 6)
d <- 1

sim <- 100

for(nn in 2){
  n <- n.seq[nn]
  
  for(pp in 1:3){
    p <- p.seq[pp]
    
    for(ll in 1:2){
      
      q <- c(2, 3)[ll]
      KK <- c(1, 100)[ll]
      
      common.out <- 
        array(NA, dim = c(sim, 4, 3)) # (x1, x2, xi)
      
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
          
          G.seq <- sort(round(seq(2 * p, n / min(5, n/(p * 2.5)), length.out = 4)))
          
          for(gg in 1:length(G.seq)){
            G <- G.seq[gg]
            cts <- common.two.step(xx, G, thr = Inf, ll, ceiling(log(n)), 'avg')
            common.out[ii, gg, jj] <- max(cts$stat)
          }
        }
        save(common.out, file = paste('common2_n', n, 'p', p, 'q', q, 'K', KK, '.RData', sep = ''))
      }
    }
  }
}


n.seq <- 250 * c(4, 8, 16)
p.seq <- 25 * c(2, 4, 6)

# model

n <- n.seq[1]
p <- p.seq[1]
print(G.seq <- sort(round(seq(2 * p, n / min(5, n/(p * 2.5)), length.out = 4))))
for(ll in 1:2){
  q <- c(2, 3)[ll]
  KK <- c(1, 100)[ll]
  load(file = paste('~/downloads/sim/sim_lanc/common1_n', n, 'p', p, 'q', q, 'K', KK, '.RData', sep = ''))
  
  dimnames(common.out)[[3]] <- c('ma', 'ar', 'none')
  dimnames(common.out)[[2]] <- G.seq
  
  if(ll == 1) out0 <- common.out
  if(ll == 2) out1 <- common.out
}

# model

n <- n.seq[1]
p <- p.seq[1]
print(G.seq <- round(n * c(1/10, 1/8, 1/6, 1/4)))
ll <- 1
q <- c(2, 3)[ll]
KK <- c(1, 100)[ll]

for(ll in 1:2){
  if(ll == 1) load(file = paste('~/downloads/sim/sim_lanc/common1_n', n, 'p', p, 'q', q, 'K', KK, '.RData', sep = ''))
  if(ll == 2) load(file = paste('~/downloads/sim/sim_lanc/common2_n', n, 'p', p, 'q', q, 'K', KK, '.RData', sep = ''))
  dimnames(common.out)[[3]] <- c('ma', 'ar', 'none')
  dimnames(common.out)[[2]] <- G.seq
  
  if(ll == 1) out0 <- common.out
  if(ll == 2) out1 <- common.out
}

######################

qu <- .9

# out0/out1
apply(out0, c(2, 3), quantile, qu, TRUE) /
  apply(out1, c(2, 3), quantile, qu, TRUE)

# ma/ar
apply(out0[,, 1], 2, quantile, qu, TRUE)/apply(out0[,, 2], 2, quantile, qu, TRUE)
apply(out1[,, 1], 2, quantile, qu, TRUE)/apply(out1[,, 2], 2, quantile, qu, TRUE)

apply(out0, c(2, 3), quantile, qu, TRUE)
apply(out1, c(2, 3), quantile, qu, TRUE)

par(mfcol = c(2, 3))
for(ll in 1:3){
  boxplot(out0[,, ll], main = paste('common0_n', n, 'p', p, sep = ''), ylim = c(1, 4))
  boxplot(out1[,, ll], main = paste('common1_n', n, 'p', p, sep = ''), ylim = c(1, 4))
}

######################

n.seq <- 250 * c(4, 8, 16)
p.seq <- 25 * c(2, 4, 6)

qu <- c(.8, .9, .95)
yy <- xx <- c()
for(nn in 1:3){
  n <- n.seq[nn]
  for(pp in 1:3){
    p <- p.seq[pp]
    G.seq <- round(n * c(1/10, 1/8, 1/6, 1/4))
    
    for(ll in 1:2){
      q <- c(2, 3)[ll]
      KK <- c(1, 100)[ll]
      for(gg in 1:4){
        G <- G.seq[gg]
        xx <- rbind(xx, c(n, p, G, q))
      }
      
      load(file = paste('~/downloads/sim/sim_lanc/common1_n', n, 'p', p, 'q', q, 'K', KK, '.RData', sep = ''))
      out <- rbind(common.out[,, 1], common.out[,, 2], common.out[,, 3])
      load(file = paste('~/downloads/sim/sim_lanc/common2_n', n, 'p', p, 'q', q, 'K', KK, '.RData', sep = ''))
      out <- rbind(out, rbind(common.out[,, 1], common.out[,, 2], common.out[,, 3]))
      tmp <- c()
      for(jj in 1:3) tmp <- cbind(tmp, apply(out, 2, quantile, qu[jj], TRUE))
      yy <- rbind(yy, tmp)
    }
  }
}

df <- data.frame(n = xx[, 1], p = xx[, 2], G = xx[, 3], y80 = yy[, 1], y90 = yy[, 2], y95 = yy[, 3])

plot(xx[, 3], yy[, 1])

y <- df$y90
fit <- lm(log(y) ~ 0 + I(log(log(n))) + I(log(G)), data = df)
fit <- lm(log(y) ~ I(log(log(n)) - log(G)) + I(log(p)), data = df)
fit <- lm(y ~ I((log(n)/G)^(2/3)) + I(1/p), data = df)

summary(fit)

plot(exp(fitted(fit)), y); abline(a = 0, b = 1, col = 1)
plot(fitted(fit), df$y90); abline(a = 0, b = 1, col = 1)

common.fit.list <- list()
for(jj in 1:3){
  if(jj == 1) y <- df$y80
  if(jj == 2) y <- df$y90
  if(jj == 3) y <- df$y95
  fit <- lm(log(y) ~ 0 + I(log(log(n))) + I(log(G)), data = df)
  common.fit.list[[jj]] <- fit
  
  if(jj == 1){
    plot(exp(fitted(fit)), y); abline(a = 0, b = 1, col = 1)
  } else points(exp(fitted(fit)), y, col = jj)
}

save(common.fit.list, file = 'common_fit.RData')

n <- 4000
G <- 500
p <- 150
exp(predict(fit, list(n = n, p = p, G = G)))

######################

library(quantreg)

n.seq <- 250 * c(4, 8, 16)
p.seq <- 25 * c(2, 4, 6)

yx <- c()
for(nn in 1:3){
  n <- n.seq[nn]
  for(pp in 1:3){
    p <- p.seq[pp]
    G.seq <- round(n * c(1/10, 1/8, 1/6, 1/4))
    
    for(ll in 1:2){
      q <- c(2, 3)[ll]
      KK <- c(1, 100)[ll]  
      load(file = paste('archive/common1_n', n, 'p', p, 'q', q, 'K', KK, '.RData', sep = ''))
      for(gg in 1:4){
        G <- G.seq[gg]
        tmp <- common.out[, gg, 1]
        tmp <- c(tmp, common.out[, gg, 2])
        
        tmp <- cbind(tmp, n)
        tmp <- cbind(tmp, p)
        tmp <- cbind(tmp, G)
        
        yx <- rbind(yx, tmp)
      }
      
      load(file = paste('archive/common2_n', n, 'p', p, 'q', q, 'K', KK, '.RData', sep = ''))
      for(gg in 1:4){
        G <- G.seq[gg]
        tmp <- common.out[, gg, 1]
        tmp <- c(tmp, common.out[, gg, 2])
        
        tmp <- cbind(tmp, n)
        tmp <- cbind(tmp, p)
        tmp <- cbind(tmp, G)
        
        yx <- rbind(yx, tmp)
      }
    }
  }
}

df <- data.frame(y = yx[, 1], n = yx[, 2], p = yx[, 3], G = yx[, 4])

fit <- rq(log(y) ~ 0 + I(log(log(n))) + I(log(G)), tau = .9, data = df)
summary(fit)

rho <- function(u, tau = .5) u*(tau - (u < 0))
1 - sum(rho(fit$resid, fit$tau))/ sum(rho(log(df$y), fit$tau))

common.fit.list <- list()
qu <- c(.8, .9, .95)
for(jj in 1:3){
  fit <- rq(log(y) ~ 0 + I(log(log(n))) + I(log(G)), tau = qu[jj], data = df)
  common.fit.list[[jj]] <- fit
}

save(common.fit.list, file = 'new_common_fit.RData')

n <- 2000
p <- 100
G.seq <- round(n * c(1/10, 1/8, 1/6, 1/4))

exp(predict(fit, list(n = n, p = p, G = G.seq[2])))
exp(predict(common.fit.list[[2]], list(n = n, p = p, G = G.seq[2])))

load(file = paste('archive/common1_n', n, 'p', p, 'q', q, 'K', KK, '.RData', sep = ''))
tmp <- rbind(common.out[,, 1], common.out[,, 2])
load(file = paste('archive/common2_n', n, 'p', p, 'q', q, 'K', KK, '.RData', sep = ''))
tmp <- rbind(tmp, rbind(common.out[,, 1], common.out[,, 2]))
apply(tmp, 2, quantile, .9)


