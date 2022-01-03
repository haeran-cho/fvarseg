source("~/Documents/GitHub/favar.segment/common_seg.R")
source("~/Documents/GitHub/favar.segment/idio_seg.R")
source("~/Documents/GitHub/favar.segment/misc.R")

n <- 2000
p <- 100
q <- 3
tc <- 2.5

cp.common <- round(n * c(1/3, 2/3))
cp.common <- c()

cp.idio <- round(n * (1:3)/4)
cp.idio <- c()

ss <- sim.data(n = n, p = p, q = q, 
              cp.common = cp.common, den.common = .75, type.common = c('ma', 'ar')[2], ma.order = 2,
              cp.idio = cp.idio, size.idio = 1, burnin = 100)
x <- ss$x
x <- ss$xi

dp <- common.spec.est(t(scale(t(x), scale = FALSE)), q = NULL, ic.op = 5, max(1, floor(200^(1/3))))
dp$hl$q.hat
dev.off()

G.seq <- c(200, 300, 400)

mean.x <- apply(x, 1, mean)
xx <- x - mean.x

##

common.est.cp.list <- list()
common.stat.list <- list()

thr <- 1.5

for(ii in 1:length(G.seq)){
  G <- G.seq[ii]
  ll <- max(1, floor(G^(1/3)))
  # thr <- thr.const * p * max(sqrt(ll * log(n)/G), 1/ll, 1/p)
  
  common.stat.list[[ii]] <- 
    cts <- common.two.step(xx, G, thr = thr, ll, ceiling(log(n)), 'm', 'avg')
  common.stat.list[[ii]]$G <- G
  common.stat.list[[ii]]$ll <- ll
  common.stat.list[[ii]]$thr <- thr
  
  est.cp <- common.search.cp(cts, thr = thr, G, .5, 'max')
  est.cp <- common.check(xx, G, est.cp, thr, ll, NULL, 5, 'm', 'avg')
  common.est.cp.list[[ii]] <- est.cp
  
  matplot(cts$norm.stat, type = 'l'); abline(v = cp.common, lty = 2, col = 2, lwd = 2); abline(v = cp.idio, lty = 3, col = 8); abline(v = est.cp, col = 4, lty = 3)
  abline(h = thr, col = 3); lines(cts$stat, col = 4, lwd = 2)
  
}

print(est.cp <- bottom.up(common.est.cp.list, G.seq, .5))

##

est.cp.common <- c()
est.cp.common <- est.cp[, 1]

K <- length(est.cp.common)
if(K > 0) est.cp.common <- sort(est.cp.common)
brks <- c(0, est.cp.common, n)
idx <- rep(c(1:(length(brks) - 1)), diff(brks))

ll <- max(1, floor(min(G.seq)^(1/3)))

pcfa <- post.cp.fa(x, est.cp.common, NULL, 5, ll)
pcfa$q.seq
Gamma_c <- pcfa$Gamma_c[,, 1:(ll + 1),, drop = FALSE]

idio.est.cp.list <- list()
idio.stat.list <- list()

for(ii in 1:length(G.seq)){
  
  G <- G.seq[ii]
  thr <- 2 # thr.const * max(sqrt(ll * log(n * p) / G), 1/ll, 1/sqrt(p))
  vv <- G
  stat <- rep(0, n)
  check.cp <- est.cp <- c()
  
  while(vv <= n - G){
    
    int <- (vv - G + 1):vv
    icv <- idio.cv(xx = xx[, int, drop = FALSE], Gamma_c = Gamma_c, idx = idx, var.order = d, 
                   path.length = 10, n.folds = 1)  
    tb <- tabulate(idx[int], nbins = K + 1)
    acv <- acv.x(xx[, int, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1), drop = FALSE]
    for(kk in 1:(K + 1)) acv <- acv - tb[kk] / G * Gamma_c[,,, kk]
    mg <- make.gg(acv, d)
    beta <- idio.beta(mg$GG, mg$gg, icv$lambda)$beta
    mean(beta != 0); image(beta, col = RColorBrewer::brewer.pal(11, 'RdBu'), breaks = seq(-max(abs(beta)), max(abs(beta)), length.out = 12))
    
    diff.Gamma_x <- acv.x(xx[, int, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1), drop = FALSE] -
      acv.x(xx[, int + G, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1), drop = FALSE]
    diff.Gamma_c <- diff.Gamma_x * 0
    if(K > 0){
      common.weights <- tabulate(idx[int], nbins = K + 1) - 
        tabulate(idx[int + G], nbins = K + 1)
      for(kk in 1:(K + 1)) diff.Gamma_c <- diff.Gamma_c + 
          common.weights[kk] / G * Gamma_c[,,, kk]
    }
    
    mg <- make.gg(diff.Gamma_x - diff.Gamma_c, d)
    null.norm <- max(abs(mg$GG %*% beta - mg$gg))
    stat[vv] <- 1
    
    first <- TRUE
    tt <- vv + 1
    check.theta <- tt.max <- n - G
    while(tt <= tt.max){
      for(h in 0:ll){
        diff.Gamma_x[,, h + 1] <- diff.Gamma_x[,, h + 1] -
          xx[, tt - G, drop = FALSE] %*% t(xx[, tt - G + h, drop = FALSE]) / G +
          xx[, tt - h, drop = FALSE] %*% t(xx[, tt, drop = FALSE]) / G +
          xx[, tt, drop = FALSE] %*% t(xx[, tt + h, drop = FALSE]) / G -
          xx[, tt + G - h, drop = FALSE] %*% t(xx[, tt + G, drop = FALSE]) / G
      }
      diff.Gamma_c <- diff.Gamma_x * 0
      if(K > 0){
        common.weights <- tabulate(idx[(tt - G + 1):tt], nbins = K + 1) -
          tabulate(idx[(tt + 1):(tt + G)], nbins = K + 1)
        for(kk in 1:(K + 1)) diff.Gamma_c <- diff.Gamma_c + common.weights[kk] / G * Gamma_c[,,, kk]
      }
      mg <- make.gg(diff.Gamma_x - diff.Gamma_c, d)
      stat[tt] <- max(abs(mg$GG %*% beta - mg$gg)) / null.norm
      
      if(first & stat[tt] > thr){
        check.theta <- tt
        tt.max <- min(tt.max, check.theta + G - 1)
        check.cp <- c(check.cp, check.theta)
        first <- FALSE
      }
      tt <- tt + 1
    }  
    ts.plot(stat); abline(h = thr, col = 3); abline(v = check.cp, col = 6); abline(v = cp.idio, col = 2, lty = 3)
    
    if(check.theta < tt.max){
      hat.theta <- (check.theta:tt.max)[which.max(stat[check.theta:tt.max])]
      est.cp <- c(est.cp, hat.theta)
      vv <- hat.theta + G
      
      abline(v = hat.theta, col = 4, lty = 2)
    } else break
    
  }
  
  if(length(est.cp) > 0){
    idio.est.cp.list[[ii]] <- est.cp
  } else idio.est.cp.list[[ii]] <- NA
  idio.stat.list[[ii]] <- list(stat = stat, G = G, thr = thr, check.cp = check.cp)
  
}

bottom.up(idio.est.cp.list, G.seq, .5)

################

thr.const <- .22
norm.type <- c('m', 'f', '2')[1]
agg.over.freq <- 'avg'
tt.by <- ceiling(log(n))

G <- 200

ll <- max(1, floor(G^(1/3)))
thr <- 1.5 # thr.const * p * max(sqrt(ll * log(n)/G), 1/ll, 1/p)

w <- bartlett.weights(((-ll):ll)/ll)
len <- 2 * ll
thetas <- 2 * pi * (0:len)/(len + 1)

null <- x1 <- array(0, dim = c(p, p, ll + 1))
for(h in 0:ll){
  ind.l <- (1 + h):G
  ind.r <- (n - G + 1 + h):n
  x1[,, h + 1] <- x[, ind.l - h] %*% t(x[, ind.l]) + x[, ind.r - h] %*% t(x[, ind.r])
  x1[,, h + 1] <- x1[,, h + 1]/(2 * G)
  for(ii in 1:(ll + 1)) {
    null[,, ii] <- null[,, ii] + 
      x1[,, h + 1] * w[h + ll + 1] * exp(-(0 + 1i) * h * thetas[ii]) +
      (h > 0) * t(x1[,, h + 1]) * w[h + ll + 1] * exp((0 + 1i) * h * thetas[ii])
  }
}
null <- null/(2 * pi)
null.norm <- apply(abs(null), 3, norm, type = norm.type)
if(norm.type %in% c('f', '2')) null.norm <- null.norm / sqrt(p)

tt.seq <- round(seq(G, n - G, by = tt.by))
stat0 <- common.stat0(x, G, thr, ll, tt.seq)

null.norm <- stat0[G, ]

norm.stat <- abs(t(t(stat0) / null.norm))
if(agg.over.freq == 'avg') stat <- apply(norm.stat, 1, mean)
if(agg.over.freq == 'max') stat <- apply(norm.stat, 1, max)

matplot(tt.seq, norm.stat[tt.seq, ], type = 'l'); abline(v = cp.common, lty = 2, col = 2); abline(h = thr, col = 3); lines(tt.seq, stat[tt.seq], col = 4, lwd = 2)

tt.list <- common.tt.list(stat, G, thr, tt.seq, tt.by)

if(length(tt.list) > 0){
  for(ii in 1:length(tt.list)){
    s <- min(tt.list[[ii]]); e <- max(tt.list[[ii]])
    stat0 <- common.stat1(x, G, thr, ll, s, e, stat0)
  }
  norm.stat <- abs(t(t(stat0) / null.norm))
  if(norm.type == 'f') norm.stat <- norm.stat * sqrt(p)
  
  if(agg.over.freq == 'avg') stat <- apply(norm.stat, 1, mean)
  if(agg.over.freq == 'max') stat <- apply(norm.stat, 1, max)
} 

matplot(norm.stat, type = 'l'); abline(v = cp.common, lty = 2, col = 2); abline(h = thr, col = 3); lines(stat, col = 4, lwd = 2)

cts <- list(norm.stat = norm.stat, stat = stat, null.norm = null.norm)
est.cp <- common.search.cp(cts, thr, G, eta = .5, rule = c('max', 'over')[1])
est.cp

common.check(x, G, est.cp, thr, ll, NULL, ic.op = 5, norm.type, agg.over.freq, cts$null.norm)

##

est.cp.common <- cp.common

d <- 1
K <- length(est.cp.common)
brks <- c(0, est.cp.common, n)
idx <- rep(c(1:(length(brks) - 1)), diff(brks))

mean.x <- apply(x, 1, mean)
xx <- x - mean.x
ll <- max(1, floor(min(G.seq)^(1/3)))

pcfa <- post.cp.factor.analysis(xx, est.cp.common, q = NULL, ic.op = 5, ll = ll)
pcfa$q.seq
Gamma_c <- pcfa$Gamma_c[,, 1:(ll + 1),, drop = FALSE]

thr <- 1.5 # 2.2 * max(sqrt(ll * log(n * p) / G), 1/ll, 1/sqrt(p))
vv <- G
stat <- rep(0, n)
check.cp <- est.cp <- c()
  
while(vv <= n - G){
  
  int <- (vv - G + 1):vv
  icv <- idio.cv(xx = xx[, int, drop = FALSE], Gamma_c = Gamma_c, idx = idx, var.order = d, 
                 path.length = 10, n.folds = 1)  
  tb <- tabulate(idx[int], nbins = K + 1)
  acv <- acv.x(xx[, int, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1)]
  for(kk in 1:(K + 1)) acv <- acv - tb[kk] / G * Gamma_c[,,, kk]
  mg <- make.gg(acv, d)
  beta <- idio.beta(mg$GG, mg$gg, icv$lambda)$beta
  mean(beta != 0); fields::imagePlot(beta, col = RColorBrewer::brewer.pal(11, 'RdBu'), breaks = seq(-max(abs(beta)), max(abs(beta)), length.out = 12))
  
  # null.norm <- max(abs(acv[,, 1]))
  # norm(acv[,, 1], 'f')/sqrt(p)
  # norm(acv[,, 1], 'm')
  
  diff.Gamma_x <- acv.x(xx[, int, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1)] -
    acv.x(xx[, int + G, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1)]
  diff.Gamma_c <- diff.Gamma_x * 0
  if(K > 0){
    common.weights <- tabulate(idx[int], nbins = K + 1) - 
      tabulate(idx[int + G], nbins = K + 1)
    for(kk in 1:(K + 1)) diff.Gamma_c <- diff.Gamma_c + 
        common.weights[kk] / G * Gamma_c[,,, kk]
  }
  
  mg <- make.gg(diff.Gamma_x - diff.Gamma_c, d)
  null.norm <- max(abs(mg$GG %*% beta - mg$gg))
  stat[vv] <- 1
  
  first <- TRUE
  tt <- vv + 1
  check.theta <- tt.max <- n - G
  while(tt <= tt.max){
    for(h in 0:ll){
      diff.Gamma_x[,, h + 1] <- diff.Gamma_x[,, h + 1] -
        xx[, tt - G, drop = FALSE] %*% t(xx[, tt - G + h, drop = FALSE]) / G +
        xx[, tt - h, drop = FALSE] %*% t(xx[, tt, drop = FALSE]) / G +
        xx[, tt, drop = FALSE] %*% t(xx[, tt + h, drop = FALSE]) / G -
        xx[, tt + G - h, drop = FALSE] %*% t(xx[, tt + G, drop = FALSE]) / G
    }
    diff.Gamma_c <- diff.Gamma_x * 0
    if(K > 0){
      common.weights <- tabulate(idx[(tt - G + 1):tt], nbins = K + 1) -
        tabulate(idx[(tt + 1):(tt + G)], nbins = K + 1)
      for(kk in 1:(K + 1)) diff.Gamma_c <- diff.Gamma_c + common.weights[kk] / G * Gamma_c[,,, kk]
    }
    mg <- make.gg(diff.Gamma_x - diff.Gamma_c, d)
    stat[tt] <- max(abs(mg$GG %*% beta - mg$gg)) / null.norm
    
    if(first & stat[tt] > thr){
      check.theta <- tt
      tt.max <- min(tt.max, check.theta + G - 1)
      check.cp <- c(check.cp, check.theta)
      first <- FALSE
    }
    tt <- tt + 1
    
  }  
  
  ts.plot(stat); abline(h = thr, col = 3); abline(v = cp.idio, col = 2, lty = 3)
  
  if(check.theta < tt.max){
    hat.theta <- (check.theta:tt.max)[which.max(stat[check.theta:tt.max])]
    est.cp <- c(est.cp, hat.theta)
    vv <- hat.theta + G
  } else break
  
}
