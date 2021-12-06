source("~/Dropbox/v/half baked/var/cp/code/hy/gdfm_functions.R")

n <- 2000
p <- 100
q <- 3
tc <- 2.5

cp.common <- round(n * c(1/3, 2/3))
cp.idio <- round(n * (1:3)/4)

x <- sim.data(n = n, p = p, q = q, 
              cp.common = cp.common, den.common = .75, type.common = c('ma', 'ar')[2], ma.order = 2,
              cp.idio = cp.idio, size.idio = 1, burnin = 100, seed = 12334)
dp <- common.spec.est(t(scale(t(x), scale = FALSE)), q = NULL, ic.op = 5, max(1, floor(200^(1/3))))
dp$hl$q.hat
dev.off()

##

thr.const <- .22
norm.type <- c('m', 'f', '2')[1]
agg.over.freq <- 'avg'
tt.by <- ceiling(log(n))

G <- 200

ll <- max(1, floor(G^(1/3)))
thr <- thr.const * p * max(sqrt(ll * log(n)/G), 1/ll, 1/p)

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
norm.stat <- abs(t(t(stat0) / null.norm))
if(agg.over.freq == 'avg') stat <- apply(norm.stat, 1, mean)
if(agg.over.freq == 'max') stat <- apply(norm.stat, 1, max)

# matplot(norm.stat, type = 'l'); abline(v = cp.common, lty = 2, col = 2); abline(h = thr, col = 3); lines(stat, col = 4, lwd = 2)

tt.list <- common.tt.list(norm.stat, G, thr, tt.seq, tt.by)

if(length(tt.list) > 0){
  for(ii in 1:length(tt.list)){
    s <- min(tt.list[[ii]]); e <- max(tt.list[[ii]])
    stat0 <- common.stat1(x, G, thr, ll, s, e, stat0)
  }
  norm.stat <- abs(t(t(stat0)/apply(abs(null), 3, norm, type = norm.type)))
  if(norm.type == 'f') norm.stat <- norm.stat * sqrt(p)
  
  if(agg.over.freq == 'avg') stat <- apply(norm.stat, 1, mean)
  if(agg.over.freq == 'max') stat <- apply(norm.stat, 1, max)
} 

# matplot(norm.stat, type = 'l'); abline(v = cp.common, lty = 2, col = 2); abline(h = thr, col = 3); lines(stat, col = 4, lwd = 2)

cts <- list(norm.stat = norm.stat, stat = stat, null.stat = null.stat)
est.cp <- common.search.cp(cts, thr, G, eta = .05, rule = c('max', 'over')[2])
est.cp

common.check(x, G, est.cp, thr, ll, NULL, ic.op = 5, norm.type, agg.over.freq, cts$null.stat)

##

est.cp.common <- cp.common

K <- length(est.cp.common)
brks <- c(0, est.cp.common, n)
idx <- rep(c(1:(length(brks) - 1)), diff(brks))

mean.x <- apply(x, 1, mean)
xx <- x - mean.x
ll <- max(1, floor(min(G.seq)^(1/3)))

pcfa <- post.cp.factor.analysis(xx, est.cp.common, q, ic.op, ll)
Gamma_c <- pcfa$Gamma_c[,, 1:(ll + 1), ]

thr <- thr.const * max(sqrt(ll * log(n * p) / G), 1/ll, 1/sqrt(p))
vv <- G
stat <- rep(0, n)
check.cp <- est.cp <- c()
  
while(vv <= n - G){
  
  int <- (vv - G + 1):vv
  icv <- idio.cv(zz = xx[, int, drop = FALSE], Gamma_c = Gamma_c, idx = idx, var.order = d, 
                 path.length = path.length, n.folds = n.folds)  
  tb <- tabulate(idx[int], nbins = K + 1)
  acv <- acv.x(xx[, int, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1)]
  for(kk in 1:(K + 1)) acv <- acv - tb[kk] / G * Gamma_c[,,, kk]
  mg <- make.gg(acv, d)
  beta <- idio.beta(mg$GG, mg$gg, icv$lambda)$beta
  # mean(beta != 0); image(beta)
  
  diff.Gamma_x <- acv.x(xx[, int, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1)] -
    acv.x(xx[, int + G, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1)]
  diff.Gamma_c <- diff.Gamma_x * 0
  if(K > 0){
    common.weights <- tabulate(idx[int], nbins = K + 1) - 
      tabulate(idx[int + G], nbins = K + 1)
    for(kk in 1:(K + 1)) diff.Gamma_c <- diff.Gamma_c + 
        common.weights[kk] / G * Gamma_c[,,, kk]
  }
  
  tt <- vv
  check.theta <- tt.max <- n - G - 1
  while(tt <= tt.max){
    mg <- make.gg(diff.Gamma_x - diff.Gamma_c, d)
    stat[tt] <- max(abs(mg$GG %*% beta - mg$gg))
    
    if(stat[tt] > thr){
      check.theta <- tt; tt.max <- min(tt.max, check.theta + G - 1)
      check.cp <- c(check.cp, check.theta)
    }
    
    tt <- tt + 1
    
    for(h in 0:ll){
      diff.Gamma_x[,, h + 1] <- diff.Gamma_x[,, h + 1] * G -
        xx[, tt - G, drop = FALSE] %*% t(xx[, tt - G + h, drop = FALSE]) +
        xx[, tt - h, drop = FALSE] %*% t(xx[, tt, drop = FALSE]) +
        xx[, tt, drop = FALSE] %*% t(xx[, tt + h, drop = FALSE]) -
        xx[, tt + G - h, drop = FALSE] %*% t(xx[, tt + G, drop = FALSE])
    }
    diff.Gamma_x <- diff.Gamma_x / G
    diff.Gamma_c <- diff.Gamma_x * 0
    if(K > 0){
      common.weights <- tabulate(idx[(tt - G + 1):tt], nbins = K + 1) - 
        tabulate(idx[(tt + 1):(tt + G)], nbins = K + 1)
      for(kk in 1:(K + 1)) diff.Gamma_c <- diff.Gamma_c + common.weights[kk] / G * Gamma_c[,,, kk]
    }
  }  
  
  ts.plot(stat); abline(h = thr, col = 3); abline(v = cp.idio, col = 2, lty = 3)
  
  # do something with stat
  
  hat.theta <- (check.theta:tt.max)[which.max(stat[check.theta:tt.max])]
  est.cp <- c(est.cp, hat.theta)
  vv <- hat.theta + G
  
}
