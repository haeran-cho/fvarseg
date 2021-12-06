source("~/Dropbox/v/half baked/var/cp/code/hy/gdfm_functions.R")

source("~/Dropbox/v/half baked/var/station/code/haeran/fnets.R")
source("~/Dropbox/v/half baked/var/station/code/haeran/common.R")
source("~/Dropbox/v/half baked/var/station/code/haeran/idio.R")
source("~/Dropbox/v/half baked/var/station/code/haeran/omega.R")

n <- 2000
p <- 100
q <- 3
tc <- 2.5

cp.common <- round(n * c(1/3, 2/3))
cp.idio <- round(n * (1:3)/4)

x <- sim.data(n = n, p = p, q = q, 
              cp.common = cp.common, den.common = .5, type.common = c('ma', 'ar')[2], ma.order = 2,
              cp.idio = cp.idio, size.idio = 1, burnin = 100, seed = 1511)
dp <- common.spec.est(t(scale(t(x), scale = FALSE)), q = NULL, ic.op = 5, max(1, floor(200^(1/3))))
dp$hl$q.hat
dev.off()

##

thr.const <- .22
norm.type <- c('m', 'f', '2')[2]
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

est.cp <- common.search.cp(list(norm.stat = norm.stat, stat = stat), thr, G, eta = .5)
matplot(norm.stat, type = 'l'); abline(v = cp.common, lty = 2, col = 2, lwd = 2); abline(v = cp.idio, lty = 3, col = 6); abline(v = est.cp, col = 4, lty = 3); abline(h = thr, col = 3); lines(stat, col = 4, lwd = 2)

est.cp <- common.check(x, G, est.cp, thr, ll, q, ic.op, norm.type, agg.over.freq)

matplot(cts$norm.stat, type = 'l'); abline(v = cp.common, lty = 2, col = 2, lwd = 2); abline(v = cp.idio, lty = 3, col = 6); abline(v = est.cp, col = 4, lty = 3); abline(h = thr, col = 3); lines(cts$stat, col = 4, lwd = 2)

