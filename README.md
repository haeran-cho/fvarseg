# favar.segment
Change point detection in factor-adjusted VAR processes

sd <- sim.data(n = 2000, p = 100, q = 0, type.common = 'ma', den.common = .5, d = 1)
x <- sd$x

fs <- fvar.seg(x)

cs <- fs$common.out
cs$est.cp

par(mar = rep(1, 4), mfrow = c(2, 2))
for(rr in 1:length(cs$G.seq)){
  matplot(cs$est.cp.list[[rr]]$norm.stat, type = 'l'); lines(cs$est.cp.list[[rr]]$stat, col = 4, lwd = 2)
  abline(v = sd$cp.common, col = 2, lty = 3); abline(v = cs$est.cp.list[[rr]]$cp, col = 4, lty = 3); abline(h = cs$est.cp.list[[rr]]$thr, col = 3)
}

is <- fs$idio.out

is$est.cp  
par(mar = rep(1, 4), mfrow = c(2, 2))
for(rr in 1:length(is$G.seq)){
  ts.plot(is$est.cp.list[[rr]]$stat); abline(h = is$est.cp.list[[rr]]$thr, col = 4)
  abline(v = sd$cp.idio, col = 2, lty = 3); abline(v = is$est.cp.list[[rr]]$cp, col = 4, lty = 2)
}

is$G.seq
is$thr

