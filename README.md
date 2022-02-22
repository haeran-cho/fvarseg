# favar.segment
Change point detection in factor-adjusted VAR processes

Generate a piecewise stationary factor-adjusted VAR process.
```
sd <- sim.data(n = 2000, p = 100, q = 2, 
  cp.common = round(n * 1:3/4), den.common = .5, type.common = 'ma', 
  cp.idio = round(n * c(3, 5)/8), d = 1, seed = 123)
x <- sd$x
````

Apply FVARseg with default settings.
```
fs <- fvar.seg(x)
```

Change points detected from the common component.
```
cs <- fs$common.out
cs$est.cp

par(mar = rep(1, 4), mfrow = c(2, 2))
for(rr in 1:length(cs$G.seq)){
  matplot(cs$est.cp.list[[rr]]$norm.stat, type = 'l'); lines(cs$est.cp.list[[rr]]$stat, col = 4, lwd = 2)
  abline(v = sd$cp.common, col = 2, lty = 3); abline(v = cs$est.cp.list[[rr]]$cp, col = 4, lty = 3); abline(h = cs$est.cp.list[[rr]]$thr, col = 3)
}
```

Change points detected from the idiosyncratic component.
```
is <- fs$idio.out
is$est.cp  

par(mar = rep(1, 4), mfrow = c(2, 2))
for(rr in 1:length(is$G.seq)){
  ts.plot(is$est.cp.list[[rr]]$stat); abline(h = is$est.cp.list[[rr]]$thr, col = 4)
  abline(v = sd$cp.idio, col = 2, lty = 3); abline(v = is$est.cp.list[[rr]]$cp, col = 4, lty = 2); abline(v = is$est.cp.list[[rr]]$check.cp, col = 6, lty = 2)
}
```

