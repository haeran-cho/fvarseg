# favar.segment
Contains a method for high-dimensional time series segmentation under a piecewise stationary factor-adjusted vector autoregressive model. See 

> High-dimensional time series segmentation via factor-adjusted vector autoregressive modelling

by Haeran Cho, Idris Eckley, Paul Fearnhead and Hyeyoung Maeng [arXiv:](https://arxiv.org/abs/) for full details.

## Installation

To install `fvarseg` from GitHub:

```
devtools::install_github("https://github.com/haeran-cho/fvarseg")
```

## Usage

We can generate an example dataset used in the above paper for simulation studies as
```
out <- sim.data(n = 2000, p = 100, q = 2, d = 1,
  cp.common = 1:3/4, den.common = .5, type.common = 'ma', 
  cp.idio = c(3, 5)/8, seed = 123)
x <- out$x
````

Apply FVARseg with default settings.
```
fs <- fvar.seg(x)
```

Change points detected from the common component.
```
cs <- fs$common.out
cs$est.cp
```

Visualise the results.
```
par(mar = rep(2, 4), mfrow = c(2, 2))
for(rr in 1:length(cs$G.seq)){
  matplot(cs$est.cp.list[[rr]]$norm.stat, type = 'l', xlab = 'time', ylab = '', main = paste('G = ', cs$est.cp.list[[rr]]$G, sep = '')) # change point detector from each frequency
  lines(cs$est.cp.list[[rr]]$stat, col = 4, lwd = 2) # aggregation
  abline(v = sd$cp.common, col = 2, lty = 3) # true change points 
  abline(v = cs$est.cp.list[[rr]]$cp, col = 4, lty = 3) # change point estimators 
  abline(h = cs$est.cp.list[[rr]]$thr, col = 3) # threshold
}
```

Change points detected from the idiosyncratic component.
```
is <- fs$idio.out
is$est.cp  
```

Visualise the results.
```
par(mar = rep(2, 4), mfrow = c(2, 2))
for(rr in 1:length(is$G.seq)){
  plot(is$est.cp.list[[rr]]$stat, type = 'l', xlab = 'time', ylab = '', main = paste('G = ', is$est.cp.list[[rr]]$G, sep = '')) # change point detector 
  abline(v = sd$cp.idio, col = 2, lty = 3) # true change points
  abline(v = is$est.cp.list[[rr]]$cp, col = 4, lty = 2) # change point estimators
  abline(h = is$est.cp.list[[rr]]$thr, col = 4) # threshold
}
```

