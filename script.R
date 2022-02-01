devtools::install_github("https://github.com/Dom-Owens-UoB/fnets")

setwd("~/Documents/GitHub/fnets.segment")
source('idio_seg.R')
source('common_seg.R')
source('misc.R')

n <- 2000
p <- 100
q <- 2
cp.common <- round(n * 1:2/3)
cp.idio <- round(n * 1:3/4)
d <- 2

ss <- sim.data2(n, p, q, cp.common, 1, 'ma', cp.idio, 1, d, seed = 111)

x <- ss$x

x <- ss$xi

# if common != 0
cs <- common.seg(x, G.seq = NULL, thr = NULL, tt.by = round(log(dim(x)[2])), 
                 demean = TRUE, agg.over.freq = c('avg', 'max')[1], 
                 rule = c('eta', 'epsilon')[1], eta = .5, epsilon = .1, do.check = FALSE, do.plot = FALSE)
cs$est.cp

# if we know common = 0, to compare against VARDetect
cs <- list(est.cp = c(), 
           ll.seq = max(1, floor(min((n/10)^(1/3), n/(2 * log(n)))))) 


is <- idio.seg(x, common.seg.out = cs, G.seq = NULL, thr = NULL, d = d, demean = TRUE,
               cv.args = list(path.length = 10, n.folds = 1, do.cv = FALSE), 
               rule = c('eta', 'epsilon')[2], eta = .5, epsilon = .1)
is$est.cp  

rr <- 3
ts.plot(is$est.cp.list[[rr]]$norm.stat); abline(h = is$est.cp.list[[rr]]$thr, col = 4)
abline(v = cp.idio, col = 2, lty = 3); abline(v = is$est.cp.list[[rr]]$cp, col = 4, lty = 2)

