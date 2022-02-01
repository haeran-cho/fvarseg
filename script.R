devtools::install_github("https://github.com/Dom-Owens-UoB/fnets")

source('~/Documents/GitHub/fnets.segment/idio_seg.R')
source('~/Documents/GitHub/fnets.segment/common_seg.R')
source('~/Documents/GitHub/fnets.segment/misc.R')

n <- 2000
p <- 50
q <- 2
cp.common <- round(n * 1:2/3)
cp.idio <- round(n * 1:3/4)
d <- 1

ss <- sim.data2(n, p, q, cp.common, 1, 'ma', cp.idio, 1, d, seed = 1)

x <- ss$x

x <- ss$xi

cs <- common.seg(x, G.seq = NULL, thr = NULL, tt.by = round(log(dim(x)[2])), 
                 demean = TRUE, agg.over.freq = c('avg', 'max')[1], 
                 rule = c('eta', 'epsilon')[1], eta = .5, epsilon = .1, do.check = FALSE)
cs$est.cp

cs <- list(est.cp = matrix(cp.common, ncol = 1)) # if we know common = 0 with cp.common = c()

is <- idio.seg(x, common.seg.out = cs, G.seq = NULL, thr = NULL, d = d, demean = TRUE,
               cv.args = list(path.length = 10, n.folds = 1, do.cv = FALSE), 
               rule = c('eta', 'epsilon')[2], eta = .5, epsilon = .1)

is$est.cp  
  

