devtools::install_github("https://github.com/Dom-Owens-UoB/fnets")

library(VARDetect)

setwd("~/Documents/GitHub/fnets.segment")
source('idio_seg.R')
source('common_seg.R')
source('misc.R')

n <- 2000
p <- 100
q <- 0
d <- 2

cp.common <- c()
cp.common <- round(n * 1:3/4)

cp.idio <- c()
cp.idio <- round(n * 1:3/4)

ss <- sim.data2(n, p, q,
                cp.common = cp.common, den.common = .7, type.common = c('ma', 'ar')[2],
                cp.idio = cp.idio, size.idio = .8, 
                d = d, 
                do.scale = !FALSE, seed = NULL)
x <- ss$x

x <- ss$xi

# pcf0 <- post.cp.fa((ss$x - 0 * ss$xi)[, 1:(2*cs$G.seq[2])], cs$G.seq[2], q = q, ic.op = 5, cs$ll.seq[1])
# pcf <- post.cp.fa((ss$x - 0 * ss$xi), cp.common, q = q, ic.op = 5, cs$ll.seq[1])
# norm(pcf$Sigma_c[,, 2, 1] - pcf$Sigma_c[,, 2, 2], '2') / norm(pcf0$Sigma_c[,, 2, 1] - pcf0$Sigma_c[,, 2, 2], '2')
# norm(pcf$Sigma_c[,, 2, 2] - pcf$Sigma_c[,, 2, 3], '2') / norm(pcf0$Sigma_c[,, 2, 1] - pcf0$Sigma_c[,, 2, 2], '2')

# fnets:::dyn.pca(x)$hl$q.hat

# if common != 0
cs <- common.seg(x, G.seq = NULL, thr = NULL, tt.by = round(log(dim(x)[2]) * 2), 
                 demean = TRUE, agg.over.freq = c('avg', 'max')[1], 
                 rule = c('eta', 'epsilon')[1], eta = .5, epsilon = .1, do.check = FALSE, do.plot = !FALSE)
cs$est.cp

# cs0 <- common.seg(ss$x - ss$xi, G.seq = NULL, thr = NULL, tt.by = round(log(dim(x)[2])^2), 
#                  demean = TRUE, agg.over.freq = c('avg', 'max')[1], 
#                  rule = c('eta', 'epsilon')[2], eta = .5, epsilon = .1, do.check = FALSE, do.plot = !FALSE)
# cs0$est.cp
# cs0$est.cp.list[[1]]$null.norm
# 

# if we know common = 0, to compare against VARDetect
cs <- list(est.cp = c(), ll.seq = max(1, floor(min((n/10)^(1/3), n/(2 * log(n))))))

is <- idio.seg(x, common.seg.out = cs, 
               G.seq = round(seq(2 * p, n / min(5, n/(2 * p)), length.out = 4)), 
               thr = rep(1, 4), d = d, demean = TRUE,
               cv.args = list(path.length = 10, n.folds = 1, do.cv = FALSE, do.plot = !FALSE), 
               rule = c('eta', 'epsilon')[1], eta = .5, epsilon = 5 / (2.5 * p))

is <- idio.seg(x, common.seg.out = cs, d = d, demean = TRUE,
               cv.args = list(path.length = 10, n.folds = 1, do.cv = FALSE, do.plot = !FALSE), 
               rule = c('eta', 'epsilon')[1], eta = .5, epsilon = 5 / (2.5 * p))

is$est.cp  
par(mfrow = c(2, 2))
for(rr in 1:4){
  ts.plot(is$est.cp.list[[rr]]$norm.stat); abline(h = is$est.cp.list[[rr]]$thr, col = 4)
  abline(v = cp.idio, col = 2, lty = 3); abline(v = is$est.cp.list[[rr]]$cp, col = 4, lty = 2)
}


block <- VARDetect::tbss(t(ss$xi), method = 'sparse', q = d, use.BIC = FALSE)
block$cp

