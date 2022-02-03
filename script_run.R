devtools::install_github("https://github.com/Dom-Owens-UoB/fnets")

setwd("~/Documents/GitHub/fnets.segment")
#setwd("~/gdfm_new")

source('idio_seg.R')
source('common_seg.R')
source('misc.R')

N <- 10
n <- 2000
d <- 1
q <- 2

##### models 
sim.list = expand.grid(structure(list( cp.idio = list(c(NULL), c(round(n * 1:2/3))), p = c(50, 100, 150), cp.common= list(c(round(n * 1:2/3)), c(round(n * 1:3/4))), 
                                       type.common = c("ma", "ar"))))
cp.idio = vector("list", 24); for(i in c(seq(2, 6, by=2), seq(14, 18, by=2))){ cp.idio[[i]] <- c(round(n * 1:2/3)); cp.idio[[i+6]] <- c(round(n * c(3,5)/8))}
sim.list$cp.idio = cp.idio
sim.list

##### simulations
for(S in 1:dim(sim.list)[1]){
  
  cp.idio = c(sim.list[S,][1], recursive=T, use.names=F)
  p = c(sim.list[S,][2], recursive=T, use.names=F)
  cp.common = c(sim.list[S,][3], recursive=T, use.names=F)  
  type.common = unlist(sim.list[S,][4], use.names = F)
  
  est.cp.com <- vector("list", length = N)
  est.cp.idio <- vector("list", length = N)
  
  for(i in 1:N){
    
    ### 1) data generation 
    ss <- sim.data2(n, p, q,
                    cp.common = cp.common, den.common = .5, type.common = type.common,
                    cp.idio = cp.idio, size.idio = 1, d = d, do.scale = !FALSE, seed = i)
    x <- ss$x
    #x <- ss$xi
    
    ### 2) common component 
    cs <- common.seg(x, G.seq = NULL, thr = NULL, tt.by = round(log(dim(x)[2])^2), 
                     demean = TRUE, agg.over.freq = c('avg', 'max')[1], 
                     rule = c('eta', 'epsilon')[1], eta = .5, epsilon = .1, do.check = FALSE, do.plot = !FALSE)
    est.cp.com[[i]] <- cs$est.cp[, 1]
  
    ### 3) idiosyncratic component 
    is <- idio.seg(x, common.seg.out = cs, G.seq = NULL, thr = NULL, d = d, demean = TRUE,
                   cv.args = list(path.length = 10, n.folds = 1, do.cv = FALSE, do.plot = !FALSE), 
                   rule = c('eta', 'epsilon')[2], eta = .5, epsilon = .1)
    est.cp.idio[[i]] <- is$est.cp[, 1]  
    
    for(rr in 1:4){
      ts.plot(is$est.cp.list[[rr]]$norm.stat); abline(h = is$est.cp.list[[rr]]$thr, col = 4)
      abline(v = cp.idio, col = 2, lty = 3); abline(v = is$est.cp.list[[rr]]$cp, col = 4, lty = 2)
    }
    
    cat("model", S, " / run", i, " / est.cp.com: ", sapply(est.cp.com[[i]], paste, collapse=' '), " / est.cp.idio: ", sapply(est.cp.idio[[i]], paste, collapse=' ') )
    cat("\n")
  }
  
  result <- list(est.cp.com = est.cp.com, est.cp.idio = est.cp.idio, cp.common = cp.common, cp.idio = cp.idio, n = n)
  save(result, file=paste0("S", S, "_N", N, ".RData", sep = ""))
  
}
