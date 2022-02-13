rm(list=ls())

#devtools::install_github("https://github.com/Dom-Owens-UoB/fnets")

#setwd("~/Documents/GitHub/fnets.segment")
setwd("/Volumes/Seagate Expansion Drive/Current.projects/with_Haeran/Rcode/new_Haeran_ver4")
#setwd("~/gdfm_new")
#setwd("/Volumes/Extreme SSD/new_Haeran_ver4")

source('idio_seg.R')
source('common_seg.R')
source('misc.R')
#source('package.R') ## contains the codes from 
#Rcpp::sourceCpp('input.cpp')

library(VARDetect)

N <- 1
n <- 2000
d <- 1

q <- 2
ma.order <- 2




###################################################################################

### models when chi exists without changes
sim.list = expand.grid(structure(list( p = c(50, 100, 150), 
                                       cp.idio = list(c(NULL), c(round(n * c(3,5)/8)), c(round(n * 1:2/3))), 
                                       cp.common= list(c(NULL)))))
sim.list


##### simulations
for(S in 1:dim(sim.list)[1]){
  
  p = c(sim.list[S,][1], recursive=T, use.names=F)
  cp.idio = c(sim.list[S,][2], recursive=T, use.names=F)
  cp.common = c(sim.list[S,][3], recursive=T, use.names=F)  

  est.cp.idio <- vector("list", length = N)
  est.cp.VARD <- vector("list", length = N)
  time.all <- matrix(NA, ncol = 2, nrow = N)
  
  for(i in 1:N){
    
    ### 1) data generation 
    ss <- sim.data0(n, p, q,  
                    cp.common = cp.common, den.common = 0,
                    cp.idio = cp.idio, size.idio = 1, d = d, do.scale = TRUE, seed = i)
    #ss <- sim.data(n, p, q,  
    #               cp.common = cp.common, den.common = 1, type.common = type.common, ma.order = 2,
    #               cp.idio = cp.idio, size.idio = 1.1, do.scale = TRUE, seed = i)
    #ss <- sim.data2(n, p, q,
    #                cp.common = cp.common, den.common = .5, type.common = type.common,
    #                cp.idio = cp.idio, size.idio = 1.1, d = d, do.scale = !FALSE, seed = i)
    
    #x <- ss$x
    x <- ss$xi
    
    ##################################
    ### 2) idiosyncratic component 
    cs <- list(est.cp = c(), ll.seq = max(1, floor(min((n/10)^(1/3), n/(2 * log(n))))))
    
    tic <- proc.time()
    is <- idio.seg(x, common.seg.out = cs, G.seq = NULL, thr = NULL, d = d, demean = TRUE,
                   cv.args = list(path.length = 10, n.folds = 1, do.cv = FALSE, do.plot = !FALSE), 
                   rule = c('eta', 'epsilon')[1], eta = .5, epsilon = .1)
    toc <- proc.time()
    time.all[i, 1] <- (toc-tic)[3]
    
    est.cp.idio[[i]] <- is$est.cp[, 1]  
    
    for(rr in 1:4){
      ts.plot(is$est.cp.list[[rr]]$norm.stat, ylim=range(c(is$est.cp.list[[rr]]$norm.stat, is$est.cp.list[[rr]]$thr))); abline(h = is$est.cp.list[[rr]]$thr, col = 4)
      abline(v = cp.idio, col = 2, lty = 3); abline(v = is$est.cp.list[[rr]]$cp, col = 4, lty = 2)
    }
    
    ##################################
    ### 3) VARDetect 
    tic <- proc.time()
    block <- VARDetect::tbss(t(x), method = 'sparse', q = d, use.BIC = F)
    toc <- proc.time()
    time.all[i, 2] <- (toc-tic)[3]
    
    est.cp.VARD[[i]] <-  block$cp
    
  
    cat("model", S, " / run", i, 
        " / est.cp.idio: ", sapply(est.cp.idio[[i]], paste, collapse=' '),
        " / est.cp.VARD: ", sapply(est.cp.VARD[[i]], paste, collapse=' '),
        " / time(idio): ",  time.all[i, 1],
        " / time(VARD): ",  time.all[i, 2])
    cat("\n")
    
  }
  
  result <- list(est.cp.idio = est.cp.idio,  est.cp.VARD =  est.cp.VARD,
                 cp.idio = cp.idio, time.all = time.all,  n = n)
  #save(result, file=paste0("S", S, "_VARD_sizeidio1_d1_N", N, ".RData", sep = ""))
  #save(result, file=paste0("S", S, "_simdat0_den05_sizeidio1_etaeta_cpIDX1_N", N, ".RData", sep = ""))
  
}








