source("~/Documents/GitHub/favar.segment/common_seg.R")
source("~/Documents/GitHub/favar.segment/idio_seg.R")
source("~/Documents/GitHub/favar.segment/misc.R")

n <- 2000
p <- 50
q <- 3

cp.common <- round(n * c(1/3, 2/3))
#cp.common <- c()

cp.idio <- round(n * (1:3)/4)
#cp.idio <- c()

ss <- sim.data(n = n, p = p, q = q, 
              cp.common = cp.common, den.common = 1, type.common = c('ma', 'ar')[1], ma.order = 2,
              cp.idio = cp.idio, size.idio = 1, do.scale = !FALSE, burnin = 100)

x <- ss$x
#x <- ss$xi

mean.x <- apply(x, 1, mean)
xx <- x - mean.x

dp <- common.spec.est(t(scale(t(x), scale = FALSE)), q = NULL, ic.op = 5, max(1, floor(200^(1/3))))
dp$hl$q.hat

G.seq <- c(200, 300, 400)
#G.seq <- c(250, 350, 450)

##

common.list <- list()

thr <- 1.5

ll.seq <- c()

for(ii in 1:length(G.seq)){
  G <- G.seq[ii]
  ll <- max(1, floor(min(G^(1/3), n/(2 * log(n)))))
  ll.seq <- c(ll.seq, ll)
  
  #thr <- thr.const * p * max(sqrt(ll * log(n)/G), 1/ll, 1/p)
  
  common.list[[ii]] <- 
    #cts <- common.two.step(xx, G, thr = thr, ll, ceiling(log(n)), 'm', 'avg')
    cts <- common.two.step(xx, G, thr = thr, ll, ceiling(log(n)), 'avg') ### added by HYEYOUNG
  common.list[[ii]]$G <- G
  common.list[[ii]]$ll <- ll
  common.list[[ii]]$thr <- thr
  
  matplot(cts$norm.stat, type = 'l'); abline(v = cp.common, lty = 2, col = 2, lwd = 2); abline(v = cp.idio, lty = 3, col = 8); 
  abline(h = thr, col = 3); lines(cts$stat, col = 4, lwd = 2)
  
  # est.cp <- common.search.cp(cts, thr = thr, G, 'max')
  #est.cp <- common.search.cp(cts, thr = thr, G, 'over')
  est.cp <- common.search.cp(cts, thr = thr, G, 'eta') ### added by HYEYOUNG
  #est.cp <- common.search.cp(cts, thr = thr, G, 'epsilon') ### added by HYEYOUNG
  
  # est.cp <- common.check(xx, G, est.cp, thr, ll, NULL, 5, 'm', 'avg')
  est.cp <- common.check(xx, G, est.cp, thr, ll, NULL, 5, 'avg') ### added by HYEYOUNG
  common.list[[ii]]$cp <- est.cp
  
  abline(v = est.cp, col = 4, lty = 3)
  
}

print(est.cp <- bottom.up(common.list, .5))

common.seg.out <- list(est.cp = est.cp, G.seq = G.seq, ll.seq = ll.seq,
                       est.cp.list = common.list, mean.x = mean.x)


#####

est.cp.common <- c()
est.cp.common <- common.seg.out$est.cp[, 1]

ll <- min(common.seg.out$ll.seq)

K <- length(est.cp.common)
if(K > 0) est.cp.common <- sort(est.cp.common)
brks <- c(0, est.cp.common, n)
idx <- rep(c(1:(length(brks) - 1)), diff(brks))

pcfa <- post.cp.fa(x, est.cp.common, NULL, 5, 
                   max(1, 4 * floor((min(diff(brks))/log(min(diff(brks))))^(1/3))))
pcfa$q.seq
Gamma_c <- pcfa$Gamma_c[,, 1:(ll + 1),, drop = FALSE]

if(FALSE){
  dp <- fnets:::dyn.pca(xx[, (brks[1] + 1):brks[2]], q = NULL, ic.op = 5, kern.const = 4)
  dp$hl$q.hat
  rr <- 1
  fields:::imagePlot(dp$acv$Gamma_c[,, rr] - Gamma_c[,, rr, 1])
  fields:::imagePlot(drop(acv.x(ss$x[, (brks[1] + 1):brks[2]] - ss$xi[, (brks[1] + 1):brks[2]], 1)$Gamma_x[,, rr]) - Gamma_c[,, rr, 1])
}

idio.list <- list()
ii <- d <- 1
rule <- c('max', 'over')[2]
epsilon <- .1

for(ii in 1:length(G.seq)){
  
  G <- G.seq[ii]
  thr <- 4 # thr.const * max(sqrt(ll * log(n * p) / G), 1/ll, 1/sqrt(p))
  vv <- G
  stat <- rep(0, n)
  check.cp <- est.cp <- c()
  do.beta <- TRUE
  
  while(vv <= n - G){
    
    int <- (vv - G + 1):vv
    
    if(do.beta){
      icv <- idio.cv(xx = xx[, int, drop = FALSE], Gamma_c = Gamma_c, idx = idx[int], var.order = d, 
                     path.length = 10, n.folds = 1)  
      tb <- tabulate(idx[int], nbins = K + 1)
      acv <- acv.x(xx[, int, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1), drop = FALSE]
      for(kk in 1:(K + 1)) acv <- acv - tb[kk] / G * Gamma_c[,,, kk]
      mg <- make.gg(acv, d)
      beta <- idio.beta(mg$GG, mg$gg, icv$lambda)$beta
      print(null.norm <- max(abs(mg$GG %*% beta - mg$gg)))
    }
    
    if(FALSE){
      A_hat0 <- t(beta); A0 <- ss$A.list[[1]]
      mean(beta != 0); image(beta, col = RColorBrewer::brewer.pal(11, 'RdBu'), breaks = seq(-max(abs(beta)), max(abs(beta)), length.out = 12))
      norm(A_hat0 - A0, 'F')/norm(A0, 'F')
      norm(A_hat0 - A0, '2')/norm(A0, '2')
      fpr <- tpr <- rep(0, 100)
      for(jj in 1:length(fpr)){
        A_hat <- A_hat0
        A_hat[abs(A_hat) < seq(min(abs(A_hat0)) * .999, max(abs(A_hat0)) * .999, length.out = length(fpr))[jj]] <- 0
        fpr[jj] <- sum(A_hat[A0 == 0] != 0)/sum(A0 == 0)
        tpr[jj] <- sum(A_hat[A0 != 0] != 0)/sum(A0 != 0)
      }
      plot(fpr, tpr, type = 'l', xlim = c(0, 1), ylim = c(0, 1))
    }
    
    if(FALSE){
      Gamma_c0 <- Gamma_c * 0
      brks <- c(0, est.cp.common, n)
      for(kk in 1:(K + 1)) Gamma_c0[,,, kk] <- acv.x((ss$x - ss$xi)[, (brks[kk] + 1):brks[kk + 1]], (dim(Gamma_c)[3] - 1)/2)$Gamma_x
      icv0 <- idio.cv(xx = x[, int], Gamma_c = Gamma_c0, idx = idx[int], var.order = d, 
                      path.length = 10, n.folds = 1)
      tb <- tabulate(idx[int], nbins = K + 1)
      acv <- acv.x(xx[, int, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1), drop = FALSE]
      for(kk in 1:(K + 1)) acv <- acv - tb[kk] / G * Gamma_c[,,, kk]
      mg <- make.gg(acv, d)
      beta <- idio.beta(mg$GG, mg$gg, icv$lambda)$beta
      ####
      icv0 <- idio.cv(xx = xi[, int, drop = FALSE], Gamma_c = Gamma_c * 0, idx = idx[int], var.order = d, 
                      path.length = 10, n.folds = 1)
      acv <- acv.x(xi[, int, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1), drop = FALSE]
      mg0 <- make.gg(acv, d)
      beta <- idio.beta(mg0$GG, mg0$gg, icv0$lambda)$beta
      ####
    }
    
    # beta <- t(A[[1]])
    
    diff.Gamma_x <- acv.x(xx[, int, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1), drop = FALSE] -
      acv.x(xx[, int + G, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1), drop = FALSE]
    diff.Gamma_c <- diff.Gamma_x * 0
    if(K > 0){
      common.weights <- tabulate(idx[int], nbins = K + 1) -
        tabulate(idx[int + G], nbins = K + 1)
      for(kk in 1:(K + 1)) diff.Gamma_c <- diff.Gamma_c + common.weights[kk] / G * Gamma_c[,,, kk]
    }
    
    mgd <- make.gg(diff.Gamma_x - diff.Gamma_c, d)
    # print(null.norm <- max(abs(mg$GG %*% beta - mg$gg)))
    stat[vv] <- max(abs(mgd$GG %*% beta - mgd$gg)) / null.norm
    
    if(FALSE){
      ####
      xi <- ss$xi
      diff.Gamma_xi <- acv.x(xi[, int, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1), drop = FALSE] -
        acv.x(xi[, int + G, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1), drop = FALSE]
      mg0 <- make.gg(diff.Gamma_xi, d)
      max(abs(mg0$GG %*% t(A[[1]]) - mg0$gg))
      max(abs(mg0$GG %*% beta - mg0$gg))
      ####
    }
    
    first <- TRUE
    tt <- vv + 1
    check.theta <- tt.max <- n - G
    while(tt <= tt.max){
      for(h in 0:ll){
        diff.Gamma_x[,, h + 1] <- diff.Gamma_x[,, h + 1] -
          xx[, tt - G, drop = FALSE] %*% t(xx[, tt - G + h, drop = FALSE]) / G +
          xx[, tt - h, drop = FALSE] %*% t(xx[, tt, drop = FALSE]) / G +
          xx[, tt, drop = FALSE] %*% t(xx[, tt + h, drop = FALSE]) / G -
          xx[, tt + G - h, drop = FALSE] %*% t(xx[, tt + G, drop = FALSE]) / G
      }
      diff.Gamma_c <- diff.Gamma_x * 0
      if(K > 0){
        common.weights <- tabulate(idx[(tt - G + 1):tt], nbins = K + 1) -
          tabulate(idx[(tt + 1):(tt + G)], nbins = K + 1)
        for(kk in 1:(K + 1)) diff.Gamma_c <- diff.Gamma_c + common.weights[kk] / G * Gamma_c[,,, kk]
      }
      mgd <- make.gg(diff.Gamma_x - diff.Gamma_c, d)
      stat[tt] <- max(abs(mgd$GG %*% beta - mgd$gg)) / null.norm
      
      if(first & stat[tt] > thr){
        check.theta <- tt
        tt.max <- min(tt.max, check.theta + G - 1)
        check.cp <- c(check.cp, check.theta)
        first <- FALSE
      }
      tt <- tt + 1
    }  
    ts.plot(stat, ylim=range(c(stat, thr))); abline(h = thr, col = 3); abline(v = check.cp, col = 6); abline(v = cp.idio, col = 2, lty = 3)
    
    if(check.theta < tt.max){
      hat.theta <- (check.theta:tt.max)[which.max(stat[check.theta:tt.max])]
      if(rule == 'max'){
        est.cp <- c(est.cp, hat.theta)
        vv <- hat.theta + G
        do.beta <- TRUE
      } else if(rule == 'over'){
        int <- max(1, hat.theta - round(epsilon * G) + 1):min(hat.theta + round(epsilon * G), n)
        if(sum(stat[int] < thr) == 0){
          est.cp <- c(est.cp, hat.theta)
          vv <- hat.theta + G
          do.beta <- TRUE
        } else{
          vv <- tt.max + 1
          do.beta <- FALSE
        }
      }
      abline(v = hat.theta, col = 4, lty = 2)
    } else break
    
  }
  
  idio.list[[ii]] <- list(norm.stat = stat, G = G, thr = thr, check.cp = check.cp)
  if(length(est.cp) > 0){
    idio.list[[ii]]$cp <- est.cp
  } else idio.list[[ii]]$cp <- NULL
  
}

bottom.up(common.list, .5)
bottom.up(idio.list, .5)
