library(lpSolve)
library(foreach)
library(doParallel)

## idio
# if est.cp.common = c() and q = 0, it becomes var segmentation
idio.seg <- function(x, G.seq, d = 1, thr.const, demean = TRUE,
                     common.seg.out, q = NULL, ic.op = 5,
                     path.length = 10, n.folds = 1){
  
  p <- dim(x)[1]
  n <- dim(x)[2]
  
  est.cp.common <- common.seg.out$est.cp
  K <- length(est.cp.common)
  if(K > 0) est.cp.common <- sort(est.cp.common)
  brks <- c(0, est.cp.common, n)
  idx <- rep(c(1:(length(brks) - 1)), diff(brks))
  
  ll <- common.seg.out$ll.seq[1]
  
  if(demean) mean.x <- apply(x, 1, mean) else mean.x <- rep(0, p)
  xx <- x - mean.x
  
  pcfa <- post.cp.fa(xx, est.cp.common, q, ic.op, ll)
  Gamma_c <- pcfa$Gamma_c[,, 1:(ll + 1),, drop = FALSE]

  idio.list <- list()
  
  for(ii in 1:length(G.seq)){
    
    G <- G.seq[ii]
    thr <- thr.const * max(sqrt(ll * log(n * p) / G), 1/ll, 1/sqrt(p))
    vv <- G
    stat <- rep(0, n)
    check.cp <- est.cp <- c()
    
    while(vv <= n - G){
      
      int <- (vv - G + 1):vv
      icv <- idio.cv(xx = xx[, int, drop = FALSE], Gamma_c = Gamma_c, idx = idx, var.order = d, 
                     path.length = path.length, n.folds = n.folds)  
      tb <- tabulate(idx[int], nbins = K + 1)
      acv <- acv.x(xx[, int, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1), drop = FALSE]
      for(kk in 1:(K + 1)) acv <- acv - tb[kk] / G * Gamma_c[,,, kk]
      mg <- make.gg(acv, d)
      beta <- idio.beta(mg$GG, mg$gg, icv$lambda)$beta
      # mean(beta != 0); fields::imagePlot(beta, col = RColorBrewer::brewer.pal(11, 'RdBu'), breaks = seq(-max(abs(beta)), max(abs(beta)), length.out = 12))
      
      diff.Gamma_x <- acv.x(xx[, int, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1), drop = FALSE] - 
        acv.x(xx[, int + G, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1), drop = FALSE]
      diff.Gamma_c <- diff.Gamma_x * 0
      if(K > 0){
        common.weights <- tabulate(idx[int], nbins = K + 1) - 
          tabulate(idx[int + G], nbins = K + 1)
        for(kk in 1:(K + 1)) diff.Gamma_c <- diff.Gamma_c + 
            common.weights[kk] / G * Gamma_c[,,, kk]
      }
      
      mg <- make.gg(diff.Gamma_x - diff.Gamma_c, d)
      null.norm <- max(abs(mg$GG %*% beta - mg$gg))
      stat[vv] <- 1
      
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
        mg <- make.gg(diff.Gamma_x - diff.Gamma_c, d)
        stat[tt] <- max(abs(mg$GG %*% beta - mg$gg)) / null.norm
        
        if(first & stat[tt] > thr){
          check.theta <- tt
          tt.max <- min(tt.max, check.theta + G - 1)
          check.cp <- c(check.cp, check.theta)
          first <- FALSE
        }
        tt <- tt + 1
      }  
      # ts.plot(stat); abline(h = thr, col = 3); abline(v = check.cp, col = 6); abline(v = cp.idio, col = 2, lty = 3)
      
      if(check.theta < tt.max){
        hat.theta <- (check.theta:tt.max)[which.max(stat[check.theta:tt.max])]
        est.cp <- c(est.cp, hat.theta)
        vv <- hat.theta + G
      } else break
      
    }
    
    idio.list[[ii]]$cp <- est.cp
    idio.list[[ii]]$stat <- stat
    idio.list[[ii]]$G <- G
    idio.list[[ii]]$thr <- thr
    idio.list[[ii]]$check.cp <- check.cp
  }
  
  est.cp <- bottom.up(idio.list, eta)
  
  out <- list(est.cp = est.cp, est.cp.list = idio.list, mean.x = mean.x)
  return(out)
  
}

#' @keywords internal
idio.beta <- function(GG, gg, lambda, n.cores = min(parallel::detectCores() - 1, 3)){
  
  p <- dim(gg)[2]
  d <- dim(gg)[1]/dim(gg)[2]
  beta <- gg * 0
  
  f.obj <- rep(1, 2 * p * d)
  f.con <- rbind(-GG, GG)
  f.con <- cbind(f.con,-f.con)
  f.dir <- rep('<=', 2 * p * d)
  
  cl <- parallel::makePSOCKcluster(n.cores)
  doParallel::registerDoParallel(cl)
  
  beta <- foreach::foreach(ii = 1:p, .combine = 'cbind', .multicombine = TRUE, .export = c('lp')) %dopar% {
    b1 <- rep(lambda, p * d) - gg[, ii]
    b2 <- rep(lambda, p * d) + gg[, ii]
    f.rhs <- c(b1, b2)
    lpout <- lp('min', f.obj, f.con, f.dir, f.rhs)
    lpout$solution[1:(p * d)] - lpout$solution[-(1:(p * d))]
  }
  parallel::stopCluster(cl)
  
  out <- list(beta = beta, lambda = lambda, Gamma = Gamma)
  return(out)
  
}

#' @keywords internal
idio.cv <- function(xx, Gamma_c, idx, lambda.max = NULL, var.order = 1, 
                    path.length = 10, n.folds = 1){
  
  nn <- ncol(xx)
  p <- nrow(xx)
  dd <- max(1, min(max(var.order), dim(Gamma_c)[3] - 1))
  K <- dim(Gamma_c)[4] - 1
  
  if(is.null(lambda.max)) lambda.max <- max(abs(xx %*% t(xx)/nn)) * 1
  lambda.path <- round(exp(seq(log(lambda.max), log(lambda.max * .005), length.out = path.length)), digits = 10)
  
  cv.err.mat <- matrix(0, nrow = path.length, ncol = length(var.order))
  ind.list <- split(1:nn, ceiling(n.folds*(1:nn)/nn)) 
  
  for(fold in 1:n.folds){ 
    ind <- 1:ceiling(length(ind.list[[fold]]) * .5)
    tx <- xx[, ind.list[[fold]][ind]]
    tb <- tabulate(idx[ind.list[[fold]][ind]], nbins = K + 1)
    train.acv <- acv.x(tx, dd)$Gamma_x[,, 1:(dd + 1), drop = FALSE]
    for(ii in 1:(K + 1)) train.acv <- train.acv - tb[ii]/length(ind) * Gamma_c[,, 1:(dd + 1), ii]
    
    ind <- setdiff(1:length(ind.list[[fold]]), ind)
    tx <- xx[, ind.list[[fold]][ind]]
    tb <- tabulate(idx[ind.list[[fold]][ind]], nbins = K + 1)
    test.acv <- acv.x(tx, dd)$Gamma_x[,, 1:(dd + 1), drop = FALSE]
    for(ii in 1:(K + 1)) test.acv <- test.acv - tb[ii]/length(ind) * Gamma_c[,, 1:(dd + 1), ii]
    
    for(jj in 1:length(var.order)){
      if(var.order[jj] >= dim(train.acv)[3]){
        cv.err.mat[, jj] <- Inf
        next
      }
      mg <- make.gg(train.acv, var.order[jj])
      gg <- mg$gg; GG <- mg$GG
      mg <- make.gg(test.acv, var.order[jj])
      test.gg <- mg$gg; test.GG <- mg$GG
      for(ii in 1:path.length){
        train.beta <- idio.beta(GG, gg, lambda = lambda.path[ii])$beta
        beta.gg <- t(train.beta) %*% test.gg
        cv.err.mat[ii, jj] <- cv.err.mat[ii, jj] + 
          sum(diag(test.acv[,, 1] - beta.gg - t(beta.gg) + t(train.beta) %*% test.GG %*% train.beta))
      }
    }
  }
  
  cv.err.mat[cv.err.mat < 0] <- Inf
  lambda.min <- min(lambda.path[apply(cv.err.mat, 1, min) == min(apply(cv.err.mat, 1, min))])
  order.min <- min(var.order[apply(cv.err.mat, 2, min) == min(apply(cv.err.mat, 2, min))])
  
  matplot(lambda.path, cv.err.mat, type = 'b', col = 2:(length(var.order) + 1), pch = 2:(length(var.order) + 1), log = 'x', xlab = 'λ (log scale)', ylab = 'CV error')
  abline(v = lambda.min)
  legend('topleft', legend = var.order, col = 2:(length(var.order) + 1), pch = 2:(length(var.order) + 1), lty = 1)
  
  out <- list(lambda = lambda.min, var.order = order.min,  cv.error = cv.err.mat, lambda.path = lambda.path)
  return(out) 
  
}

#' @keywords internal
post.cp.fa <- function(xx, est.cp.common, q = NULL, ic.op = 5, ll){
  
  p <- dim(xx)[1]
  n <- dim(xx)[2]
  
  brks <- c(0, est.cp.common, n)
  len <- 2 * ll
  thetas <- 2 * pi * (0:len)/(len + 1)
  w <- Bartlett.weights(((-ll):ll)/ll)
  
  Gamma_c <- array(0, dim = c(p, p, 2 * ll + 1, length(brks) - 1))
  q.seq <- rep(0, length(brks) - 1)
  
  for(jj in 1:(length(brks) - 1)){
    int <- (brks[jj] + 1):brks[jj + 1]
    nn <- length(int)
    if(!is.null(q)){
      qq <- as.integer(q)
      ax <- acv.x(xx, ll, w)
      Gamma_x <- ax$Gamma_x
      Gamma_xw <- ax$Gamma_xw
      Sigma_x <- aperm(apply(Gamma_xw, c(1, 2), fft), c(2, 3, 1)) / (2 * pi)  
      sv <- list(1:(ll + 1))
      if(qq > 0) for(ii in 1:(ll + 1)) sv[[ii]] <- svd(Sigma_x[, , ii], nu = qq, nv = 0) 
    }
    if(is.null(q)){
      q.max <- min(50, floor(sqrt(min(nn - 1, p))))
      fne <- hl.factor.number(xx[, int, drop = FALSE], q.max, ll, w, do.plot = FALSE, center = FALSE)
      qq <- fne$q.hat[ic.op]
      Gamma_x <- fne$Gamma_x
      Sigma_x <- fne$Sigma_x
      sv <- fne$sv
    }
    q.seq[jj] <- qq
    
    if(qq >= 1){
      Sigma_c <- Sigma_x * 0
      for(ii in 1:(ll + 1)){
        Sigma_c[,, ii] <- sv[[ii]]$u[, 1:qq, drop = FALSE] %*% diag(sv[[ii]]$d[1:qq], qq) %*% Conj(t(sv[[ii]]$u[, 1:qq, drop = FALSE]))
        if(ii > 1){
          Sigma_c[,, 2 * ll + 1 - (ii - 1) + 1] <- Conj(Sigma_c[,, ii])
        }
      }
      Gamma_c[,,, jj] <- aperm(apply(Sigma_c, c(1, 2), fft, inverse = TRUE), c(2, 3, 1)) * (2 * pi) / (2 * ll + 1)
    } 
  }
  ls <- list(Gamma_c = Re(Gamma_c), q.seq = q.seq)
  
  return(ls)
  
}

