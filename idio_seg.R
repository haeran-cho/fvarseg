library(lpSolve)
library(foreach)
library(doParallel)




## idio
# if est.cp.common = c() and q = 0, it becomes var segmentation
idio.seg <- function(x, G.seq, est.cp.common = c(), thr.const, 
                     d = 1, q = NULL, ic.op = 5){
  
  p <- dim(x)[1]
  n <- dim(x)[2]
  if(length(est.cp.common) > 0) est.cp.common <- sort(est.cp.common)
  brks <- c(0, est.cp.common, n)
  
  mean.x <- apply(x, 1, mean)
  xx <- x - mean.x
  
  ll <- max(1, floor(min(G.seq)^(1/3)))
  
  pcfa <- post.cp.factor.analysis(xx, est.cp.common, q, ic.op, ll)
  Gamma_c <- pcfa$Gamma_c[,, 1:(ll + 1), ]
  ## QQ
  # if(sum(pcfa$q.seq) == 0) if there's no factor can use smaller threshold?
    
  common.est.cp.list <- list()
  common.stat.list <- list()
  
  for(ii in 1:length(G.seq)){
    
    G <- G.seq[ii]
    thr <- thr.const * sqrt(ll * log(n * p) / G)
    
    idx <- rep(c(1:(length(brks) - 1)), diff(brks))
    
    diff.Gamma_x <- array(0, dim = c(p, p, ll + 1))
    Gamma_i <- array(0, dim = c(p, p, ll + 1))
    beta <- matrix(0, nrow = p * d, ncol = p)
    
    tt <- G
    do.beta <- TRUE
    while(tt <= n - G){
      ind <- (tt - G + 1):tt
      tb <- table(idx[ind])
      tb.idx <- as.numeric(names(tb))
      
      if(do.beta){
        icv <- idio.cv(xx[, ind, drop = FALSE], Gamma_c, idx)  
      }
      
      # tmp <- acv.x(xx[, ind, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1)]
      # for(ii in tb.idx) tmp <- tmp - tb[ii]/G * Gamma_c[,,, ii]
      # Gamma_i[,, tt - G + 1, ] <- tmp
    }
    
    
    v0 <- G
    
    while(v0 < n-G){
      
      ###################################################
      ##### STEP 0-1: cv for choosing an optimal lambda 
      
      ### cv
      #c1 = seq(0.1, 3, by=0.1)
      #lambda = c1*sqrt(2*log(p)/G)
      # c1 = seq(0.03, 2, length.out = 10) 
      # lambda = c1*max(sqrt(ll*log(n*p)/G), 1/ll, 1/sqrt(p))
      
      ## haeran change
      lambda.max <- max(abs(Gamma_xi[,, v0, 1]))
      lambda <- round(exp(seq(log(lambda.max), log(lambda.max * .0001), length.out = 10)), digits = 10)
      
      #cv.err <- rep(NA, length(lambda))
      #b.nzr <- rep(NA, length(lambda))
      
      ###############################################################################
      ### estimating local auto cov matrix of idiosyncratic part
      #n2 <- 10 # number of replicate for computing cv error
      #n1 <- G - n2 + 1 # length of train set
      n1.train <- ceiling(G*0.5)
      n1.val <- G - n1.train
      
      ###############################################################################
      ### train sets: estimating local auto cov matrix of idiosyncratic part
      
      #Gamma_xi_train <- array(0, dim = c(p, p, n, 2 * ll + 1))
      Gamma_xi_train <- array(0, dim = c(p, p, 2 * ll + 1))
      x2 <- array(0, dim = c(p, p, 2 * ll + 1))
      idx.est.cp.common <- rep(c(1:(length(sort.e.cp)-1)), diff(sort.e.cp))
      
      ## tt = v0 - n1.val
      tt <- v0 - n1.val
      ind.L <- (tt - n1.train + 1):tt
      tb <- table(idx.est.cp.common[ind.L])
      tb.idx <- as.numeric(names(tb))
      Gamma_chi_train <- array(0, dim = c(p, p, 2 * ll + 1))
      
      for(h in 0:(ll - 1)){
        ind <- (tt - n1.train + 1 + h):tt
        x2[,, h + 1] <- x[, ind - h] %*% t(x[, ind])
        for(i in 1:length(tb)){
          Gamma_chi_train[,, h + 1] <- Gamma_chi_train[,, h + 1] + (auto.cov[,,h + 1, tb.idx[i]]*tb[i])
          if(h != 0) Gamma_chi_train[, , 2 * ll + 1 - h + 1] <- t(Gamma_chi_train[, , h + 1])
        }
        Gamma_xi_train[,,h+1] <- (x2[,, h + 1] - Gamma_chi_train[,, h + 1]) / n1.train
        if(h != 0) Gamma_xi_train[,,2 * ll + 1 - h + 1] <- t(Gamma_xi_train[,,h+1])
      }
      
      ###############################################################################
      ### validation sets: estimating local auto cov matrix of idiosyncratic part
      
      Gamma_xi_val <- array(0, dim = c(p, p, 2 * ll + 1))
      x2 <- array(0, dim = c(p, p, 2 * ll + 1))
      idx.est.cp.common <- rep(c(1:(length(sort.e.cp)-1)), diff(sort.e.cp))
      
      ## tt = v0
      tt <- v0
      ind.L <- (tt - n1.val + 1):tt
      tb <- table(idx.est.cp.common[ind.L])
      tb.idx <- as.numeric(names(tb))
      Gamma_chi_val <- array(0, dim = c(p, p, 2 * ll + 1))
      
      for(h in 0:(ll - 1)){
        ind <- (tt - n1.val + 1 + h):tt
        x2[,, h + 1] <- x[, ind - h] %*% t(x[, ind])
        for(i in 1:length(tb)){
          Gamma_chi_val[,, h + 1] <- Gamma_chi_val[,, h + 1] + (auto.cov[,,h + 1, tb.idx[i]]*tb[i])
          if(h != 0) Gamma_chi_val[, , 2 * ll + 1 - h + 1] <- t(Gamma_chi_val[, , h + 1])
        }
        Gamma_xi_val[,, h+1] <- (x2[,, h + 1] - Gamma_chi_val[,, h + 1]) / n1.val
        if(h != 0) Gamma_xi_val[,, 2 * ll + 1 - h + 1] <- t(Gamma_xi_val[,, h+1])
      }
      
      ###############################################################################
      ### cv error
      
      ### setup parallel backend 
      cl <- parallel::makeCluster(7, setup_strategy = "sequential")
      registerDoParallel(cl)
      
      cv.err <- foreach(kk=1:length(lambda), .combine = 'c') %dopar% {
        
        #cv.err <- rep(NA, length(lambda))
        #for(k in 1:length(lambda)){
        
        ### 1) train set
        betahat <- matrix(0, ncol=p, nrow=p)
        xTx <- Gamma_xi_train[, , 1]
        xTy <- Gamma_xi_train[, , 2]
        
        for(j in 1:p){
          ### old
          #bvec=c(c(-xTy[, j], xTy[, j]) - lambda[kk]*rep(1, 2*p))
          #Amat=rbind(cbind(-xTx, xTx), -cbind(-xTx, xTx))
          #result <- linprog::solveLP( cvec=rep(1, 2*p), bvec=bvec, Amat=Amat, const.dir = rep( ">=", length(bvec) ), lpSolve=T )
          #betahat[, j] <- result$solution[1:p] - result$solution[(p+1):(2*p)] 
          #plot(c(betahat[, j]), type="b", main=paste(j))
          
          ### new
          f.rhs = c(rep(lambda[kk], p) - xTy[, j], rep(lambda[kk], p) + xTy[, j])
          f.con = cbind(rbind(-xTx, xTx), -rbind(-xTx, xTx))
          lpout <- lpSolve::lp("min", rep(1, 2*p), f.con, rep("<=",2*p), f.rhs)
          betahat[, j] <- lpout$solution[1:p]-lpout$solution[-(1:p)]
          
        }
        #sum(betahat!=0)
        #sum(A0!=0)
        #sum((A0-betahat)^2)
        #image(betahat)
        #image(A0)
        
        ### 2) validation set
        xTx <- Gamma_xi_val[, , 1]
        xTy <- Gamma_xi_val[, , 2]
        
        sqrt(sum((xTx%*%betahat - xTy)^2))
        #cv.err[k] <- sqrt(sum((xTx%*%betahat - xTy)^2)) # Frobenius norm
        
        ##cv.err[k] <- max(abs(xTx%*%betahat - xTy)) # entrywise maximum norm
        ##b.nzr[k] <- sum(betahat!=0)
      }
      parallel::stopCluster(cl)
      
      plot(lambda, cv.err, type="b", main=paste("cv error / G=", G))
      abline(v=lambda[which.min(cv.err)], col=2)
      
      ### HERE
      #opt.lambda <- lambda[which.min(cv.err)]
      opt.lambda <- max(0.2, lambda[which.min(cv.err)]) # giving a lower bound
      
      opt.lam <- c(opt.lam, opt.lambda)
      
      
      ###################################################
      ##### STEP 0-2: Estimating beta from first G-interval
      
      #c1 <- 1.2
      #lambda = c1*sqrt(2*log(p)/G)
      
      betahat <- matrix(0, ncol=p, nrow=p)
      xTx <- Gamma_xi_train[, , 1]
      xTy <- Gamma_xi_train[, , 2]
      
      for(j in 1:p){
        
        ### old
        #bvec=c(c(-xTy[, j], xTy[, j]) - opt.lambda*rep(1, 2*p))
        #Amat=rbind(cbind(-xTx, xTx), -cbind(-xTx, xTx))
        #result <- linprog::solveLP( cvec=rep(1, 2*p), bvec=bvec, Amat=Amat, const.dir = rep( ">=", length(bvec) ), lpSolve=T )
        #betahat[, j] <- result$solution[1:p] - result$solution[(p+1):(2*p)] 
        #plot(c(betahat[, j]), type="b", main=paste(j))
        
        ### new
        f.rhs = c(rep(opt.lambda, p) - xTy[, j], rep(opt.lambda, p) + xTy[, j])
        f.con = cbind(rbind(-xTx, xTx), -rbind(-xTx, xTx))
        lpout <- lpSolve::lp("min", rep(1, 2*p), f.con, rep("<=",2*p), f.rhs)
        betahat[, j] <- lpout$solution[1:p]-lpout$solution[-(1:p)]
        
      }
      
      
      ###################################################
      ##### STEP 1: T2 TEST STATISTICS
      
      #setup parallel backend 
      #cl <- parallel::makeCluster(4, setup_strategy = "sequential")
      #registerDoParallel(cl)
      
      ### HERE
      thr.idio = thr.const.i*max(sqrt(ll*log(n*p)/G), 1/ll, 1/sqrt(p)) 
      #thr.idio = opt.lambda * 3.0
      
      TST <- rep(0, n)
      
      for(tt in (v0+1):(n-G)){
        
        xTx.1 <- Gamma_xi[, , tt, 1]
        xTy.1 <- Gamma_xi[, , tt, 2]
        
        xTx.2 <- Gamma_xi[, , tt+G, 1]
        xTy.2 <- Gamma_xi[, , tt+G, 2]
        
        D <- (xTx.1%*%betahat-xTy.1) - (xTx.2%*%betahat-xTy.2)
        TST[tt] <- max(abs(D))
        
      }
      
      
      plot(TST, type="l", main=paste("T2 / G=", G), ylab="", xlab="t", ylim=range(c(TST, thr.idio)))
      abline(h=thr.idio, col=3)
      if(!is.null(cp.idio)){
        abline(v=cp.idio, col=2, lty=1, lwd=1) 
      }
      if(length(which(TST > thr.idio))>0){
        abline(v=min(which(TST > thr.idio)), col=4, lty=2, lwd=2) 
      }
      #abline(v=examine.range[1], col=4, lty=2, lwd=2)
      #abline(v=examine.range[length(examine.range)], col=4, lty=2, lwd=2)
      
      ###################################################
      ##### finding a maximum
      survived.idx <- which(TST > thr.idio)
      
      if(length(survived.idx)>0){
        theta.check <- min(survived.idx)
        examine.range <- c(max(theta.check-G, v0):min(theta.check+G, n-G))
        #examine.range <- c(max(theta.check-G, v0+G):min(theta.check+G, n-G))
        
        abline(v=examine.range[1], col=4, lty=1, lwd=2)
        abline(v=examine.range[length(examine.range)], col=4, lty=1, lwd=2)
        
        if(theta.check %in% examine.range){
          
          th.new <- which(TST==max(TST[examine.range]))
          
          est.cp.idio.G <- c(est.cp.idio.G, th.new)
          #v0 <- min(theta.check + 2*G, th.new + 1.5*G)
          v0 <- min(theta.check + 2*G, th.new + 1.2*G)
          abline(v=th.new, col=2, lty=2)
          abline(v=v0, col=3)
          
        }else{
          break
        }
      }else{
        break
      }
      
    } #  while
    
    est.cp.idio.Gs[[k]] <- est.cp.idio.G
  }
  
  return(est.cp.idio.Gs)
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
  ll <- dim(Gamma_c)[3] - 1
  
  if(is.null(lambda.max)) lambda.max <- max(abs(xx %*% t(xx)/nn)) * 1
  lambda.path <- round(exp(seq(log(lambda.max), log(lambda.max * .001), length.out = path.length)), digits = 10)
  
  cv.err.mat <- matrix(0, nrow = path.length, ncol = length(var.order))
  ind.list <- split(1:nn, ceiling(n.folds*(1:nn)/nn)) 
  
  for(fold in 1:n.folds){ 
    ind <- 1:ceiling(length(ind.list[[fold]]) * .5)
    
    tx <- xx[, ind.list[[fold]][ind]]
    tb <- table(idx[ind])
    tb.idx <- as.numeric(names(tb))
    train.acv <- acv.x(tx, ll)$Gamma_x[,, 1:(ll + 1)]
    for(ii in tb.idx) train.acv <- train.acv - tb[ii]/length(ind) * Gamma_c[,,, ii]
    
    
    ind <- setdiff(ind.list[[fold]], ind)
    tx <- xx[, ind.list[[fold]][ind]]
    tb <- table(idx[ind])
    tb.idx <- as.numeric(names(tb))
    test.acv <- acv.x(tx, ll)$Gamma_x[,, 1:(ll + 1)]
    for(ii in tb.idx) test.acv <- test.acv - tb[ii]/length(ind) * Gamma_c[,,, ii]
    
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
        cv.err.mat[ii, jj] <- cv.err.mat[ii, jj] + sum(diag(test.acv[,, 1] - beta.gg - t(beta.gg) + t(train.beta) %*% test.GG %*% train.beta))
      }
    }
  }
  
  cv.err.mat
  lambda.min <- lambda.path[which.min(apply(cv.err.mat, 1, min))]
  order.min <- var.order[which.min(apply(cv.err.mat, 2, min))]
  
  matplot(lambda.path, cv.err.mat, type = 'b', col = 2:(length(var.order) + 1), pch = 2:(length(var.order) + 1), log = 'x', xlab = 'λ (log scale)', ylab = 'CV error')
  abline(v = lambda.min)
  legend('topleft', legend = var.order, col = 2:(length(var.order) + 1), pch = 2:(length(var.order) + 1), lty = 1)
  
  out <- list(lambda = lambda.min, var.order = order.min,  cv.error = cv.err.mat, lambda.path = lambda.path)
  return(out) 
  
}

#' @keywords internal
post.cp.factor.analysis <- function(xx, est.cp.common, q = NULL, ic.op = 5, ll){
  
  p <- dim(xx)[1]
  n <- dim(xx)[2]
  
  brks <- c(0, est.cp.common, n)
  len <- 2 * ll
  thetas <- 2 * pi * (0:len)/(len + 1)
  w <- bartlett.weights(((-ll):ll)/ll)
  
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
      fne <- factor.number.est(xx[, int, drop = FALSE], q.max, ll, w)
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
      Gamma_c[,,, jj] <-  (aperm(apply(Sigma_c, c(1, 2), fft, inverse = TRUE), c(2, 3, 1)) ) * (2 * pi) / (2 * ll + 1)
      Gamma_c[,,, jj] <- Re(Gamma_c[,,, jj])
    } 
  }
  ls <- list(Gamma_c = Gamma_c, q.seq = q.seq)
  
  return(ls)
  
}

## misc

