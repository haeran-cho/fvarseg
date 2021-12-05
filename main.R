sim.data <- function(n, p, q,  
                     cp.common = c(), den.common = 1, type.common = c('ma', 'ar')[1], ma.order = 0,
                     cp.idio = c(), size.idio = 1, burnin = 100, seed = NULL){
  
  if(is.null(seed)) set.seed(seed)
  
  # common component
  if(q > 0){
    brks <- c(0, c(cp.common, n) + burnin)
    
    if(type.common == 'ma'){
      u <- sqrt(rep(c(1, .5, 1.5), ceiling(q/3))[1:q]) * matrix(rnorm(q * (n + burnin)), nrow = q)
      chi <- matrix(0, nrow = p, ncol = n + burnin)
      for(ii in 1:(ma.order + 1)){
        bb <- matrix(rnorm(p * q), nrow = p) 
        chi[, ii:brks[2]] <- chi[, ii:brks[2]] + bb %*% u[, ii:brks[2] - ii + 1]
        
        if(length(cp.common) >= 1){  
          for(j in 1:length(cp.common)){
            cp.ind <- sample(p, floor(den.common * p))
            bb[cp.ind, ] <- matrix(rnorm(length(cp.ind) * q), ncol = q)
            chi[, (brks[j + 1] + 1):brks[j + 2]] <- chi[, (brks[j + 1] + 1):brks[j + 2]] + 
              bb %*% u[, (brks[j + 1] + 1):brks[j + 2] - ii + 1]
          }
        }
      }
    }
    
    if(type.common == 'ar'){ 
      u <- matrix(rnorm(q * (n + burnin)), nrow = q)
      chi <- matrix(0, nrow = p, ncol = n + burnin)
      bb <- matrix(runif(p * q, -1, 1), nrow = p)
      
      for(k in 0:length(cp.common)){
        tmp <- matrix(0, nrow = p, ncol = n + burnin)
        if(k == 0) ar.coef <- matrix(runif(p * q, -.8, .8), nrow = p)
        if(k >= 1){
          cp.ind <- sample(p, floor(den.common * p))
          ar.coef[cp.ind, ] <- -ar.coef[cp.ind, ]
        }
        for(ii in 1:p){
          uu <- u
          for(j in 1:q){
            for(tt in 2:(n + burnin)) {
              uu[j, tt] <- ar.coef[ii, j] * uu[j, tt - 1] + uu[j, tt]
            }
            tmp[ii, ] <- tmp[ii, ] + bb[ii, j] * uu[j, ]
          }
        }
        chi[, (brks[k + 1] + 1):brks[k + 2]] <- tmp[, (brks[k + 1] + 1):brks[k + 2]]
      }
    }
    chi <- chi[, -(1:burnin)] 
  } else chi <- matrix(0, nrow = p, ncol = n)
  
  ## idio commonponent
  
  sig_vep <- toeplitz(.3^(1:p - 1)) 
  xi <- vep <- t(mvtnorm::rmvnorm(n + burnin, sigma = sig_vep))
  
  A <- matrix(0, nrow = p, ncol = p)
  c.A <- .4
  diag(A) <- c.A
  A[row(A) + 1 == col(A)] <- - c.A
  
  brks <- c(0, c(cp.idio, n) + burnin)
  
  for(tt in (brks[1] + 1 + 1):brks[2]) xi[, tt] <- vep[, tt] + A %*% xi[, tt - 1]
  if(length(cp.idio) >= 1){
    for(j in 1:length(cp.idio)){
      A[row(A) + 1 == col(A)] <- - size.idio * A[row(A) + 1 == col(A)]
      for(tt in (brks[j + 1] + 1):brks[j + 2]) xi[, tt] <- vep[, tt] + A %*% xi[, tt - 1]
    }
  }
  x <- xi <- xi[, -(1:burnin)]
  
  if(q > 0) x <- x + chi/apply(chi, 1, sd) * apply(xi, 1, sd) 

  return(x)
  
}

## common

common.seg <- function(x, G.seq, thr.const, tt.by = round(log(dim(x)[2])), q = NULL, ic.op = 5,
                       norm.type = c('f', '2', 'm'), agg.over.freq = c('avg', 'max'), eta = .5, do.check = FALSE){
  
  p <- dim(x)[1]
  n <- dim(x)[2]
  
  common.est.cp.list <- list()
  common.stat.list <- list()
  
  for(ii in 1:length(G.seq)){
    G <- G.seq[ii]
    ll <- max(1, floor(G^(1/3)))
    thr <- thr.const * p * max(sqrt(ll * log(n)/G), 1/ll, 1/p)
    
    common.stat.list[[ii]] <- cts <- common.two.step(x, G, thr, ll, tt.by, norm.type, agg.over.freq)
    common.stat.list[[ii]]$G <- G
    common.stat.list[[ii]]$ll <- ll
    common.stat.list[[ii]]$thr <- thr
    
    est.cp <- common.search.cp(cts, thr, G, eta)
    if(do.check) est.cp <- common.check(x, G, est.cp, thr, ll, q, ic.op, norm.type, agg.over.freq)
    common.est.cp.list[[ii]] <- est.cp
    # matplot(cts$norm.stat, type = 'l'); abline(v = cp.common, lty = 2, col = 2, lwd = 2); abline(v = cp.idio, lty = 3, col = 6); abline(v = est.cp, col = 4, lty = 3); abline(h = thr, col = 3); lines(cts$stat, col = 4, lwd = 2)
    
  }
  
  out <- list(est.cp.list  = common.est.cp.list, stat.list = common.stat.list)
  return(out)
  
}

#' @keywords internal
common.two.step <- function(x, G, thr, ll, tt.by, norm.type, agg.over.freq){  
  
  p <- dim(x)[1]
  n <- dim(x)[2]
  
  w <- weights.Bartlett(((-ll):ll)/ll)
  len <- 2 * ll
  thetas <- 2 * pi * (0:len)/(len + 1)
  
  null <- x1 <- array(0, dim = c(p, p, ll + 1))
  for(h in 0:ll){
    ind.l <- (1 + h):G
    ind.r <- (n - G + 1 + h):n
    x1[,, h + 1] <- x[, ind.l - h] %*% t(x[, ind.l]) + x[, ind.r - h] %*% t(x[, ind.r])
    x1[,, h + 1] <- x1[,, h + 1]/(2 * G)
    for(ii in 1:(ll + 1)) {
      null[,, ii] <- null[,, ii] + 
        x1[,, h + 1] * w[h + ll + 1] * exp(-(0 + 1i) * h * thetas[ii]) +
        (h > 0) * t(x1[,, h + 1]) * w[h + ll + 1] * exp((0 + 1i) * h * thetas[ii])
    }
  }
  null <- null/(2 * pi)
  null.norm <- apply(abs(null), 3, norm, type = norm.type)
  if(norm.type %in% c('f', '2')) null.norm <- null.norm / sqrt(p)
  
  tt.seq <- round(seq(G, n - G, by = tt.by))
  stat0 <- common.stat0(x, G, thr, ll, tt.seq)
  norm.stat <- abs(t(t(stat0) / null.norm))
  if(agg.over.freq == 'avg') stat <- apply(norm.stat, 1, mean)
  if(agg.over.freq == 'max') stat <- apply(norm.stat, 1, max)
  
  # matplot(norm.stat, type = 'l'); abline(v = cp.common, lty = 2, col = 2); abline(h = thr, col = 3); lines(stat, col = 4, lwd = 2)

  tt.list <- common.tt.list(norm.stat, G, thr, tt.seq, tt.by)
  
  if(length(tt.list) > 0){
    for(ii in 1:length(tt.list)){
      s <- min(tt.list[[ii]]); e <- max(tt.list[[ii]])
      stat0 <- common.stat1(x, G, thr, ll, s, e, stat0)
    }
    norm.stat <- abs(t(t(stat0)/apply(abs(null), 3, norm, type = norm.type)))
    if(norm.type == 'f') norm.stat <- norm.stat * sqrt(p)

    if(agg.over.freq == 'avg') stat <- apply(norm.stat, 1, mean)
    if(agg.over.freq == 'max') stat <- apply(norm.stat, 1, max)
  } 
  
  ls <- list(norm.stat = norm.stat, stat = stat)
  return(ls)
  
}

#' @keywords internal
common.search.cp <- function(cts, thr, G, eta = .5){
  
  n <- length(cts$stat)
  est.cp <- c()
  new.norm.stat <- cts$norm.stat * (cts$norm.stat > thr)
  survived <- apply(new.norm.stat > 0, 1, max) == 1
  new.stat <- cts$stat
  
  while(sum(survived) > 0){
    agg.stat <- apply(new.norm.stat, 1, max)
    mv <- max(agg.stat)
    theta.hat0 <- min(which(agg.stat == mv))
    lstar <- which.max(new.norm.stat[theta.hat0, ])
    int <- max(1, theta.hat0 - round(eta * G) + 1):min(theta.hat0 + round(eta * G), n)
    if(new.norm.stat[theta.hat0, lstar] >= max(new.norm.stat[int, lstar])){
      theta.hat <- int[1] + which.max(new.stat[int]) - 1
      est.cp <- c(est.cp, theta.hat)
    }
    int <- max(1, theta.hat - G + 1):min(theta.hat + G, n)
    survived[int] <- FALSE
    new.norm.stat[!survived, ] <- 0
  }
   
  est.cp
  
}

#' @keywords internal
common.tt.list <- function(eval, G, thr, tt.seq, tt.by){
  
  if(!is.matrix(eval)) eval <- matrix(eval, ncol = 1)
  n <- dim(eval)[1]
  which.survived <- which(apply(eval[tt.seq, ] > thr, 1, max) == 1)
  tt.list <- list()
  if(length(which.survived) > 0){
    a <- cbind(c(1, which(diff(which.survived) > 2) + 1), 
               c(which(diff(which.survived) > 2), length(which.survived)))
    
    for(ii in 1:(dim(a)[1])){
      l <- max(tt.seq[which.survived[a[ii, 1]]] - tt.by + 2, G)
      r <- min(tt.seq[which.survived[a[ii, 2]]] + tt.by - 1, n - G)
      tt.list[[ii]] <- l:r
    }
  } 
  
  return(tt.list)
  
}

#' @keywords internal
common.stat0 <- function(x, G, thr, ll, tt.seq){  
  
  p <- dim(x)[1]
  
  w <- weights.Bartlett(((-ll):ll)/ll)
  len <- 2 * ll
  thetas <- 2 * pi * (0:len)/(len + 1)
  
  x2 <- array(0, dim = c(p, p, ll + 1))  
  stat0 <- matrix(0, nrow = n, ncol = ll + 1)
  
  for(ttt in 1:length(tt.seq)){
    eval <- array(0, dim = c(p, p, ll + 1))  
    tt <- tt.seq[ttt]
    for(h in 0:ll){
      ind.l <- (tt - G + 1 + h):tt 
      ind.r <- (tt + 1 + h):(tt + G)
      x2[,, h + 1] <- x[, ind.l - h] %*% t(x[, ind.l]) - x[, ind.r - h] %*% t(x[, ind.r])
      
      for(ii in 1:(ll + 1)){    
        eval[,, ii] <- eval[,, ii] + 
          x2[,, h + 1] * w[h + ll + 1] * exp(-(0 + 1i) * h * thetas[ii]) +
          (h > 0) * t(x2[,, h + 1]) * w[h + ll + 1] * exp((0 + 1i) * h * thetas[ii]) 
      }
    }
    eval <- eval/(2 * pi * G)
    for(ii in 1:(ll + 1)) stat0[tt, ii] <- norm(eval[,, ii], type = '2')
  }
  
  stat0
  
}

#' @keywords internal
common.stat1 <- function(x, G, thr, ll, s, e, stat0){  
  
  p <- dim(x)[1]
  
  w <- weights.Bartlett(((-ll):ll)/ll)
  len <- 2 * ll
  thetas <- 2 * pi * (0:len)/(len + 1)
  
  x2 <- eval <- array(0, dim = c(p, p, ll + 1)) 
  ttt <- 1; tt <- s
  for(h in 0:ll){
    ind.l <- (tt - G + 1 + h):tt 
    ind.r <- (tt + 1 + h):(tt + G) 
    x2[,, h + 1] <- x[, ind.l - h] %*% t(x[, ind.l]) - x[, ind.r - h] %*% t(x[, ind.r])
    for(ii in 1:(ll + 1)){    
      eval[,, ii] <- eval[,, ii] + 
        x2[,, h + 1] * w[h + ll + 1] * exp(-(0 + 1i) * h * thetas[ii]) +
        (h > 0) * t(x2[,, h + 1]) * w[h + ll + 1] * exp((0 + 1i) * h * thetas[ii])
    }
  }
  eval <- eval/(2 * pi * G)
  for(ii in 1:(ll + 1)) stat0[tt, ii] <- norm(eval[,, ii], type = '2')
  
  for(ttt in 2:(e - s + 1)){
    tt <- (s:e)[ttt]
    eval <- eval * 0
    for(h in 0:ll){
      x2[,, h + 1] <- x2[,, h + 1] - 
        x[, tt - G, drop = FALSE]%*%t(x[, tt - G + h, drop = FALSE]) +
        x[, tt - h, drop = FALSE]%*%t(x[, tt, drop = FALSE]) +
        x[, tt, drop = FALSE]%*%t(x[, tt + h, drop = FALSE]) -
        x[, tt + G - h, drop = FALSE]%*%t(x[, tt + G, drop = FALSE])
      
      for(ii in 1:(ll + 1)){    
        eval[,, ii] <- eval[,, ii] + 
          x2[,, h + 1] * w[h + ll + 1] * exp(-(0 + 1i) * h * thetas[ii]) +
          (h > 0) * t(x2[,, h + 1]) * w[h + ll + 1] * exp((0 + 1i) * h * thetas[ii])
      }
    }
    eval <- eval/(2 * pi * G)
    for(ii in 1:(ll + 1)) stat0[tt, ii] <- norm(eval[,, ii], type = '2')
  }
  
  stat0
  
}

#' @keywords internal
common.spec.est <- function(xx, q = NULL, ic.op = 5, ll){
  
  p <- dim(xx)[1]
  n <- dim(xx)[2]
  
  len <- 2 * ll
  thetas <- 2 * pi * (0:len)/(len + 1)
  w <- weights.Bartlett(((-ll):ll)/ll)
  
  if(!is.null(q)){
    qq <- NA
    q <- as.integer(q)
    Gamma_x <- Gamma_xw <- array(0, dim = c(p, p, 2 * ll + 1))
    for(h in 0:ll){
      Gamma_x[, , h + 1] <- xx[, 1:(n - h)] %*% t(xx[, 1:(n - h) + h])/n 
      Gamma_xw[, , h + 1] <- Gamma_x[, , h + 1] * w[h + ll + 1]
      if(h != 0){
        Gamma_x[, , 2 * ll + 1 - h + 1] <- t(Gamma_x[, , h + 1])
        Gamma_xw[, , 2 * ll + 1 - h + 1] <- t(Gamma_xw[, , h + 1])
      }
    }
    Sigma_x <- aperm(apply(Gamma_xw, c(1, 2), fft), c(2, 3, 1)) / (2 * pi)  
    sv <- list(1:(ll + 1))
    if(q > 0) for(ii in 1:(ll + 1)) sv[[ii]] <- svd(Sigma_x[, , ii], nu = q, nv = 0) 
  }
  if(is.null(q)){
    q.max <- min(50, floor(sqrt(min(n - 1, p))))
    qq <- factor.number.est(xx, q.max, ll, w)
    q <- qq$q.hat[ic.op]
    Sigma_x <- qq$Sigma_x
    sv <- qq$sv
  }
  
  Sigma_c <- Sigma_x * 0
  if(q >= 1){
    for(ii in 1:(ll + 1)){
      Sigma_c[,, ii] <- sv[[ii]]$u[, 1:q, drop = FALSE] %*% diag(sv[[ii]]$d[1:q], q) %*% Conj(t(sv[[ii]]$u[, 1:q, drop = FALSE]))
      if(ii > 1) Sigma_c[,, 2 * ll + 1 - (ii - 1) + 1] <- Conj(Sigma_c[,, ii])
    }
  }
  out <- list(q = q, hl = qq, Sigma_c = Sigma_c, Sigma_x = Sigma_x)
  return(out)
  
}

#' @keywords internal
common.check <- function(x, G, est.cp, thr, ll, q = NULL, ic.op = 5, norm.type, agg.over.freq){
  
  p <- dim(x)[1]; n <- dim(x)[2]
  mean.x <- apply(x, 1, mean)
  xx <- x - mean.x
  
  check <- rep(FALSE, length(est.cp))
  if(length(est.cp) > 0){
    null1 <- common.spec.est(xx[, 1:G], q, ic.op, ll)$Sigma_c[,, 1:(1 + ll)]
    null2 <- common.spec.est(xx[, 1:G], q, ic.op, ll)$Sigma_c[,, 1:(1 + ll)]
    null <- (null1 + null2)/2
    null.norm <- apply(abs(null), 3, norm, type = norm.type)
    if(norm.type %in% c('2', 'f')) null.norm <- null.norm / sqrt(p)
    
    for(ii in 1:length(est.cp)){
      cse1 <- common.spec.est(xx[, max(1, est.cp[ii] - G):est.cp[ii]], q, ic.op, ll)
      cse2 <- common.spec.est(xx[, (est.cp[ii] + 1):min(n, est.cp[ii] + G)], q, ic.op, ll)
      if(cse1$q > 0 & cse2$q > 0){
        stat0 <- apply(cse1$Sigma_c[,, 1:(1 + ll)] - cse2$Sigma_c[,, 1:(1 + ll)], 3, norm, type = '2')
        norm.stat <- stat0 / null.norm
        if(agg.over.freq == 'avg') stat <- mean(norm.stat)
        if(agg.over.freq == 'max') stat <- max(norm.stat)
        if(stat > thr) check[ii] <- TRUE
      }
    }
    return(est.cp[check])
    
  } else return(est.cp)
  
}

## idio

idio.seg <- function(x, G.seq, est.cp.common = est.cp.common, thr.const, q = NULL, ic.op = 5){
  
  p <- dim(x)[1]
  n <- dim(x)[2]
  
  mean.x <- apply(x, 1, mean)
  xx <- x - mean.x
  Gamma_i <- cp.factor.analysis(x, est.cp.common, q, ic.op)
  ll <- (dim(Gamma_i)[3] - 1)/2
    
  common.est.cp.list <- list()
  common.stat.list <- list()
  
  for(ii in 1:length(G.seq)){
    
    G <- G.seq[ii]
    thr <- thr.const * sqrt(ll * log(n * p) / G)
    
    
    ### estimating local auto cov matrix of idiosyncratic part
    
    Gamma_xi <- array(0, dim = c(p, p, n, 2 * ll + 1))
    #Gamma_xi <- array(0, dim = c(p, p, n, 2))
    x2 <- array(0, dim = c(p, p, 2 * ll + 1))
    idx.est.cp.common <- rep(c(1:(length(sort.e.cp)-1)), diff(sort.e.cp))
    
    ## tt=G
    tt <- G
    ind.L <- (tt - G + 1):tt
    tb <- table(idx.est.cp.common[ind.L])
    tb.idx <- as.numeric(names(tb))
    Gamma_chi <- array(0, dim = c(p, p, 2 * ll + 1))
    
    for(h in 0:(ll - 1)){
      #for(h in 0:1){
      ind <- (tt - G + 1 + h):tt
      x2[,, h + 1] <- x[, ind - h] %*% t(x[, ind])
      for(i in 1:length(tb)){
        Gamma_chi[,, h + 1] <- Gamma_chi[,, h + 1] + (auto.cov[,,h + 1, tb.idx[i]]*tb[i])
        if(h != 0) Gamma_chi[, , 2 * ll + 1 - h + 1] <- t(Gamma_chi[, , h + 1])
      }
      Gamma_xi[,, tt, h+1] <- (x2[,, h + 1] - Gamma_chi[,, h + 1])/G
      if(h != 0) Gamma_xi[,, tt, 2 * ll + 1 - h + 1] <- t(Gamma_xi[,, tt, h+1])
    }
    
    ## tt=(G + 1):(n - G)
    for(tt in (G + 1):(n)){
      
      ind.L <- (tt - G + 1):tt
      tb <- table(idx.est.cp.common[ind.L])
      tb.idx <- as.numeric(names(tb))
      Gamma_chi <- array(0, dim = c(p, p, 2 * ll + 1))
      
      for(h in 0:(ll - 1)){
        #for(h in 0:1){
        x2[,, h + 1] <- x2[,, h + 1] - 
          x[, tt - G, drop = FALSE]%*%t(x[, tt - G + h, drop = FALSE]) +
          x[, tt - h, drop = FALSE]%*%t(x[, tt, drop = FALSE]) 
        for(i in 1:length(tb)){
          Gamma_chi[,, h + 1] <- Gamma_chi[,, h + 1] + (auto.cov[,,h + 1, tb.idx[i]]*tb[i])
          if(h != 0) Gamma_chi[, , 2 * ll + 1 - h + 1] <- t(Gamma_chi[, , h + 1])
        }
        Gamma_xi[,, tt, h+1] <- (x2[,, h + 1] - Gamma_chi[,, h + 1])/G
        if(h != 0) Gamma_xi[,, tt, 2 * ll + 1 - h + 1] <- t(Gamma_xi[,, tt, h+1])
      }
    }
    
    ########################################################################
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
cp.factor.analysis <- function(x, est.cp.common, q = NULL, ic.op = 5){
  
  p <- dim(x)[1]
  n <- dim(x)[2]
  mean.x <- apply(x, 1, mean)
  xx <- x - mean.x
  
  brks <- c(0, sort(est.cp.common), n)
  ll <- max(1, floor(min(diff(brks))^(1/3)))
  len <- 2 * ll
  thetas <- 2 * pi * (0:len)/(len + 1)
  w <- weights.Bartlett(((-ll):ll)/ll)
  
  Gamma_i <- array(0, dim = c(p, p, 2 * ll + 1, length(brks) - 1))
  
  for(jj in 1:(length(brks) - 1)){
    int <- (brks[jj] + 1):brks[jj + 1]
    nn <- length(int)
    if(!is.null(q)){
      qq <- as.integer(q)
      Gamma_x <- Gamma_xw <- array(0, dim = c(p, p, 2 * ll + 1))
      for(h in 0:ll){
        Gamma_x[, , h + 1] <- xx[, (brks[ii] + 1):(brks[ii + 1] - h)] %*% t(xx[, (brks[ii] + 1):(brks[ii + 1] - h) + h])/nn 
        Gamma_xw[, , h + 1] <- Gamma_x[, , h + 1] * w[h + ll + 1]
        if(h != 0){
          Gamma_x[, , 2 * ll + 1 - h + 1] <- t(Gamma_x[, , h + 1])
          Gamma_xw[, , 2 * ll + 1 - h + 1] <- t(Gamma_xw[, , h + 1])
        }
      } 
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
    
    Gamma_c <- Sigma_c <- Sigma_x * 0
    if(qq >= 1){
      for(ii in 1:(ll + 1)){
        Sigma_c[,, ii] <- sv[[ii]]$u[, 1:qq, drop = FALSE] %*% diag(sv[[ii]]$d[1:qq], qq) %*% Conj(t(sv[[ii]]$u[, 1:qq, drop = FALSE]))
        if(ii > 1){
          Sigma_c[,, 2 * ll + 1 - (ii - 1) + 1] <- Conj(Sigma_c[,, ii])
        }
      }
      Gamma_c <-  (aperm(apply(Sigma_c, c(1, 2), fft, inverse = TRUE), c(2, 3, 1)) ) * (2 * pi) / (2 * ll + 1)
      Gamma_c <- Re(Gamma_c)
    } 
    Gamma_i[,,, jj] <- Gamma_x - Gamma_c
  }
  
  return(Gamma_i)
  
}

## misc

#' @keywords internal
weights.Bartlett <- function(x) 1 - abs(x)

#' @keywords internal
factor.number.est <- function(xx, q.max, ll, w){
  
  p <- dim(xx)[1]; n <- dim(xx)[2]
  
  p.seq <- floor(3*p/4 + (1:10) * p/40) 
  n.seq <- n - (9:0) * floor(n/20)
  const.seq <- seq(.001, 2, by = .01)
  IC <- array(0, dim = c(q.max + 1, length(const.seq), 10, 2 * 3))
  
  Gamma_x <- Gamma_xw <- array(0, dim = c(p, p, 2 * ll + 1))
  
  for(kk in 1:10){
    nn <- n.seq[kk]
    pp <- p.seq[kk]
    pen <- c((1/ll^2 + sqrt(ll/nn) + 1/pp) * log(min(pp, ll^2, sqrt(nn/ll))), 
             1/sqrt(min(pp, ll^2, sqrt(nn/ll))), 
             1/min(pp, ll^2, sqrt(nn/ll)) * log(min(pp, ll^2, sqrt(nn/ll))))
    
    for(h in 0:(ll - 1)){
      Gamma_x[, , h + 1] <- xx[, 1:(nn - h)] %*% t(xx[, 1:(nn - h) + h])/nn 
      Gamma_xw[, , h + 1] <- Gamma_x[, , h + 1] * w[h + ll + 1]
      if(h != 0){
        Gamma_x[, , 2 * ll + 1 - h + 1] <- t(Gamma_x[, , h + 1])
        Gamma_xw[, , 2 * ll + 1 - h + 1] <- t(Gamma_xw[, , h + 1])
      }
    } 
    Sigma_x <- aperm(apply(Gamma_xw, c(1, 2), fft), c(2, 3, 1)) / (2 * pi)  
    sv <- list(1:(ll + 1))
    
    tmp <- rep(0, q.max + 1)
    for(ii in 1:(ll + 1)){
      if(kk == length(n.seq)) nu <- q.max else nu <- 0
      sv[[ii]] <- svd(Sigma_x[1:pp, 1:pp, ii], nu = nu, nv = 0) 
      dd <- sum(sv[[ii]]$d)
      tmp[1] <- tmp[1] + dd/pp/(2 * ll + 1)
      for(jj in 1:q.max) {
        dd <- dd - sv[[ii]]$d[jj]
        tmp[jj + 1] <- tmp[jj + 1] + dd/pp/(2 * ll + 1)
      }
      for(jj in 1:length(const.seq)){
        for(pen.op in 1:3){
          IC[, jj, kk, 3 * 0 + pen.op] <- tmp + (0:q.max) * const.seq[jj] * pen[pen.op]
          IC[, jj, kk, 3 * 1 + pen.op] <- log(tmp) + (0:q.max) * const.seq[jj] * pen[pen.op]
        }
      }
    }
  }
  
  q.mat <- apply(IC, c(2, 3, 4), which.min)
  Sc <- apply(q.mat, c(1, 3), var)
  q.hat <- rep(0, 6)
  for(ii in 1:6){
    ss <- Sc[, ii]
    if(min(ss) > 0){
      q.hat[ii] <- q.mat[which.min(ss), ii] - 1
    } else{
      q.hat[ii] <- q.mat[which(ss[-length(const.seq)] != 0 & ss[-1] == 0)[1] + 1, 10, ii] - 1
    }
  }
  
  par(mfrow = c(2, 3))
  for(ii in 1:6){
    plot(const.seq, q.mat[, 10, ii] - 1, type = 'b', pch = 1, col = 2, bty = 'n', axes = FALSE, xlab = 'constant', ylab = '', main = paste('IC ', ii))
    box()
    axis(1, at = pretty(range(const.seq)))
    axis(2, at = pretty(range(q.mat[, 10, ii] - 1)), col = 2, col.ticks = 2, col.axis = 2)
    par(new = TRUE)
    plot(const.seq, Sc[, ii], col = 4, pch = 2, type = 'b', bty = 'n', axes = FALSE, xlab = '', ylab = '')
    axis(4, at = pretty(range(Sc[, ii])), col = 4, col.ticks = 4, col.axis = 4)
    legend('topright', legend = c('q', 'Sc'), col = c(2, 4), lty = c(1, 1), pch = c(1, 2), bty = 'n')
  }
  
  ls <- list(q.hat = q.hat, Gamma_x = Gamma_x, Sigma_x = Sigma_x, sv = sv)
  return(ls)
  
}

#' @keywords internal
bottom.up <- function(est.cp.list, G.seq, eta){
  
  est.cp <- c()
  if(length(est.cp.list) > 0){
    
    Gs.giving.cp <- which(unlist(lapply(est.cp, length))>0)
    
    if(length(Gs.giving.cp)==0){
      est.cp.com <- c()
    }else{
      first.G.giving.cp <- Gs.giving.cp[1]
      est.cp.com <- c(est.cp[[first.G.giving.cp]])
      if(length(Gs.giving.cp)>1){
        
        for(i in Gs.giving.cp[-1]){
          cp.idx <- rep(NA, length(est.cp[[i]]))
          for(j in 1:length(cp.idx)){
            cp.idx[j] <- ifelse( min(abs(est.cp[[i]][j] - est.cp.com)) <= round(G.seq[i]* eta), FALSE, TRUE)
          }
          est.cp.com <- c(est.cp.com, est.cp[[i]][cp.idx])
        }
      }
    }
  }
  
  return(est.cp.com)
}

#' @keywords internal
finding.dH <- function(chp, true.chp, T){
  
  est.pnts <- sort(unique(c(0, chp, T)))
  true.pnts <- sort(unique(c(0, true.chp, T)))
  
  d <- abs(matrix(rep(est.pnts, length(true.pnts)), nrow=length(est.pnts))-
             matrix(rep(true.pnts, length(est.pnts)), nrow=length(est.pnts), byrow=T))
  
  D1 <- max(apply(d, 2, min)) * 100 / T
  D2 <- max(apply(d, 1, min)) * 100 / T
  
  dH <- mean((abs(D1-D2) + D1 + D2)/2)
  
  return(dH)
  
}



