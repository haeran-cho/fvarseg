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
  
  len <- 2 * ll
  thetas <- 2 * pi * (0:len)/(len + 1)
  w <- bartlett.weights(((-ll):ll)/ll)
  
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
  
  w <- bartlett.weights(((-ll):ll)/ll)
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
  
  w <- bartlett.weights(((-ll):ll)/ll)
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
  w <- bartlett.weights(((-ll):ll)/ll)
  
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