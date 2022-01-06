common.seg <- function(x, G.seq, thr.const, tt.by = round(log(dim(x)[2])), q = NULL, ic.op = 5, demean = TRUE,
                       norm.type = c('f', '2', 'm'), ## possibly redundant if we're to do self-normalisation
                       agg.over.freq = c('avg', 'max'), eta = .5, do.check = FALSE){
  
  p <- dim(x)[1]
  n <- dim(x)[2]
  
  if(demean) mean.x <- apply(x, 1, mean) else mean.x <- rep(0, p)
  xx <- x - mean.x
  
  common.list <- list()
  ll.seq <- c()

  for(ii in 1:length(G.seq)){
    G <- G.seq[ii]
    ll <- max(1, min(floor(4 * (G/log(G))^(1/3)), floor(n/(2 * log(n)))))
    ll.seq <- c(ll.seq, ll)
    thr <- thr.const * p * max(sqrt(ll * log(n)/G), 1/ll, 1/p)
    
    common.list[[ii]] <- cts <- common.two.step(xx, G, thr, ll, tt.by, norm.type, agg.over.freq)
    common.list[[ii]]$G <- G
    common.list[[ii]]$ll <- ll
    common.list[[ii]]$thr <- thr
    
    est.cp <- common.search.cp(cts, thr, G, eta)
    if(do.check) est.cp <- common.check(xx, G, est.cp, thr, ll, q, ic.op, norm.type, agg.over.freq)
    common.list[[ii]]$cp <- est.cp
    # matplot(cts$norm.stat, type = 'l'); abline(v = cp.common, lty = 2, col = 2, lwd = 2); abline(v = cp.idio, lty = 3, col = 6); abline(v = est.cp, col = 4, lty = 3); abline(h = thr, col = 3); lines(cts$stat, col = 4, lwd = 2)
    
  }
  
  est.cp <- bottom.up(common.list, eta)
  
  out <- list(est.cp = est.cp, G.seq = G.seq, ll.seq = ll.seq,
              est.cp.list = common.list, mean.x = mean.x)
  return(out)
  
}

#' @keywords internal
common.two.step <- function(xx, G, thr, ll, tt.by, norm.type, agg.over.freq){  
  
  p <- dim(xx)[1]
  n <- dim(xx)[2]
  
  len <- 2 * ll
  thetas <- 2 * pi * (0:len)/(len + 1)
  w <- Bartlett.weights(((-ll):ll)/ll)
  
  # null <- x1 <- array(0, dim = c(p, p, ll + 1))
  # for(h in 0:ll){
  #   ind.l <- (1 + h):G
  #   ind.r <- (n - G + 1 + h):n
  #   x1[,, h + 1] <- xx[, ind.l - h] %*% t(xx[, ind.l]) + xx[, ind.r - h] %*% t(xx[, ind.r])
  #   x1[,, h + 1] <- x1[,, h + 1]/(2 * G)
  #   for(ii in 1:(ll + 1)) {
  #     null[,, ii] <- null[,, ii] + 
  #       x1[,, h + 1] * w[h + ll + 1] * exp(-(0 + 1i) * h * thetas[ii]) +
  #       (h > 0) * t(x1[,, h + 1]) * w[h + ll + 1] * exp((0 + 1i) * h * thetas[ii])
  #   }
  # }
  # null <- null/(2 * pi)
  # null.norm <- apply(abs(null), 3, norm, type = norm.type)
  # if(norm.type %in% c('f', '2')) null.norm <- null.norm / sqrt(p)
  
  tt.seq <- round(seq(G, n - G, by = tt.by))
  stat0 <- common.stat.by(xx, G, ll, tt.seq)
  
  null.norm <- stat0[tt.seq[1], ]
  
  norm.stat <- abs(t(t(stat0) / null.norm))
  if(agg.over.freq == 'avg') stat <- apply(norm.stat, 1, mean)
  if(agg.over.freq == 'max') stat <- apply(norm.stat, 1, max)
  
  # matplot(tt.seq, norm.stat[tt.seq, ], type = 'l'); abline(v = cp.common, lty = 2, col = 2); abline(h = thr, col = 3); lines(tt.seq, stat[tt.seq], col = 4, lwd = 2)
  
  # tt.list <- common.tt.list(norm.stat, G, thr, tt.seq, tt.by)
  if(tt.by > 1){
    tt.list <- common.tt.list(stat, G, thr, tt.seq, tt.by)
    if(length(tt.list) > 0){
      for(ii in 1:length(tt.list)){
        s <- min(tt.list[[ii]]); e <- max(tt.list[[ii]])
        stat0 <- common.stat.seq(xx, G, ll, s, e, stat0)
      }
      norm.stat <- abs(t(t(stat0) / null.norm))

      if(agg.over.freq == 'avg') stat <- apply(norm.stat, 1, mean)
      if(agg.over.freq == 'max') stat <- apply(norm.stat, 1, max)
    } 
  }
  
  ls <- list(norm.stat = norm.stat, stat = stat, null.norm = null.norm)
  return(ls)
  
}

#' @keywords internal
common.search.cp <- function(cts, thr, G, eta = .5, rule = c('max', 'over')[1]){
  
  n <- length(cts$stat)
  est.cp <- c()
  new.stat <- cts$stat * (cts$stat > thr)
  survived <- new.stat > 0
  
  while(sum(survived) > 0){
    mv <- max(new.stat)
    theta.hat <- min(which(new.stat == mv))
    int <- max(1, theta.hat - round(eta * G) + 1):min(theta.hat + round(eta * G), n)
    if(rule == 'max') if(new.stat[theta.hat] >= max(new.stat[int])) est.cp <- c(est.cp, theta.hat)
    if(rule == 'over') if(sum(new.stat[int] < thr) == 0) est.cp <- c(est.cp, theta.hat)
    int <- max(1, theta.hat - G + 1):min(theta.hat + G, n)
    survived[int] <- FALSE
    new.stat[!survived] <- 0
  }
  
  est.cp
  
}

#' #' @keywords internal
#' literal implementation of alg 1
#' common.search.cp <- function(cts, thr, G, eta = .5){
#'   
#'   n <- length(cts$stat)
#'   est.cp <- c()
#'   new.norm.stat <- cts$norm.stat * (cts$norm.stat > thr)
#'   survived <- apply(new.norm.stat > 0, 1, max) == 1
#'   new.stat <- cts$stat
#'   
#'   while(sum(survived) > 0){
#'     agg.stat <- apply(new.norm.stat, 1, max)
#'     mv <- max(agg.stat)
#'     theta.hat0 <- min(which(agg.stat == mv))
#'     lstar <- which.max(new.norm.stat[theta.hat0, ])
#'     int <- max(1, theta.hat0 - round(eta * G) + 1):min(theta.hat0 + round(eta * G), n)
#'     if(new.norm.stat[theta.hat0, lstar] >= max(new.norm.stat[int, lstar])){
#'       theta.hat <- int[1] + which.max(new.stat[int]) - 1
#'       est.cp <- c(est.cp, theta.hat)
#'     }
#'     int <- max(1, theta.hat - G + 1):min(theta.hat + G, n)
#'     survived[int] <- FALSE
#'     new.norm.stat[!survived, ] <- 0
#'   }
#'   
#'   est.cp
#'   
#' }

#' @keywords internal
common.tt.list <- function(eval, G, thr, tt.seq, tt.by){
  
  if(!is.matrix(eval)) eval <- matrix(eval, ncol = 1)
  n <- dim(eval)[1]
  which.survived <- which(apply(eval[tt.seq,, drop = FALSE] > thr, 1, max) == 1)
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
common.stat.by <- function(xx, G, ll, tt.seq){  
  
  p <- dim(xx)[1]
  
  w <- Bartlett.weights(((-ll):ll)/ll)
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
      x2[,, h + 1] <- xx[, ind.l - h] %*% t(xx[, ind.l]) - xx[, ind.r - h] %*% t(xx[, ind.r])
      
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
common.stat.seq <- function(xx, G, ll, s, e, stat0){  
  
  p <- dim(xx)[1]
  
  w <- Bartlett.weights(((-ll):ll)/ll)
  len <- 2 * ll
  thetas <- 2 * pi * (0:len)/(len + 1)
  
  x2 <- eval <- array(0, dim = c(p, p, ll + 1)) 
  ttt <- 1; tt <- s
  for(h in 0:ll){
    ind.l <- (tt - G + 1 + h):tt 
    ind.r <- (tt + 1 + h):(tt + G) 
    x2[,, h + 1] <- xx[, ind.l - h] %*% t(xx[, ind.l]) - xx[, ind.r - h] %*% t(xx[, ind.r])
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
        xx[, tt - G, drop = FALSE]%*%t(xx[, tt - G + h, drop = FALSE]) +
        xx[, tt - h, drop = FALSE]%*%t(xx[, tt, drop = FALSE]) +
        xx[, tt, drop = FALSE]%*%t(xx[, tt + h, drop = FALSE]) -
        xx[, tt + G - h, drop = FALSE]%*%t(xx[, tt + G, drop = FALSE])
      
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
  nn <- dim(xx)[2]
  
  len <- 2 * ll
  thetas <- 2 * pi * (0:len)/(len + 1)
  w <- Bartlett.weights(((-ll):ll)/ll)
  
  if(!is.null(q)){
    qq <- NA
    q <- as.integer(q)
    ax <- acv.x(xx, ll, w)
    Gamma_x <- ax$Gamma_x
    Gamma_xw <- ax$Gamma_xw
    Sigma_x <- aperm(apply(Gamma_xw, c(1, 2), fft), c(2, 3, 1)) / (2 * pi)  
    sv <- list(1:(ll + 1))
    if(q > 0) for(ii in 1:(ll + 1)) sv[[ii]] <- svd(Sigma_x[, , ii], nu = q, nv = 0) 
  }
  if(is.null(q)){
    q.max <- min(50, floor(sqrt(min(nn - 1, p))))
    qq <- hl.factor.number(xx, q.max, ll, w, do.plot = FALSE, center = FALSE)
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
common.check <- function(xx, G, est.cp, thr, ll, q = NULL, ic.op = 5, 
                         norm.type, # redundant?
                         agg.over.freq, null.norm = NULL){
  
  p <- dim(xx)[1]; n <- dim(xx)[2]
  
  check <- rep(FALSE, length(est.cp))
  if(length(est.cp) > 0){
    if(is.null(null.norm)){
      null.norm <- rep(0, ll + 1)
      null <- common.spec.est(xx[, 1:G], q, ic.op, ll)$Sigma_c[,, 1:(1 + ll)] - 
        common.spec.est(xx[, 1:G + G], q, ic.op, ll)$Sigma_c[,, 1:(1 + ll)]
      for(ii in 1:(ll + 1)) null.norm[ii] <- norm(null[,, ii], type = '2')
      # null2 <- common.spec.est(xx[, 1:G], q, ic.op, ll)$Sigma_c[,, 1:(1 + ll)]
      # # null <- (null1 + null2)/2
      # null.norm <- apply(cbind(apply(abs(null1), 3, norm, type = norm.type),
      #                        apply(abs(null2), 3, norm, type = norm.type)), 1, max)
      # if(norm.type %in% c('2', 'f')) null.norm <- null.norm / sqrt(p)
    }
    for(jj in 1:length(est.cp)){
      cse1 <- common.spec.est(xx[, max(1, est.cp[jj] - G):est.cp[jj]], q, ic.op, ll)
      cse2 <- common.spec.est(xx[, (est.cp[jj] + 1):min(n, est.cp[jj] + G)], q, ic.op, ll)
      if(cse1$q > 0 || cse2$q > 0){
        stat0 <- apply(cse1$Sigma_c[,, 1:(1 + ll)] - cse2$Sigma_c[,, 1:(1 + ll)], 3, norm, type = '2')
        norm.stat <- stat0 / null.norm
        if(agg.over.freq == 'avg') stat <- mean(norm.stat)
        if(agg.over.freq == 'max') stat <- max(norm.stat)
        if(stat > thr) check[jj] <- TRUE
      }
    }
    return(est.cp[check])
    
  } else return(est.cp)
  
}
