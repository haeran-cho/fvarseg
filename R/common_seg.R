#' @title Segment factor-driven common component
#' @details See Algorithm 1 of Cho, Eckley, Fearnhead and Maeng (2022) for further details.
#' @param x input time series matrix, with each row representing a variable
#' @param center whether to de-mean the input \code{x} row-wise
#' @param G.seq an integer vector of bandwidth; if \code{G.seq = NULL}, a default choice \code{G.seq = round(n * 1/c(10, 8, 6, 4))} is used
#' @param thr a vector of thresholds which is of the same length as \code{G.seq}; if \code{thr = NULL}, a default choice based on numerical experiments is used
#' @param alpha used when \code{thr = NULL} for default threshold selection, which relates to the quantile used in numerical experiments; \code{alpha = 0.05} and \code{alpha = 0.1} are supported
#' @param tt.by an integer specifying the grid over which the test statistic is computed, which is \code{round(seq(G, dim(x)[2] - G, by = tt.by))} for each bandwidth \code{G}
#' @param eta a constant between \code{0} and \code{1}; each local maximiser of the test statistic within its \code{eta * G}-environment for the common component is deemed as a change point estimator. Also the bottom-up merging across the multiple bandwidths \code{G.seq} depends on this parameter
#'
#' @return a list containing the following fields:
#' \item{est.cp}{ a matrix containing the change point estimators in the first column and the finest bandwidth at which each is detected in the second column }
#' \item{G.seq}{ an integer vector of bandwidths }
#' \item{thr}{ a vector of thresholds which is of the same length as \code{G.seq} }
#' \item{est.cp.list}{ a list containing the following fields:
#' \itemize{
#' \item{\code{cp}}{ change point estimators }
#' \item{\code{G}}{ bandwidth }
#' \item{\code{ll}}{ kernel window size for spectral density estimation }
#' \item{\code{norm.stat}}{ a matrix containing test statistic values at Fourier frequencies }
#' \item{\code{stat}}{ a vector containing test statistic values across multiple frequencies }
#' }}
#' \item{mean.x}{ if \code{center = TRUE}, returns a vector containing row-wise sample means of \code{x}; if \code{center = FALSE}, returns a vector of zeros}
#'
#' @examples 
#' \dontrun{
#' out <- sim.data(n = 2000, p = 50, q = 2, d = 1,
#' cp.common = 1:3/4, den.common = .5, type.common = 'ma', 
#' cp.idio = c(3, 5)/8, seed = 123)
#' cs <- common.seg(out$x)
#' cs$est.cp
#' }
#' @importFrom stats predict.lm
#' @references Cho, H., Eckley, I., Fearnhead, P. & Maeng, H. (2022) High-dimensional time series segmentation via factor-adjusted vector autoregressive modelling. arXiv preprint arXiv:2204.02724
#' @export
common.seg <- function(x, center = TRUE, G.seq = NULL, thr = NULL, alpha = .1,
                       tt.by = floor(2 * log(dim(x)[2])), eta = .5){
  
  p <- dim(x)[1]
  n <- dim(x)[2]
  
  if(alpha == .1) COMMON_INDEX <- 2 else if(alpha == .05) COMMON_INDEX <- 3

  if(is.null(G.seq)) G.seq <- round(n * c(1/10, 1/8, 1/6, 1/4)) else G.seq <- round(sort(G.seq, decreasing = FALSE))
  if(is.null(thr) | length(thr) != length(G.seq)){
    thr <- c()
    for(ii in 1:length(G.seq)) thr <- c(thr, exp(stats::predict.lm(common.fit.list[[COMMON_INDEX]], list(n = n, p = p, G = G.seq[ii]))))
  }
  
  # rule <- match.arg(rule, c('eta', 'epsilon'))
  # agg.over.freq <- match.arg(agg.over.freq, c('avg', 'max'))
  
  rule <- 'eta'
  agg.over.freq <- 'avg'
  do.check <- FALSE
  epsilon <- 0

  if(center) mean.x <- apply(x, 1, mean) else mean.x <- rep(0, p)
  xx <- x - mean.x
  
  common.list <- list()
  
  for(ii in 1:length(G.seq)){
    G <- G.seq[ii]
    ll <- max(1, floor(min(G^(1/3), n/(2 * log(n)))))

    common.list[[ii]] <- cts <- common.two.step(xx, G, thr[ii], ll, tt.by, agg.over.freq)
    common.list[[ii]]$G <- G
    common.list[[ii]]$ll <- ll
    common.list[[ii]]$thr <- thr[ii]
    
    est.cp <- common.search.cp(cts, thr[ii], G, rule, eta, epsilon)
    if(do.check) est.cp <- common.check(xx, G, est.cp, thr[ii], ll, q = NULL, ic.op = 5, agg.over.freq)
    common.list[[ii]]$cp <- est.cp
  }
  
  est.cp <- bottom.up(common.list, eta)
  
  out <- list(est.cp = est.cp, G.seq = G.seq, thr = thr,
              est.cp.list = common.list, mean.x = mean.x)
  return(out)
  
}

#' @keywords internal
common.two.step <- function(xx, G, thr, ll, tt.by, agg.over.freq = 'avg'){  
  
  p <- dim(xx)[1]
  n <- dim(xx)[2]
  
  len <- 2 * ll
  thetas <- 2 * pi * (0:len)/(len + 1)
  w <- Bartlett.weights(((-ll):ll)/ll)

  tt.seq <- round(seq(G, n - G, by = tt.by))
  stat0 <- common.stat.by(xx, G, ll, tt.seq)
  
  null.norm <- stat0[tt.seq[1], ]
  
  norm.stat <- abs(t(t(stat0) / null.norm))
  if(agg.over.freq == 'avg') stat <- apply(norm.stat, 1, mean)
  if(agg.over.freq == 'max') stat <- apply(norm.stat, 1, max)
  
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
  
  ls <- list(norm.stat = norm.stat, stat = stat)
  return(ls)
  
}

#' @keywords internal
common.search.cp <- function(cts, thr, G, rule, eta = .5, epsilon = .1){
  
  n <- length(cts$stat)
  est.cp <- c()
  new.stat <- cts$stat * (cts$stat > thr)
  survived <- new.stat > 0
  
  while(sum(survived) > 0){
    mv <- max(new.stat)
    hat.theta <- min(which(new.stat == mv))
    if(rule == 'eta'){
      int <- max(1, hat.theta - round(eta * G) + 1):min(hat.theta + round(eta * G), n)
      if(new.stat[hat.theta] >= max(new.stat[int])) est.cp <- c(est.cp, hat.theta)
    }
    if(rule == 'epsilon'){
      int <- max(1, G, hat.theta - round(epsilon * G) + 1):min(hat.theta + round(epsilon * G), n - G, n)
      if(sum(new.stat[int] < thr) == 0) est.cp <- c(est.cp, hat.theta)
    }
    int <- max(1, hat.theta - G + 1):min(hat.theta + G, n)
    survived[int] <- FALSE
    new.stat[!survived] <- 0
  }
  
  est.cp
  
}

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
  n <- dim(xx)[2]
  
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
    qq <- hl.factor.number(xx, q.max, ll, w, center = FALSE)
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
common.check <- function(xx, G, est.cp, thr, ll, q = NULL, ic.op = 5, agg.over.freq = 'avg'){
  
  p <- dim(xx)[1]; n <- dim(xx)[2]
  
  check <- rep(FALSE, length(est.cp))
  if(length(est.cp) > 0){
    null.norm <- rep(0, ll + 1)
    null <- common.spec.est(xx[, 1:G], q, ic.op, ll)$Sigma_c[,, 1:(1 + ll)] - 
      common.spec.est(xx[, 1:G + G], q, ic.op, ll)$Sigma_c[,, 1:(1 + ll)]
    for(ii in 1:(ll + 1)) null.norm[ii] <- norm(null[,, ii], type = '2')
    
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
