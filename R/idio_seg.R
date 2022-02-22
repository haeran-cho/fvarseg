#' @title Segment idiosyncratic VAR process
#' @description 
#' @details See Algorithm 2 of Cho, Eckley, Fearnhead and Maeng (2022) for further details.
#' @param x input time series matrix, with each row representing a variable
#' @param center whether to de-mean the input \code{x} row-wise
#' @param common.out output from \link[fvarseg]{common.seg}; if \code{common.out = NULL}, \code{x} is regarded as a piecewise stationary VAR process
#' @param q an integer specifying the number of factors. If \code{q = NULL}, the factor number is estimated by an information criterion-based approach of Hallin and Liška (2007) for each segment
#' @param d an integer specifying the VAR order
#' @param G.seq a sequence of integers [TO DO]
#' @param thr 
#' @param eta 
#' @param cv.args a list specifying the tuning parameters required for Dantzig selector tuning parameter selection via cross-validation. It contains:
#' \itemize{
#'    \item{\code{n.folds}}{ number of folds}
#'    \item{\code{path.length}}{ number of regularisation parameter values to consider; a sequence is generated automatically based in this value}
#'    \item{\code{do.cv}}{ if \code{do.cv = FALSE}, a fixed value is selected from a sequence of 10 values chosen in a data-driven way}
#' }
#' @return a list containing the following fields:
#' \item{est.cp}{ }
#' \item{G.seq}{ }
#' \item{thr}{ }
#' \item{est.cp.list}{ a list containing the following fields:
#' \itemize{
#' \item{\code{cp}}{ }
#' \item{\code{G}}{ }
#' \item{\code{stat}}{ }
#' \item{\code{check.cp}}{ }
#' }}
#' \item{mean.x}{ if \code{center = TRUE}, returns a vector containing row-wise sample means of \code{x}; if \code{center = FALSE}, returns a vector of zeros}
#'
#' @importFrom fnets dyn.pca yw.cv var.dantzig 
#' @references H. Cho, I. Eckley, P. Fearnhead and H. Maeng (2022) High-dimensional time series segmentation via factor-adjusted vector autoregressive modelling. arXiv preprint arXiv: TODO
#' @references Hallin, M. & Liška, R. (2007) Determining the number of factors in the general dynamic factor model. Journal of the American Statistical Association, 102(478), 603--617.
#' @export
idio.seg <- function(x, center = TRUE, common.out = NULL, q = NULL, d = 1, 
                     G.seq = NULL, thr = NULL, eta = .5, 
                     cv.args = list(path.length = 10, n.folds = 1, do.cv = FALSE)){
  
  p <- dim(x)[1]
  n <- dim(x)[2]
  
  IDIO_INDEX <- 3
  load("data/idio_thr.RData")

  if(center) mean.x <- apply(x, 1, mean) else mean.x <- rep(0, p)
  xx <- x - mean.x
  
  # rule <- match.arg(rule, c('eta', 'epsilon'))
  rule <- 'eta'
  
  if(!is.null(common.out)){
    est.cp.common <- common.out$est.cp[, 1]
    K <- length(est.cp.common)
    if(K >= 1) est.cp.common <- sort(est.cp.common)
    ll <- floor(min(cs$G.seq[1]^(1/3), n/(2 * log(n))))
  } else {
    est.cp.common <- c()
    ll <- floor(min(c(n/10, 2 * p)^(1/3), n/(2 * log(n))))
  }
  brks <- c(0, est.cp.common, n)
  idx <- rep(c(1:(length(brks) - 1)), diff(brks))
  lll <- max(1, 4 * floor((min(diff(brks))/log(min(diff(brks))))^(1/3)))
  ll <- max(1, d, min(ll, lll))
  pcfa <- post.cp.fa(xx, est.cp.common, q, 5, lll)
  Gamma_c <- pcfa$Gamma_c[,, 1:(ll + 1),, drop = FALSE]
  
  if(is.null(G.seq)){
    if(K >= 1 | sum(pcfa$q.seq > 0)){
      G.seq <- round(seq(2.5 * p, n / min(4, n/(3 * p)), length.out = 4))
    } else G.seq <- round(seq(2 * p, n / min(5, n/(2 * p)), length.out = 4))
  } else G.seq <- round(sort(G.seq, decreasing = FALSE))
  
  if(is.null(thr) | length(thr) != length(G.seq)){
    thr <- c()
    if(K >= 1 | sum(pcfa$q.seq > 0)){
      for(ii in 1:length(G.seq)) thr <- c(thr, exp(predict(idio.fit.list[[IDIO_INDEX]], list(n = n, p = p, G = G.seq[ii]))))
    } else thr <- rep(1, length(G.seq))
  }
  
  idio.list <- list()
  for(ii in 1:length(G.seq)){
    
    G <- G.seq[ii]
    vv <- G
    stat <- rep(0, n)
    check.cp <- est.cp <- c()
    do.beta <- TRUE
    
    while(vv <= n - G){
      int <- (vv - G + 1):vv
      if(do.beta){
        dpca <- fnets:::dyn.pca(xx[, int], q = q, ic.op = 5, mm = lll)
        if(cv.args$do.cv){
          ycv <- fnets:::yw.cv(xx[, int], method = 'ds', var.order = d, 
                               n.folds = cv.args$n.folds, path.length = cv.args$path.length,
                               q = dpca$q, do.plot = FALSE)
          lambda <- ycv$lambda
        } else{
          lambda.max <- max(abs(xx[, int] %*% t(xx[, int])/G))
          lambda <- round(exp(seq(log(lambda.max), log(lambda.max * .0001), length.out = 10)),
                          digits = 10)[4]
        }
        mg <- make.gg(dpca$acv$Gamma_i, d)
        ive <- fnets:::var.dantzig(mg$GG, mg$gg, lambda)
        beta <- ive$beta
        
        dpca.l <- fnets:::dyn.pca(xx[, int[1:round(G/2)]], q = dpca$q, ic.op = 5, mm = ll)
        dpca.r <- fnets:::dyn.pca(xx[, int[-(1:round(G/2))]], q = dpca$q, ic.op = 5, mm = ll)
        dgi <- dpca.l$acv$Gamma_i[,, 1:(ll + 1)] - dpca.r$acv$Gamma_i[,, 1:(ll + 1)]
        null.norm <- max(abs(dgi))
     }
      
      diff.Gamma_x <- acv.x(xx[, int, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1), drop = FALSE] -
        acv.x(xx[, int + G, drop = FALSE], ll)$Gamma_x[,, 1:(ll + 1), drop = FALSE]
      diff.Gamma_c <- diff.Gamma_x * 0
      if(K > 0){
        common.weights <- tabulate(idx[int], nbins = K + 1) -
          tabulate(idx[int + G], nbins = K + 1)
        for(kk in 1:(K + 1)) diff.Gamma_c <- diff.Gamma_c + common.weights[kk] / G * Gamma_c[,,, kk]
      }
      
      mgd <- make.gg(diff.Gamma_x - diff.Gamma_c, d)
      stat[vv] <- max(abs(mgd$GG %*% beta - mgd$gg)) / null.norm
      
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
        
        if(first & stat[tt] > thr[ii]){
          check.theta <- tt
          tt.max <- min(tt.max, check.theta + G - 1)
          check.cp <- c(check.cp, check.theta)
          first <- FALSE
        }
        
        tt <- tt + 1
      }  

      if(check.theta < tt.max){
        hat.theta <- (check.theta:tt.max)[which.max(stat[check.theta:tt.max])]
        
        int <- max(1, vv, hat.theta - 2 + 1):min(hat.theta + 2, tt.max, n)
        if(sum(stat[int] < thr[ii]) == 0){
          est.cp <- c(est.cp, hat.theta)
          vv <- hat.theta + G
          do.beta <- TRUE
        } else{
          vv <- tt.max + 1
          check.cp <- setdiff(check.cp, check.theta)
          do.beta <- FALSE
        }
      } else break
    }
    idio.list[[ii]] <- list(cp = est.cp, stat = stat, G = G, thr = thr[ii], check.cp = check.cp)
  }
  
  est.cp <- bottom.up(idio.list, eta)
  
  out <- list(est.cp = est.cp, G.seq = G.seq, thr = thr, est.cp.list = idio.list, mean.x = mean.x)
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
  
  Sigma_c <- Gamma_c <- array(0, dim = c(p, p, 2 * ll + 1, length(brks) - 1))
  q.seq <- rep(0, length(brks) - 1)
  
  for(jj in 1:(length(brks) - 1)){
    int <- (brks[jj] + 1):brks[jj + 1]
    nn <- length(int)
    if(!is.null(q)){
      qq <- as.integer(q)
      ax <- acv.x(xx[, int, drop = FALSE], ll, w)
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
      for(ii in 1:(ll + 1)){
        Sigma_c[,, ii, jj] <- sv[[ii]]$u[, 1:qq, drop = FALSE] %*% diag(sv[[ii]]$d[1:qq], qq) %*% Conj(t(sv[[ii]]$u[, 1:qq, drop = FALSE]))
        if(ii > 1){
          Sigma_c[,, 2 * ll + 1 - (ii - 1) + 1, jj] <- Conj(Sigma_c[,, ii, jj])
        }
      }
      Gamma_c[,,, jj] <- aperm(apply(Sigma_c[,,, jj], c(1, 2), fft, inverse = TRUE), c(2, 3, 1)) * (2 * pi) / (2 * ll + 1)
    } 
  }
  ls <- list(Gamma_c = Re(Gamma_c), q.seq = q.seq)
  
  return(ls)
  
}
