#' @title Simulate a piecewise stationary factor-adjusted VAR process
#' @description Generate time series used in the simulation studies of Cho, Eckley, Fearnhead and Maeng (2022) for further details.
#' 
#' @param n sample size
#' @param p number of variables
#' @param q number of dynamic factors
#' @param d VAR order
#' @param cp.common a vector specifying the re-scaled locations of the change points between \code{0} and \code{1} in the common component; possible to set \code{cp.common = c()} (no change point) 
#' @param den.common a value between \code{0} and \code{1} specifying the cross-sectional density of each change point
#' @param type.common if \code{type.common = 'ma'}, factors are loaded as innovations of a moving average process with order \code{ma.order}; if \code{type.common = 'ar'}, factors are loaded as innovations of an autoregressive process of order \code{1}
#' @param ma.order order of the factor-driven moving average process; used when \code{type.common = 'ma'}
#' @param cp.idio a vector specifying the re-scaled locations of the change points between \code{0} and \code{1} in the idiosyncratic component; possible to set \code{cp.idio = c()} (no change point) 
#' @param size.idio at each change point, each of VAR parameter matrices has its sign changed and is multiplied by \code{size.idio}
#' @param do.scale if \code{do.scale = TRUE}, each variable of the common component is scaled to have the same sample variance as the corresponding idiosyncratic variable
#' @param seed an integer setting the seed of the random number generator
#' @return a list containing
#' \item{x}{ generated piecewise stationary factor-adjusted vector autoregressive process }
#' \item{xi}{ generated piecewise stationary vector autoregressive process }
#' \item{A.list}{ a list containing the VAR parameter matrices over the segments }
#' \item{cp.common}{ input parameter }
#' \item{cp.idio}{ input parameter }
#' @examples
#' out <- sim.data(n = 2000, p = 50, q = 2, d = 1,
#' cp.common = 1:3/4, den.common = .5, type.common = 'ma', 
#' cp.idio = c(3, 5)/8, seed = 123)
#' @references Cho, H., Eckley, I., Fearnhead, P. & Maeng, H. (2022) High-dimensional time series segmentation via factor-adjusted vector autoregressive modelling. arXiv preprint arXiv: TODO
#' @importFrom stats rnorm runif sd
#' @export
sim.data <- function(n, p, q = 2, d = 1,
                      cp.common = c(1:3) / 4, den.common = .5, type.common = c('ma', 'ar'), ma.order = 2,
                      cp.idio = c(3, 5) / 8, size.idio = 1, 
                      do.scale = TRUE, seed = NULL){
  
  cp.common <- round(cp.common * n)
  cp.idio <- round(cp.idio * n)
  type.common <- match.arg(type.common, c('ma', 'ar'))
  burnin <- 100
  if(!is.null(seed)) set.seed(seed)
  
  # common component
  if(q > 0){
    if(type.common == 'ma'){
      
      brks <- c(0, c(cp.common, n) + burnin)
      u <- sqrt(rep(c(1, .5, 1.5), ceiling(q / 3))[1:q]) * matrix(stats::rnorm(q * (n + burnin)), nrow = q)
      chi <- matrix(0, nrow = p, ncol = n + burnin)
      for(jj in 1:(ma.order + 1)){
        bb <- matrix(stats::rnorm(p * q), nrow = p) 
        chi[, jj:brks[2]] <- chi[, jj:brks[2]] + bb %*% u[, jj:brks[2] - jj + 1]
        
        if(length(cp.common) >= 1){  
          for(kk in 1:length(cp.common)){
            cp.ind <- sample(p, floor(den.common * p))
            bb[cp.ind, ] <- matrix(stats::rnorm(length(cp.ind) * q), ncol = q)
            chi[, (brks[kk + 1] + 1):brks[kk + 2]] <- chi[, (brks[kk + 1] + 1):brks[kk + 2]] + 
              bb %*% u[, (brks[kk + 1] + 1):brks[kk + 2] - jj + 1]
          }
        }
      }
      chi <- chi[, -(1:burnin)] 
    }
    
    if(type.common == 'ar'){  
      
      brks <- c(0, cp.common, n)
      trunc.lags <- min(20, round(n/log(n)))
      u <- matrix(stats::rnorm((n + trunc.lags) * q), ncol = q)
      
      chi <- matrix(0, nrow = p, ncol = n)
      for(kk in 0:length(cp.common)){
        tmp <- matrix(0, nrow = p, ncol = n)
        if(kk == 0){
          a <- matrix(stats::runif(p * q, -1, 1), ncol = q)
          alpha <- matrix(stats::runif(p * q, -.8, .8), ncol = q)
        } else if(kk >= 1){
          cp.ind <- sample(p, floor(den.common * p))
          a[cp.ind, ] <- - a[cp.ind, ]
          alpha[cp.ind, ] <- - alpha[cp.ind, ]
        }
        for(ii in 1:p){
          for(jj in 1:q){
            coeffs <- alpha[ii, jj] * as.numeric(var.to.vma(as.matrix(a[ii, jj]), trunc.lags))
            for(tt in 1:n) tmp[ii, tt] <- tmp[ii, tt] + coeffs %*% u[(tt + trunc.lags):tt, jj]
          }
        }
        chi[, (brks[kk + 1] + 1):brks[kk + 2]] <- tmp[, (brks[kk + 1] + 1):brks[kk + 2]]
      }
    }
    
  } else chi <- matrix(0, nrow = p, ncol = n)
  
  ## idio component
  
  prob <- 1/p
  vep <- matrix(stats::rnorm((n + burnin) * p), nrow = p)
  
  A.list <- list()
  brks <- c(0, cp.idio, n)
  xi <- matrix(0, nrow = p, ncol = n)
  for(kk in 0:length(cp.idio)){
    tmp <- vep
    if(kk == 0){
      A1 <- A2 <- matrix(0, p, p)
      index <- sample(c(0, 1), p^2, TRUE, prob = c(1 - prob, prob))
      A1[which(index == 1)] <- .4
      A1 <- A1 / norm(A1, "2")
      if(d == 2){
        A1 <- A1 * .5
        index <- sample(c(0, 1), p^2, TRUE, prob = c(1 - prob, prob))
        A2[which(index == 1)] <- .4
        A2 <- A2 / norm(A2, "2") * .5
      }
    } else if(kk >= 1){
      A1 <- - size.idio * A1; A2 <- - size.idio * A2
    }
    for(tt in 3:(n + burnin)) tmp[, tt] <- tmp[, tt] + A1 %*% tmp[, tt - 1] + A2 %*% tmp[, tt - 2]
    tmp <- tmp[, -(1:burnin)]
    xi[, (brks[kk + 1] + 1):brks[kk + 2]] <- tmp[, (brks[kk + 1] + 1):brks[kk + 2]]
    if(d == 1) A <- A1 else if(d == 2) A <- cbind(A1, A2)
    A.list[[kk + 1]] <- A
  }
  
  if(q > 0){
    if(do.scale) chi <- chi/apply(chi, 1, stats::sd) * apply(xi, 1, stats::sd) 
    x <- chi + xi
  } else x <- xi
  
  out <- list(x = x, xi = xi, A.list = A.list, cp.common = cp.common, cp.idio = cp.idio)
  return(out)
}

#' misc functions

#' @title Bartlett weights
#' @keywords internal
Bartlett.weights <- function(x) 1 - abs(x)

#' @keywords internal 
acv.x <- function(xx, ll, w = NULL){
  
  p <- dim(xx)[1]
  nn <- dim(xx)[2]
  Gamma_x <- array(0, dim = c(p, p, 2 * ll + 1))
  if(!is.null(w)) Gamma_xw <- Gamma_x
  for(h in 0:ll){
    Gamma_x[, , h + 1] <- xx[, 1:(nn - h)] %*% t(xx[, 1:(nn - h) + h])/nn 
    if(!is.null(w)) Gamma_xw[, , h + 1] <- Gamma_x[, , h + 1] * w[h + ll + 1]
    if(h != 0){
      Gamma_x[, , 2 * ll + 1 - h + 1] <- t(Gamma_x[, , h + 1])
      if(!is.null(w)) Gamma_xw[, , 2 * ll + 1 - h + 1] <- t(Gamma_xw[, , h + 1])
    }
  } 
  if(is.null(w)) ls <- list(Gamma_x = Gamma_x)
  if(!is.null(w)) ls <- list(Gamma_x = Gamma_x, Gamma_xw = Gamma_xw)
  
  return(ls)
  
}

#' @keywords internal
make.gg <- function(acv, d){
  
  p <- dim(acv)[1]
  gg <- matrix(0, nrow = p * d, ncol = p)
  GG <- matrix(0, p * d, p * d)
  for(ll in 1:d){
    gg[(ll - 1) * p + 1:p, ] <- acv[,, ll + 1]
    for(lll in ll:d){
      GG[(ll - 1) * p + 1:p, (lll - 1) * p + 1:p] <- t(acv[,, 1 + lll - ll])
      GG[(lll - 1) * p + 1:p, (ll - 1) * p + 1:p] <- acv[,, 1 + lll - ll]
    }
  }
  out <- list(gg = gg, GG = GG)
  return(out)
  
}

#' @keywords internal
var.to.vma <- function(A, trunc.lags){
  
  d <- dim(A)[1]; s <- dim(A)[2]
  l <- s/d
  B <- array(0, dim = c(d, d, trunc.lags + 1))
  B[,, 1] <- diag(1, d)
  for(ii in 1:trunc.lags){
    for(jj in 1:min(ii, l)) B[,, ii + 1] <- B[,, ii + 1] + B[,, ii - jj + 1] %*% A[, (jj - 1) * d + 1:d]
  }
  B
  
}

#' @title Dynamic PCA
#' @importFrom stats fft
#' @keywords internal
dyn.pca <- function(xx, q = NULL, ic.op = 5, kern.const = 4, mm = NULL){
  
  p <- dim(xx)[1]
  n <- dim(xx)[2]
  
  if(is.null(mm)) mm <- min(max(1, kern.const * floor((n/log(n))^(1/3))), floor(n/4) - 1) else mm <- min(max(mm, 1, kern.const * floor((n/log(n))^(1/3))), floor(n/4) - 1)
  len <- 2 * mm
  w <- Bartlett.weights(((-mm):mm)/mm)
  
  # dynamic pca
  
  if(!is.null(q)){
    q <- as.integer(q); hl <- NA
    Gamma_x <- Gamma_xw <- array(0, dim = c(p, p, 2 * mm + 1))
    for(h in 0:(mm - 1)){
      Gamma_x[, , h + 1] <- xx[, 1:(n - h)] %*% t(xx[, 1:(n - h) + h])/n
      Gamma_xw[, , h + 1] <- Gamma_x[, , h + 1] * w[h + mm + 1]
      if(h != 0){
        Gamma_x[, , 2 * mm + 1 - h + 1] <- t(Gamma_x[, , h + 1])
        Gamma_xw[, , 2 * mm + 1 - h + 1] <- t(Gamma_xw[, , h + 1])
      }
    }
    Sigma_x <- aperm(apply(Gamma_xw, c(1, 2), fft), c(2, 3, 1)) / (2 * pi)
    sv <- list(1:(mm + 1))
    if(q > 0) for(ii in 1:(mm + 1)) sv[[ii]] <- svd(Sigma_x[, , ii], nu = q, nv = 0)
  }
  if(is.null(q)){
    q.max <- min(50, floor(sqrt(min(n - 1, p))))
    hl <- hl.factor.number(xx, q.max, mm, w, center = FALSE)
    q <- hl$q.hat[ic.op]
    Gamma_x <- hl$Gamma_x
    Sigma_x <- hl$Sigma_x
    sv <- hl$sv
  }
  
  Gamma_c <- Gamma_i <- Sigma_c <- Sigma_i <- Sigma_x * 0
  if(q >= 1){
    for(ii in 1:(mm + 1)){
      Sigma_c[,, ii] <- sv[[ii]]$u[, 1:q, drop = FALSE] %*% diag(sv[[ii]]$d[1:q], q) %*% Conj(t(sv[[ii]]$u[, 1:q, drop = FALSE]))
      if(ii > 1){
        Sigma_c[,, 2 * mm + 1 - (ii - 1) + 1] <- Conj(Sigma_c[,, ii])
      }
    }
    Gamma_c <- aperm(apply(Sigma_c, c(1, 2), fft, inverse = TRUE), c(2, 3, 1)) * (2 * pi) / (2 * mm + 1)
    Gamma_c <- Re(Gamma_c)
  }
  Sigma_i <- Sigma_x - Sigma_c
  Gamma_i <- Gamma_x - Gamma_c
  
  spec <- list(Sigma_x = Sigma_x, Sigma_c = Sigma_c, Sigma_i = Sigma_i)
  acv <- list(Gamma_x = Gamma_x, Gamma_c = Re(Gamma_c), Gamma_i = Re(Gamma_i))
  
  out <- list(q = q, hl = hl, spec = spec, acv = acv, kern.const = kern.const)
  return(out)
  
}

#' @keywords internal
#' @importFrom stats var
hl.factor.number <- function(x, q.max, ll, w, center = TRUE){

  p <- dim(x)[1]; n <- dim(x)[2]
  if(center) mean.x <- apply(x, 1, mean) else mean.x <- rep(0, p)
  xx <- x - mean.x

  p.seq <- floor(3*p/4 + (1:10) * p/40)
  n.seq <- n - (9:0) * floor(n/20)
  const.seq <- seq(.001, 2, by = .01)
  IC <- array(0, dim = c(q.max + 1, length(const.seq), 10, 2 * 3))

  for(kk in 1:10){
    nn <- n.seq[kk]
    pp <- p.seq[kk]
    pen <- c((1/ll^2 + sqrt(ll/nn) + 1/pp) * log(min(pp, ll^2, sqrt(nn/ll))),
             1/sqrt(min(pp, ll^2, sqrt(nn/ll))),
             1/min(pp, ll^2, sqrt(nn/ll)) * log(min(pp, ll^2, sqrt(nn/ll))))

    ax <- acv.x(xx[, 1:nn, drop =FALSE], ll, w)
    Gamma_x <- ax$Gamma_x
    Gamma_xw <- ax$Gamma_xw
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
  Sc <- apply(q.mat, c(1, 3), stats::var)
  q.hat <- rep(0, 6)
  for(ii in 1:6){
    ss <- Sc[, ii]
    if(min(ss) > 0){
      q.hat[ii] <- q.mat[which.min(ss), ii] - 1
    } else{
      q.hat[ii] <- q.mat[which(ss[-length(const.seq)] != 0 & ss[-1] == 0)[1] + 1, 10, ii] - 1
    }
  }

  ls <- list(q.hat = q.hat, Gamma_x = Gamma_x, Sigma_x = Sigma_x, sv = sv)
  return(ls)

}

#' @keywords internal
bottom.up <- function(ls, eta){
  
  est.cp <- matrix(NA, nrow = 0, ncol = 2)
  dimnames(est.cp)[[2]] <- c('cp', 'G')
  
  kk.seq <- seq(1, length(ls))
  while(length(kk.seq) > 0){
    kk <- kk.seq[1]
    cand.cp <- ls[[kk]]$cp
    if(length(cand.cp) > 0){
      if(dim(est.cp)[1] == 0){
        est.cp <- rbind(est.cp, cbind(cand.cp, ls[[kk]]$G))
      } else{
        for(ii in 1:length(cand.cp)){
          if(min(abs(cand.cp[ii] - est.cp[, 1])) > ls[[kk]]$G * eta) est.cp <- 
              rbind(est.cp, c(cand.cp[ii], ls[[kk]]$G))
        }
      }
    }
    kk.seq <- setdiff(kk.seq, kk)
  }

  return(est.cp)
}
