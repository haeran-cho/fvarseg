#' @title Simulate a piecewise stationary factor-adjusted VAR process
#' @importFrom fnets var.to.vma
#' @export
sim.data <- function(n, p, q = 2,  
                      cp.common = round(n * 1:3 / 4), den.common = .5, type.common = c('ma', 'ar'), ma.order = 2,
                      cp.idio = round(n * c(3, 5) / 8), d = 1, size.idio = 1, 
                      do.scale = TRUE, seed = NULL){
  
  type.common <- match.arg(type.common, c('ma', 'ar'))
  burnin <- 100
  if(!is.null(seed)) set.seed(seed)
  
  # common component
  if(q > 0){
    if(type.common == 'ma'){
      
      brks <- c(0, c(cp.common, n) + burnin)
      u <- sqrt(rep(c(1, .5, 1.5), ceiling(q / 3))[1:q]) * matrix(rnorm(q * (n + burnin)), nrow = q)
      chi <- matrix(0, nrow = p, ncol = n + burnin)
      for(jj in 1:(ma.order + 1)){
        bb <- matrix(rnorm(p * q), nrow = p) 
        chi[, jj:brks[2]] <- chi[, jj:brks[2]] + bb %*% u[, jj:brks[2] - jj + 1]
        
        if(length(cp.common) >= 1){  
          for(kk in 1:length(cp.common)){
            cp.ind <- sample(p, floor(den.common * p))
            bb[cp.ind, ] <- matrix(rnorm(length(cp.ind) * q), ncol = q)
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
      u <- matrix(rnorm((n + trunc.lags) * q), ncol = q)
      
      chi <- matrix(0, nrow = p, ncol = n)
      for(kk in 0:length(cp.common)){
        tmp <- matrix(0, nrow = p, ncol = n)
        if(kk == 0){
          a <- matrix(runif(p * q, -1, 1), ncol = q)
          alpha <- matrix(runif(p * q, -.8, .8), ncol = q)
        } else if(kk >= 1){
          cp.ind <- sample(p, floor(den.common * p))
          a[cp.ind, ] <- - a[cp.ind, ]
          alpha[cp.ind, ] <- - alpha[cp.ind, ]
        }
        for(ii in 1:p){
          for(jj in 1:q){
            coeffs <- alpha[ii, jj] * as.numeric(fnets:::var.to.vma(as.matrix(a[ii, jj]), trunc.lags))
            for(tt in 1:n) tmp[ii, tt] <- tmp[ii, tt] + coeffs %*% u[(tt + trunc.lags):tt, jj]
          }
        }
        chi[, (brks[kk + 1] + 1):brks[kk + 2]] <- tmp[, (brks[kk + 1] + 1):brks[kk + 2]]
      }
    }
    
  } else chi <- matrix(0, nrow = p, ncol = n)
  
  ## idio component
  
  prob <- 1/p
  vep <- matrix(rnorm((n + burnin) * p), nrow = p)
  
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
    if(do.scale) chi <- chi/apply(chi, 1, sd) * apply(xi, 1, sd) 
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
hl.factor.number <- function(x, q.max, ll, w, do.plot = FALSE, center = TRUE){

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

  if(do.plot){
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
