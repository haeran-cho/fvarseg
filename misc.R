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
  
  # sig_vep <- toeplitz(.3^(1:p - 1)) 
  # xi <- vep <- t(mvtnorm::rmvnorm(n + burnin, sigma = sig_vep))
  # 
  # A <- matrix(0, nrow = p, ncol = p)
  # c.A <- .4
  # diag(A) <- c.A
  # A[row(A) + 1 == col(A)] <- - c.A
  # 
  # brks <- c(0, c(cp.idio, n) + burnin)
  # 
  # for(tt in (brks[1] + 1 + 1):brks[2]) xi[, tt] <- vep[, tt] + A %*% xi[, tt - 1]
  # if(length(cp.idio) >= 1){
  #   for(j in 1:length(cp.idio)){
  #     A[row(A) + 1 == col(A)] <- - size.idio * A[row(A) + 1 == col(A)]
  #     for(tt in (brks[j + 1] + 1):brks[j + 2]) xi[, tt] <- vep[, tt] + A %*% xi[, tt - 1]
  #   }
  # }
  # x <- xi <- xi[, -(1:burnin)]
  
  A.list <- list()
  brks <- c(0, cp.idio, n)
  xi <- c()
  for(tt in 1:(length(brks) - 1)){
    tmp <- fnets::sim.var(brks[tt + 1] - brks[tt], p, Gamma = diag(1, p), heavy = FALSE)
    xi <- cbind(xi, tmp$data)
    A.list[[tt]] <- tmp$A
  }
  x <- xi
  
  if(q > 0){
    chi <- chi/apply(chi, 1, sd) * apply(xi, 1, sd) 
    x <- x + chi
  }
  
  out <- list(x = x, xi = xi, A.list = A.list)
  return(out)
  
}

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
