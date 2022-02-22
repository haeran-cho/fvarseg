# TODO: remove load

#' @title Segment factor-adjusted VAR process
#' @description 
#' @details See Cho, Eckley, Fearnhead and Maeng (2022) for further details.
#' @param x input time series matrix, with each row representing a variable
#' @param demean whether to de-mean the input \code{x} row-wise
#' @param q an integer specifying the number of factors. If \code{q = NULL}, the factor number is estimated by an information criterion-based approach of Hallin and Liška (2007) for each segment
#' @param d an integer specifying the VAR order
#' @param eta 
#' @param common.args a list specifying the tuning parameters required for segmenting the factor-driven common component, see also \link[fvarseg]{common.seg}. It contains 
#' \itemize{
#'    \item{\code{G.seg}}{ }
#'    \item{\code{thr}}{ }
#'    \item{\code{tt.by}}{ }
#' }
#' @param idio.args a list specifying the tuning parameters required for segmenting the idiosyncratic VAR process, see also \link[fvarseg]{idio.seg}. It contains 
#' \itemize{
#'    \item{\code{G.seg}}{ }
#'    \item{\code{thr}}{ }
#' }
#' @param cv.args a list specifying the tuning parameters required for Dantzig selector tuning parameter selection via cross-validation. It contains:
#' \itemize{
#'    \item{\code{n.folds}}{ number of folds}
#'    \item{\code{path.length}}{ number of regularisation parameter values to consider; a sequence is generated automatically based in this value}
#'    \item{\code{do.cv}}{ if \code{do.cv = FALSE}, a fixed value is selected from a sequence of 10 values chosen in a data-driven way}
#' }
#' @return a list containing the following fields:
#' \item{common.out}{ output from \link[fvarseg]{common.seg}
#' \itemize{
#' \item{\code{est.cp}}{ }
#' \item{\code{G.seq}}{ }
#' \item{\code{thr}}{ }
#' \item{\code{est.cp.list}}{ See \link[fvarseg]{common.seg} for further details }
#' }}
#' \item{idio.out}{ output from \link[fvarseg]{idio.seg}
#' \itemize{
#' \item{\code{est.cp}}{ }
#' \item{\code{G.seq}}{ }
#' \item{\code{thr}}{ }
#' \item{\code{est.cp.list}}{ See \link[fvarseg]{common.seg} for further details }
#' }}
#' \item{mean.x}{ if \code{center = TRUE}, returns a vector containing row-wise sample means of \code{x}; if \code{center = FALSE}, returns a vector of zeros}
#'
#' @importFrom quantreg predict.rq
#' @importFrom fnets dyn.pca yw.cv var.dantzig 
#' @references H. Cho, I. Eckley, P. Fearnhead and H. Maeng (2022) High-dimensional time series segmentation via factor-adjusted vector autoregressive modelling. arXiv preprint arXiv: TODO
#' @references Hallin, M. & Liška, R. (2007) Determining the number of factors in the general dynamic factor model. Journal of the American Statistical Association, 102(478), 603--617.
#' @export
fvar.seg <- function(x, center = TRUE, q = NULL, d = 1, eta = .5,
                     common.args = list(G.seq = NULL, thr = NULL, tt.by = floor(2 * log(dim(x)[2]))),
                     idio.args = list(G.seq = NULL, thr = NULL),
                     cv.args = list(path.length = 10, n.folds = 1, do.cv = FALSE)){
  
  p <- dim(x)[1]
  n <- dim(x)[2]
  
  if(center) mean.x <- apply(x, 1, mean) else mean.x <- rep(0, p)
  xx <- x - mean.x
  
  cs <- common.seg(xx, center = FALSE, G.seq = common.args$G.seq, thr = common.args$thr, tt.by = common.args$tt.by, eta = eta)
  
  is <- idio.seg(xx, center = FALSE, common.out = cs, q = q, d = d, 
                 G.seq = idio.args$G.seq, thr = idio.args$thr, eta = eta,
                 cv.args = cv.args)
  
  cs$mean.x <- NULL
  is$mean.x <- NULL
  
  out <- list(common.out = cs, idio.out = is, mean.x = mean.x)
  return(out)
  
}