#' @title Segment factor-adjusted VAR process
#' @description Segment high-dimensional time series using the two-stage segmentation method proposed in Cho, Eckley, Fearnhead and Maeng (2022).
#' It first detects change points from the factor-driven common component, then from the idiosyncratic VAR process.
#' @details See Cho, Eckley, Fearnhead and Maeng (2022) for further details.
#' @param x input time series matrix, with each row representing a variable
#' @param center whether to de-mean the input \code{x} row-wise
#' @param q an integer specifying the number of factors. If \code{q = NULL}, the factor number is estimated by an information criterion-based approach of Hallin and Liška (2007) for each segment
#' @param d an integer specifying the VAR order
#' @param eta a constant between \code{0} and \code{1}; each local maximiser of the test statistic within its \code{eta * G}-environment for the common component is deemed as a change point estimator. Also the bottom-up merging across the multiple bandwidths \code{G.seq} depends on this parameter
#' @param common.args a list specifying the tuning parameters required for segmenting the factor-driven common component, see also \link[fvarseg]{common.seg}. It contains 
#' \itemize{
#'    \item{\code{G.seg}}{ an integer vector of bandwidth; see \code{fvarseg}[common.seg] for the default choice when \code{G.seq = NULL} }
#'    \item{\code{thr}}{ a vector of thresholds which is of the same length as \code{G.seq}; if \code{thr = NULL}, a default choice based on simulations is used }
#'    \item{\code{tt.by}}{ an integer specifying the grid over which the test statistic is computed, which is \code{round(seq(G, dim(x)[2] - G, by = tt.by))} for each bandwidth \code{G} }
#' }
#' @param idio.args a list specifying the tuning parameters required for segmenting the idiosyncratic VAR process, see also \link[fvarseg]{idio.seg}. It contains 
#' \itemize{
#'    \item{\code{G.seg}}{ an integer vector of bandwidth; see \code{fvarseg}[idio.seg] for the default choice when \code{G.seq = NULL} }
#'    \item{\code{thr}}{ a vector of thresholds which is of the same length as \code{G.seq}; if \code{thr = NULL}, a default choice based on simulations is used }
#' }
#' @param cv.args a list specifying the tuning parameters required for Dantzig selector tuning parameter selection via cross-validation. It contains:
#' \itemize{
#'    \item{\code{n.folds}}{ number of folds}
#'    \item{\code{path.length}}{ number of regularisation parameter values to consider; a sequence is generated automatically based in this value}
#'    \item{\code{do.cv}}{ if \code{do.cv = FALSE}, a fixed value is selected from a sequence of 10 values chosen in a data-driven way}
#' }
#' @return a list containing the following fields:
#' \item{common.out, idio.out }{ output from \link[fvarseg]{common.seg} and \link[fvarseg]{idio.seg}
#' \itemize{
#' \item{\code{est.cp}}{ a matrix containing the change point estimators in the first column and the finest bandwidth at which each is detected in the second column }
#' \item{\code{G.seq}}{ an integer vector of bandwidths }
#' \item{\code{thr}}{ a vector of thresholds which is of the same length as \code{G.seq} }
#' \item{\code{est.cp.list}}{ a list containing various quantities related to the segmentation; see \link[fvarseg]{common.seg} and \link[fvarseg]{idio.seg} for further details }
#' }}
#' \item{mean.x}{ if \code{center = TRUE}, returns a vector containing row-wise sample means of \code{x}; if \code{center = FALSE}, returns a vector of zeros}
#'
#' @examples
#' \dontrun{
#' out <- sim.data(n = 2000, p = 50, q = 2, d = 1,
#' cp.common = 1:3/4, den.common = .5, type.common = 'ma', 
#' cp.idio = c(3, 5)/8, seed = 123)
#' fs <- fvar.seg(out$x, q = NULL, d = 1)
#' fs$common.out$est.cp
#' fs$idio.out$est.cp
#' }
#' @references Cho, H., Eckley, I., Fearnhead, P. & Maeng, H. (2022) High-dimensional time series segmentation via factor-adjusted vector autoregressive modelling. arXiv preprint arXiv:2204.02724
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