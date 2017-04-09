#' Calculate and plot Venn diagrams in 2D and 3D.
#'
#' @name vennplot-package
#' @docType package
#' @import stringr
#' @import rgl
#' @importFrom grDevices rainbow
#' @importFrom graphics par plot.new plot.window polygon text title
#' @importFrom stats runif
#' @importFrom utils combn
#' @importFrom Rcpp evalCpp
#' @useDynLib vennplot
#' @examples
#' # arbitrary intersection sizes
#' combinations = c(A=1.8, B=0.9,C=1.3, D = 1.3,E = 1.6, AC=0.3,
#'                  AD= 0.3,BE = 0.3, AE = 0.4, f = 0.7,g =0.8,
#'                  h = 0.5, gf = 0.2, Bh = 0.1,i = 1,j = 0.4,
#'                  k=0.7,l = 1.4,kl = 0.2,m = 0.5,lm = 0.2,
#'                  o = 0.8,p = 0.9, op = 0.3)
#' ve = vennplot(combinations)
#'
#' # binary dataset
#' combinations = sharks[,c(1,3:5,8)]
#' vennplot(combinations = combinations)
NULL
