#' Matrix Variate Distributions
#'
#' \pkg{matdist} provides tools for computing densities and generating samples
#' from matrix variate distributions, including matrix normal, Wishart, matrix t,
#' matrix Dirichlet, matrix beta distributions. Special cases and automatic
#' rank adjustments are not yet supported. For detailed discussion, we refer to
#' a \href{https://www.crcpress.com/Matrix-Variate-Distributions/Gupta-Nagar/p/book/9781584880462}{book}
#' with thorough description by Gupta and Nagar (1999).
#'
#' @section Composition:
#' Current version involves following distributions,
#' \describe{
#'   \item{\link{mxDir}}{matrix Dirichlet distribution}
#'   \item{\link{mxWishart}}{Wishart distribution}
#'   \item{\link{mxbeta}}{matrix Beta distribution}
#'   \item{\link{mxnorm}}{matrix Normal distribution}
#'   \item{\link{mxt}}{matrix t-distribution}
#' }
#' For each function, currently we support prefix \code{d} and \code{r}
#' for \emph{density evaluation} and \emph{random generation}, respectively.
#'
#' @docType package
#' @name matdist
#' @import Rlinsolve
#' @importFrom stats rWishart rnorm
#' @importFrom utils installed.packages
#' @importFrom RcppZiggurat zrnorm
#' @importFrom Rcpp evalCpp
#' @useDynLib matdist, .registration=TRUE
NULL


