#' Wishart Distribution
#'
#' Functions \code{dmxWishart} and \code{rmxWishart} are for evaluating
#' densities and generating random samples from corresponding Wishart distribution.
#'
#' @param X a \code{(p-by-p)} matrix whose density be computed.
#' @param df degrees of freedom
#' @param S a \code{(p-by-p)} scale/covariance matrix.
#' @param n the number of samples to be generated.
#'
#' @examples
#' ## generate 10 samples with S=diag(3) with df=10
#' smW = rmxWishart(10, df=10, S=diag(3))
#'
#' ## evaluate each sample's density
#' t(apply(smW, 3, dmxWishart, df=10, S=diag(3)))
#'
#' @name mxWishart
#' @rdname mxWishart
NULL

#' @rdname mxWishart
#' @export
dmxWishart <- function(X, df, S=diag(nrow(X))){
  ## Preprocessing
  #   1. should be a matrix
  if (!is.matrix(X)){
    stop("* dmxWishart : an input X should be a matrix.")
  } else {
    p = ncol(X)
  }
  if (nrow(X)!=p){
    stop("* dmxWishart : an input X should be a square matrix.")
  }
  #   2. S
  if (!missing(S)){
    if (!is.matrix(S)){
      stop("* dmxWishart : scale S should be a matrix.")
    }
    if (!all(dim(S)==p)){
      stop("* dmxWishart : scale S does not have corresponding size.")
    }
    if (!isSymmetric(S, tol=sqrt(.Machine$double.eps))){
      stop("* dmxWishart : scale S must be a symmetric matrix.")
    }
  }
  evals_S = tryCatch(eigen(S, only.values = TRUE)$values, error=function(e)e)
  if ((inherits(evals_S, "error"))||(any(evals_S<=0))){
    stop("* dmxWishart : S is not a valid scale parameter.")
  }


  #   3. degree of freedom
  n = df
  if ((n-round(n)>sqrt(.Machine$double.eps))||(is.na(n))||(is.infinite(n))||(n<p)){
    stop("* dmxWishart : 'df' (degree of freedom) is not a proper integer >= p.")
  }
  n = round(n)

  ## Main Computation
  #   1. cholesky decomposition on both matrices
  cholX = tryCatch(chol(X), error=function(e)e)
  cholS = tryCatch(chol(S), error=function(e)e)
  if (inherits(cholX, "error")){
    stop("* dmxWishart : cannot compute chol(X).")
  }
  if (inherits(cholS, "error")){
    stop("* dmxWishart : cannot compute chol(S).")
  }
  #   2. relevant compotations
  #   2-1. materials for numerator
  detX = choldet(cholX)
  solveSX = cholsolve(S,X)
  val_num = (detX^((n-p-1)/2))*exp(-0.5*sum(diag(solveSX)))
  #   2-2. materials for denominator
  detS = choldet(cholS)
  val_den = (2^(n*p/2))*(detS^(n/2))*multigamma(p,(n/2))

  ## Return output
  output = val_num/val_den
  return(output)
}



#' @rdname mxWishart
#' @export
rmxWishart <- function(n, df, S){
  ## Preprocessing
  #   1. number of replicates
  if ((n-round(n)>sqrt(.Machine$double.eps))||(is.na(n))||(is.infinite(n))||(n<1)){
    stop("* rmxWishart : number of samples 'n' is not a proper integer >= 1.")
  }
  n = as.integer(n)
  ##  2. S
  if (!is.matrix(S)){
    stop("* rmxWishart : scale S should be a matrix.")
  } else {
    p = ncol(S)
  }
  evals_S = tryCatch(eigen(S, only.values = TRUE)$values, error=function(e)e)
  if ((inherits(evals_S, "error"))||(any(evals_S<=0))){
    stop("* rmxWishart : S is not a valid scale parameter.")
  }
  if (!isSymmetric(S, tol=sqrt(.Machine$double.eps))){
    stop("* rmxWishart : scale S must be a symmetric matrix.")
  }
  ##  3. degree of freedom
  if ((df-round(df)>sqrt(.Machine$double.eps))||(is.na(df))||(is.infinite(df))||(df<p)){
    stop("* rmxWishart : 'df' (degree of freedom) is not a proper integer >= ncol(S).")
  }
  dfval = round(df)

  ## Simply generate

  if (n==1){
    tmp    = rWishart(1, df=dfval, Sigma=S)
    output = tmp[,,1]
  } else {
    output = rWishart(n, df=dfval, Sigma=S)
  }

  ## Return
  return(output)
}



