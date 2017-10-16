#' Matrix Variate Normal Distribution
#'
#' Functions \code{dmxnorm} and \code{rmxnorm} are for evaluating
#' densities and generating random samples from corresponding matrix
#' variate normal distribution.
#'
#' @param X a \code{(p-by-n)} matrix whose density be computed.
#' @param n the number of samples to be generated.
#' @param M a \code{(p-by-n)} mean matrix.
#' @param U a \code{(p-by-p)} \emph{left} scale matrix.
#' @param V an \code{(n-by-n)} \emph{right} scale matrix.
#' @param log a logical; \code{TRUE} to return log density, \code{FALSE} otherwise.
#'
#' @examples
#' ## generate 1000 samples from {M,U,V = diag()}
#' samples1000= rmxnorm(1000, diag(3))
#'
#' ## test for LLN : taking average of 1000 samples
#' average1000 = apply(samples1000, c(1,2), mean)
#'
#' ## evaluate the density for average matrix
#' density1000 = dmxnorm(average1000, diag(3))
#'
#' @name mxnorm
#' @rdname mxnorm
NULL

#' @rdname mxnorm
#' @export
dmxnorm <- function(X, M=array(0,c(p,n)), U=diag(nrow(M)), V=diag(ncol(M)), log=FALSE){
  ## Preprocessing
  #   1. should be a matrix
  if (!is.matrix(X)){
    stop("* dmxnorm : an input X should be a matrix.")
  } else {
    p = nrow(X)
    n = ncol(X)
  }
  #   2. other U, V
  if (!missing(U)){
    if (!all(dim(U)==p)){
      stop("* dmxnorm : left scale U does not have corresponding size.")
    }
    if (!isSymmetric(U, tol=sqrt(.Machine$double.eps))){
      stop("* dmxnorm : left scale U must be a symmetric matrix.")
    }
  }
  if (!missing(V)){
    if (!all(dim(V)==n)){
      stop("* dmxnorm : right scale V does not have corresponding size.")
    }
    if (!isSymmetric(V, tol=sqrt(.Machine$double.eps))){
      stop("* dmxnorm : right scale V must be a symmetric matrix.")
    }
  }
  evals_U = tryCatch(eigen(U, only.values = TRUE)$values, error=function(e)e)
  if ((inherits(evals_U, "error"))||(any(evals_U<=0))){
    stop("* dmxnorm : U is not a valid scale parameter.")
  }
  evals_V = tryCatch(eigen(V, only.values = TRUE)$values, error=function(e)e)
  if ((inherits(evals_V, "error"))||(any(evals_V<=0))){
    stop("* dmxnorm : S is not a valid scale parameter.")
  }
  #   3. no NA or Infs
  if (any(is.infinite(X))||any(is.infinite(M))||any(is.infinite(U))||any(is.infinite(V))){
    stop("* dmxnorm : no infinite values allowed.")
  }
  if (any(is.na(X))||any(is.na(M))||any(is.na(U))||any(is.na(V))){
    stop("* dmxnorm : no NA values allowed.")
  }

  ## Main Computation
  #   1. cholesky decomposition
  cholU = tryCatch(chol(U), error=function(e)e)
  cholV = tryCatch(chol(V), error=function(e)e)
  if (inherits(cholU, "error")){
    stop("* dmxnorm : cannot compute chol(U).")
  }
  if (inherits(cholV, "error")){
    stop("* dmxnorm : cannot compute chol(V).")
  }
  detU = (prod(diag(cholU))^2)
  detV = (prod(diag(cholV))^2)
  #   2. Numerator
  tmp = cholsolve(cholU,X-M)
  tmp = t(X-M)%*%tmp
  res = cholsolve(cholV,tmp)
  val_num = exp(-0.5*sum(diag(res)))
  #   3. Denominator
  den1 = ((2*pi)^((n*p)/2))
  den2 = (detV^(p/2))
  den3 = (detU^(n/2))
  val_den = den1*den2*den3

  ## Return results
  val_final = val_num/val_den
  if (log==TRUE){
    val_final = log(val_final)
    return(val_final)
  } else if (log==FALSE){
    return(val_final)
  } else {
    stop("* dmxnorm : 'log' is a logical variable.")
  }
}

#' @rdname mxnorm
#' @export
rmxnorm <- function(n, M, U=diag(nrow(M)), V=diag(ncol(M))){
  ## Preprocessing
  #   1. should be a matrix
  if (!is.matrix(M)){
    stop("* rmxnorm : M should be a matrix.")
  }
  #   2. other U, V
  if (!missing(U)){
    if (!all(dim(U)==nrow(M))){
      stop("* rmxnorm : left scale U does not have corresponding size.")
    }
    if (!isSymmetric(U, tol=sqrt(.Machine$double.eps))){
      stop("* rmxnorm : left scale U must be a symmetric matrix.")
    }
  }
  if (!missing(V)){
    if (!all(dim(V)==ncol(M))){
      stop("* rmxnorm : right scale V does not have corresponding size.")
    }
    if (!isSymmetric(V, tol=sqrt(.Machine$double.eps))){
      stop("* rmxnorm : right scale V must be a symmetric matrix.")
    }
  }
  evals_U = tryCatch(eigen(U, only.values = TRUE)$values, error=function(e)e)
  if ((inherits(evals_U, "error"))||(any(evals_U<=0))){
    stop("* rmxnorm : U is not a valid scale parameter.")
  }
  evals_V = tryCatch(eigen(V, only.values = TRUE)$values, error=function(e)e)
  if ((inherits(evals_V, "error"))||(any(evals_V<=0))){
    stop("* rmxnorm : S is not a valid scale parameter.")
  }
  #   3. no NA or Infs
  if (any(is.infinite(M))||any(is.infinite(U))||any(is.infinite(V))){
    stop("* rmxnorm : no infinite values allowed.")
  }
  if (any(is.na(M))||any(is.na(U))||any(is.na(V))){
    stop("* rmxnorm : no NA values allowed.")
  }
  #   4. n : integer
  if ((n-round(n)>sqrt(.Machine$double.eps))||(is.na(n))||(is.infinite(n))||(n<1)){
    stop("* rmxnorm : n is not a proper integer >= 1.")
  }
  n = as.integer(n)

  ## Main Computation
  #   1. cholesky decomposition
  cholU = tryCatch(chol(U), error=function(e)e)
  cholV = tryCatch(chol(V), error=function(e)e)
  if (inherits(cholU, "error")){
    stop("* rmxnorm : cannot compute chol(U).")
  }
  if (inherits(cholV, "error")){
    stop("* rmxnorm : cannot compute chol(V).")
  }
  #   2. use of RcppZiggurat
  if ("RcppZiggurat" %in% installed.packages()[,1]){
    ZiggFlag = TRUE
  } else {
    ZiggFlag = FALSE
  }
  #   3. case branching
  m = nrow(M)
  p = ncol(M)
  if (n==1){
    output = rmxnorm.single(m,p,M,cholU,cholV,ZiggFlag)
  } else {
    output = array(0,c(m,p,n))
    for (i in 1:n){
      tmpout = rmxnorm.single(m,p,M,cholU,cholV,ZiggFlag)
      output[,,i] = tmpout
    }
  }

  ## Return output
  return(output);
}



# single rnd with Ziggurat ------------------------------------------------
#' @keywords internal
#' @noRd
rmxnorm.single <- function(m,p,M,cholU,cholV,ZiggFlag){
  # random number generation
  if (ZiggFlag){
    X = matrix(zrnorm(m*p),nrow=m)
  } else {
    X = matrix(rnorm(m*p), nrow=m)
  }
  # transform
  output = M + (t(cholU)%*%X%*%cholV);
  return(output);
}


#  ------------------------------------------------------------------------
## CHOLESKY WITH PIVOTING ALLOWED.
#
# R <- chol(sigma, pivot = TRUE)
# R[, order(attr(R, "pivot"))]
