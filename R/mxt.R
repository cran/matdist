#' Matrix Variate t-Distribution
#'
#' Functions \code{dmxt} and \code{rmxt} are for evaluating
#' densities and generating random samples from corresponding matrix
#' variate t-distribution.
#'
#'
#' @param X a \code{(p-by-m)} matrix whose density be computed.
#' @param df degrees of freedom.
#' @param n the number of samples to be generated.
#' @param M a \code{(p-by-n)} mean matrix.
#' @param U a \code{(p-by-p)} \emph{left} spread matrix.
#' @param V an \code{(n-by-n)} \emph{right} spread matrix.
#'
#' @examples
#' ## generate 5 samples from {M,U,V = diag()} of df = 10
#' st5 = rmxt(5, df=10, M=diag(3))
#'
#' ## evaluate density for the 3rd sample
#' dmxt(st5[,,3], df=10, M=diag(3))
#'
#' @name mxt
#' @rdname mxt
NULL

#' @rdname mxt
#' @export
dmxt <- function(X, df, M=array(0,c(p,m)), U=diag(nrow(M)), V=diag(ncol(M))){
  ## Preprocessing
  #   1. should be a matrix
  if (!is.matrix(X)){
    stop("* dmxt : an input X should be a matrix.")
  } else {
    p = nrow(X)
    m = ncol(X)
  }
  #   2. degree of freedom
  n = df
  if ((n-round(n)>sqrt(.Machine$double.eps))||(is.na(n))||(is.infinite(n))){
    stop("* dmxt : 'df' (degree of freedom) is not a proper integer.")
  }
  n = round(n)
  #   3. U and V
  if (!missing(U)){
    if (!all(dim(U)==p)){
      stop("* dmxt : left scale U does not have corresponding size.")
    }
    if (!isSymmetric(U, tol=sqrt(.Machine$double.eps))){
      stop("* dmxt : left scale U must be a symmetric matrix.")
    }
  }
  if (!missing(V)){
    if (!all(dim(V)==m)){
      stop("* dmxt : right scale V does not have corresponding size.")
    }
    if (!isSymmetric(V, tol=sqrt(.Machine$double.eps))){
      stop("* dmxt : right scale V must be a symmetric matrix.")
    }
  }
  evals_U = tryCatch(eigen(U, only.values = TRUE)$values, error=function(e)e)
  if ((inherits(evals_U, "error"))||(any(evals_U<=0))){
    stop("* dmxt : U is not a valid scale parameter.")
  }
  evals_V = tryCatch(eigen(V, only.values = TRUE)$values, error=function(e)e)
  if ((inherits(evals_V, "error"))||(any(evals_V<=0))){
    stop("* dmxt : V is not a valid scale parameter.")
  }

    ## Main Computation
  #   1. cholesky decomposition
  cholS = tryCatch(chol(U), error=function(e)e)
  cholO = tryCatch(chol(V), error=function(e)e)
  if (inherits(cholS, "error")){stop("* dmxt : cannot compute chol(U).")}
  if (inherits(cholO, "error")){stop("* dmxt : cannot cmopute chol(V).")}
  #   2. term-by-term computation
  #   2-1. term1 : multigamma and constants
  t1_num = multigamma(p,(n+m+p-1)/2)
  t1_den = (pi^(m*p/2))*multigamma(p,(n+p-1)/2)
  t1     =  (t1_num)/(t1_den)
  #   2-2. determinants
  detS = choldet(cholS)
  detO = choldet(cholO)
  t2   = (detS^(-m/2))*(detO^(-p/2))
  #   2-3. last term
  tmp1 = cholsolve(cholO,t(X-M))
  tmp1 = (X-M)%*%tmp1
  tmp2 = cholsolve(cholS,tmp1)+diag(p)
  t3   = ((base::det(tmp2))^(-(n+m+p-1)/2))

  ## Return output
  output = t1*t2*t3
  return(output)
}

#' @rdname mxt
#' @export
rmxt <- function(n, df, M, U=diag(nrow(M)), V=diag(ncol(M))){
  ## Basic setting
  if (!is.matrix(M)){
    stop("* rmxt : M should be a matrix of proper size.")
  } else {
    p = nrow(M)
    m = ncol(M)
  }

  evals_U = tryCatch(eigen(U, only.values = TRUE)$values, error=function(e)e)
  if ((inherits(evals_U, "error"))||(any(evals_U<=0))){
    stop("* rmxt : U is not a valid scale parameter.")
  }
  evals_V = tryCatch(eigen(V, only.values = TRUE)$values, error=function(e)e)
  if ((inherits(evals_V, "error"))||(any(evals_V<=0))){
    stop("* rmxt : V is not a valid scale parameter.")
  }

  ## Use of others
  invU = tryCatch(chol2inv(U), error=function(e)e)
  if (inherits(invU, "error")){stop("* rmxt : cannot compute chol(U).")}
  output = array(0,c(p,m,n))
  for (i in 1:n){
    S = rmxWishart(1,df+p-1,invU);
    cholS = tryCatch(chol(S), error=function(e)e);
    if (inherits(cholS, "error")){stop("* rmxt : cannot compute chol(S).")}
    X = rmxnorm(1,array(0,c(p,m)),diag(p),V)
    invS2 = forwardsolve(t(cholS),diag(p))
    output[,,i] = (t(invS2)%*%X+M)
  }
  if (n==1){
    output = output[,,1]
  }
  ## Return output
  return(output)
}
