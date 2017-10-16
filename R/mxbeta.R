#' Matrix Variate Beta Distributions.
#'
#' We have 2 types of Beta distributions, namely, Type 1 and Type 2.
#' Like other functions, \code{dmxbeta1} and \code{dmxbeta2} are
#' to evaluate densities, while \code{rmxbeta1} and \code{rmxbeta2} generate
#' random samples accordingly to the model provided.
#'
#' @param X a \code{(p-by-p)} matrix whose density be computed.
#' @param a positive real number \code{> (p-1)/2}.
#' @param b positive real number \code{> (p-1)/2}.
#' @param n the number of samples to be generated.
#' @param p the dimension for sample matrix.
#' @param S a \code{(p-by-p)} covariance matrix required for Type 1 Beta sampling.
#'
#' @examples
#' ## Sample Generations
#' sbeta1 = rmxbeta1(10, p=4, a=6, b=2, S=diag(4))
#' sbeta2 = rmxbeta2(10, p=4, a=6, b=2)
#'
#' @name mxbeta
#' @rdname mxbeta
NULL

#' @rdname mxbeta
#' @export
dmxbeta1 <- function(X, a, b){
  ## Preprocessing
  #   1. should be a matrix
  if (!is.matrix(X)){
    stop("* dmxbeta1 : an input X should be a matrix.")
  } else {
    p = nrow(X)
  }
  if (p!=ncol(X)){
    stop("* dmxbeta1 : X should be a square matrix.")
  }
  #   2. symmetric
  if (!isSymmetric(X, tol=sqrt(.Machine$double.eps))){
    stop("* dmxbeta1 : X should be a symmetric matrix.")
  }
  #   3. Parameters a,b
  if (a<=(p-1)/2){
    stop("* dmxbeta1 : parameter 'a' is invalid.")
  }
  if (b<=(p-1)/2){
    stop("* dmxbeta1 : parameter 'b' is invalid.")
  }

  ## Main Computation
  #   1. eigenvalue and determinants
  eigvalsX = tryCatch(eigen(X, only.values=TRUE)$values, error=function(e)e)
  if (inherits(eigvalsX, "error")){
    stop("* dmxbeta1 : eigendecomposition of X fails.")
  }
  if ((any(eigvalsX<=0))||(any(eigvalsX>=1))){
    stop("* dmxbeta1 : the condition {0 < X < diag(p)} is not met.")
  }
  detX = prod(eigvalsX)

  eigvalsImX = tryCatch(eigen((diag(p)-X), only.values=TRUE)$values, error=function(e)e)
  if (inherits(eigvalsImX, "error")){
    stop("* dmxbeta1 : eigendecomposition of diag(p)-X fails.")
  }
  if (any(eigvalsImX<=0)){
    stop("* dmxbeta1 : the condition {diag(p)-X > 0} is not met.")
  }
  detImX = prod(eigvalsImX)

  #   2. term decompositions
  tm1 = 1/multibeta(p,a,b)
  tm2 = (detX^(a-((p+1)/2)))
  tm3 = (detImX^(b-((p+1)/2)))
  output = tm1*tm2*tm3

  ## Return output
  return(output)
}

#' @rdname mxbeta
#' @export
dmxbeta2 <- function(X, a, b){
  ## Preprocessing
  #   1. should be a matrix
  if (!is.matrix(X)){
    stop("* dmxbeta2 : an input X should be a matrix.")
  } else {
    p = nrow(X)
  }
  if (p!=ncol(X)){
    stop("* dmxbeta2 : X should be a square matrix.")
  }
  #   2. symmetric
  if (!isSymmetric(X, tol=sqrt(.Machine$double.eps))){
    stop("* dmxbeta2 : X should be a symmetric matrix.")
  }
  #   3. Parameters a,b
  if (a<=(p-1)/2){
    stop("* dmxbeta2 : parameter 'a' is invalid.")
  }
  if (b<=(p-1)/2){
    stop("* dmxbeta2 : parameter 'b' is invalid.")
  }

  ## Main Computation
  #   1. cholesky decomposition
  eigvalsX = tryCatch(eigen(X, only.values=TRUE)$values, error=function(e)e)
  if (inherits(eigvalsX, "error")){
    stop("* dmxbeta2 : eigendecomposition of X fails.")
  }
  if (any(eigvalsX<=0)){
    stop("* dmxbeta2 : the condition {0 < X} is not met.")
  }
  detX = prod(eigvalsX)

  eigvalsIpX = tryCatch(eigen(diag(p)+X, only.values=TRUE)$values, error=function(e)e)
  if (inherits(eigvalsIpX, "error")){
    stop("* dmxbeta2 : eigendecomposition of diag(p)+X fails.")
  }
  detIpX = prod(eigvalsIpX)

  #   2. term decompositions
  tm1 = 1/multibeta(p,a,b)
  tm2 = (detX^(a-((p+1)/2)))
  tm3 = (detIpX^(-(a+b)))
  output = tm1*tm2*tm3

  ## Return output
  return(output)
}


#' @rdname mxbeta
#' @export
rmxbeta1 <- function(n, p, a, b, S){
  ## Preprocessing
  ## For this case, leave p, a, b to the latter single parts
  if ((n-round(n)>sqrt(.Machine$double.eps))||(is.na(n))||(is.infinite(n))||(n<1)){
    stop("* rmxbeta1 : n is not a proper integer >= 1.")
  }
  n = as.integer(n)
  if ((p-round(p)>sqrt(.Machine$double.eps))||(is.na(p))||(is.infinite(p))||(p<1)){
    stop("* rmxbeta1 : p is not a proper integer >= 1.")
  }
  p = as.integer(p)
  ## Iteration
  if (n==1){
    output = rmxbeta1.single(p,a,b,S)
  } else {
    output = array(0,c(p,p,n))
    for (i in 1:n){
      output[,,i] = rmxbeta1.single(p,a,b,S)
    }
  }
  ## Return
  return(output)
}


#' @rdname mxbeta
#' @export
rmxbeta2 <- function(n, p, a, b){
  ## Preprocessing
  ## For this case, leave p, a, b to the latter single parts
  if ((n-round(n)>sqrt(.Machine$double.eps))||(is.na(n))||(is.infinite(n))||(n<1)){
    stop("* rmxbeta2 : n is not a proper integer >= 1.")
  }
  n = as.integer(n)
  if ((p-round(p)>sqrt(.Machine$double.eps))||(is.na(p))||(is.infinite(p))||(p<1)){
    stop("* rmxbeta2 : p is not a proper integer >= 1.")
  }
  p = as.integer(p)
  ## Iteration
  if (n==1){
    output = rmxbeta2.single(p,a,b)
  } else {
    output = array(0,c(p,p,n))
    for (i in 1:n){
      output[,,i] = rmxbeta2.single(p,a,b)
    }
  }
  ## Return
  return(output)
}




# single functions --------------------------------------------------------
#' @keywords internal
#' @noRd
rmxbeta1.single <- function(p,a,b,S){
  ## parameters
  if (2*a>=p){
    n1 = 2*a
    n2 = 2*b
    if (missing(S)){
      S = diag(p)
    }
    X = tryCatch(rmxnorm(1,array(0,c(p,n1)),U=S,V=diag(n1)), error=function(e)e)
    S2= tryCatch(rmxWishart(1,n2,S), error=function(e)e)
    if ((inherits(X, "error"))||(inherits(S2, "error"))){
      stop("* rmxbeta1 : generation from rmxnorm and rmxWishart failed.")
    }
    S = S2+X%*%t(X)
    Shalf=tryCatch(matrixsqrt(S), error=function(e)e)
    if (inherits(Shalf, "error")){
      stop("* rmxbeta1 : acquiring matrix square root failed.")
    }
    Z = tryCatch(solve(Shalf,X), error=function(e)e)
    if (inherits(Z, "error")){
      stop("* rmxbeta1 : solving linear system failed.")
    }
    output = Z%*%t(Z)
  } else {
    p_new = 2*a
    n1_new = p
    n2_new = 2*(a+b)-p
    if (missing(S)){
      S = diag(p_new)
    }

    X = tryCatch(rmxnorm(1,array(0,c(p_new,n1_new)),U=S), error=function(e)e)
    S2= tryCatch(rmxWishart(1,n2_new,S), error=function(e)e)
    if ((inherits(X, "error"))||(inherits(S2, "error"))){
      stop("* rmxbeta1 : generation from rmxnorm and rmxWishart failed.")
    }
    S = S2+X%*%t(X)
    Shalf=tryCatch(matrixsqrt(S), error=function(e)e)
    if (inherits(Shalf, "error")){
      stop("* rmxbeta1 : acquiring matrix square root failed.")
    }
    Z = tryCatch(solve(Shalf,X), error=function(e)e)
    if (inherits(Z, "error")){
      stop("* rmxbeta1 : solving linear system failed.")
    }
    output = t(Z)%*%Z
  }
  return(output)
}

#' @keywords internal
#' @noRd
rmxbeta2.single <- function(p,a,b){
  ## parameters
  n1 = 2*a
  n2 = 2*b
  ## generate from Wishart
  S1 = tryCatch(rmxWishart(1,n1,diag(p)), error=function(e)e)
  S2 = tryCatch(rmxWishart(1,n2,diag(p)), error=function(e)e)
  if ((inherits(S1, "error"))||(inherits(S2, "error"))){
    stop("* rmxbeta2 : generation using rmxWishart failed.")
  }
  S2halfinv = tryCatch(solve((matrixsqrt(S2))), error=function(e)e)
  if (inherits(S2halfinv, "error")){
    stop("* rmxbeta2 : acquiring inverse of matrix symmetric square root failed.")
  }
  V = S2halfinv%*%S1%*%S2halfinv
  return(V)
}
