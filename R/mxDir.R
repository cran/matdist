#' Matrix Variate Dirichlet Distributions
#'
#' We have 2 types of Dirichlet distributions, namely, Type 1 and Type 2.
#' Like other functions, \code{dmxDir1} and \code{dmxDir2} are
#' to evaluate densities, while \code{rmxDir1} and \code{rmxDir2} generate
#' random samples accordingly to the model provided.
#'
#' @param X a \code{(p-by-p-by-r)} array whose density be computed.
#' @param a weight vector of length \code{r+1}.
#' @param p the dimension for sample matrix.
#' @param S a \code{(p-by-p)} covariance matrix required for sampling from Wishart.
#'
#' @examples
#' ## Sample Generations
#' sDir1 = rmxDir1(4,a=c(6,7,8,9),S=diag(4))
#' sDir2 = rmxDir2(4,a=c(9,8,7,6))
#'
#' @name mxDir
#' @rdname mxDir
NULL

#' @rdname mxDir
#' @export
dmxDir1 <- function(X, a){
  ## Preprocessing
  ##  1. if X is 3d array
  if (length(dim(X))==2){
    r = 1
  } else if (length(dim(X))==3){
    r = dim(X)[3]
  }
  ##  2. square property & turn it into 3d array
  if (dim(X)[1]!=dim(X)[2]){
    stop("* dmxDir : input 'X' should be square matrix/matrices.")
  } else {
    p = dim(X)[1]
  }
  if (r==1){
    tmpX = X
    X    = array(0,c(p,p,1))
    X[,,1]= tmpX
  }
  ##  3. vector 'a' : length and value
  if (length(a)!=(r+1)){
    stop("* dmxDir : input vector 'a' is not matching its length.")
  }
  if (any(a<=((p-1)/2))){
    stop("* dmxDir : every element in 'a' should not be larger than (p-1)/2.")
  }
  ##  4. symmetricity check
  for (i in 1:r){
    # 4-1. check if symmetric
    tgt = X[,,i]
    if (!isSymmetric(tgt, tol=sqrt(.Machine$double.eps))){
      stop("* dmxDir : X should be a symmetric matrix.")
    }
  }
  ##  5. determinant check : THIS SHOULD BE DIFFERENT FOR TYPE2 DIRICHLET DISTRIBUTIONS
  ##                         LAST ELEMENT IS FOR I-SUM(U_i) OR I+SUM(V_i)
  sumX3    = sum3d(X)
  vec_detX = rep(0,r+1)
  for (i in 1:r){
    eigvals = tryCatch(eigen(X[,,i], only.values=TRUE)$values, error=function(e)e)
    if (inherits(eigvals, "error")){
      stop("* dmxDir : eigen decomposition is invalid.")
    }
    if (any(eigvals<=0)){
      stop("* dmxDir : X is not positive definite.")
    }
    if (any(eigvals>=1)){
      stop("* dmxDir : X has larger eigenvalues than 1.")
    }
    vec_detX[i] = prod(eigvals)
  }
  eigvals = tryCatch(eigen(diag(p)-sumX3, only.values=TRUE)$values, error=function(e)e)
  if (inherits(eigvals, "error")){
    stop("* dmxDir : eigen decomposition for summation of matrices failed.")
  }
  if (any(eigvals<=0)){
    stop("* dmxDir : sum of X is not positive definite.")
  }
  if (any(eigvals>=1)){
    stop("* dmxDir : sum of X has larger eigenvalues than 1.")
  }
  vec_detX[r+1] = prod(eigvals)

  ## Main Computation
  ## Output
  output = prod(vec_detX^(a-((p+1)/2)))/multibeta_vec(p,a)
  return(output)
}

#' @rdname mxDir
#' @export
dmxDir2 <- function(X, a){
  ## Preprocessing
  ##  1. if X is 3d array
  if (length(dim(X))==2){
    r = 1
  } else if (length(dim(X))==3){
    r = dim(X)[3]
  }
  ##  2. square property & turn it into 3d array
  if (dim(X)[1]!=dim(X)[2]){
    stop("* dmxDir : input 'X' should be square matrix/matrices.")
  } else {
    p = dim(X)[1]
  }
  if (r==1){
    tmpX = X
    X    = array(0,c(p,p,1))
    X[,,1]= tmpX
  }
  ##  3. vector 'a' : length and value
  if (length(a)!=(r+1)){
    stop("* dmxDir : input vector 'a' is not matching its length.")
  }
  if (any(a<=((p-1)/2))){
    stop("* dmxDir : every element in 'a' should not be larger than (p-1)/2.")
  }
  ##  4. symmetricity check
  for (i in 1:r){
    # 4-1. check if symmetric
    tgt = X[,,i]
    if (!isSymmetric(tgt, tol=sqrt(.Machine$double.eps))){
      stop("* dmxDir : X should be a symmetric matrix.")
    }
  }
  ##  5. determinant check : THIS SHOULD BE DIFFERENT FOR TYPE2 DIRICHLET DISTRIBUTIONS
  ##                         LAST ELEMENT IS FOR I-SUM(U_i) OR I+SUM(V_i)
  sumX3    = sum3d(X)
  vec_detX = rep(0,r+1)
  for (i in 1:r){
    eigvals = tryCatch(eigen(X[,,i], only.values=TRUE)$values, error=function(e)e)
    if (inherits(eigvals, "error")){
      stop("* dmxDir : eigen decomposition is invalid.")
    }
    if (any(eigvals<=0)){
      stop("* dmxDir : X is not positive definite.")
    }
    vec_detX[i] = prod(eigvals)
  }
  eigvals = tryCatch(eigen(diag(p)+sumX3, only.values=TRUE)$values, error=function(e)e)
  if (inherits(eigvals, "error")){
    stop("* dmxDir : eigen decomposition for summation of matrices failed.")
  }
  if (any(eigvals<=0)){
    stop("* dmxDir : sum of X is not positive definite.")
  }
  vec_detX[r+1] = prod(eigvals)
  vec_exponent = rep(0,r+1)
  vec_exponent[1:r] = a[1:r]-((p+1)/2)
  vec_exponent[r+1] = -sum(a)

  ## Main Computation
  ## Output
  output = prod(vec_detX^vec_exponent)/multibeta_vec(p,a)
  return(output)
}



# Type 1 : Theorem 6.2.1 from Gupta
#' @rdname mxDir
#' @export
rmxDir1 <- function(p, a, S){
  if ((p<=1)||(abs(p-round(p))>sqrt(.Machine$double.eps))){
    stop("* rmxDir1 : 'p' should be an integer > 1.")
  }
  ## S should be of size (p-by-p), symmetric, and PD
  if (!missing(S)){
    if ((length(dim(S))!=2)||(dim(S)[1]!=p)){
      stop("* rmxDir1 : an input S should be a square matrix")
    }
    if (!checkSPD(S)){
      stop("* rmxDir1 : an input S should be symmetric and positive definite.")
    }
  } else {
    S = diag(p)
  }
  ## about a
  if ((any(is.infinite(a)))||(any(is.na(a)))||(any(a<=0))||(!is.vector(a))){
    stop("* rmxDir1 : an input 'a' is not valid.")
  }
  ## Main Computation
  r = length(a)-1
  nvals = 2*a[1:r]
  mval  = 2*a[r+1]
  #   1. generate B
  B = tryCatch(rmxWishart(1,mval,S), error=function(e)e)
  if (inherits(B, "error")){
    stop("* rmxDir1 : generation from Wishart for B failed.")
  }
  #   2. generate Si's
  aggregateS = array(0,c(p,p,r))
  for (i in 1:r){
    tmpSi = tryCatch(rmxWishart(1,nvals[i],S), error=function(e)e)
    if (inherits(tmpSi, "error")){
      stop("* rmxDir1 : generation from Wishart for S_i failed.")
    }
    aggregateS[,,i] = tmpSi
  }
  #   3. other key ones
  sumS = sum3d(aggregateS)+B
  Shalf = tryCatch(t(chol(sumS)), error=function(e)e)
  if (inherits(Shalf, "error")){
    stop("* rmxDir1 : cholesky on {sum S_i + B} failed.")
  }
  Shalfinv = solve(Shalf)
  #   4. computing the output
  output = array(0,c(p,p,r))
  for (i in 1:r){
    output[,,i] = Shalfinv%*%aggregateS[,,i]%*%t(Shalfinv)
  }
  if (r==1){
    output = output[,,1]
  }
  return(output)
}


# Type 2 : Theorem 6.2.2 from Gupta
#' @rdname mxDir
#' @export
rmxDir2 <- function(p, a){
  ## p as integer
  if ((p<=1)||(abs(p-round(p))>sqrt(.Machine$double.eps))){
    stop("* rmxDir2 : 'p' should be an integer > 1.")
  }
  ## a
  r = length(a)-1
  nvals = 2*a[1:r]
  mval  = 2*a[r+1]

  ## generate B
  B = tryCatch(rmxWishart(1, mval, diag(p)), error=function(e)e)
  if (inherits(B, "error")){
    stop("* rmxDir2 : generating B from Wishart failed.")
  }
  Bsqrt = tryCatch(matrixsqrt(B), error=function(e)e)
  if (inherits(Bsqrt, "error")){
    stop("* rmxDir2 : finding matrix square root for B failed.")
  }
  Bsqrtinv = tryCatch(solve(Bsqrt), error=function(e)e)
  if (inherits(Bsqrtinv, "error")){
    stop("* rmxDir2 : inverting matrix square root failed.")
  }
  ## generate Si and output directly
  output = array(0,c(p,p,r))
  for (i in 1:r){
    tmpSi = tryCatch(rmxWishart(1,nvals[i],S=diag(p)), error=function(e)e)
    if (inherits(tmpSi, "error")){
      stop("* rmxDir2 : generating S_i from Wishart failed.")
    }
    output[,,i] = Bsqrtinv%*%tmpSi%*%Bsqrtinv
  }
  if (r==1){
    output = output[,,1]
  }
  return(output)
}
