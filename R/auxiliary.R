# Applications of Cholesky Decompositions ---------------------------------
# 1. cholsolve : solve(A,B) using cholesky decomposition 'cholA' of A
#' @noRd
#' @keywords internal
cholsolve <- function(cholA,B){
  tmpx  = forwardsolve(t(cholA),B)
  x     = backsolve(cholA,tmpx)
  return(x)
}

# 2. choldet : compute det(A) based on 'cholA'
#' @noRd
#' @keywords internal
choldet <- function(cholU){
  detU = (prod(diag(cholU))^2)
  return(detU)
}


# multivariate distributions ----------------------------------------------
# 1. multigamma : evaluate Gamma_p(a)
#    Refer to Gupta (THM1.4.1)
#' @noRd
#' @keywords internal
multigamma <- function(p,a){
  output = pi^(p*(p-1)/4)
  for (j in 1:p){
    output = output*gamma(a+((1-j)/2))
  }
  return(output)
}

# 2. multibeta : evaluate Beta_p(a,b)
#    Refer to Gupta (1.4.8); page 20
#' @noRd
#' @keywords internal
multibeta <- function(p,a,b){
  output = multigamma(p,a)*multigamma(p,b)/multigamma(p,(a+b))
  return(output)
}

# 3. multibeta_vec : where we have 'veca'
#' @noRd
#' @keywords internal
multibeta_vec <- function(pval, veca){
  tm1 = prod(unlist(lapply(veca, multigamma, p=pval)))
  tm2 = multigamma(pval, sum(veca))
  output = tm1/tm2
  return(output)
}


# Miscellany --------------------------------------------------------------
# 1. square root of a matrix
# https://stackoverflow.com/questions/28227111/square-root-of-a-singular-matrix-in-r
#' @noRd
#' @keywords internal
matrixsqrt <- function(X){
  Y = eigen(X)
  T = Y$vectors
  Tinv = solve(T)
  Jsqrt = diag(sqrt(Y$values))
  Xsqrt = T%*%Jsqrt%*%Tinv
  return(Xsqrt)
}

# 2. summation of 3d array in 3rd dimension
#' @noRd
#' @keywords internal
sum3d <- function(X){
  if (length(dim(X))!=3){
    stop("* sum3d : not suitable.")
  }
  m = dim(X)[1]
  n = dim(X)[2]
  output = array(0,c(m,n))
  for (i in 1:m){
    for (j in 1:n){
      output[i,j] = sum(X[i,j,])
    }
  }
  return(output)
}

# 3. check if Symmetric and Positive Definite
#' @noRd
#' @keywords internal
checkSPD <- function(X){
  # 3-1. size equal
  p = nrow(X)
  if (ncol(X)!=p){
    stop("* checker : the given matrix is not square sized.")
    return(FALSE)
  }
  # 3-2. symmetric
  if (!isSymmetric(X, tol=sqrt(.Machine$double.eps))){
    stop("* checker : the given matrix is not symmetric up to sqrt{machine epsilon}.")
    return(FALSE)
  }
  # 3-3. positive definite
  eigvals <- tryCatch(eigen(X, only.values = TRUE)$values, error = function(e)e)
  if (inherits(eigvals, "error")){
    stop("* checker : eigendecomposition of an input matrix failed.")
    return(FALSE)
  }
  if (any(eigvals<=0)){
    stop("* checker : an input is not positive definite.")
    return(FALSE)
  }
  return(TRUE)
}
