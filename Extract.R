# Sample input

set.seed(1)
x <- matrix(sample.int(15, size = 10*5, replace = TRUE), nrow = 10, ncol = 5)
k <- 3

# The model is:
#     
# Factor equation:
#    | y[t] |   |   I        0   |   | y[t] |   |   0  |
#    |      | = |                | x |      | + |      |
#    | x[t] |   | L[y,t]  L[f,t] |   | f[t] |   | e[t] |
#  
# VAR equation: 
#    | y[t] |            | y[t-1] |   
#    |      | = B[t-1] x |        | + u[t]
#    | f[t] |            | f[t-1] |   
#       
# where L[t] = (L[y,t] ; L[f,t]) and B[t] are coefficients, f[t] are factors, 
# e[t] ~ N(0,V[t]) and u[t] ~ N(0,Q[t]).


# Logic: factor loadings(lam) multiplied by the factors(fac) + error terms(e[t]) should give back x[t].
# (lam' %*% fac)' - x[t] = e(t)
# t(Extract(x,k)[[2]] %*% t(Extract(x,k)[[1]])) - x[t] = e[t]

Extract <- function(x,k){
        
  # This function extracts the first k principal components from   
  # t*n matrix x. The loadings are normalized so that lam'lam/n = I.
  
  nrow <- dim(x)[1]
  ncol <- dim(x)[2]
  
  xx <- t(x) %*% x
  e <- eigen(t(xx))
  
  # e$values <- round(e$values,2)
  # e$vectors <- round(e$vectors,10)
  
  eigvmat.dec <- diag(x = e$values, ncol, ncol)
  
  # Check : xx %*% e$vectors - e$vectors %*% diag(x = e$values, ncol, ncol)
  # inc.order <- sort(e$values)
  # eigvmat <- diag(x = inc.order, ncol, ncol)
  
  n <- dim(xx)[1]
  evc <- matrix(0, n, n)
  
  for (i in 1:n){
    evc[ ,i] <- e$vectors[ ,i]
  }
  
  lam <- sqrt(n) * evc[ , c(1:k)]
  
  # Check : round(t(lam) %*% lam/n) == diag(1,n,n)
  
  fac <- x %*% lam/n
  list <- list(fac, lam)
  
  return(list)
}
