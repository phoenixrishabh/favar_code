source("StandardizeColumns.R")
source("Extract.R")
source("Olssvd.R")
source("Mlag2.R")

library(pracma)

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
# e[t]~N(0,V[t]) and u[t]~N(0,Q[t]).

set.seed(1)
x <- matrix(sample.int(15, size = 10*5, replace = TRUE), nrow = 10, ncol = 5)

set.seed(2)
y <- matrix(sample.int(15, size = 10*5, replace = TRUE), nrow = 10, ncol = 2)


nfac <- 1    # number of factors to be extracted
nlag <- 2    # number of lags of factors
p <- 2       # number of macro variables
y_true <- 1  # include y[t] in the factor equation if y_true = 1, if it is NOT
             # included, then the coefficient/loading L[y,t] is zero for all periods

Olspcdfm <- function(x, y, nfac, nlag, p, y_true){
        
  # This function computes the loadings matrix L(t) and the errors of the
  # factor equation and the VAR equation.
        
  n <- ncol(x)
  t <- nrow(y)        
  r <- nfac + p

  X_st <- StandardizeColumns(x)
  FPC2 <- as.matrix(unlist(Extract(X_st, nfac)[1]))
  LPC <- as.matrix(unlist(Extract(X_st, nfac)[2]))
  FPC <- cbind(y,FPC2)
  YX <- cbind(y,X_st)

  L = t(Olssvd(YX, FPC))
  YF = FPC
  Lf = LPC
  
  # Obtain L (the loadings matrix)
  if (y_true == 1) {
     L <- t(Olssvd(YX, FPC))
  } else if (y_true == 0) {
     L <- rbind(cbind(diag(p), matrix(0, p, nfac)), cbind(matrix(0, n, p), LPC))
  } else 
  print("Enter either 0 or 1")
  
  # Obtain the errors of the factor equation
  e <- YX - YF %*% t(L)
  sigma2 <- diag(diag((t(e)%*%e)/t))
  
  # Obtain the errors of the VAR equation
  yy <- YF[c((nlag+1):nrow(YF)),]
  xx <- Mlag2(YF, nlag)
  xx <- xx[c((nlag+1):nrow(xx)),]
  beta_OLS <- solve(t(xx) %*% xx) %*% (t(xx) %*% yy)
  sigmaf <- t(yy - xx %*% beta_OLS) %*% (yy - xx %*% beta_OLS)/(t-nlag-1)
  beta_var <- kron(sigmaf, solve(t(xx) %*% xx))

  for (i in 1:nlag){
       a<- (i-1)*r+1
       b <- i*r
       g <- beta_OLS[c(a:b),1:r]
       bb <- as.vector(t(g))
  }
 
  output <- list(L, bb, beta_OLS, sigma2, sigmaf)
  return(output)
 
} 

