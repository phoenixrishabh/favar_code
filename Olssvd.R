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

# Sample input
set.seed(2)
y <- matrix(sample.int(15, size = 10*5, replace = TRUE), nrow = 10, ncol = 2)

p <- 1

set.seed(3)
lagy <- matrix(rexp(nrow(y) *((ncol(y)*p) + 1)), nrow(y), ncol(y)*p + 1)

Olssvd <- function(y,lagy){
 
  # This function calculates the OLS coefficient estimates for a matrix y   
  # regressed on its lagged values. 
  
  t <- nrow(y)
  
  # Check: round((lagy),4) == round(((svd(lagy)$u) %*% diag(svd(lagy)$d) %*% t(svd(lagy)$v)),4)            
  # The above check corresponds to matrix = (u)(sigma)(v')        
 
  beta <- (svd(lagy)$v * matrix(rep(1/t(svd(lagy)$d), each = nrow(svd(lagy)$v)), 
           nrow = nrow(svd(lagy)$v))) %*% (t(svd(lagy)$u) %*% y)
 
  return(beta)  
 
}


