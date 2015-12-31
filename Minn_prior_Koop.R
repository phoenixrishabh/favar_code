                                                                                                                                
Minn_prior_KOOP <- function(gamma,M,p,K){

# This is the version of the Minnesota prior with no dependence on the
# standard deviations of the univariate regressions. This prior allows
# online estimation and forecasting of the large TVP-VAR.

# 1. Minnesota Mean on VAR regression coefficients

A_prior <- t(cbind(0.9*diag(M), matrix(0,(p-1)*M,M)))
m <- c(A_prior)
a_prior <- matrix(m,length(m),1)

# 2. Minnesota Variance on VAR regression coefficients

# Create an array of dimensions K x M, which will contain the K diagonal   
# elements of the covariance matrix, in each of the M equations.
V_i <- matrix(0, K/M, M)

for (i in 1:M){  # for each i-th equation
    for (j in 1:K/M){   # for each j-th RHS variable        
        V_i[j,i] <- gamma/((ceiling(j/M))^2) # variance on own lags           
        # Note: the "ceil((j-1)/M^2)" command finds the associated lag 
        # number for each parameter.         
    }
}

# Now V (MINNESOTA VARIANCE) is a diagonal matrix with diagonal elements the V_i'  
V_i_T = t(V_i)
V_prior <- diag(c(V_i_T))  # this is the prior variance of the vector alpha

output <- list(a_prior,V_prior)
return(output)

}
