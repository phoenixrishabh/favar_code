factor_0.mean <- matrix(0,k,1)
factor_0.var <- diag(10,k)

# Initialize matrices
factor_0_prmean = factor_0.mean
factor_0_prvar = factor_0.var


factor_pred <- matrix(0,k,t)
factor_update <- matrix(0,k,t)

Rf_t <- array(0, dim=c(k,k,t))
Sf_t <- array(0, dim=c(k,k,t))

x_t_predf <- matrix(0,t,q)
ef_t <- matrix(0,q,t)



