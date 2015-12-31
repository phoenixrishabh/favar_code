source("Minn_prior_Koop.R")

# Initial condition on lambda_t
lambda_0.mean <- matrix(0,q,r)
lambda_0.var <- diag(1,r)

# Initial condition on beta_t

# Obtain a Minnesota-type prior
b_prior <- Minn_prior_KOOP (0.1, r, nlag, m)[1] 
Vb_prior <- Minn_prior_KOOP (0.1, r, nlag, m)[2] 

beta_0.mean = b_prior
beta_0.var = Vb_prior

# Initial condition on the covariance matrices
V_0 <- 0.1 * diag(q) 
V_0[c(1:p),c(1:p)] = 0
Q_0 <- 0.1 * diag(r)


# Initialize matrices
lambda_0_prmean = lambda_0.mean
lambda_0_prvar = lambda_0.var

beta_0_prmean = beta_0.mean
beta_0_prvar = beta_0.var

lambda_pred <- array(0,c(q,r,t))
lambda_update <- array(0,c(q,r,t))

for (j in 1:t){
lambda_pred[c(1:r),c(1:r),j] <- diag(r)
lambda_update[c(1:r),c(1:r),j] <- diag(r)
}

beta_pred = matrix(0, m, t)
beta_update = matrix(0, m, t)
------------------------------------------------------------------------------------
        
        
        
Rl_t = array(r,r,q,t);
Sl_t = zeros(r,r,q,t);
Rb_t = zeros(m,m,t);
Sb_t = zeros(m,m,t);

x_t_pred = zeros(t,q);
e_t = zeros(q,t);
lambda_t = zeros(q,r,t);
beta_t = zeros(k,k,t);

Q_t = zeros(r,r,t);
V_t = zeros(q,q,t);

# Decay and forgetting factors
l_1 = l(1); l_2 = l(2); l_3 = l(3); l_4 = l(4);

# Define lags of the factors to be used in the state (VAR) equation         
yy = FPC(nlag+1:t,:);      
xx = mlag2(FPC,nlag); xx = xx(nlag+1:t,:);
templag = mlag2(FPC,nlag); templag = templag(nlag+1:t,:);      
[Flagtemp,m] = create_RHS_NI(templag,r,nlag,t);  
Flag = [zeros(k,m); Flagtemp];