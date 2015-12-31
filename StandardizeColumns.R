# Sample input
x <- matrix(c(NaN,1,2,3,NaN,NaN,1,2,3,NaN),nrow=5,ncol=2)

StandardizeColumns <- function(x){
        
  # This function transforms the columns to have mean 0 and variance 1. 
  # Only non-missing values are standardized.
  # Missing values are retained to NaN.
        
  # Args:
  # x: Data matrix with row and column numbers specified (allows NaN values).        
  # Returns: Matrix with standardized columns. 
  # Find number of variables in the input matrix.
        
  n <- dim(x)[2]
  z <- matrix(0, nrow=dim(x)[1], ncol=dim(x)[2])
  y <- x*z
        
  for (i in 1:n){
    f <- which (1 - is.nan (x[,i]) != 0)
  for (q in f[1]:f[length(f)]){
    y[q,i] <- ((x[q,i] - apply(x, 2, mean, na.rm = TRUE)[i]) /
                        apply(x, 2, sd, na.rm = TRUE)[i])
  }
  }
  
  return (y)
}


