p <- 2
set.seed(1)
x <- matrix(sample.int(15, size = 10*5, replace = TRUE), nrow = 10, ncol = 5)

Mlag2 <- function(x,p){
  # This function creates VAR lags of the order p.
        
  Traw <- nrow(x)
  N <- ncol(x)
  xlag <- matrix(0,Traw, N*p)
 
  for (i in 1:p){
       a <- N*(i-1)+1
       b <- N*i
       xlag[c((p+1):Traw), c(a:b)] <- x[c((p+1-i):(Traw-i)), c(1:N)]
  }
 
  return(xlag)
 
}
