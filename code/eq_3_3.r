Fi <- function(Bi, x){
#evaluate model at given parameters B (globally), and abundances x
x = flatten_dbl(list(x))
Tx <-  x%o%x%o%x
return(sum(Bi*Tx))
}

model <- function(x){
#evaluate model at a given parameters B, and abundances x
return(c(Fi(B[1,,,],x),Fi(B[2,,,],x),Fi(B[3,,,],x)))
}

random_tensor <- function(d, n){
#sample B
B <- array(rnorm(n^d), dim = c( n,n,n,n ))
return(B)
}