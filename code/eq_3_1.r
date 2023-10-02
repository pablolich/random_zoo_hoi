Fi <- function(Bi, x){
#evaluate model at given parameters B (globally), and abundances x
x = flatten_dbl(list(x))
xmod = c(x, 1)
Tx <-  xmod
return(sum(Bi*Tx))
}

model <- function(x){
#evaluate model at a given parameters B, and abundances x
return(c(Fi(B[1,],x),Fi(B[2,],x),Fi(B[3,],x)))
}

random_tensor <- function(d, n){
#sample B
B <- array(rnorm((n+1)^(d+1)), dim = c( n,n+1 ))
return(B)
}