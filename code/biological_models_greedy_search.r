# input
# vector r
# matrix A
# diag(A) < 0

build_circulant <-function(x) {
  n <- length(x)
  suppressWarnings(
    matrix(x[matrix(1:n,n+1,n+1,byrow=T)[c(1,n:2),1:n]],n,n))
}

dxdt <- function(x){
  # simply evaluate the equations
  ## \dot{x}_i = x_i(r_i + \sum_j aij xj / dj)
  ## where if A_{ij} >0, i != j, dj = 1 + \sum{k} A_{ik} x_k over the Aik > 0
  ## if A_{ij} < 0, i != j, dj = 1 + \sum{k} A_{kj} x_k over the Akj > 0
  ## if i = j, then dj = 1
  # build matrix of denominators
  denom <- 1+ as.vector(Ap %*% abs(x))
  denom_row <- denom %o% rep(1, n)
  denom_col <- rep(1, n) %o% denom
  denom <- (Ap > 0) * denom_row + (An < 0) * denom_col + diag(rep(1, length(r)))
  denom[denom == 0] <- 1 # if no prey
  # equations
  dx <- x * (r + as.vector((A / denom) %*% x))
  return(dx)
}

find_root <- function(){
  n <- length(r)
  # sample initial condition
  x0 <- rexp(n + 1, 1)
  x0 <- x0[1:n] / x0[n+1]
  tmp <- rootSolve::multiroot(dxdt, start = x0, rtol = 10^-16, atol = 10^-16)
  return(tmp)
}

has_feasible <- function(ntry = 1000){
  for (i in 1:ntry){
    tmp <- find_root()
    if (all(tmp$root > 10^-9) & (tmp$estim.precis < 10^-16) & all(abs(tmp$f.root) < 10^-16)){
      return(1)
    } 
  }
  return(0)
}

library(tidyverse)

results <- tibble()
for (nn in 2:2){
  for (kk in 1:5000){
    n <- nn
    r <- rnorm(n) + sign(rnorm(1))
    A <- matrix(rnorm(n * n), n, n)
    diag(A) <- -abs(diag(A)) - 1
    Anod <- A # without diagonal
    diag(Anod) <- 0
    #A[upper.tri(A)] <- 0
    #A[lower.tri(A)] <- 0
    # Split the matrix
    Ap <- Anod * (A > 0) # positive
    An <- Anod * (A < 0) # negative
    Ad <- diag(diag(A)) # diagonal
    results <- bind_rows(results, tibble(s = n ,  kk = kk, feas = has_feasible()))
    #print(results[nrow(results),] %>% as.numeric())
    if (kk %% 10 == 0){
      print(results %>% group_by(s, feas) %>% count() %>% as.matrix())
    }
  }
}

results %>% ggplot(aes(x = as.factor(s), fill = as.factor(feas))) + geom_bar()
# print prob feas
results %>% group_by(s) %>% summarise(pf = sum(feas) / length(feas))