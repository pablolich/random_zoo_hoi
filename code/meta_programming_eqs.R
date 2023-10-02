library(pracma)
library(tidyverse)

write_kss_model <- function(n, d){
  sink(paste0("eq_", n, "_", d, ".r"))
  #write function to evaluate equation i of the model
  cat("Fi <- function(Bi, x){\n")
  #meta-comment
  cat("#evaluate model at given parameters B (globally), and abundances x\n")
  cat("x = flatten_dbl(list(x))\n")
  cat("xmod = c(x, 1)\n")
  outerx <- paste0("xmod",paste0(rep("%o%xmod",d), collapse = ""))
  cat("Tx <- ", outerx)
  cat("\nreturn(sum(Bi*Tx))\n}\n\n")
  
  
  #write function to outputing evaluations of model
  #signature of function
  cat("model <- function(x){\n")
  #meta-comment
  cat("#evaluate model at a given parameters B, and abundances x\n")
  #return a vector of evaluations for each equation
  cat("return(c(")
  comas <- paste0(rep(",",d), collapse="")
  cat(paste0("Fi(B[",1,comas,"],x)"))
  if (n>1){
    for (i in 2:n){
      cat(paste0(",Fi(B[",i,comas,"],x)"))
    }
  }
  cat("))\n}\n\n")
  
  #write a function to sample a tensor of appropriate dimensions
  #signature of function
  cat("random_tensor <- function(d, n){\n")
  #meta-comment
  cat("#sample B\n")
  ns <- paste0("n", paste0(rep(",n+1", d), collapse=""))
  cat("B <- array(rnorm((n+1)^(d+1)), dim = c(", ns,"))\n")
  cat("return(B)\n}")
  
  sink()
}

#test
maxd = 3#6
maxn = 3#8
results = expand.grid(2:maxn, 1:maxd)
#add extra column
results[,3] = rep(0,nrow(results))
colnames(results)<-c("n", "d", "pf")
for (nd in 1:nrow(results)){
  n <- results[nd,1]
  d <- results[nd,2]
  #write model functions
  write_kss_model(n, d)
  #source model functions
  source(paste0("eq_", n, "_", d, ".r"))
  #find solutions
  search_sol <- function(x0){
    tmp <- fsolve(model, x0)
    return(-sum(tmp$x[tmp$x < 0]))
  }
  
  # for each system, try finding a feasible solution
  prob <- 0
  for (j in 1:500){
    B <<- random_tensor(d, n)
    tmp <- list(value = 10)
    for(i in 1:10){
      if(tmp$value > 0) tmp <- optim(par = rnorm(n), fn = search_sol)
    }
    if (tmp$value == 0) {
      print("success")
    }
    prob <- prob + (tmp$value == 0)
  }
  results[nd,3] = prob/500
}
