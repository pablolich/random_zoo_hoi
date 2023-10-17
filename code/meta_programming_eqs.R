library(pracma)
library(tidyverse)
library(ggplot2)
library(lhs)
library(rootSolve)

write_kss_model <- function(n, d){
  sink(paste0("eq_", n, "_", d, ".r"))
  #write function to evaluate equation i of the model
  cat("Fi <- function(Bi, x){\n")
  #meta-comment
  cat("#evaluate model at given parameters B (globally), and abundances x\n")
  cat("x = flatten_dbl(list(x))\n")
  cat("xmod = c(x, 1)\n")
  outerx = "xmod"
  if (d>1){
    outerx <- paste0(outerx,paste0(rep("%o%xmod",d-1), collapse = ""))
  }
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


#script to time multiroot vs fsolve algorithms
#test
maxd = 6#6
maxn = 8
tol = 1e-6
nsim = 5000
nattemptsmax = nsim*2
results = expand.grid(2:maxn, 4:4)
#add extra column
results[,3] = rep(0,nrow(results))
colnames(results)<-c("n", "d", "pf")
ptm <- proc.time()
for (nd in 1:nrow(results)){
  n <- results[nd,1]
  d <- results[nd,2]
  #write model functions
  write_kss_model(n, d)
  #source model functions
  source(paste0("eq_", n, "_", d, ".r"))
  #find solutions
  search_sol <- function(x0){ #time different algorithms
    tmp <- rootSolve::multiroot(model, x0, rtol = tol, maxiter=niterations)
    #check value of function at the supposed root
    if (sum(tmp$f.root^2) > tol){
      return(1)
    }
    else return(-sum(tmp$root[tmp$root < 0]))
  }
  # for each system, try finding a feasible solution
  prob <- 0
  sim = 0
  nattempts = 0
  niterations <<- 100
  ninit = 40*n
  while (sim<nsim){
    B <<- random_tensor(d, n)
    tmp <- list(value = 10)
    #create a grid of initial conditions uniformly distributed in the simplex
    init_conds_simplex = randomLHS(ninit, n+1)
    #project back to Rn
    init_conds = init_conds_simplex/init_conds_simplex[,n+1]
    init_conds = init_conds[,1:n]
    convergence_history = rep(0, ninit)
    for(i in 1:(ninit)){
      #tmp <- optim(par = init_conds[i,], fn = search_sol)
      tmp = search_sol(init_conds[i,])
      convergence_history[i] = tmp
      #convergence_history[i] = tmp$value
      if (tmp == 0){
        print(paste0("Diversity: ", n, " Interaction order: ", d,  
                     " attempt #: ", i, " in simulation ", sim, 
                     " was succesful"))
        break
      }
      else if (tmp != 1){
        print(paste0("Diversity: ", n, " Interaction order: ", d,  
                     " attempt #: ", i, " in simulation ", sim, 
                     " was not succesful"))
      }
    }
    #only count if at least one of the searches has converged
    if (any(convergence_history != 1)){
      sim = sim + 1
      prob <- prob + (tmp == 0)
      nattempts = 0 #reset number of attempts, since convergence has been reached
    }
    else print(paste0("Minimizer did not converge, trial n: ", nattempts))
    nattempts = nattempts + 1
    if (nattempts>nattemptsmax){
      niterations <<- niterations*2
      print(paste0("multiroot didn't converge; increasing number of iterations to: ", 
                   niterations))
    }
  }
  results[nd,3] = prob/nsim
  #plot
  p<-ggplot(results, aes(x= n, y=pf))+
    geom_line(aes(group = as.factor(d),
                  color = as.factor(d)))+
    geom_point(aes(color = as.factor(d)))
  print(p)
}
tfsolve = proc.time() - ptm
