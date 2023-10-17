library(pracma)
library(tidyverse)
library(ggplot2)
library(lhs)
library(profvis)
library(rootSolve)



#find solutions
search_sol_multi <- function(x0){ #time different algorithms
  tmp <- rootSolve::multiroot(model, x0, rtol = 1e-6)
  #check value of function at the supposed root
  if (sum(tmp$f.root^2) > tol){
    return(1)
  }
  else return(-sum(tmp$root[tmp$root < 0]))
}

search_sol_fsolve <- function(x0){ #time different algorithms
  tmp <- fsolve(model, x0, tol = 1e-6)
  #check value of function at the supposed root
  if (sum(tmp$fval^2) > tol){
    return(1)
  }
  else return(-sum(tmp$x[tmp$x < 0]))
}


#script to time multiroot vs fsolve algorithms
#test
maxd = 2#6
maxn = 4
tol = 1e-6
nsearches = 100
nattemptsmax = nsearches*2
ninit = 10
resultsfsolve = expand.grid(2:maxn, 1:maxd)
resultsmultiroot = expand.grid(2:maxn, 1:maxd)
#add extra column
resultsfsolve[,3] = rep(0,nrow(resultsfsolve))
resultsfmultiroot[,3] = rep(0,nrow(resultsmultiroot))
colnames(resultsfsolve)<-c("n", "d", "pf")
colnames(resultsmultiroot)<-c("n", "d", "pf")



ptm <- proc.time()
for (nd in 1:nrow(resultsfsolve)){
  n <- resultsfsolve[nd,1]
  d <- resultsfsolve[nd,2]
  #write model functions
  write_kss_model(n, d)
  #source model functions
  source(paste0("eq_", n, "_", d, ".r"))
  # for each system, try finding a feasible solution
  prob <- 0
  nsearched = 0
  nattempts = 0
  niterations = 100
  while (nsearched<nsearches){
    B <<- random_tensor(d, n)
    tmp <- list(value = 10)
    #create a grid of initial conditions uniformly distributed in the simplex
    init_conds_simplex = randomLHS(ninit*n, n+1)
    #project back to Rn
    init_conds = init_conds_simplex/init_conds_simplex[,n+1]
    init_conds = init_conds[,1:n]
    convergence_history = rep(0, ninit*n)
    for(i in 1:(ninit*n)){
      tmp <- optim(par = init_conds[i,], fn = search_sol_multi)
      convergence_history[i] = tmp$value
      if (tmp$value == 0){
        print(paste0("Diversity: ", n, " Interaction order: ", d,  
                     " Search #: ", nsearched, " was succesful after ", 
                     nattempts, " attempts."))
        break
      }
    }
    #only count if at least one of the searches has converged
    if (any(convergence_history != 1)){
      nsearched = nsearched + 1
      prob <- prob + (tmp$value == 0)
      nattempts = 0 #reset number of attempts, since convergence has been reached
    }
    nattempts = nattempts + 1
    if (nattempts>nattemptsmax){
      niterations = niterations*2
      print(paste0("Minimizer is not able to converge; increasing number of iterations to: ", 
            niterations))
    }
  }
  resultsfsolve[nd,3] = prob/nsearches
}
tfsolve = proc.time() - ptm


ggplot(resultsfsolve, aes(x= n, y=pf))+
  geom_point(aes(color = as.factor(d)))


ptm <- proc.time()
for (nd in 1:nrow(resultsmultiroot)){
  n <- resultsmultiroot[nd,1]
  d <- resultsmultiroot[nd,2]
  #write model functions
  write_kss_model(n, d)
  #source model functions
  source(paste0("eq_", n, "_", d, ".r"))
  # for each system, try finding a feasible solution
  prob <- 0
  nsearches = 10
  ninit = 10
  nsearched = 0
  while (nsearched<nsearches){
    B <<- random_tensor(d, n)
    tmp <- list(value = 10)
    #create a grid of initial conditions uniformly distributed in the simplex
    init_conds_simplex = randomLHS(ninit, n+1)
    #project back to Rn
    init_conds = init_conds_simplex/init_conds_simplex[,n+1]
    init_conds = init_conds[,1:n]
    convergence_history = rep(0, ninit)
    for(i in 1:(ninit)){
      tmp <- optim(par = init_conds[i,], fn = search_sol_)
      convergence_history[i] = tmp$value
      if (tmp$value == 0){
        print(paste0("Diversity: ", n, " Interaction order: ", d,  
                     " Search #: ", nsearched, " was succesful"))
        break
      }
    }
    #only count if at least one of the searches has converged
    if (any(convergence_history != 1)){
      nsearched = nsearched + 1
      prob <- prob + (tmp$value == 0)
    }
  }
  resultsmultiroot[nd,3] = prob/nsearches
}
tmultiroot = proc.time() - ptm


ggplot(resultsmultiroot, aes(x= n, y=pf))+
  geom_point(aes(color = as.factor(d)))
