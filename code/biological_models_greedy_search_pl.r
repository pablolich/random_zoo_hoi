library(lhs)

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
  return(r + as.vector((A / denom) %*% x))
}

search_finely <- function(x0, tol){
  n = length(x0)
  resolution = 10*n
  perturbations = randomLHS(resolution, n)
  x0_grid = x0*(1 + 0.1*perturbations)
  for (i in 1:resolution){
    print(paste0("Searching finely...  ", i, " of ", resolution))
    x0_pert = x0_grid[i,]
    tmp = rootSolve::multiroot(dxdt, x0, rtol = tol)
    if (tmp$estim.precis<tol){
      return(tmp)
    }
  }
  return(tmp)
}

search_sol <- function(x0, tol, niterations){
  tmp <- rootSolve::multiroot(dxdt, x0, rtol = tol)
  # #expand grid around problematic point to search exhaustively around hard areas
  # if (tmp$estim.precis>tol){
  #   tmp = search_finely(x0, tol)
  # }
  #did the search converge?
  if (tmp$estim.precis > tol){
    return(-1)
  }
  else{
    #is the solution feasible?
    if (all(tmp$root > 0)){ 
      #run search with higher precision to make sure that solution is positive
      tmp  = rootSolve::multiroot(dxdt, tmp$root, rtol = 1e-6*tol)
      if (all(tmp$root > 0)){
        return(1)
      }
      else{return(0)}
      }
    else{ return(0) }
  }
}

sample_x0 <- function(n){
  #create a grid of initial conditions uniformly distributed in the simplex
  init_conds_simplex = randomLHS(1, n+1)
  #project back to Rn
  init_conds = init_conds_simplex/init_conds_simplex[,n+1]
  init_conds = init_conds[,1:n]
  return(init_conds)
}

get_max_nneigh <- function(n){
  return((n-1)%/%2)
}

circ <- function(x) {
  n <- length(x)
  suppressWarnings(
    matrix(x[matrix(1:n,n+1,n+1,byrow=T)[c(1,n:2),1:n]],n,n))
}

build_struct_A <- function(n, nneigh){
  #get vector
  struct_vec = c(-1, rep(1, nneigh), rep(0, n-1-2*nneigh), rep(-1, nneigh))
  A_struct = circ(struct_vec)
  return(A_struct)
}

sample_parameters <- function(n, nneigh){
  r <- rnorm(n) + sign(rnorm(1))
  Avals <- matrix(abs(rnorm(n * n)), n, n)
  Asigns = build_struct_A(n, nneigh)
  A = Avals*Asigns
  return(list(r, A))
}

#set maximum diversity
maxn = 5
#set tolerance for determining if a number is 0
tol = 1e-6
#set number of simulations for each diversity value
nsim = 1000
#set random seed
set.seed(1)
#initialize data.frame for storing
results = data.frame("n"=seq(5, maxn), "pf"=rep(0, maxn-1))
colnames(results)<-c("n", "pf")

#find roots for each n
for (i in 1:nrow(results)){
  #set number of species
  n = results[i,1]
  #set structure to be fully connected
  nneigh = get_max_nneigh(n)
  #set maximum number of successful trials depending on n
  ntrialsmax = 40*n
  #set maximum number of failed searches before trying a different parameter set
  nfailsmax = ntrialsmax/10
  #set number of searches resulting in feasible solutions to 0
  nfeasible = 0
  #reset simulation number back to 0
  sim=0
  #number of initial conditions
  ninit = 40*n
  
  #perform many searches for a given n, and many parameter sets, 
  #to get a good estimate of Pf
  while (sim < nsim){
    #sample new set of parameters for each simulation
    pars = sample_parameters(n, nneigh)
    #parse out each of them
    r <<- pars[[1]]
    A <<- pars[[2]]
    diag(A) <- -abs(diag(A)) - 1
    Anod <- A # without diagonal
    diag(Anod) <- 0
    #A[upper.tri(A)] <- 0
    #A[lower.tri(A)] <- 0
    # Split the matrix
    Ap <<- Anod * (A > 0) # positive
    An <<- Anod * (A < 0) # negative
    Ad <<- diag(diag(A)) # diagonal
    #set feasiblility flag for this parameter set to 0
    feasible = 0
    #reset number of, trials, and failed trials to 0
    ntrials = 0
    nfails = 0
    
    #attempt to find root from many initial conditions
    while (ntrials < ntrialsmax){
      #sample initial condition
      x0 = sample_x0(n)
      #return whether solution is feasible, not feasible, or not solution
      tmp = search_sol(x0, tol)
      if (tmp == 1){
        #solution converged, and negative part is small enough, stop searching
        feasible = 1
        print(paste0("Diversity: ", n,
                     " attempt #: ", ntrials, " in simulation ", sim, 
                     " is feasible!"))
        break
      }
      else if (tmp == 0){
        #solution converged, but with negative components
        ntrials = ntrials + 1
        print(paste0("Diversity: ", n,
                     " attempt #: ", ntrials, " in simulation ", sim, 
                     " is not feasible"))
      }
      else {
        #solution did not converged, try again
        nfails = nfails + 1
        print(paste0("Diversity: ", n,
                     " attempt #: ", ntrials, " in simulation ", sim, 
                     " did not converge"))
      }
      if (nfails > nfailsmax){
        #when enough runs have not converged, try a different parameterization
        sim = sim - 1
        print(paste0("Diversity: ", n,  
                     " attempt #: ", ntrials, " in simulation ", sim, 
                     " maximum number of fails reached, try a different",
                     "parameter set"))
        break
      }
    }
    #once search is finished, add the result to running total
    nfeasible = nfeasible + feasible
    #increase simulation counter
    sim = sim + 1
  }
  #once all simulations are completed for a given n, estimate probability of 
  #feasibility numerically
  results[i,2] = nfeasible/nsim
}
