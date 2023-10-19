model <- function(x, pars){
  return(pars[[2]] + pars[[1]]%*%x)
}

sample_parameters <- function(n, nneigh){
  r <- rnorm(n)# + sign(rnorm(1))
  A = matrix(rnorm(n*n), n, n)
  return(list(A, r))
}

#set maximum diversity
maxn = 2
#set model
model = dxdt_type1
#name of parametrization
parameters_label = "randomzoo"
#name of file to save
output_name = paste0("glv", parameters_label)
#set tolerance for determining if a number is 0
tol = 1e-6
#set number of simulations for each diversity value
nsim = 2000
#set random seed
set.seed(1)

#initialize data.frame for storing
results = data.frame("n"=seq(2, maxn), "pf"=0)#rep(0, maxn-1))
colnames(results)<-c("n", "pf")

source("search_feasible_greedy.r")