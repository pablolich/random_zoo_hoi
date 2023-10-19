#set maximum diversity
maxn = 2
#set model
model = dxdt_type1
#name of parametrization
parameters_label = "standard"
#name of file to save
output_name = paste0("dxdt_type1", parameters_label)
#set tolerance for determining if a number is 0
tol = 1e-6
#set number of simulations for each diversity value
nsim = 2000
#set random seed
set.seed(1)
#set structure to be fully connected
nneigh = get_max_nneigh(n)
#nneigh =0
#initialize data.frame for storing
results = data.frame("n"=seq(2, maxn), "pf"=0)#rep(0, maxn-1))
colnames(results)<-c("n", "pf")

dxdt_type1 <- function(x, pars){
   # simply evaluate the equations
  ## \dot{x}_i = x_i(r_i + \sum_j aij xj / dj)
  ## where if A_{ij} >0, i != j, dj = 1 + \sum{k} A_{ik} x_k over the Aik > 0
  ## if A_{ij} < 0, i != j, dj = 1 + \sum{k} A_{kj} x_k over the Akj > 0
  ## if i = j, then dj = 1
  #make convenient parameters out each of them
  diag(pars[[1]]) <- -abs(diag(pars[[1]])) - 1
  Anod <- pars[[1]] # without diagonal
  diag(Anod) <- 0
  #A[upper.tri(A)] <- 0
  #A[lower.tri(A)] <- 0
  # Split the matrix
  Ap <- Anod * (pars[[1]] > 0) # positive
  An <- Anod * (pars[[1]] < 0) # negative
  Ad <- diag(diag(pars[[1]])) # diagonal
  n=length(x)
  # build matrix of denominators
  denom <- 1 + as.vector(Ap %*% abs(x))
  denom_row <- denom %o% rep(1, n)
  denom_col <- rep(1, n) %o% denom
  denom <- (Ap > 0) * denom_row + (An < 0) * denom_col + diag(rep(1, n))
  denom[denom == 0] <- 1 # if no prey
  # equations
  return(pars[[2]] + as.vector((pars[[1]] / denom) %*% x))
}

