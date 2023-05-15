# Get Bertini from here
# https://bertini.nd.edu/download.html
# move it to /usr/local/bin and make runnable
# From CL:
# sudo cp bertini /usr/local/bin
# sudo chmod 777 /usr/local/bin/bertini


# Install algstat stuff:
# if (!requireNamespace("devtools")) install.packages("devtools")
# devtools::install_github("dkahle/mpoly")
# devtools::install_github("coneill-math/m2r")
# devtools::install_github("dkahle/latte")
# devtools::install_github("dkahle/bertini")
# devtools::install_github("dkahle/algstat")

# Load bertini library
library(bertini)
# Set the path to the program
set_bertini_path("/usr/local/bin/")

next_indices <- function(indices, n){
  if (all(indices == n)) return(NULL)
  indices[1] <- indices[1] + 1
  stop <- FALSE
  tocheck <- 1
  while(!stop){
    if (indices[tocheck] > n){
      # if possible, add to the next
      if (tocheck < length(indices)){
        indices[tocheck + 1] <- indices[tocheck + 1] + 1
        indices[tocheck] <- 1
        tocheck <- tocheck + 1
      } else {
        return(indices)
      }
    } else {
      return(indices)
    }
  }
  return(-1)
}

build_poly_dg2 <- function(Bi, X){
  pp <- character(0)
  dims <- dim(Bi)
  n <- dims[1]
  indices <- rep(1, length(dims))
  pp <- paste(as.numeric(R.utils::extract.array(Bi,indices =  indices)), paste(X[indices], collapse = " "))
  while(!is.null(indices)){
    indices <- next_indices(indices, n)
    if (!is.null(indices)){
      pp <- paste(pp, "+", as.numeric(R.utils::extract.array(Bi,indices =  indices)), paste(X[indices], collapse = " "))  
    }
  }
  return(pp)
}

build_linear <- function(Bi, X){
  return(paste(paste(Bi, X), collapse = " + "))
}

# n is the actual number of species
# the tensor has size n+1 in each of the d dimensions
count_solutions <- function(n,d){
  # If one wants to truncate the random numbers
  B <- array(round(rnorm((n+1)^d), 2), dim = rep((n+1), d))
  # Full number of digits
  #B <- array(rnorm((n+1)^d), dim = rep((n+1), d))
  X <- paste("x", 1:(n+1), sep = "")
  X[n + 1] <- 1
  p <- list()
  for (i in 1:(n)) {
    if (d == 2) p[[i]] <- build_linear(B[i,], X)
    if (d == 3) p[[i]] <- build_poly_dg2(B[i,,], X)
    if (d == 4) p[[i]] <- build_poly_dg2(B[i,,,], X)
    if (d == 5) p[[i]] <- build_poly_dg2(B[i,,,,], X)
    if (d == 6) p[[i]] <- build_poly_dg2(B[i,,,,,], X)
  }
  out <- bertini(bertini_input(mp(p)))$real_finite_solutions
  #print(out)
  if (all(out == FALSE)){
    nsol <- 0
    npos <- 0
  } else {
    nsol <- nrow(out)
    npos <- 0
    if (nsol == 1) npos <- all(out > 0) * 1
    if (nsol > 1) npos <- sum(apply(out, 1, function(x) all(x > 0)))
  }
  return(c(n = n, d = d, nsol = nsol, npos = npos ))
}
