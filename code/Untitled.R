nvariables = function(n){
  return(n*(n+1))
}

nequations = function(n){
  neqs = 0
  for (i in 1:n){
    neqs = neqs + choose(i+n-1, n-1)
  }
  return(neqs)
}

check = function(n){
  return(2*factorial(2*n-1)/(factorial(n)*factorial(n-1))-1)
}

check2 = function(n){
  return(choose(2*n, n))
}

nmax = 7
nvar = c()
neq = c()
guess = c()

for (i in 1:nmax){
  nvari = nvariables(i)
  neqi = nequations(i)
  nvar = c(nvar, nvari)
  neq = c(neq, neqi)
  guess = c(guess, i^i)
}

plot(1:nmax, nvar)
lines(1:nmax, neq)
lines(1:nmax, guess)