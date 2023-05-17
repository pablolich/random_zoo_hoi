library("tidyverse")

even_binomial = function(x, size, prob){
  if (x %% 2 == 0){
    #x is even
    return(dbinom(x, size, prob))
  }
  else {
    return(0)
  }
}


#load stefanos data
results <- tibble()

for (n in c(3,4,5)){
  for (d in c(2,3,4,5,6)){
    f1 <-  paste0("n_", n, "_d_", d, ".csv")
    if (file.exists(f1)){
      dt <- read_csv(f1)
      n <- dt$n[1]
      d <- dt$d[1]
      obs <- mean(dt$npos)
      expected <- 1/(2^n) * (d-1)^(n / 2)
      exppfeas <- 1 - (1 - 1/(2^n))^((d-1)^(n / 2))
      exppois <- 1- (exp(-2^(-n) * (d-1)^(n/2)))
      expbinom <- 1 - ((1 + (-2^n + 4^n*(-1 + d)^(n/2))^(-1))^(-1 + d)^n*(1 - 1/(2^n*(-1 + d)^(n/2)))^(-1 + d)^n)
      #vector with possible number of solutions
      n_sols = seq(0, (d-1)^n)
      pfeas <- mean(dt$npos > 0)
      results <- rbind(results, tibble(
        n = n, 
        d = d, 
        observed = obs,
        expected = expected,
        pfeas = pfeas,
        exppfeas = exppfeas,
        exppois = exppois,
        expbinom = expbinom
      ))
    }
  }
}
