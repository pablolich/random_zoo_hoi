results = tibble()
for (n in seq(1:8)){
  for (d in seq(1:6)){
    results = rbind(results, 
                    tibble(n=n,
                           d=d,
                           comp = d^n))
  }
}