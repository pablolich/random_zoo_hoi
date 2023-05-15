rm(list = ls())
library(tidyverse)
source("count_bertini.R")

time_and_result <- function(n, d, nsim){
  z <- microbenchmark::microbenchmark(
    A <- t(replicate(n = nsim, count_solutions(n, d))),
    times = 1
  )
  return(list(A, z$time / 10^9))
}



success <- 0
for (d in c(2,3,4,5,6)){
  for (n in c(3,4,5)){
    res <- tibble()
    timing <- tibble()
    while(nrow(res) < 5000){
      tmp <- NULL
      try({tmp <- time_and_result(n, d, 250)})
      if (!is.null(tmp)){
        res <- rbind(res, tmp[[1]])
        timing <- rbind(timing, tibble(n = n, d = d, time = tmp[[2]]))
      }
      print(nrow(res))
    }
    write_csv(res, file = paste0("n_", n, "_d_", d, ".csv"))
    write_csv(timing, file = paste0("n_", n, "_d_", d, "_times.csv"))
  }
}

