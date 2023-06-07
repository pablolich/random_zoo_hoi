library(tidyverse)
# all files 
allf <- list.files("../../data/merged_files/", pattern = "n")

pfeas_binom <- function(n, d){
  return(1 - (1 - (2^-n) * d^(-n/2) )^(d^n))
}

results <- tibble()
for (ff in allf){
  dt <- read_delim(paste0("../../data/merged_files/", ff))
  n <- dt$n[1]
  d <- dt$d[1] - 1
  pfeas <- mean(dt$npos > 0)
  p_binom <- pfeas_binom(n, d)
  results <- rbind(results, tibble(n = n, d = d, pfeas = pfeas, p_binom = p_binom))
}