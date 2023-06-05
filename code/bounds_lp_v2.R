library(lpSolve)

get_bounds <- function(n, d){
  max_num_pos_roots <- d^n
  mean_num_pos_roots <- d^(n/2) / 2^n
  # if there is the right file, compute var
  var <- NULL
  pfeas <- -1
  path = "../data/merged_files/"
  fn <- paste0(path, "n_", n, "_d_", d + 1, ".csv")
  if (file.exists(fn)){
    dt <- read.csv(fn, sep = " ")
    mean_num_pos_roots <- mean(dt$npos)
    var <- mean((dt$npos - mean_num_pos_roots)^2)
    third <- mean((dt$npos - mean_num_pos_roots)^3)
    fourth <- mean((dt$npos - mean_num_pos_roots)^4)
    pfeas <- mean(dt$npos > 0)
  }
  # coeff function
  coeff_func <- c(0, rep(1, max_num_pos_roots))
  # coeff for simplex
  coeff_simplex <- rep(1, (max_num_pos_roots + 1))
  # coefficients for mean
  coeff_mean <- 0:max_num_pos_roots
  # maximize/minimize coefficients^T p 
  # constraints
  # 1^T p = 1
  f.obj <- coeff_func
  f.con <- matrix(c(coeff_simplex, coeff_mean), nrow=2, byrow=TRUE)
  f.dir <- c("==", "==")
  f.rhs <- c(1, mean_num_pos_roots)
  out_max <- lp("max", f.obj, f.con, f.dir, f.rhs)
  out_min <- lp("min", f.obj, f.con, f.dir, f.rhs)
  out_max2 <- list(objval = -1)
  out_min2 <- list(objval = -1)
  if (!is.null(var)){
    coeff_var <- (coeff_mean - mean_num_pos_roots)^2
    coeff_skew <- (coeff_mean - mean_num_pos_roots)^3
    coeff_kurt <- (coeff_mean - mean_num_pos_roots)^4
    f.obj <- coeff_func
    f.con <- matrix(c(coeff_simplex, 
                      coeff_mean, 
                      coeff_var, 
                      coeff_skew,
                      coeff_kurt
                      ), nrow=5, byrow=TRUE)
    f.dir <- c("==", "==", "==", "==", "==")
    f.rhs <- c(1, mean_num_pos_roots, 
               var, 
               third,
               fourth
               )
    out_max2 <- lp("max", f.obj, f.con, f.dir, f.rhs)
    out_min2 <- lp("min", f.obj, f.con, f.dir, f.rhs)
  } 
  return(c(pfeas, out_max$objval, out_min$objval, out_max2$objval, out_min2$objval))  
}

library(tidyverse)
results <- tibble()
for (n in 1:8){
  for(d in 1:6){
    tmp <- get_bounds(n, d)
    results <- rbind(results, tibble(n = n, 
                                     d = d, 
                                     p_feas = tmp[1], 
                                     max_feas1 = tmp[2],
                                     max_feas2 = tmp[4],
                                     min_feas1 = tmp[3],
                                     min_feas2 = tmp[5]))
    print(results[nrow(results),])
  }
}
results <- results %>% mutate(p_feas = ifelse(p_feas > 0, p_feas, NA),
                              #max_feas2 = ifelse(max_feas2 <= 1, max_feas2, NA),
                              min_feas2 = ifelse(min_feas2 < 0, NA, min_feas2))
write_csv(results, file = "bounds_lp_v2.csv")

results <- results %>% drop_na()

# plot
toplot <- results %>% pivot_longer(names_to = "bound", values_to = "pfeas",cols = -c(n, d))

show(ggplot(toplot, aes(x = n, y = pfeas, colour = bound)) + geom_point() + geom_line() + facet_wrap(~d))

     