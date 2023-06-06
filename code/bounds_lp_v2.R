  library(lpSolve)
  library(tidyverse)
  
  lp_optimization<-function(max_num_pos_roots, mean_num_pos_roots, var, opt_dir, dt){
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
    out_opt <- lp(opt_dir, f.obj, f.con, f.dir, f.rhs)
    if (!is.null(var)){
      coeff_var <- (coeff_mean - mean_num_pos_roots)^2
      #coeff_skew <- (coeff_mean - mean_num_pos_roots)^3
      #coeff_kurt <- (coeff_mean - mean_num_pos_roots)^4
      f.obj <- coeff_func
      f.con <- matrix(c(coeff_simplex,
                        coeff_mean,
                        coeff_var
                        ),
                      nrow=3, byrow=TRUE)
      f.dir <- c("==", "==", "==")
      f.rhs <- c(1, mean_num_pos_roots,
                 var)
      out_opt <- lp(opt_dir, f.obj, f.con, f.dir, f.rhs)
    }
    return(out_opt)
  }
  
  get_ps = function(n, d, dt){
    max_num_pos_roots <- d^n
    mean_num_pos_roots <- d^(n/2) / 2^n
    # if there is the right file, compute var
    dt <- read.csv(fn, sep = " ")
    mean_num_pos_roots <- mean(dt$npos)
    var <- mean((dt$npos - mean_num_pos_roots)^2)
    pfeas <- mean(dt$npos > 0)
    out_max = lp_optimization(max_num_pos_roots, mean_num_pos_roots, var, "min")
    return(rbind(out_max$solution, out_max$constraints[1:(d^n+1),2]))
  }
  
  get_bounds <- function(n, d, dt){
    max_num_pos_roots <- d^n
    mean_num_pos_roots <- d^(n/2) / 2^n
    dt <- read.csv(fn, sep = " ")
    mean_num_pos_roots <- mean(dt$npos)
    var <- mean((dt$npos - mean_num_pos_roots)^2)
    pfeas <- mean(dt$npos > 0)
    out_max = lp_optimization(max_num_pos_roots, mean_num_pos_roots, var, "max")
    out_min = lp_optimization(max_num_pos_roots, mean_num_pos_roots, var, "min")
    return(c(pfeas, out_max$objval, out_min$objval))  
  }
  
  library(tidyverse)
  #preallocate results for p
  results <- tibble()
  p_mat = tibble()
  for (n in 1:8){
    for(d in 1:6){
      path = "../data/merged_files/"
      fn <- paste0(path, "n_", n, "_d_", d + 1, ".csv")
      if (file.exists(fn)){
        dt =  read.csv(fn, sep = " ")
        tmp <- get_bounds(n, d, dt)
        results <- rbind(results, tibble(n = n, 
                                         d = d, 
                                         p_feas = tmp[1], 
                                         max_feas = tmp[2],
                                         min_feas = tmp[3]))
        sol = get_ps(n, d, dt)
        p_vec = sol[1,]
        n_sols = sol[2,]
        #append 0s
        #p_vec_long = c(p_vec, rep(0, dmax^nmax-length(p_vec)))
        p_mat = rbind(p_mat, 
                      tibble(n=n, 
                             d=d, 
                             npos = n_sols,
                             p_vec = p_vec))
        print(results[nrow(results),])
      } else{
        next
      }
    }
  }
  results <- results %>% mutate(p_feas = ifelse(p_feas > 0, p_feas, NA),
                                min_feas = ifelse(min_feas < 0, NA, min_feas))
  write_csv(results, file = "bounds_lp_v2.csv")
  
  results <- results %>% drop_na()
  
  # plot
  toplot <- results %>% pivot_longer(names_to = "bound", values_to = "pfeas",cols = -c(n, d))
  
  show(ggplot(toplot, aes(x = n, y = pfeas, colour = bound)) + geom_point() + geom_line() + facet_wrap(~d))
  
  ggplot(p_mat %>% filter(p_vec>1e-9))+
    geom_col(aes(x = npos, y=p_vec)) +
    facet_wrap(n~d, scales="free")
