library("tidyverse")
library("utils")

#setwd("/Users/pablolechon/Desktop/phd/random_zoo_hoi/code/Rcode") #mac

even_odd_binomial = function(n, p, k){
  if (n %% 2 == 0){
    #even distributiv = on
    prob_even = ((2*(1 - p)^(-2*k + n)*p^(2*k)*choose(n, 2*k))/(1 + (1 - 2*p)^n))
    prob = as.vector(rbind(prob_even, 0))
  } else{
    #odd distribution
    prob_odd = ((-2*(1 - p)^(-1 - 2*k + n)*p^(1 + 2*k)*choose(n, 2*k + 1))/(-1 + (1 - 2*p)^n))
    prob = as.vector(rbind(0, prob_odd))
  }
  return(prob)
}

even_odd_binomial_2 = function(n, p){
  if (n %% 2 == 0){
    #sequece of points to evaluate
    k = 2*seq(0, n/2)
    #even distributiv = on
    prob_even = dbinom(k, n, p)
    normalization = sum(prob_even)
    prob = as.vector(rbind(prob_even/normalization, 0))
    #eliminate last extra zero
    prob = head(prob, -1)
  } else{
    #sequence of points to evaluate
    k = 2*seq(0, (n-1)/2)+1
    #odd distribution
    prob_odd = dbinom(k, n, p)
    normalization = sum(prob_odd)
    prob = as.vector(rbind(0, prob_odd/normalization))
  }
  return(prob)
}

expected_value = function(probabilities, k){
  return(sum(probabilities*k)) 
}

fun = function(p, n, d, k){
  p_vec = even_odd_binomial_2(d^n, p)
  E = expected_value(p_vec, k)
  return(E - sqrt(d)^n)
}

f_x = function(n, i){
  return((1-1/2^n)^i)
}

f2_x = function(n, i){
  return(log(1-1/2^n)^2*(1-1/2^n)^(i))
}

f3_x = function(n, i){
  return((1 - 2^(-n))^i*log(1 - 2^(-n))^3)
}

expected_value_expansion = function(n, d, var, skw, order){
  mean_x = sqrt(d^n)
  #evaluate function and second derivative at mean
  f_mean = f_x(n, mean_x)
  f2_mean = f2_x(n, mean_x)
  f3_mean = f2_x(n, mean_x)
  
  if (order == 1){
    approx = f_mean 
  } else if (order == 2){
    approx = f_mean + 1/2*f2_mean*var 
    
  } else{
    approx = f_mean + 1/2*f2_mean*var + 1/6*f3_mean*skw
  }
  return(approx)
} 

var_curve = function(v_inf, n, d){
  return(v_inf*d^(n/2))
}

merged = T
#load stefanos data
results <- tibble()

for (n in c(3,4,5)){
  for (d in c(2,3,4,5,6)){
  #for (d in c(3,5)){
    f1 <-  paste0("n_", n, "_d_", d, ".csv")
    if (file.exists(f1)){
      dt <- read_csv(f1, col_types = cols())
      n <- dt$n[1]
      d <- dt$d[1]
      if (n == 5 & d == 4){
        if (merged){
          #load more simulations
          more_sim = read.csv("../../data/expected_n_roots_dim_4_div_5_s_5000_normal_stabtrue.csv",
                              sep = "\t")
          #clean for merge
          names(more_sim) = c('d', 'n', 'sim', 'nsol', 'npos', 'max_eig')
          cleaned = more_sim %>% group_by(d, n, sim) %>% 
            slice_min(max_eig) %>% 
            ungroup() %>% 
            select(c(d, n, nsol, npos)) %>% 
            relocate(n, .before = d)
          #merge
          dt = rbind(dt, cleaned)
        } else{
          dt <- read_csv(f1, col_types = cols())
        }
      }
      if ((d-1) %% 2  == 0){
        dt_test = dt %>% filter(nsol %% 2 == 0)
      } else{
        dt = dt %>% filter(nsol %% 2 != 0)
      }
      obs <- mean(dt$npos)
      obs_real = mean(dt$nsol)
      dt$npos = dt$npos/obs_real
      #calculate experimental variance of the simulations
      var_pos_t = mean((dt$npos)^2) - (mean(dt$npos))^2
      var_pos = sd(dt$npos)^2
      var_real = mean((dt$nsol)^2) - (mean(dt$nsol))^2
      if (var_real == 0){
        skw_real = 0
      } else{
        skw_real = mean(((dt$nsol - obs_real)/sqrt(var_real))^3)
      }
      print(skw_real)
      expected <- 1/(2^n) * (d-1)^(n / 2)
      exppfeas <- 1 - (1 - 1/(2^n))^((d-1)^(n / 2))
      exppois <- 1- exp(-((-1 + d)^(n/2)/2^n))
      expbinom <- 1 - (1 + (-1 + 2^n)/(2^n*(-1 + (-1 + d)^(n/2))))^(-1 + d)^n*(1 - (-1 + d)^(-1/2*n))^(-1 + d)^n
      #vector of quantiles
      k_vec = seq(0, (d-1)^n)
      #get p for even binomial matching expected value of real roots 
      if (d > 2){
        p_res = uniroot(fun, c(0.00001, 0.99999), tol = 0.00001, n = n, 
                        d = d-1, k = k_vec)
        p = p_res$root
      } else {
        p = 1
      }
      probs = even_odd_binomial_2((d-1)^n, p)
      expbinom_even_odd = 1 - sum(probs*(1-1/(2^n))^(k_vec))
      #calculate by taylor expanding the expected value
      exptaylor = 1 - expected_value_expansion(n, d-1, var_real, skw_real, 3)
      var_pred = var_curve(0.2832072, n, d-1) #hard code Vinf
      v_inf = var_real/sqrt(d-1)^n
      pfeas <- mean(dt$npos > 0)
      results <- rbind(results, tibble(
        n = n, 
        d = d, 
        obs_mean_pos = obs,
        obs_mean_real = obs_real,
        obs_variance_pos_t = var_pos_t,
        obs_variance_pos = var_pos,
        obs_variance_real = var_real,
        expected = expected,
        pfeas = pfeas,
        error = abs(pfeas-exptaylor),
        exppfeas = exppfeas,
        exppois = exppois,
        expbinom = expbinom,
        expbinom_even = expbinom_even_odd,
        exptaylor = exptaylor,
        var_pred = var_pred,
        v_inf = v_inf
      ))
    }
  }
}

par(mfcol = c(1, 2))
plot(results$obs_mean_pos, results$expected, main = "# equil"); abline(c(0,1))
plot(results$pfeas, results$exppfeas, main = "p feas 10000", pch = 20); abline(c(0,1))
points(results$pfeas, results$expbinom_even, col = "red", pch = 20)
points(results$pfeas, results$exppois, col = "blue", pch = 4)
points(results$pfeas, results$expbinom, col = "green", pch = 0)
points(results$pfeas, results$exptaylor, col = "purple", pch = 3)
legend(.01, .8, legend=c("Jensen", "Poisson", "Binomial", "Binomial (even or odd)",
                         "Taylor"),
       col=c("black", "blue", "green", "red", "purple"), pch = c(20, 4, 0, 20, 3), 
       cex = 0.8)



###quick check of histograms 
# 
# n = 5
# d = 4
# par(mfcol = c(1, 2))
# merged = T
# f1 <-  paste0("n_", n, "_d_", d, ".csv")
# if (merged){
#   #load more simulations
#   more_sim = read.csv("../../data/expected_n_roots_dim_4_div_5_s_5000_normal_stabtrue.csv",
#                       sep = "\t")
#   #clean for merge
#   names(more_sim) = c('d', 'n', 'sim', 'nsol', 'npos', 'max_eig')
#   cleaned = more_sim %>% group_by(d, n, sim) %>% 
#     slice_min(max_eig) %>% 
#     ungroup() %>% 
#     select(c(d, n, nsol, npos)) %>% 
#     relocate(n, .before = d)
#   #merge
#   dt = rbind(dt, cleaned)
# } else{
#   dt <- read_csv(f1, col_types = cols())
# }
# #get probability of real roots
# k_vec = seq(0, (d-1)^n)
# #get p for even binomial matching expected value of real roots 
# if (d > 2){
#   p_res = uniroot(fun, c(0.00001, 0.99999), tol = 0.00001, n = n, 
#                   d = d-1, k = k_vec)
#   p = p_res$root
# } else {
#   p = 1
# }
# probs = even_odd_binomial_2((d-1)^n, p)
# points = probs*nrow(dt)
# hist(dt$nsol, breaks = seq(max(dt$nsol)), 
#      main = "p_real 10000")
# points(seq(max(dt$nsol))+0.5, points[seq(max(dt$nsol))])
# 
# obs_mean = mean(dt$nsol)
# dist_mean = sum(k_vec*probs)
