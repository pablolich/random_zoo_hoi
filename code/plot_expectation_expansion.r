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

results <- tibble()
for (n in c(1,2,3,4,5,6,7)){
  for (d in c(2,3,4,5,6)){
    #for (d in c(3,5)){
    f1 <-  paste0("n_", n, "_d_", d, ".csv")
    if (file.exists(f1)){
      dt <- read.table(f1, sep = "\t", header = T)
      if (length(dt) == 1){
        dt <- read.table(f1, sep = ",", header = T)
      }
      colnames(dt) = c("n", "d", "nsol", "npos")
      n <- dt$n[1]
      d <- dt$d[1]
      if ((d-1) %% 2  == 0){
        dt_test = dt %>% filter(nsol %% 2 == 0)
      } else{
        dt = dt %>% filter(nsol %% 2 != 0)
      }
      var_real = mean((dt$nsol)^2) - (mean(dt$nsol))^2
      if (var_real == 0){
        skw_real = 0
      } else{
        skw_real = mean(((dt$nsol - mean(dt$nsol))/sqrt(var_real))^3)
      }
      #jensen
      #calculate by taylor expanding the expected value
      exptaylor1 = 1 - expected_value_expansion(n, d-1, var_real, skw_real, 1)
      exptaylor2 = 1 - expected_value_expansion(n, d-1, var_real, skw_real, 2)
      exptaylor3 = 1 - expected_value_expansion(n, d-1, var_real, skw_real, 3)
      pfeas <- mean(dt$npos>0)
      results <- rbind(results, tibble(
        n = n, 
        d = d, 
        pfeas = pfeas,
        pfeas_exp1 = exptaylor1,
        pfeas_exp2 = exptaylor2,
        pfeas_exp3 = exptaylor3,
      ))
    }
  }
}

dat0 = results %>% select(c(pfeas, n, d))
dat1 = results %>% select(c(pfeas, pfeas_exp1, d))
dat2 = results %>% select(c(pfeas, pfeas_exp2, d))
dat3 = results %>% select(c(pfeas, pfeas_exp3, d))
dat4 = results %>% select(c(pfeas_exp2, n, d))
dat5 = results %>% select(c(n, d))

#save data
write.table(dat0, '../../data/pfeas_sims.dat', sep = " ", row.names = F,
            quote = F)
write.table(dat1, '../../data/pfeas_approx1.dat', sep = " ", row.names = F,
            quote=F)
write.table(dat2, '../../data/pfeas_approx2.dat', sep = " ", row.names = F,
            ,quote=F)
write.table(dat3, '../../data/pfeas_approx3.dat', sep = " ", row.names = F,
            quote=F)
write.table(dat4, '../../data/pfeas_approx2_n.dat', sep = " ", row.names = F,
            quote=F)
