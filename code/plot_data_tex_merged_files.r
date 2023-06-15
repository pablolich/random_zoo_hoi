expected_value_expansion = function(n, d, var, skw, order){
  if (order == 1){
    approx = (1 - 1/2^n)^(d^(n/2))
  } else if (order == 2){
    approx = (1 - 1/2^n)^(d^(n/2))*(1+1/2*log(1-1/2^n)^2*var)
    
  } else{
    approx = (1 - 1/2^n)^(d^(n/2))*(1+1/2*log(1-1/2^n)^2*var+1/6*log(1-1/2^n)^3*skw)
  }
  return(approx)
} 

results <- tibble()
for (n in c(1,2,3,4,5,6,7,8)){
  for (d in c(2,3,4,5,6,7)){
    f1 <-  paste0("../data/growing_file/n_", n, "_d_", d, ".csv")
    if (file.exists(f1)){
      dt <- read.table(f1, sep = " ", header = T)
      #add names to columns
      colnames(dt) = c("n", "d", "nsol", "npos")
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
      #calculate by taylor expanding the expected value
      exptaylor1 = 1 - expected_value_expansion(n, d-1, var_real, skw_real, 1)
      exptaylor2 = 1 - expected_value_expansion(n, d-1, var_real, skw_real, 2)
      exptaylor3 = 1 - expected_value_expansion(n, d-1, var_real, skw_real, 3)
      pfeas <- mean(dt$npos>0)
      results <- rbind(results, tibble(
        n = n, 
        d = d, 
        nsol_av = mean(dt$nsol),
        npos_av = mean(dt$npos),
        pfeas = pfeas,
        pfeas_exp1 = exptaylor1,
        pfeas_exp2 = exptaylor2,
        pfeas_exp3 = exptaylor3,
        var_real = var_real,
        skw_real = skw_real
      ))
    }
  }
}

dat0 = results %>% select(c(pfeas, n, d))
dat1 = results %>% select(c(pfeas, pfeas_exp1, d))
dat2 = results %>% select(c(pfeas, pfeas_exp2, d))
dat3 = results %>% select(c(pfeas, pfeas_exp3, d))
dat4 = results %>% select(c(pfeas_exp2, n, d))
dat5 = results %>% mutate(dd = d-1) %>% select(c(n, nsol_av, dd)) %>% rename(d = dd)
dat6 = results %>% mutate(dd = d-1) %>% select(c(n, npos_av, dd)) %>% rename(d = dd)

#save data
write.table(dat0, '../data/pfeas_sims.dat', sep = " ", row.names = F,
            quote = F)
write.table(dat1, '../data/pfeas_approx1.dat', sep = " ", row.names = F,
            quote=F)
write.table(dat2, '../data/pfeas_approx2.dat', sep = " ", row.names = F,
            quote=F)
write.table(dat3, '../data/pfeas_approx3.dat', sep = " ", row.names = F,
            quote=F)
write.table(dat4, '../data/pfeas_approx2_n.dat', sep = " ", row.names = F,
            quote=F)
write.table(dat5, '../data/exp_eq.dat', sep = " ", row.names = F,
            quote=F)
write.table(dat6, '../data/exp_eq_pos.dat', sep = " ", row.names = F,
            quote=F)