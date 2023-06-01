setwd("~/Desktop/random_zoo_hoi/code")

big_file = tibble()
nsim = 1000
for (n in c(1,2,3,4,5,6,7,8)){
  for (d in c(2, 3, 4, 5, 6)){
    f1 <- paste0("../data/check_kss_variance/n_", n, "_d_", d, ".csv")
    if (file.exists(f1)){
      dt <- read.table(f1, sep = "\t", header = T)
      ###this has change for next versions \t
      colnames(dt) = c("n", "d", "nsol", "npos")
      if ((d-1) %% 2  == 0){
        dt_test = dt %>% filter(nsol %% 2 == 0)
      } else{
        dt = dt %>% filter(nsol %% 2 != 0)
      }
      big_file = rbind(big_file, dt)
    }
  }
}
  
data = big_file %>% group_by(n, d) %>% count(nsol, name = "tot_nsol") %>% 
  mutate(mean_nsol = sum(nsol*tot_nsol)/nsim)
