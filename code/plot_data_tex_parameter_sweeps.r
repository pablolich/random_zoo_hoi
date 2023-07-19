library(tidyverse)

parameter_sweeps_results <- read.table("../data/parameter_sweeps.csv", 
                                       sep = " ", 
                                       header = F)
#add names to columns
colnames(parameter_sweeps_results) = c("d", "n", "nsol", "npos")
pfeas_results = parameter_sweeps_results %>% 
  group_by(d, n) %>% 
  summarise(pfeas = mean(npos>0)) %>% 
  select(c(pfeas, n, d)) %>% 
  mutate(pfeas_jensen = 1-(1-1/2^n)^(sqrt(d^n)))



#save data
write.table(pfeas_results %>% select(pfeas, n, d), '../data/pfeas_sims.dat', sep = " ", row.names = F,
            quote = F)
write.table(pfeas_results %>% select(pfeas, pfeas_jensen, n, d), '../data/pfeas_jensen_sims.dat', sep = " ", row.names = F,
            quote = F)
