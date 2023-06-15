parameter_sweeps_results <- read.table("../data/parameter_sweeps.csv", 
                                       sep = " ", 
                                       header = T)
#add names to columns
colnames(parameter_sweeps_results) = c("d", "n", "nsol", "npos")
pfeas_results = parameter_sweeps_results %>% 
  group_by(d, n) %>% 
  summarise(pfeas = mean(npos>0)) %>% 
  select(c(pfeas, n, d))
#save data
write.table(dat0, '../data/pfeas_sims.dat', sep = " ", row.names = F,
            quote = F)