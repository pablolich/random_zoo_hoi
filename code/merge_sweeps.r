library(tidyverse)

#merge all sweeps into one file
system("cat *.csv > ../merged_sweeps.csv")
#load the file

parameter_sweeps_results_merged <- read.table("../data/merged_sweeps.csv", 
                                               sep = " ", 
                                               header = F)

colnames(parameter_sweeps_results_merged) = c("d", "n", "nsol", "npos")
pfeas_results = parameter_sweeps_results_merged %>% 
  group_by(d, n) %>% 
  summarise(pfeas = mean(npos>0)) %>% 
  select(c(pfeas, n, d))


#save data
write.table(pfeas_results %>% select(pfeas, n, d), '../data/pfeas_sims_merged.dat', sep = " ", row.names = F,
            quote = F)
