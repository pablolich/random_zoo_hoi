data = read.table("../data/cr_feasibility.csv", header = T)
colnames(data) = c("n", "m", "nsol", "npos")

dt_feas = data %>% 
  group_by(n, m) %>% 
  summarise(pfeas = mean(npos>0))
