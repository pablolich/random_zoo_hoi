data =  read.table("../data/p_feas_min.csv", header = F)
colnames(data) = c("n", "d", "pfeas_min")
data = data %>% relocate(pfeas_min, .before = d)

results_merged = merge(results, data, by = c("n", "d"))

ggplot(results_merged)+
  geom_point(aes(x = n, y = pfeas_min), 
             shape = 24,
             fill = "black")+
  geom_line(aes(x = n, y = pfeas_min, 
                group = as.factor(d),
                color = as.factor(d)),
            linetype = "dashed")+
  geom_point(aes(x = n, y = pfeas_exp1),
             shape = 25,
             fill = "black")+
  geom_line(aes(x = n, y = pfeas_exp1, 
                group = as.factor(d),
                color = as.factor(d)),
            linetype = "solid")+
  geom_point(aes(x = n, y = pfeas, 
                 color =as.factor(d)))+
  geom_line(aes(x = n, y = pfeas, 
                group = as.factor(d),
                color = as.factor(d)))

write.table(data, '../data/pfeas_min.dat', sep = " ", row.names = F,
            quote=F)
