library(tidyverse)

parameter_sweeps_glv_func_results <- read.table("../data/parameter_sweeps_glvfunc.csv", 
                                       sep = " ", 
                                       header = F)
colnames(parameter_sweeps_glv_func_results) = c("n", "c", "nsol", "npos", "nonlin")

#data wrang
toplot = parameter_sweeps_glv_func_results %>% 
  group_by(n, c, nonlin) %>% 
  slice_max(npos, with_ties = F) %>% 
  summarise(nsol_av = mean(nsol),
            npos_av = mean(npos))

#save
write.table(toplot %>% ungroup() %>% filter(nonlin==1) %>% select(npos_av, n), '../data/n_feasib_sols.dat', sep = " ", row.names = F,
            quote = F)
write.table(toplot %>% ungroup() %>% filter(nonlin==0) %>% select(npos_av, n), '../data/n_feasib_sols_control.dat', sep = " ", row.names = F,
            quote = F)
#plot
ggplot(toplot) +
  geom_point(aes(x = n, y=npos_av,
                 color = as.factor(nonlin)))+
  geom_line(aes(x = n, y = npos_av, 
                group = nonlin,
                color = as.factor(nonlin)))



##############################################################################

#Now look at eigenvlaues
e
eigenvalues <- read.table("../data/eigenvalues.csv", sep = " ", header = T)
colnames(eigenvalues) = c("n", "c", "num", "real", "imaginary")

eigenvals_big = eigenvalues %>% filter(n==15, c<3)


ggplot(eigenvals_big)+
  geom_point(aes(x=real, y = imaginary,
                 color = as.factor(num)),
             alpha = 0.2)+
  facet_wrap(~n+c, scales = "free")

 