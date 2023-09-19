library(tidyverse)

integration_results <- read.table("../data/kss_simulations.csv", 
                                       sep = " ", 
                                       header = F)

colnames(integration_results) = c("n", "d", "nf")

dat = integration_results %>%
  group_by(n, d) %>% 
  summarize(nf_av = mean(nf))

ggplot(data = dat,
       aes(x = n, y = nf_av))+
  geom_point()