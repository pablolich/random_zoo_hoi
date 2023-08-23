#!/usr/bin/env Rscript

library(tidyverse)


dat = read.table("../data/times.csv")
colnames(dat) = c("n", "d", "t_tensor", "t_poly")

dat_av = dat %>% group_by(n, d) %>%
	      summarize(t_tensor_av = mean(t_tensor),
		     t_poly_av = mean(t_poly)) %>% 
  	      pivot_longer(!c(n, d), names_to = "method",
			   values_to = "times")%>%
	      filter(n>1)

pdf("../data/times.pdf")
p = ggplot(dat_av, aes(x = n, y = times))+
    geom_point(aes(shape = as.factor(method), 
		   color = as.factor(d)),
    	       size = 2)+
    scale_y_continuous(trans = "log10")+
    scale_shape_manual(values=c(16,3))
plot(p)
dev.off()


