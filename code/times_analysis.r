#!/usr/bin/env Rscript

library(tidyverse)


dat = read.table("../data/times_deformed.csv")
colnames(dat) = c("n", "d", "t_deformed", "t_total")

dat_av = dat %>% group_by(n, d) %>%
	      summarize(t_deformed_av = mean(t_deformed),
		     t_total_av = mean(t_total)) %>% 
  	      pivot_longer(!c(n, d), names_to = "method",
			   values_to = "times")%>%
	      filter(n>1)

pdf("../data/times_deformed.pdf")
p = ggplot(dat_av, aes(x = n, y = times))+
    geom_point(aes(shape = as.factor(method), 
		   color = as.factor(d)),
    	       size = 2)+
    theme(aspect.ratio=0.7)+
    scale_y_continuous(trans = "log10")+
    scale_shape_manual(values=c(16,3))
plot(p)
dev.off()


