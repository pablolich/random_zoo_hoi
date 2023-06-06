library(cowplot)
results <- tibble()
big_file = tibble()
for (n in c(1,2,3,4,5,6,7,8)){
  for (d in c(1,2,3,4,5,6)){
    f1 <-  paste0("../data/merged_files/n_", n, "_d_", d+1, ".csv")
    if (file.exists(f1)){
      dt <- read.table(f1, sep = " ", header = T)
      if (d %% 2  == 0){
        dt_test = dt %>% filter(nsol %% 2 == 0)
      } else{
        dt = dt %>% filter(nsol %% 2 != 0)
      }
      dt$d = d
      p0pos <- mean(dt$npos==0) #probability of no positive roots
      nsolav=mean(dt$nsol)
      results <- rbind(results, tibble(
        n = n, 
        d = d, 
        nsolav=nsolav,
        p0pos = p0pos,
        indep_assump = (1-1/2^n)^((d)^(n/2))
      ))
      big_file = rbind(big_file, dt)
    }
  }
}

ggplot(results, aes(x=n, y = log(p0pos/indep_assump), 
                    color = as.factor(d)))+
  geom_point()+
  geom_line()+
  facet_wrap(~d)

ggplot(big_file, aes(x=nsol, y=npos))+
  scale_fill_brewer(palette = "Spectral", direction = -1)+
  geom_point(aes(fill = as.factor(n)), 
             shape=21, color="grey")+
  facet_wrap(~d)
