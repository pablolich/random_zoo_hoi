#for d = 4, and different n
d = 5
big_file = tibble()
for (n in c(0,1,2,3,4,5,6,7,8)){
  f1 <-  paste0("n_", n, "_d_", d, ".csv")
  if (file.exists(f1)){
    dt <- read.table(f1, sep = "\t", header = T)
    if (length(dt) == 1){
      dt <- read.table(f1, sep = ",", header = T)
    }
    colnames(dt) = c("n", "d", "nsol", "npos")
    if ((d-1) %% 2  == 0){
      dt_test = dt %>% filter(nsol %% 2 == 0)
    } else{
      dt = dt %>% filter(nsol %% 2 != 0)
    }
    big_file = rbind(big_file, dt)
  }
}

count_dat = big_file %>%group_by(d, n) %>% 
  count(npos, name = "count_pos")

ggplot(data = count_dat,
       aes(x = npos, y = count_pos)) +
  geom_bar(stat="identity")+
  scale_x_discrete(limit =c("0", "1", "2", "3", "4", "5", "6", "7"),
                   breaks = seq(0, 7))+
  facet_wrap(~n, nrow = 1)
  