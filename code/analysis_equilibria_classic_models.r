library(tidyverse)
data = read.table("../data/neq_classic_models.csv", sep = ",", header = T)

data %>% 
  group_by(n) %>% 
  summarise(avnsol = mean(nsol), avnpos = mean(npos)) %>%
  pivot_longer(c(avnsol, avnpos), 
               names_to = "solution",
               values_to = "average") %>% 
  ggplot(aes(x=n, y=average))+
  geom_line(aes(group = solution), 
            color = "grey")+
  geom_point(aes(color = solution, group = solution))+
  theme(aspect.ratio = 1)


data_jac = read.table("../data/neq_classic_models_jac.csv", sep = ",", header = T)

data_jac %>% 
  group_by(n) %>% 
  summarise(avnsol = mean(nsol), avnpos = mean(npos)) %>%
  pivot_longer(c(avnsol, avnpos), 
               names_to = "solution",
               values_to = "average") %>% 
  ggplot(aes(x=n, y=average))+
  geom_line(aes(group = solution), 
            color = "grey")+
  geom_point(aes(color = solution, group = solution))+
  theme(aspect.ratio = 1)

data_rand = read.table("../data/neq_classic_models_rand.csv", sep = ",", header = T)

data_rand %>% 
  group_by(n) %>% 
  summarise(avnsol = mean(nsol), avnpos = mean(npos)) %>%
  pivot_longer(c(avnsol, avnpos), 
               names_to = "solution",
               values_to = "average") %>% 
  ggplot(aes(x=n, y=average))+
  geom_line(aes(group = solution), 
            color = "grey")+
  geom_point(aes(color = solution, group = solution))+
  theme(aspect.ratio = 1)

