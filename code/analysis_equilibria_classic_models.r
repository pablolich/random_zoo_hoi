library(tidyverse)
dev.off()

data = read.table("../data/neq_glv_func_resp.csv", sep = ",", header = T)

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

data_cr = read.table("../data/neq_cr_.csv", sep = ",", header = T)

data_cr %>% 
  group_by(m, n) %>% 
  summarise(avnsol = mean(nsol), 
            avnpos = mean(npos)) %>% 
  pivot_longer(c(avnsol, avnpos), 
               names_to = "solution",
               values_to = "average") %>% 
  ggplot(aes(x=n, y=average))+
  geom_line(aes(group =  as.factor(m),
                color =  as.factor(m)))+
  geom_point(aes(color =  as.factor(m)))+
  theme(aspect.ratio = 1)+
  facet_wrap(~solution)

data_constr = read.table("../data/neq_glv_func_resp_constr.csv", sep = ",", header = T)

data_constr %>% 
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

data_sym = read.table("../data/neq_glv_func_resp_symm.csv", sep = ",", header = T)

data_constr %>% 
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



