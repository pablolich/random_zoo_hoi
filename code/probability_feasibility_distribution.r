library("tidyverse")

#load data
data = read.csv("../data/expected_n_roots_dim_6_div_6_s_1000_normal.csv", sep = "", 
                header = T )

names(data) = c('d', 'n', 'sim', 'nsol', 'npos', 'stable')

n_sim = 1000

df = data %>% group_by(d, n) %>% 
  count(positive, name = 'count') %>% 
  mutate(total_eq = sum(positive*count),
         mean_eq = total_eq/n_sim,
         var_eq = sum(count*(positive - mean_eq)^2)/n_sim)

ggplot(df, aes(x = n, y = mean_eq))+
  geom_point(aes(color = as.factor(d)))+
  geom_function(fun = expected_eq, args = list(d = 2),
                lty = 2)+
  geom_function(fun = expected_eq, args = list(d = 3))+
  geom_function(fun = expected_eq, args = list(d = 4))+
  geom_function(fun = expected_eq, args = list(d = 5))+
  geom_function(fun = expected_eq, args = list(d = 6))+
  theme(aspect.ratio = 1)
